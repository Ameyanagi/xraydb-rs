//! Regression tests for bugs found in the 0.1.2 review.
//!
//! Each test here corresponds to a defect that shipped in 0.1.2; they exist to keep
//! those specific behaviours from coming back.

use xraydb::{CrossSectionKind, XrayDb};

/// The materials table contains an `air` entry whose formula the 0.1.2 parser
/// (`chemical-formula`) rejected, so `material_mu_named("air", …)` always failed.
#[test]
fn air_is_usable_as_a_material() {
    let db = XrayDb::new();
    let mu = db
        .material_mu_named("air", &[10_000.0], CrossSectionKind::Total, None)
        .expect("air must be usable");
    assert!(mu[0].is_finite() && mu[0] > 0.0, "mu = {}", mu[0]);
}

/// The ion-chamber path is the main consumer of the `air` entry, and inherited the
/// same failure.
#[test]
fn ionchamber_accepts_air() {
    let db = XrayDb::new();
    let fluxes = db
        .ionchamber_fluxes(&[("air", 1.0)], 1.0, 10.0, 10_000.0, 1e-6, true, true)
        .expect("air ion chamber must work");
    assert!(fluxes.incident > 0.0);
    assert!(fluxes.transmitted < fluxes.incident);
}

/// Deuterium was accepted by `xray_delta_beta` but rejected by `material_mu`.
#[test]
fn deuterium_works_in_every_formula_path() {
    let db = XrayDb::new();
    assert!(
        db.material_mu("D2O", 1.1, &[10_000.0], CrossSectionKind::Total)
            .is_ok()
    );
    assert!(db.xray_delta_beta("D2O", 1.1, 10_000.0).is_ok());

    // D is counted as H, so heavy water and light water have the same composition.
    let heavy = db.mass_fractions("D2O").unwrap();
    let light = db.mass_fractions("H2O").unwrap();
    assert_eq!(heavy, light);
}

/// Scientific-notation stoichiometry was accepted by one parser and not the other.
#[test]
fn scientific_notation_stoichiometry_works_in_every_path() {
    let db = XrayDb::new();
    assert!(
        db.material_mu("Zn1.e-5Fe3O4", 5.2, &[10_000.0], CrossSectionKind::Total)
            .is_ok()
    );
    assert!(db.xray_delta_beta("Zn1.e-5Fe3O4", 5.2, 10_000.0).is_ok());
}

/// The two public formula entry points must accept exactly the same language.
#[test]
fn material_mu_and_delta_beta_accept_the_same_formulas() {
    let db = XrayDb::new();
    let formulas: Vec<&str> = db.materials().iter().map(|m| m.formula).collect();

    for formula in
        formulas
            .iter()
            .copied()
            .chain(["D2O", "Zn1.e-5Fe3O4", "Fe.7Mg.3O", "Ru1wt%SiO2"])
    {
        let mu = db.material_mu(formula, 1.0, &[10_000.0], CrossSectionKind::Total);
        let db_ = db.xray_delta_beta(formula, 1.0, 10_000.0);
        assert_eq!(
            mu.is_ok(),
            db_.is_ok(),
            "'{formula}': material_mu ok={}, xray_delta_beta ok={}",
            mu.is_ok(),
            db_.is_ok()
        );
    }
}

/// `multilayer_reflectivity` computed `total_layers - 1` before checking that any
/// layers existed, underflowing on an empty stackup or zero periods.
#[cfg(feature = "optics")]
#[test]
fn multilayer_rejects_degenerate_stacks_instead_of_panicking() {
    use xraydb::MultilayerParams;

    let db = XrayDb::new();
    let theta = [0.005_f64];

    let empty = db.multilayer_reflectivity(MultilayerParams {
        stackup: &[],
        thickness: &[],
        density: &[],
        theta: &theta,
        ..Default::default()
    });
    assert!(empty.is_err(), "empty stackup must be an error");

    let zero_periods = db.multilayer_reflectivity(MultilayerParams {
        stackup: &["Pt"],
        thickness: &[100.0],
        density: &[21.45],
        theta: &theta,
        n_periods: 0,
        ..Default::default()
    });
    assert!(zero_periods.is_err(), "n_periods = 0 must be an error");
}

/// Unpolarized light has no Parratt recursion here; it used to reach an
/// `unreachable!()` guard only after doing work.
#[cfg(feature = "optics")]
#[test]
fn multilayer_rejects_unpolarized_up_front() {
    use xraydb::{MultilayerParams, Polarization};

    let db = XrayDb::new();
    let err = db.multilayer_reflectivity(MultilayerParams {
        stackup: &["Pt"],
        thickness: &[100.0],
        density: &[21.45],
        theta: &[0.005],
        polarization: Polarization::Unpolarized,
        ..Default::default()
    });
    assert!(err.is_err());
}

/// `darwin_width` used to return `Some` with empty arrays and zeroed widths when the
/// zeta grid degenerated, instead of signalling "no result" with `None`.
#[cfg(feature = "optics")]
#[test]
fn darwin_width_signals_degenerate_cases_with_none() {
    use xraydb::DarwinParams;

    let db = XrayDb::new();

    // Bragg condition unsatisfiable at 100 eV.
    let low = db
        .darwin_width(DarwinParams {
            energy: 100.0,
            ..Default::default()
        })
        .unwrap();
    assert!(low.is_none());

    // Any Some(_) result must carry a non-empty grid.
    let ok = db
        .darwin_width(DarwinParams {
            energy: 10_000.0,
            ..Default::default()
        })
        .unwrap()
        .expect("Si(111) at 10 keV diffracts");
    assert!(!ok.zeta.is_empty());
    assert!(!ok.intensity.is_empty());
    assert_eq!(ok.zeta.len(), ok.intensity.len());
    assert_eq!(ok.zeta.len(), ok.rocking_curve.len());
    assert!(ok.energy_width > 0.0);
}

/// Invalid hkl and unsupported crystals must be reported, not silently accepted.
#[cfg(feature = "optics")]
#[test]
fn darwin_width_validates_its_inputs() {
    use xraydb::DarwinParams;

    let db = XrayDb::new();
    for bad in [
        DarwinParams {
            hkl: (1, 1, 0),
            ..Default::default()
        }, // mixed parity
        DarwinParams {
            hkl: (0, 0, 0),
            ..Default::default()
        },
        DarwinParams {
            crystal: "NaCl",
            ..Default::default()
        },
        DarwinParams {
            energy: -1.0,
            ..Default::default()
        },
        DarwinParams {
            order: 0,
            ..Default::default()
        },
    ] {
        assert!(db.darwin_width(bad).is_err(), "expected error for {bad:?}");
    }
}

/// `guess_edge` used to answer confidently for any energy; a tolerance lets callers
/// reject nonsense.
#[test]
fn guess_edge_respects_its_tolerance() {
    let db = XrayDb::new();

    let exact = db.guess_edge(7112.0, None, Some(1.0)).unwrap();
    assert_eq!(exact.element, "Fe");
    assert_eq!(exact.edge, "K");
    assert_eq!(exact.difference, 0.0);

    // Far above every tabulated K edge.
    assert!(db.guess_edge(200_000.0, None, Some(500.0)).is_none());
    assert!(db.guess_edge(200_000.0, None, None).is_some());
}

/// The binary-search `guess_edge` must agree with an exhaustive scan.
#[test]
fn guess_edge_matches_a_brute_force_scan() {
    let db = XrayDb::new();
    let filter = ["K", "L3", "L2", "L1", "M5"];

    for query in [50.0, 500.0, 3_000.0, 7_112.0, 12_345.0, 60_000.0, 99_999.0] {
        let got = db.guess_edge(query, None, None).expect("some edge exists");

        let brute = db
            .raw()
            .xray_levels
            .iter()
            .filter(|l| l.absorption_edge > 0.0 && filter.contains(&l.iupac_symbol.as_str()))
            .min_by(|a, b| {
                (a.absorption_edge - query)
                    .abs()
                    .total_cmp(&(b.absorption_edge - query).abs())
            })
            .expect("brute force finds something");

        assert_eq!(
            got.energy, brute.absorption_edge,
            "query {query}: indexed {got:?} vs brute {}/{} @ {}",
            brute.element, brute.iupac_symbol, brute.absorption_edge
        );
    }
}

/// The precomputed log tables must produce the same numbers as transforming the
/// tables per call, which is what 0.1.2 did.
#[test]
fn chantler_precomputed_logs_match_the_naive_transform() {
    let db = XrayDb::new();

    for element in ["Fe", "Si", "Au", "H", "U"] {
        let row = db
            .raw()
            .chantler
            .iter()
            .find(|c| c.element == element)
            .expect("element has Chantler data");

        let (emin, emax) = db.chantler_energy_range(element).unwrap();
        let log_e: Vec<f64> = row.energy.iter().map(|v| v.ln()).collect();
        let log_f2: Vec<f64> = row
            .f2
            .iter()
            .map(|&v| if v.abs() < 1e-99 { 1e-99 } else { v }.ln())
            .collect();

        for energy in [1_000.0_f64, 7_112.0, 10_000.0, 25_000.0, 80_000.0] {
            let e = energy.clamp(emin, emax);
            let got = db.f2_chantler_at(element, energy).unwrap();

            // Reproduce the 0.1.2 path: transform the tables, then interpolate.
            let idx = log_e.partition_point(|&v| v < e.ln());
            let expected = if idx == 0 {
                log_f2[0].exp()
            } else if idx >= log_e.len() {
                log_f2[log_f2.len() - 1].exp()
            } else {
                let t = (e.ln() - log_e[idx - 1]) / (log_e[idx] - log_e[idx - 1]);
                (log_f2[idx - 1] + t * (log_f2[idx] - log_f2[idx - 1])).exp()
            };

            let rel = (got - expected).abs() / expected.max(1e-30);
            assert!(rel < 1e-9, "{element} @ {energy}: {got} vs {expected}");
        }
    }
}

/// Scalar `_at` methods must agree with the batch forms exactly.
#[test]
fn scalar_and_batch_apis_agree() {
    let db = XrayDb::new();
    let energies = [1_000.0, 7_112.0, 10_000.0, 55_555.0];

    for kind in [
        CrossSectionKind::Total,
        CrossSectionKind::Photo,
        CrossSectionKind::Coherent,
        CrossSectionKind::Incoherent,
    ] {
        let batch = db.mu_elam("Fe", &energies, kind).unwrap();
        for (i, &e) in energies.iter().enumerate() {
            let scalar = db.mu_elam_at("Fe", e, kind).unwrap();
            assert_eq!(scalar, batch[i], "mu_elam {kind:?} @ {e}");
        }
    }

    let batch = db.f1_chantler("Fe", &energies).unwrap();
    for (i, &e) in energies.iter().enumerate() {
        assert_eq!(db.f1_chantler_at("Fe", e).unwrap(), batch[i]);
    }

    let batch = db.f2_chantler("Fe", &energies).unwrap();
    for (i, &e) in energies.iter().enumerate() {
        assert_eq!(db.f2_chantler_at("Fe", e).unwrap(), batch[i]);
    }

    let batch = db
        .material_mu("SiO2", 2.2, &energies, CrossSectionKind::Total)
        .unwrap();
    for (i, &e) in energies.iter().enumerate() {
        let scalar = db
            .material_mu_at("SiO2", 2.2, e, CrossSectionKind::Total)
            .unwrap();
        assert!(
            (scalar - batch[i]).abs() <= 1e-12 * batch[i].abs(),
            "material_mu @ {e}: {scalar} vs {}",
            batch[i]
        );
    }
}

/// Parsing must be deterministic: the same formula must give bit-identical results
/// across repeated calls, which a `HashMap`-ordered summation did not guarantee.
#[test]
fn delta_beta_is_bit_reproducible() {
    let db = XrayDb::new();
    let first = db.xray_delta_beta("Mn(SO4)2(H2O)7", 2.1, 12_000.0).unwrap();
    for _ in 0..50 {
        let again = db.xray_delta_beta("Mn(SO4)2(H2O)7", 2.1, 12_000.0).unwrap();
        assert_eq!(first.delta.to_bits(), again.delta.to_bits());
        assert_eq!(first.beta.to_bits(), again.beta.to_bits());
    }
}

/// Every element with Chantler data must survive a scalar query at both ends of its
/// range without producing NaN.
#[test]
fn chantler_endpoints_are_finite_for_every_element() {
    let db = XrayDb::new();
    for row in &db.raw().chantler {
        let (emin, emax) = db.chantler_energy_range(&row.element).unwrap();
        for e in [emin, (emin + emax) / 2.0, emax] {
            let f2 = db.f2_chantler_at(&row.element, e).unwrap();
            assert!(f2.is_finite(), "{} f2 @ {e} = {f2}", row.element);
            let f1 = db.f1_chantler_at(&row.element, e).unwrap();
            assert!(f1.is_finite(), "{} f1 @ {e} = {f1}", row.element);
        }
    }
}

/// Every element in the database must produce finite Elam cross-sections across the
/// whole supported range.
#[test]
fn elam_is_finite_for_every_element_with_data() {
    let db = XrayDb::new();
    let energies: Vec<f64> = (0..40).map(|i| 100.0 * 1.3_f64.powi(i)).collect();

    for row in &db.raw().photoabsorption {
        let mu = db
            .mu_elam(&row.element, &energies, CrossSectionKind::Total)
            .unwrap();
        for (e, v) in energies.iter().zip(mu.iter()) {
            assert!(v.is_finite() && *v > 0.0, "{} @ {e} = {v}", row.element);
        }
    }
}
