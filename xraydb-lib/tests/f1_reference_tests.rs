//! Validates `f1_chantler` against reference values from upstream XrayDB.
//!
//! The fixture in `tests/data/f1_chantler_reference.csv` was produced once with
//! upstream's own Python implementation (`uv run --with xraydb`); nothing here needs
//! Python. See the fixture header for how the reference values are defined.

use xraydb::XrayDb;

struct Case {
    element: String,
    energy: f64,
    f1: f64,
}

fn load_reference() -> Vec<Case> {
    let raw = include_str!("data/f1_chantler_reference.csv");
    raw.lines()
        .filter(|l| !l.trim_start().starts_with('#') && !l.trim().is_empty())
        .filter_map(|l| {
            let mut parts = l.split(',');
            let element = parts.next()?.to_string();
            let energy = parts.next()?.parse().ok()?;
            let f1 = parts.next()?.parse().ok()?;
            Some(Case {
                element,
                energy,
                f1,
            })
        })
        .collect()
}

/// The headline claim: the global cubic spline reproduces upstream's wide-window f1.
#[test]
fn f1_matches_upstream_across_the_tabulated_range() {
    let db = XrayDb::new();
    let cases = load_reference();
    assert!(
        cases.len() > 1500,
        "fixture looks truncated: {}",
        cases.len()
    );

    let mut worst_abs = 0.0_f64;
    let mut worst: Option<&Case> = None;

    for case in &cases {
        let got = db
            .f1_chantler_at(&case.element, case.energy)
            .unwrap_or_else(|e| panic!("{} @ {}: {e}", case.element, case.energy));
        assert!(
            got.is_finite(),
            "{} @ {}: got {got}",
            case.element,
            case.energy
        );
        let diff = (got - case.f1).abs();
        if diff > worst_abs {
            worst_abs = diff;
            worst = Some(case);
        }
    }

    if let Some(c) = worst {
        println!(
            "worst f1 deviation: {worst_abs:.3e} at {} @ {} eV (reference {:+.6})",
            c.element, c.energy, c.f1
        );
    }
    // Upstream's own values are only defined to ~1e-10 in the fixture's text form.
    assert!(
        worst_abs < 1e-6,
        "worst deviation {worst_abs:.3e} exceeds 1e-6"
    );
}

/// The improvement this change delivers, stated as a test so it cannot silently regress.
///
/// Linear interpolation on the same grid is reconstructed here and must be dramatically
/// worse than the spline, especially at absorption edges.
#[test]
fn the_spline_beats_linear_interpolation_by_a_wide_margin() {
    let db = XrayDb::new();
    let cases = load_reference();

    let mut spline_worst = 0.0_f64;
    let mut linear_worst = 0.0_f64;
    let mut spline_sum = 0.0_f64;
    let mut linear_sum = 0.0_f64;

    for case in &cases {
        let row = db
            .raw()
            .chantler
            .iter()
            .find(|c| c.element == case.element)
            .expect("element has Chantler data");

        // Reconstruct the old linear interpolation.
        let e = case.energy;
        let idx = row
            .energy
            .partition_point(|&v| v < e)
            .clamp(1, row.energy.len() - 1);
        let (lo, hi) = (idx - 1, idx);
        let span = row.energy[hi] - row.energy[lo];
        let linear = if span > 0.0 {
            let t = (e - row.energy[lo]) / span;
            row.f1[lo] + t * (row.f1[hi] - row.f1[lo])
        } else {
            row.f1[hi]
        };

        let spline = db.f1_chantler_at(&case.element, e).expect("f1");

        spline_worst = spline_worst.max((spline - case.f1).abs());
        linear_worst = linear_worst.max((linear - case.f1).abs());
        spline_sum += (spline - case.f1).abs();
        linear_sum += (linear - case.f1).abs();
    }

    let n = cases.len() as f64;
    println!(
        "mean |error|  spline {:.3e}  linear {:.3e}  ({:.0}x better)",
        spline_sum / n,
        linear_sum / n,
        (linear_sum / spline_sum.max(f64::MIN_POSITIVE)).min(1e9)
    );
    println!("worst |error| spline {spline_worst:.3e}  linear {linear_worst:.3e}");

    assert!(
        linear_worst > 0.1,
        "expected linear interpolation to be visibly wrong somewhere, got {linear_worst:.3e}"
    );
    assert!(
        spline_worst < linear_worst / 1000.0,
        "spline {spline_worst:.3e} should be far better than linear {linear_worst:.3e}"
    );
}

/// Absorption edges are where linear interpolation was worst, and where f1 matters most
/// (anomalous scattering, MAD phasing). Pin a few explicitly.
#[test]
fn f1_is_accurate_at_absorption_edges() {
    let db = XrayDb::new();
    // (element, edge energy in eV, upstream wide-window f1)
    let edges: &[(&str, f64, f64)] = &[
        ("Fe", 7112.0, -9.186206),
        ("Cu", 8979.0, -8.892794),
        ("Au", 11919.0, -17.769546),
        ("Zn", 9659.0, -8.785750),
        ("U", 17166.0, -13.993460),
    ];
    for &(el, energy, want) in edges {
        let got = db.f1_chantler_at(el, energy).expect("f1");
        assert!(
            (got - want).abs() < 1e-5,
            "{el} @ {energy}: got {got:+.6}, want {want:+.6}"
        );
    }
}

/// Caesium's grid repeats 11.4 eV, which makes upstream raise
/// `ValueError: x must be strictly increasing`. We fit the spline to the strictly
/// increasing subsequence instead, so every element stays queryable.
#[test]
fn caesium_is_queryable_despite_its_duplicated_grid_point() {
    let db = XrayDb::new();
    let (emin, emax) = db.chantler_energy_range("Cs").expect("Cs range");

    for e in [emin, 11.4, 100.0, 5_000.0, 50_000.0, emax] {
        let got = db.f1_chantler_at("Cs", e).expect("Cs f1");
        assert!(got.is_finite(), "Cs @ {e}: got {got}");
    }

    // Away from the duplicate, upstream does work and we should agree with it.
    let got = db.f1_chantler_at("Cs", 50_000.0).expect("Cs f1");
    assert!(
        (got - (-0.352786)).abs() < 1e-3,
        "Cs @ 50 keV: got {got:+.6}"
    );
}

/// Every element must produce finite f1 across its whole tabulated range.
#[test]
fn f1_is_finite_for_every_element() {
    let db = XrayDb::new();
    for row in &db.raw().chantler {
        let (emin, emax) = db.chantler_energy_range(&row.element).expect("range");
        let steps = 200;
        for i in 0..=steps {
            let t = f64::from(i) / f64::from(steps);
            let e = (emin.ln() + (emax.ln() - emin.ln()) * t).exp();
            let got = db.f1_chantler_at(&row.element, e).expect("f1");
            assert!(got.is_finite(), "{} @ {e}: got {got}", row.element);
        }
    }
}

/// The spline must pass exactly through the tabulated points.
#[test]
fn f1_interpolates_the_knots_exactly() {
    let db = XrayDb::new();
    for element in ["Fe", "Si", "Au", "U"] {
        let row = db
            .raw()
            .chantler
            .iter()
            .find(|c| c.element == element)
            .expect("has data");

        for (i, (&e, &want)) in row.energy.iter().zip(row.f1.iter()).enumerate() {
            if i == 0 || i + 1 == row.energy.len() || e > 1e6 {
                continue;
            }
            // Skip repeated abscissae, where a spline cannot represent the jump.
            if row.energy[i - 1] == e || row.energy[i + 1] == e {
                continue;
            }
            let got = db.f1_chantler_at(element, e).expect("f1");
            assert!(
                (got - want).abs() < 1e-9 * want.abs().max(1.0),
                "{element} knot {i} @ {e}: got {got:+.9}, want {want:+.9}"
            );
        }
    }
}

/// Scalar and batch forms must agree exactly.
#[test]
fn scalar_and_batch_f1_agree() {
    let db = XrayDb::new();
    let energies = [1_000.0, 7_112.0, 10_000.0, 55_555.0, 300_000.0];
    let batch = db.f1_chantler("Fe", &energies).expect("batch");
    for (i, &e) in energies.iter().enumerate() {
        assert_eq!(db.f1_chantler_at("Fe", e).expect("scalar"), batch[i]);
    }
}
