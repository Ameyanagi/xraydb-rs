use approx::assert_relative_eq;
use xraydb::{CrossSectionKind, XrayDb, XrayDbError};

#[test]
fn test_material_mu_water() {
    let db = XrayDb::new();
    // Water at density 1.0 g/cm³
    let mu = db
        .material_mu("H2O", 1.0, &[10000.0], CrossSectionKind::Total)
        .unwrap();
    assert!(mu[0] > 0.0, "mu for water at 10 keV = {}", mu[0]);
    assert!(mu[0] < 100.0, "mu for water at 10 keV = {}", mu[0]);
}

#[test]
fn test_material_mu_sio2() {
    let db = XrayDb::new();
    let mu = db
        .material_mu("SiO2", 2.65, &[10000.0], CrossSectionKind::Total)
        .unwrap();
    assert!(mu[0] > 0.0);
}

#[test]
fn test_material_mu_multiple_energies() {
    let db = XrayDb::new();
    let energies = vec![5000.0, 10000.0, 20000.0];
    let mu = db
        .material_mu("Fe2O3", 5.24, &energies, CrossSectionKind::Total)
        .unwrap();
    assert_eq!(mu.len(), 3);
    for val in &mu {
        assert!(*val > 0.0);
    }
}

#[test]
fn test_material_mu_supports_weight_percent_formula() {
    let db = XrayDb::new();
    let energies = vec![8000.0, 10000.0, 15000.0];
    let mu = db
        .material_mu("Ru1wt%SiO2", 2.65, &energies, CrossSectionKind::Total)
        .unwrap();
    assert_eq!(mu.len(), energies.len());
    assert!(mu.iter().all(|v| v.is_finite() && *v > 0.0));
}

#[test]
fn test_material_mu_fe2o3_formula_matches_mass_fraction_path() {
    let db = XrayDb::new();
    let energies = vec![5000.0, 10000.0, 20000.0];
    let density = 5.24;

    let by_formula = db
        .material_mu("Fe2O3", density, &energies, CrossSectionKind::Total)
        .unwrap();

    let fe_mass = 2.0 * db.molar_mass("Fe").unwrap();
    let o_mass = 3.0 * db.molar_mass("O").unwrap();
    let total = fe_mass + o_mass;
    let by_fractions = db
        .material_mu_from_mass_fractions(
            &[("Fe", fe_mass / total), ("O", o_mass / total)],
            density,
            &energies,
            CrossSectionKind::Total,
        )
        .unwrap();

    for (lhs, rhs) in by_formula.iter().zip(by_fractions.iter()) {
        assert_relative_eq!(lhs, rhs, epsilon = 1e-12, max_relative = 1e-12);
    }
}

#[test]
fn test_material_mu_from_mass_fractions_normalizes_input() {
    let db = XrayDb::new();
    let energies = vec![5000.0, 10000.0, 20000.0];
    let density = 5.24;

    let mu_a = db
        .material_mu_from_mass_fractions(
            &[("Fe", 2.0), ("O", 3.0)],
            density,
            &energies,
            CrossSectionKind::Total,
        )
        .unwrap();
    let mu_b = db
        .material_mu_from_mass_fractions(
            &[("Fe", 20.0), ("O", 30.0)],
            density,
            &energies,
            CrossSectionKind::Total,
        )
        .unwrap();

    for (lhs, rhs) in mu_a.iter().zip(mu_b.iter()) {
        assert_relative_eq!(lhs, rhs, epsilon = 1e-12, max_relative = 1e-12);
    }
}

#[test]
fn test_material_mu_density_validation() {
    let db = XrayDb::new();
    for density in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let err = db
            .material_mu("SiO2", density, &[10_000.0], CrossSectionKind::Total)
            .unwrap_err();
        assert!(matches!(err, XrayDbError::InvalidInput(_)), "{err}");
    }
}

#[test]
fn test_material_mu_from_mass_fractions_invalid_inputs() {
    let db = XrayDb::new();

    let err = db
        .material_mu_from_mass_fractions(&[], 2.65, &[10_000.0], CrossSectionKind::Total)
        .unwrap_err();
    assert!(matches!(err, XrayDbError::InvalidInput(_)), "{err}");

    let err = db
        .material_mu_from_mass_fractions(
            &[("Si", -0.1), ("O", 1.1)],
            2.65,
            &[10_000.0],
            CrossSectionKind::Total,
        )
        .unwrap_err();
    assert!(matches!(err, XrayDbError::InvalidInput(_)), "{err}");
}

#[test]
fn test_material_mu_all_cross_section_kinds_are_finite() {
    let db = XrayDb::new();
    let kinds = [
        CrossSectionKind::Total,
        CrossSectionKind::Photo,
        CrossSectionKind::Coherent,
        CrossSectionKind::Incoherent,
    ];

    for kind in kinds {
        let mu = db.material_mu("Fe2O3", 5.24, &[10_000.0], kind).unwrap();
        assert!(mu[0].is_finite());
        assert!(mu[0] > 0.0);
    }
}

#[test]
fn test_xray_delta_beta_si() {
    let db = XrayDb::new();
    let n = db.xray_delta_beta("Si", 2.33, 10000.0).unwrap();
    let (delta, beta, atlen) = (n.delta, n.beta, n.attenuation_length_cm);
    // delta should be small and positive
    assert!(delta > 0.0, "delta = {delta}");
    assert!(delta < 1e-3, "delta = {delta}");
    // beta should be small and positive
    assert!(beta > 0.0, "beta = {beta}");
    assert!(beta < 1e-4, "beta = {beta}");
    // attenuation length should be reasonable (in cm)
    assert!(atlen > 0.0, "atlen = {atlen}");
}

#[test]
fn test_xray_delta_beta_au() {
    let db = XrayDb::new();
    let au = db.xray_delta_beta("Au", 19.3, 10000.0).unwrap();
    let (delta, beta) = (au.delta, au.beta);
    // Gold has high Z, so delta and beta should be larger than Si
    let si = db.xray_delta_beta("Si", 2.33, 10000.0).unwrap();
    let (delta_si, beta_si) = (si.delta, si.beta);
    assert!(delta > delta_si);
    assert!(beta > beta_si);
}

#[test]
fn test_material_mu_invalid_formula() {
    let db = XrayDb::new();
    let err = db
        .material_mu("co", 1.0, &[10_000.0], CrossSectionKind::Total)
        .unwrap_err();
    assert!(matches!(err, XrayDbError::InvalidFormula { .. }));
}

#[test]
fn test_material_mu_unknown_element_symbol_is_error() {
    let db = XrayDb::new();
    let err = db
        .material_mu("SiXx", 2.3, &[10_000.0], CrossSectionKind::Total)
        .unwrap_err();
    assert!(matches!(err, XrayDbError::InvalidFormula { .. }), "{err}");
}

#[test]
fn test_xray_delta_beta_unknown_element_symbol_is_error() {
    let db = XrayDb::new();
    let err = db.xray_delta_beta("FeXx", 7.8, 10_000.0).unwrap_err();
    assert!(matches!(err, XrayDbError::InvalidFormula { .. }), "{err}");
}

#[test]
fn test_material_mu_named_requires_density_for_unknown_material() {
    let db = XrayDb::new();
    let err = db
        .material_mu_named("unobtainium", &[10_000.0], CrossSectionKind::Total, None)
        .unwrap_err();
    assert!(matches!(err, XrayDbError::InvalidInput(_)), "{err}");
}
