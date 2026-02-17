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
fn test_xray_delta_beta_si() {
    let db = XrayDb::new();
    let (delta, beta, atlen) = db.xray_delta_beta("Si", 2.33, 10000.0).unwrap();
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
    let (delta, beta, _atlen) = db.xray_delta_beta("Au", 19.3, 10000.0).unwrap();
    // Gold has high Z, so delta and beta should be larger than Si
    let (delta_si, beta_si, _) = db.xray_delta_beta("Si", 2.33, 10000.0).unwrap();
    assert!(delta > delta_si);
    assert!(beta > beta_si);
}

#[test]
fn test_material_mu_invalid_formula() {
    let db = XrayDb::new();
    let err = db
        .material_mu("co", 1.0, &[10_000.0], CrossSectionKind::Total)
        .unwrap_err();
    assert!(matches!(err, XrayDbError::InvalidFormula(_)));
}

#[test]
fn test_material_mu_unknown_element_symbol_is_error() {
    let db = XrayDb::new();
    let err = db
        .material_mu("SiUnh", 2.3, &[10_000.0], CrossSectionKind::Total)
        .unwrap_err();
    assert!(matches!(err, XrayDbError::UnknownElement(_)));
}

#[test]
fn test_xray_delta_beta_unknown_element_symbol_is_error() {
    let db = XrayDb::new();
    let err = db.xray_delta_beta("FeUnh", 7.8, 10_000.0).unwrap_err();
    assert!(matches!(err, XrayDbError::UnknownElement(_)));
}

#[test]
fn test_material_mu_named_requires_density_for_unknown_material() {
    let db = XrayDb::new();
    let err = db
        .material_mu_named("unobtainium", &[10_000.0], CrossSectionKind::Total, None)
        .unwrap_err();
    assert!(matches!(err, XrayDbError::DataError(_)));
}
