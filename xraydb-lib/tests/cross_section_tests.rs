use xraydb::{ChantlerKind, CrossSectionKind, XrayDb, XrayDbError};

#[test]
fn test_mu_elam_fe_7112() {
    let db = XrayDb::new();
    let result = db
        .mu_elam("Fe", &[7112.0], CrossSectionKind::Total)
        .unwrap();
    // Fe K-edge at 7112 eV — cross-section is a real value
    assert!(
        result[0] > 1.0,
        "mu_total for Fe at 7112 eV = {}",
        result[0]
    );
}

#[test]
fn test_mu_elam_fe_10kev() {
    let db = XrayDb::new();
    let result = db
        .mu_elam("Fe", &[10000.0], CrossSectionKind::Total)
        .unwrap();
    // Fe at 10 keV (above K-edge): total mu around 170 cm²/g
    assert!(
        result[0] > 100.0,
        "mu_total for Fe at 10 keV = {}",
        result[0]
    );
    assert!(
        result[0] < 300.0,
        "mu_total for Fe at 10 keV = {}",
        result[0]
    );
}

#[test]
fn test_mu_elam_fe_photo() {
    let db = XrayDb::new();
    let photo = db
        .mu_elam("Fe", &[10000.0], CrossSectionKind::Photo)
        .unwrap();
    let total = db
        .mu_elam("Fe", &[10000.0], CrossSectionKind::Total)
        .unwrap();
    // Photo should be the dominant contribution above K-edge
    assert!(photo[0] > 0.0);
    assert!(photo[0] < total[0]);
}

#[test]
fn test_mu_elam_decreasing_above_edge() {
    let db = XrayDb::new();
    // Use energies well above Fe K-edge (7112 eV) where mu should decrease
    let energies = vec![10000.0, 20000.0, 50000.0, 100000.0];
    let result = db
        .mu_elam("Fe", &energies, CrossSectionKind::Total)
        .unwrap();
    assert_eq!(result.len(), 4);
    for i in 1..result.len() {
        assert!(
            result[i] < result[i - 1],
            "expected decreasing mu above edge: {} vs {}",
            result[i],
            result[i - 1]
        );
    }
}

#[test]
fn test_mu_elam_various_elements() {
    let db = XrayDb::new();
    for element in &["H", "C", "O", "Si", "Fe", "Cu", "Ag", "Au", "U"] {
        let result = db
            .mu_elam(element, &[10000.0], CrossSectionKind::Total)
            .unwrap();
        assert!(
            result[0] > 0.0,
            "mu for {element} at 10 keV should be positive"
        );
    }
}

#[test]
fn test_mu_elam_energy_clamping() {
    let db = XrayDb::new();
    // Energies below 100 eV should be clamped, not error
    let result = db.mu_elam("Fe", &[10.0], CrossSectionKind::Total).unwrap();
    assert!(result[0] > 0.0);

    // Energies above 800 keV should be clamped
    let result = db
        .mu_elam("Fe", &[1_000_000.0], CrossSectionKind::Total)
        .unwrap();
    assert!(result[0] > 0.0);
}

#[test]
fn test_f1_chantler_fe() {
    let db = XrayDb::new();
    let result = db.f1_chantler("Fe", &[10000.0]).unwrap();
    // f1 in Chantler table stores f' (anomalous correction, with Z subtracted)
    // At 10 keV (well above K-edge), f' should be small (close to 0)
    assert!(result[0].abs() < 5.0, "f' for Fe at 10keV = {}", result[0]);
}

#[test]
fn test_f2_chantler_fe() {
    let db = XrayDb::new();
    let result = db.f2_chantler("Fe", &[10000.0]).unwrap();
    assert!(result[0] > 0.0, "f2 should be positive");
    assert!(result[0] < 10.0, "f2 for Fe at 10keV = {}", result[0]);
}

#[test]
fn test_mu_chantler_fe() {
    let db = XrayDb::new();
    let total = db
        .mu_chantler("Fe", &[10000.0], ChantlerKind::Total)
        .unwrap();
    let photo = db
        .mu_chantler("Fe", &[10000.0], ChantlerKind::Photo)
        .unwrap();
    assert!(total[0] > 0.0);
    assert!(photo[0] > 0.0);
    assert!(photo[0] <= total[0]);
}

#[test]
fn test_chantler_energies() {
    let db = XrayDb::new();
    let energies = db.chantler_energies("Fe", None, None).unwrap();
    assert!(energies.len() > 100);
    assert!(energies[0] > 0.5);
    assert!(*energies.last().unwrap() < 2e6);
}

#[test]
fn test_chantler_energies_filtered() {
    let db = XrayDb::new();
    let energies = db
        .chantler_energies("Fe", Some(5000.0), Some(10000.0))
        .unwrap();
    assert!(!energies.is_empty());
    assert!(energies[0] >= 5000.0);
    assert!(*energies.last().unwrap() <= 10000.0);
}

#[test]
fn test_cross_section_unknown_element_errors() {
    let db = XrayDb::new();
    assert!(matches!(
        db.mu_elam("Xx", &[10_000.0], CrossSectionKind::Total),
        Err(XrayDbError::UnknownElement(_))
    ));
    assert!(matches!(
        db.f1_chantler("Xx", &[10_000.0]),
        Err(XrayDbError::UnknownElement(_))
    ));
    assert!(matches!(
        db.f2_chantler("Xx", &[10_000.0]),
        Err(XrayDbError::UnknownElement(_))
    ));
    assert!(matches!(
        db.mu_chantler("Xx", &[10_000.0], ChantlerKind::Total),
        Err(XrayDbError::UnknownElement(_))
    ));
}

#[test]
fn test_cross_section_empty_energy_arrays() {
    let db = XrayDb::new();
    let mu = db.mu_elam("Fe", &[], CrossSectionKind::Total).unwrap();
    let f1 = db.f1_chantler("Fe", &[]).unwrap();
    let f2 = db.f2_chantler("Fe", &[]).unwrap();
    let mu_ch = db.mu_chantler("Fe", &[], ChantlerKind::Total).unwrap();

    assert!(mu.is_empty());
    assert!(f1.is_empty());
    assert!(f2.is_empty());
    assert!(mu_ch.is_empty());
}

#[test]
fn test_cross_section_single_energy_point() {
    let db = XrayDb::new();
    let energy = [7112.0];
    let mu = db.mu_elam("Fe", &energy, CrossSectionKind::Total).unwrap();
    let f1 = db.f1_chantler("Fe", &energy).unwrap();
    let f2 = db.f2_chantler("Fe", &energy).unwrap();
    let mu_ch = db.mu_chantler("Fe", &energy, ChantlerKind::Total).unwrap();

    assert_eq!(mu.len(), 1);
    assert_eq!(f1.len(), 1);
    assert_eq!(f2.len(), 1);
    assert_eq!(mu_ch.len(), 1);
    assert!(mu[0].is_finite() && mu[0] > 0.0);
}
