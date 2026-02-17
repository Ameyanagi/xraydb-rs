use std::thread;

use approx::assert_relative_eq;
use xraydb::{XrayDb, XrayDbError};

#[test]
fn test_element_count() {
    let db = XrayDb::new();
    assert_eq!(db.raw().elements.len(), 118);
}

#[test]
fn test_atomic_number_by_symbol() {
    let db = XrayDb::new();
    assert_eq!(db.atomic_number("H").unwrap(), 1);
    assert_eq!(db.atomic_number("Fe").unwrap(), 26);
    assert_eq!(db.atomic_number("Au").unwrap(), 79);
    assert_eq!(db.atomic_number("U").unwrap(), 92);
}

#[test]
fn test_atomic_number_by_name() {
    let db = XrayDb::new();
    assert_eq!(db.atomic_number("iron").unwrap(), 26);
    assert_eq!(db.atomic_number("gold").unwrap(), 79);
    assert_eq!(db.atomic_number("hydrogen").unwrap(), 1);
}

#[test]
fn test_atomic_number_by_z_string() {
    let db = XrayDb::new();
    assert_eq!(db.atomic_number("26").unwrap(), 26);
    assert_eq!(db.atomic_number("1").unwrap(), 1);
}

#[test]
fn test_atomic_number_by_invalid_z_string() {
    let db = XrayDb::new();
    assert!(matches!(
        db.atomic_number("0"),
        Err(XrayDbError::UnknownElement(_))
    ));
    assert!(matches!(
        db.atomic_number("999"),
        Err(XrayDbError::UnknownElement(_))
    ));
}

#[test]
fn test_unknown_element() {
    let db = XrayDb::new();
    assert!(db.atomic_number("Xx").is_err());
}

#[test]
fn test_symbol() {
    let db = XrayDb::new();
    assert_eq!(db.symbol("26").unwrap(), "Fe");
    assert_eq!(db.symbol("iron").unwrap(), "Fe");
    assert_eq!(db.symbol("Fe").unwrap(), "Fe");
}

#[test]
fn test_molar_mass() {
    let db = XrayDb::new();
    assert_relative_eq!(db.molar_mass("Fe").unwrap(), 55.845, epsilon = 0.01);
    assert_relative_eq!(db.molar_mass("H").unwrap(), 1.0078, epsilon = 0.001);
    assert_relative_eq!(db.molar_mass("Au").unwrap(), 196.967, epsilon = 0.01);
}

#[test]
fn test_density() {
    let db = XrayDb::new();
    assert_relative_eq!(db.density("Fe").unwrap(), 7.86, epsilon = 0.01);
    assert_relative_eq!(db.density("Au").unwrap(), 19.3, epsilon = 0.1);
}

#[test]
fn test_unknown_molar_mass_and_density() {
    let db = XrayDb::new();
    assert!(matches!(
        db.molar_mass("Xx"),
        Err(XrayDbError::UnknownElement(_))
    ));
    assert!(matches!(
        db.density("Xx"),
        Err(XrayDbError::UnknownElement(_))
    ));
}

#[test]
fn test_data_table_counts() {
    let db = XrayDb::new();
    let raw = db.raw();

    assert_eq!(raw.xray_levels.len(), 1430);
    assert_eq!(raw.xray_transitions.len(), 1807);
    assert_eq!(raw.photoabsorption.len(), 98);
    assert_eq!(raw.scattering.len(), 98);
    assert_eq!(raw.chantler.len(), 92);
    assert_eq!(raw.waasmaier.len(), 211);
    assert!(raw.compton_energies.incident.len() > 100);
    assert!(raw.keski_rahkonen_krause.len() > 1000);
}

#[test]
fn test_xray_levels_sample() {
    let db = XrayDb::new();
    let fe_k = db
        .raw()
        .xray_levels
        .iter()
        .find(|l| l.element == "Fe" && l.iupac_symbol == "K")
        .unwrap();

    // Fe K-edge: 7112 eV
    assert_relative_eq!(fe_k.absorption_edge, 7112.0, epsilon = 1.0);
    assert!(fe_k.fluorescence_yield > 0.0);
    assert!(fe_k.jump_ratio > 1.0);
}

#[test]
fn test_chantler_data_sample() {
    let db = XrayDb::new();
    let fe_chantler = db
        .raw()
        .chantler
        .iter()
        .find(|c| c.element == "Fe")
        .unwrap();

    assert!(fe_chantler.energy.len() > 100);
    assert!(fe_chantler.f1.len() == fe_chantler.energy.len());
    assert!(fe_chantler.f2.len() == fe_chantler.energy.len());
    assert!(fe_chantler.density > 7.0);
}

#[test]
fn test_waasmaier_sample() {
    let db = XrayDb::new();
    let fe = db.raw().waasmaier.iter().find(|w| w.ion == "Fe").unwrap();

    assert_eq!(fe.atomic_number, 26);
    assert_eq!(fe.scale.len(), 5);
    assert_eq!(fe.exponents.len(), 5);
}

#[test]
fn test_concurrent_db_new_and_lookup() {
    let mut handles = Vec::new();
    for _ in 0..8 {
        handles.push(thread::spawn(|| {
            let db = XrayDb::new();
            assert_eq!(db.atomic_number("Fe").unwrap(), 26);
            assert_eq!(db.symbol("79").unwrap(), "Au");
            assert!(db.molar_mass("Si").unwrap() > 28.0);
        }));
    }

    for handle in handles {
        handle.join().unwrap();
    }
}
