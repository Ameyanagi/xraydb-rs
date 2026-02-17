use approx::assert_relative_eq;
use xraydb::XrayDb;

#[test]
fn test_f0_fe_q0() {
    let db = XrayDb::new();
    let result = db.f0("Fe", &[0.0]).unwrap();
    // At q=0, f0 should equal Z (26 for Fe)
    assert_relative_eq!(result[0], 26.0, epsilon = 0.1);
}

#[test]
fn test_f0_decreasing() {
    let db = XrayDb::new();
    let result = db.f0("Fe", &[0.0, 0.5, 1.0, 2.0]).unwrap();
    // f0 decreases with increasing q
    for i in 1..result.len() {
        assert!(result[i] < result[i - 1]);
    }
}

#[test]
fn test_f0_ions_list() {
    let db = XrayDb::new();
    let ions = db.f0_ions(None).unwrap();
    assert_eq!(ions.len(), 211);
}

#[test]
fn test_f0_ions_for_element() {
    let db = XrayDb::new();
    let fe_ions = db.f0_ions(Some("Fe")).unwrap();
    assert!(fe_ions.contains(&"Fe"));
    assert!(fe_ions.len() >= 1);
}

#[test]
fn test_f0_unknown_ion() {
    let db = XrayDb::new();
    assert!(db.f0("Xx99+", &[0.0]).is_err());
}

#[test]
fn test_xray_edges_fe() {
    let db = XrayDb::new();
    let edges = db.xray_edges("Fe").unwrap();
    assert!(edges.contains_key("K"));
    assert!(edges.contains_key("L1"));
    assert!(edges.contains_key("L2"));
    assert!(edges.contains_key("L3"));
    assert_relative_eq!(edges["K"].energy, 7112.0, epsilon = 1.0);
}

#[test]
fn test_xray_edge_fe_k() {
    let db = XrayDb::new();
    let edge = db.xray_edge("Fe", "K").unwrap();
    assert_relative_eq!(edge.energy, 7112.0, epsilon = 1.0);
    assert!(edge.fluorescence_yield > 0.0);
    assert!(edge.jump_ratio > 1.0);
}

#[test]
fn test_xray_lines_fe() {
    let db = XrayDb::new();
    let lines = db.xray_lines("Fe", None, None).unwrap();
    assert!(lines.contains_key("Ka1"));
    assert!(lines.contains_key("Ka2"));
    assert!(lines.contains_key("Kb1"));
    // Ka1 energy should be around 6404 eV
    assert_relative_eq!(lines["Ka1"].energy, 6404.0, epsilon = 2.0);
}

#[test]
fn test_xray_lines_initial_level() {
    let db = XrayDb::new();
    let k_lines = db.xray_lines("Fe", Some("K"), None).unwrap();
    for (_, line) in &k_lines {
        assert_eq!(line.initial_level, "K");
    }
}

#[test]
fn test_xray_lines_excitation_energy() {
    let db = XrayDb::new();
    // Only lines excitable below 2000 eV (no K lines for Fe)
    let lines = db.xray_lines("Fe", None, Some(2000.0)).unwrap();
    for (name, _) in &lines {
        assert!(
            !name.starts_with('K'),
            "K lines should not be excitable below 2keV: {name}"
        );
    }
}

#[test]
fn test_guess_edge_fe_k() {
    let db = XrayDb::new();
    let result = db.guess_edge(7112.0, None);
    assert!(result.is_some());
    let (elem, edge) = result.unwrap();
    assert_eq!(elem, "Fe");
    assert_eq!(edge, "K");
}

#[test]
fn test_ionization_potential_argon() {
    let db = XrayDb::new();
    let pot = db.ionization_potential("argon").unwrap();
    assert_relative_eq!(pot, 26.4, epsilon = 0.1);
}

#[test]
fn test_ck_probability() {
    let db = XrayDb::new();
    // Cu has L1 -> L2 Coster-Kronig transition
    let prob = db.ck_probability("Cu", "L1", "L2", false);
    // This may or may not exist for Cu; just test it doesn't panic
    if let Ok(p) = prob {
        assert!(p >= 0.0 && p <= 1.0);
    }
}

#[test]
fn test_core_width_fe_k() {
    let db = XrayDb::new();
    let widths = db.core_width("Fe", Some("K")).unwrap();
    assert!(widths.contains_key("K"));
    let w = widths["K"];
    assert!(w > 0.0);
    assert!(w < 10.0); // K width for Fe should be a few eV
}

#[test]
fn test_core_width_all() {
    let db = XrayDb::new();
    let widths = db.core_width("Fe", None).unwrap();
    assert!(widths.len() > 1);
    assert!(widths.contains_key("K"));
}

#[test]
fn test_compton_energies() {
    let db = XrayDb::new();
    let ce = db.compton_energies(10000.0);
    assert_relative_eq!(ce.incident, 10000.0);
    // At 10 keV, Compton shift is small
    assert!(ce.xray_90deg > 0.0);
    assert!(ce.xray_90deg < 10000.0);
    assert!(ce.xray_mean > 0.0);
    assert!(ce.electron_mean > 0.0);
    assert!(ce.electron_mean < 10000.0);
}
