//! End-to-end tests that invoke the built `xraydb` binary.
//!
//! These assert known physical values, so they fail if either the CLI plumbing or
//! the underlying data regresses.

use std::process::{Command, Output};

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_xraydb")
}

fn run(args: &[&str]) -> Output {
    Command::new(bin())
        .args(args)
        .env("NO_COLOR", "1")
        .output()
        .expect("failed to execute the xraydb binary")
}

fn stdout_of(args: &[&str]) -> String {
    let out = run(args);
    assert!(
        out.status.success(),
        "`xraydb {}` failed: {}",
        args.join(" "),
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).expect("stdout is valid UTF-8")
}

fn json_of(args: &[&str]) -> serde_json::Value {
    let text = stdout_of(args);
    serde_json::from_str(&text)
        .unwrap_or_else(|e| panic!("invalid JSON from {args:?}: {e}\n{text}"))
}

#[test]
fn element_reports_known_values() {
    let v = json_of(&["--json", "element", "Fe"]);
    assert_eq!(v["z"], 26);
    assert_eq!(v["symbol"], "Fe");
    assert_eq!(v["name"], "iron");
    assert!((v["molar_mass"].as_f64().unwrap_or_default() - 55.845).abs() < 0.001);
}

#[test]
fn element_accepts_symbol_name_and_z() {
    for id in ["Fe", "fe", "iron", "IRON", "26"] {
        let v = json_of(&["--json", "element", id]);
        assert_eq!(v["z"], 26, "identifier {id}");
    }
}

#[test]
fn fe_k_edge_is_7112_ev() {
    let v = json_of(&["--json", "edge", "Fe", "K"]);
    assert_eq!(v["energy_ev"], 7112.0);

    let text = stdout_of(&["edges", "Fe"]);
    assert!(text.contains("7112"), "{text}");
}

#[test]
fn cu_kalpha1_matches_the_elam_tables() {
    let v = json_of(&["--json", "lines", "Cu", "--level", "K"]);
    let ka1 = v
        .as_array()
        .and_then(|rows| rows.iter().find(|r| r["line"] == "Ka1"))
        .expect("Cu has a Ka1 line");
    // Elam tabulates 8046.3 eV.
    assert!((ka1["energy_eV"].as_f64().unwrap_or_default() - 8046.3).abs() < 0.5);
}

#[test]
fn lines_are_sorted_by_descending_intensity() {
    let v = json_of(&["--json", "lines", "Au"]);
    let rows = v.as_array().expect("array");
    let intensities: Vec<f64> = rows
        .iter()
        .filter_map(|r| r["intensity"].as_f64())
        .collect();
    for pair in intensities.windows(2) {
        assert!(pair[0] >= pair[1], "not sorted: {intensities:?}");
    }
}

#[test]
fn water_attenuation_matches_the_known_value() {
    let v = json_of(&[
        "--json",
        "mu",
        "H2O",
        "--density",
        "1.0",
        "--energy",
        "10000",
    ]);
    let mu = v[0]["mu_per_cm"].as_f64().unwrap_or_default();
    assert!((mu - 5.33).abs() < 0.05, "got {mu}");
}

#[test]
fn energy_spec_forms_all_work() {
    // Single value
    assert_eq!(
        json_of(&["--json", "mu", "Fe", "-e", "10000"])
            .as_array()
            .map(Vec::len),
        Some(1)
    );
    // Comma list
    assert_eq!(
        json_of(&["--json", "mu", "Fe", "-e", "5000,7112,10000"])
            .as_array()
            .map(Vec::len),
        Some(3)
    );
    // start:stop:step
    assert_eq!(
        json_of(&["--json", "mu", "Fe", "-e", "1000:1400:100"])
            .as_array()
            .map(Vec::len),
        Some(5)
    );
    // start:stop/count, log-spaced
    let v = json_of(&["--json", "mu", "Fe", "-e", "100:10000/3"]);
    let rows = v.as_array().expect("array");
    assert_eq!(rows.len(), 3);
    assert!((rows[1]["energy_eV"].as_f64().unwrap_or_default() - 1000.0).abs() < 1e-6);
}

#[test]
fn material_lookup_supplies_its_own_density() {
    let v = json_of(&["--json", "material", "kapton"]);
    assert_eq!(v["formula"], "C22H10N2O5");
    assert_eq!(v["density"], 1.42);

    // `mu` on a material name needs no --density.
    let out = run(&["mu", "kapton", "-e", "10000"]);
    assert!(out.status.success());
}

#[test]
fn guess_edge_finds_fe_k() {
    let v = json_of(&["--json", "guess-edge", "--energy", "7100"]);
    assert_eq!(v["element"], "Fe");
    assert_eq!(v["edge"], "K");
    assert_eq!(v["difference_ev"], 12.0);
}

#[test]
fn guess_edge_respects_tolerance() {
    let out = run(&["guess-edge", "--energy", "200000", "--tolerance", "500"]);
    assert!(
        !out.status.success(),
        "should fail with no match in tolerance"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("within 500"), "{stderr}");
}

#[test]
fn delta_beta_is_physically_sensible() {
    let v = json_of(&[
        "--json",
        "delta-beta",
        "SiO2",
        "--density",
        "2.2",
        "--energy",
        "10000",
    ]);
    let delta = v[0]["delta"].as_f64().unwrap_or_default();
    let beta = v[0]["beta"].as_f64().unwrap_or_default();
    assert!(delta > 0.0 && delta < 1e-4, "delta = {delta}");
    assert!(beta > 0.0 && beta < delta, "beta = {beta}");
    // Critical angle for SiO2 at 10 keV is roughly 3 mrad.
    let theta_c = v[0]["crit_angle_mrad"].as_f64().unwrap_or_default();
    assert!((theta_c - 3.0).abs() < 0.5, "theta_c = {theta_c}");
}

#[test]
fn air_works_end_to_end() {
    // The 0.1.2 bug: the `air` entry's formula was unparseable.
    let out = run(&["mu", "air", "-e", "10000"]);
    assert!(
        out.status.success(),
        "air must work: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let out = run(&[
        "ionchamber",
        "--gas",
        "air",
        "--volts",
        "1.0",
        "--length",
        "30",
        "--energy",
        "10000",
        "--sensitivity",
        "1e-6",
    ]);
    assert!(
        out.status.success(),
        "air ion chamber must work: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}

#[test]
fn every_material_in_the_table_is_usable() {
    let v = json_of(&["--json", "materials"]);
    let rows = v.as_array().expect("array");
    assert!(rows.len() > 80);

    for row in rows {
        let name = row["name"].as_str().unwrap_or_default();
        let out = run(&["mu", name, "-e", "10000"]);
        assert!(
            out.status.success(),
            "material '{name}' failed: {}",
            String::from_utf8_lossy(&out.stderr)
        );
    }
}

#[test]
fn csv_output_is_well_formed() {
    let text = stdout_of(&["--csv", "mu", "Fe", "-e", "8000,10000"]);
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len(), 3);
    assert_eq!(lines[0], "energy_eV,mu_cm2_per_g,mu_per_cm,atten_len_cm");
    assert_eq!(lines[1].split(',').count(), 4);
}

#[test]
fn text_output_has_a_header_and_a_rule() {
    let text = stdout_of(&["edges", "Fe"]);
    let lines: Vec<&str> = text.lines().collect();
    assert!(lines[0].starts_with("edge"), "{:?}", lines[0]);
    assert!(
        lines[1].chars().all(|c| c == '-' || c == ' '),
        "{:?}",
        lines[1]
    );
}

#[test]
fn unknown_element_fails_with_a_readable_message() {
    let out = run(&["element", "Fx"]);
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("unknown element: Fx"), "{stderr}");
    // Not a debug-formatted enum.
    assert!(!stderr.contains("UnknownElement"), "{stderr}");
}

#[test]
fn unknown_edge_lists_the_available_ones() {
    let out = run(&["edge", "Fe", "Q9"]);
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("available:"), "{stderr}");
    assert!(stderr.contains('K'), "{stderr}");
}

#[test]
fn bare_formula_without_density_explains_itself() {
    let out = run(&["mu", "Fe2O3", "-e", "10000"]);
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("--density is required"), "{stderr}");
}

#[test]
fn malformed_energy_spec_is_rejected() {
    for bad in ["bogus", "1000:2000", "1000:2000/0"] {
        let out = run(&["mu", "Fe", "-e", bad]);
        assert!(!out.status.success(), "'{bad}' should be rejected");
    }
}

#[test]
fn info_reports_the_database_contents() {
    let v = json_of(&["--json", "info"]);
    assert_eq!(v["elements"], 118);
    assert_eq!(v["chantler_elements"], 92);
    assert!(v["materials"].as_u64().unwrap_or_default() > 80);
}

#[test]
fn help_and_version_work() {
    assert!(run(&["--help"]).status.success());
    assert!(run(&["--version"]).status.success());
    for sub in [
        "element",
        "edges",
        "edge",
        "lines",
        "mu",
        "f1f2",
        "f0",
        "delta-beta",
        "material",
        "materials",
        "guess-edge",
        "core-width",
        "ionchamber",
        "compton",
        "info",
    ] {
        let out = run(&[sub, "--help"]);
        assert!(out.status.success(), "`{sub} --help` failed");
    }
}
