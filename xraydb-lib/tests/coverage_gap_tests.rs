//! Tests aimed at the paths coverage measurement showed were unexercised.
//!
//! `cargo llvm-cov` put the library at 91.3% lines, with `transitions.rs` (77.7%) and
//! `coster_kronig.rs` (69.6%) weakest — mostly error branches and the less-travelled
//! filtering options. These fill those in.

#![cfg(feature = "embedded-data")]

use xraydb::{CrossSectionKind, XrayDb};

// ── transitions.rs ──────────────────────────────────────────────────────────

#[test]
fn xray_lines_filters_by_initial_level() {
    let db = XrayDb::new();

    let all = db.xray_lines("Au", None, None).expect("all lines");
    let k_only = db.xray_lines("Au", Some("K"), None).expect("K lines");
    let l3_only = db.xray_lines("Au", Some("L3"), None).expect("L3 lines");

    assert!(!k_only.is_empty() && !l3_only.is_empty());
    assert!(k_only.len() < all.len(), "filtering must actually narrow");
    assert!(k_only.values().all(|l| l.initial_level == "K"));
    assert!(l3_only.values().all(|l| l.initial_level == "L3"));

    // A level with no transitions yields nothing rather than erroring.
    assert!(
        db.xray_lines("Au", Some("Q9"), None)
            .expect("ok")
            .is_empty()
    );
}

#[test]
fn xray_lines_excitation_filter_drops_unreachable_levels() {
    let db = XrayDb::new();

    // Au K edge is ~80.7 keV; below it, no K lines can be excited.
    let below_k = db
        .xray_lines("Au", None, Some(20_000.0))
        .expect("below K edge");
    assert!(
        below_k.values().all(|l| l.initial_level != "K"),
        "K lines must be dropped below the K edge"
    );
    assert!(!below_k.is_empty(), "L and M lines should remain");

    let above_k = db
        .xray_lines("Au", None, Some(100_000.0))
        .expect("above K edge");
    assert!(
        above_k.values().any(|l| l.initial_level == "K"),
        "K lines must return above the K edge"
    );
    assert!(above_k.len() > below_k.len());
}

#[test]
fn xray_lines_by_level_groups_in_spectroscopic_order() {
    let db = XrayDb::new();
    let groups = db.xray_lines_by_level("Au", None).expect("grouped");

    assert!(groups.len() > 1, "Au should have several initial levels");
    assert_eq!(groups[0].0, "K", "K must sort first");

    for (level, lines) in &groups {
        assert!(!lines.is_empty(), "group {level} is empty");
        assert!(lines.iter().all(|(_, l)| &l.initial_level == level));
        // Sorted by descending intensity within the group.
        for pair in lines.windows(2) {
            assert!(pair[0].1.intensity >= pair[1].1.intensity, "group {level}");
        }
    }

    // The excitation filter propagates through the grouping.
    let filtered = db
        .xray_lines_by_level("Au", Some(20_000.0))
        .expect("filtered");
    assert!(filtered.iter().all(|(level, _)| level != "K"));
}

#[test]
fn guess_edge_honours_a_restricted_edge_list() {
    let db = XrayDb::new();

    // Restricting to L3 must not return a K edge even when one is nearer.
    let guess = db
        .guess_edge(7112.0, Some(&["L3"]), None)
        .expect("some L3 edge exists");
    assert_eq!(guess.edge, "L3");

    // An edge label that matches nothing yields None rather than erroring.
    assert!(db.guess_edge(7112.0, Some(&["Z9"]), None).is_none());

    // The signed difference points the right way on both sides.
    let below = db.guess_edge(7100.0, Some(&["K"]), None).expect("k edge");
    let above = db.guess_edge(7130.0, Some(&["K"]), None).expect("k edge");
    assert!(below.difference > 0.0, "edge above query -> positive");
    assert!(above.difference < 0.0, "edge below query -> negative");
}

#[test]
fn xray_edges_and_lines_reject_unknown_elements() {
    let db = XrayDb::new();
    assert!(db.xray_edges("Xx").is_err());
    assert!(db.xray_lines("Xx", None, None).is_err());
    assert!(db.xray_lines_by_level("Xx", None).is_err());
    assert!(db.xray_lines_by_intensity("Xx", None, None).is_err());
    assert!(db.xray_edge("Fe", "Q9").is_err());
}

// ── coster_kronig.rs ────────────────────────────────────────────────────────

#[test]
fn coster_kronig_returns_direct_and_total_probabilities() {
    let db = XrayDb::new();

    let direct = db.ck_probability("Au", "L1", "L3", false).expect("direct");
    let total = db.ck_probability("Au", "L1", "L3", true).expect("total");

    assert!(direct > 0.0 && direct <= 1.0, "direct = {direct}");
    assert!(total > 0.0 && total <= 1.0, "total = {total}");
    assert!(
        total >= direct,
        "total ({total}) includes routes via intermediate states, so it cannot be \
         smaller than direct ({direct})"
    );
}

#[test]
fn coster_kronig_errors_are_specific() {
    let db = XrayDb::new();

    // Unknown element.
    assert!(db.ck_probability("Xx", "L1", "L2", false).is_err());
    // Known element, transition that does not exist.
    let err = db
        .ck_probability("Fe", "L1", "Q9", false)
        .expect_err("no such transition");
    assert!(err.to_string().contains("L1->Q9"), "{err}");
}

// ── chantler.rs ─────────────────────────────────────────────────────────────

#[test]
fn chantler_energies_can_be_filtered() {
    let db = XrayDb::new();

    let all = db.chantler_energies("Fe", None, None).expect("all");
    let window = db
        .chantler_energies("Fe", Some(5_000.0), Some(10_000.0))
        .expect("window");

    assert!(!window.is_empty() && window.len() < all.len());
    assert!(window.iter().all(|&e| (5_000.0..=10_000.0).contains(&e)));
    assert!(window.windows(2).all(|w| w[1] >= w[0]), "must stay sorted");

    // An empty window is legal, not an error.
    assert!(
        db.chantler_energies("Fe", Some(1e8), Some(1e9))
            .expect("ok")
            .is_empty()
    );
    assert!(db.chantler_energies("Xx", None, None).is_err());
}

#[test]
fn chantler_kinds_differ_and_are_consistent() {
    use xraydb::ChantlerKind;
    let db = XrayDb::new();
    let e = 10_000.0;

    let total = db
        .mu_chantler_at("Fe", e, ChantlerKind::Total)
        .expect("total");
    let photo = db
        .mu_chantler_at("Fe", e, ChantlerKind::Photo)
        .expect("photo");
    let incoh = db
        .mu_chantler_at("Fe", e, ChantlerKind::Incoherent)
        .expect("incoh");

    assert!(total > 0.0 && photo > 0.0 && incoh > 0.0);
    assert!(photo < total, "a component cannot exceed the total");
    assert!(incoh < total);
    // At 10 keV photoabsorption dominates for iron.
    assert!(photo > incoh);

    assert!(db.mu_chantler("Xx", &[e], ChantlerKind::Total).is_err());
}

#[test]
fn chantler_energy_range_bounds_are_respected() {
    let db = XrayDb::new();
    let (emin, emax) = db.chantler_energy_range("Fe").expect("range");
    assert!(emin > 0.0 && emax > emin);

    // Queries outside the range clamp to the endpoints rather than erroring.
    let below = db.f2_chantler_at("Fe", emin / 10.0).expect("below");
    let at_min = db.f2_chantler_at("Fe", emin).expect("at min");
    assert_eq!(below, at_min);

    let above = db.f2_chantler_at("Fe", emax * 10.0).expect("above");
    let at_max = db.f2_chantler_at("Fe", emax).expect("at max");
    assert_eq!(above, at_max);

    assert!(db.chantler_energy_range("Xx").is_err());
}

// ── ionchamber.rs ───────────────────────────────────────────────────────────

#[test]
fn ionchamber_validates_its_arguments() {
    let db = XrayDb::new();
    let ok = |gases: &[(&str, f64)], len: f64| {
        db.ionchamber_fluxes(gases, 1.0, len, 10_000.0, 1e-6, true, true)
    };

    assert!(ok(&[], 30.0).is_err(), "empty gas mixture");
    assert!(
        ok(&[("nitrogen", 0.0)], 30.0).is_err(),
        "fractions sum to 0"
    );
    assert!(ok(&[("nitrogen", 1.0)], 0.0).is_err(), "zero length");
    assert!(ok(&[("nitrogen", 1.0)], -5.0).is_err(), "negative length");
    assert!(ok(&[("nitrogen", 1.0)], f64::NAN).is_err(), "NaN length");
    assert!(ok(&[("unobtainium", 1.0)], 30.0).is_err(), "unknown gas");
}

#[test]
fn ionchamber_handles_mixtures_and_gas_aliases() {
    let db = XrayDb::new();

    let pure = db
        .ionchamber_fluxes(&[("nitrogen", 1.0)], 1.0, 30.0, 10_000.0, 1e-6, true, true)
        .expect("pure N2");
    // Formula aliases resolve to the same material.
    let aliased = db
        .ionchamber_fluxes(&[("N2", 1.0)], 1.0, 30.0, 10_000.0, 1e-6, true, true)
        .expect("N2 alias");
    assert!((pure.incident - aliased.incident).abs() < 1e-6 * pure.incident.abs());

    // Fractions are normalized, so scaling them changes nothing.
    let scaled = db
        .ionchamber_fluxes(&[("nitrogen", 5.0)], 1.0, 30.0, 10_000.0, 1e-6, true, true)
        .expect("scaled");
    assert!((pure.incident - scaled.incident).abs() < 1e-6 * pure.incident.abs());

    // A mixture sits between its pure components.
    let argon = db
        .ionchamber_fluxes(&[("argon", 1.0)], 1.0, 30.0, 10_000.0, 1e-6, true, true)
        .expect("Ar");
    let mixed = db
        .ionchamber_fluxes(
            &[("nitrogen", 0.5), ("argon", 0.5)],
            1.0,
            30.0,
            10_000.0,
            1e-6,
            true,
            true,
        )
        .expect("mixture");
    let (lo, hi) = if argon.incident < pure.incident {
        (argon.incident, pure.incident)
    } else {
        (pure.incident, argon.incident)
    };
    assert!(
        mixed.incident > lo && mixed.incident < hi,
        "mixture {} should sit between {lo} and {hi}",
        mixed.incident
    );
}

#[test]
fn ionchamber_toggles_change_the_answer() {
    let db = XrayDb::new();
    let run = |compton: bool, both: bool| {
        db.ionchamber_fluxes(
            &[("nitrogen", 1.0)],
            1.0,
            30.0,
            10_000.0,
            1e-6,
            compton,
            both,
        )
        .expect("fluxes")
    };

    let base = run(true, true);
    assert!(run(false, true).incident != base.incident, "compton toggle");
    // Counting one carrier instead of two doubles the inferred flux.
    let single = run(true, false);
    assert!(
        (single.incident / base.incident - 2.0).abs() < 1e-9,
        "one carrier should double the flux: {} vs {}",
        single.incident,
        base.incident
    );
}

#[test]
fn material_mu_named_at_matches_the_batch_form() {
    let db = XrayDb::new();
    let e = 12_000.0;
    let scalar = db
        .material_mu_named_at("kapton", e, CrossSectionKind::Total, None)
        .expect("scalar");
    let batch = db
        .material_mu_named("kapton", &[e], CrossSectionKind::Total, None)
        .expect("batch");
    assert!((scalar - batch[0]).abs() <= 1e-12 * batch[0].abs());

    // Density override is honoured.
    let doubled = db
        .material_mu_named_at("kapton", e, CrossSectionKind::Total, Some(2.84))
        .expect("override");
    assert!((doubled / scalar - 2.0).abs() < 1e-9);
}

// ── waasmaier.rs ────────────────────────────────────────────────────────────

#[test]
fn f0_ion_listing_and_errors() {
    let db = XrayDb::new();

    let all = db.f0_ions(None).expect("all ions");
    let iron = db.f0_ions(Some("Fe")).expect("Fe ions");
    assert!(iron.len() < all.len() && !iron.is_empty());
    assert!(iron.iter().all(|i| i.starts_with("Fe")));

    assert!(db.f0_ions(Some("Xx")).is_err());
    assert!(db.f0("NotAnIon", &[0.0]).is_err());
    assert!(db.f0_at("NotAnIon", 0.0).is_err());
}

#[test]
fn f0_decays_with_momentum_transfer() {
    let db = XrayDb::new();
    let q: Vec<f64> = (0..=12).map(|i| f64::from(i) * 0.5).collect();
    let f0 = db.f0("Fe", &q).expect("f0");

    // At q = 0, f0 equals the electron count.
    assert!((f0[0] - 26.0).abs() < 0.1, "f0(0) = {}", f0[0]);
    for pair in f0.windows(2) {
        assert!(pair[1] <= pair[0], "f0 must decrease with q");
    }
    // Ions differ from the neutral atom.
    let fe3 = db.f0_at("Fe3+", 0.0).expect("Fe3+");
    assert!(fe3 < f0[0], "Fe3+ has fewer electrons: {fe3} vs {}", f0[0]);
}
