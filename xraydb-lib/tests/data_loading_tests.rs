//! Lifecycle of the global database handle.
//!
//! The database is a process-wide `OnceLock`, so these all live in one test function:
//! test threads share the process, and the assertions here depend on initialisation
//! order. Splitting them into separate `#[test]`s would make them race.

#![cfg(feature = "zstd")]

use xraydb::XrayDb;

/// Bytes identical to what the `embedded-data` feature compiles in.
const BLOB: &[u8] = include_bytes!("../data/xraydb.bin.zst");

#[test]
fn database_loading_lifecycle() {
    // ── Before anything is loaded ──────────────────────────────────────────
    assert!(
        !XrayDb::is_loaded(),
        "nothing should be loaded before the first constructor call"
    );

    // ── Loading from supplied bytes ────────────────────────────────────────
    let db = XrayDb::load_compressed(BLOB).expect("the shipped blob must load");
    assert!(XrayDb::is_loaded(), "is_loaded must flip after a load");
    assert_eq!(db.atomic_number("Fe").expect("Fe"), 26);
    assert_eq!(db.xray_edge("Fe", "K").expect("Fe K").energy, 7112.0);

    // ── current() returns the same database ────────────────────────────────
    let same = XrayDb::current().expect("already loaded");
    assert_eq!(same.element_count(), db.element_count());
    assert_eq!(same.data_version(), db.data_version());

    // ── First call wins ────────────────────────────────────────────────────
    // Later loads return the existing database rather than replacing it, so handles
    // obtained earlier stay valid. The bytes are not even examined, which is why
    // obvious rubbish still succeeds here.
    let again = XrayDb::load_compressed(b"not a zstd stream at all")
        .expect("returns the already-loaded database");
    assert_eq!(again.atomic_number("Fe").expect("Fe"), 26);

    let again = XrayDb::load_uncompressed(b"not postcard either")
        .expect("returns the already-loaded database");
    assert_eq!(again.atomic_number("Cu").expect("Cu"), 29);

    // ── The embedded constructor agrees with the loaded data ───────────────
    #[cfg(feature = "embedded-data")]
    {
        let embedded = XrayDb::try_new().expect("embedded");
        assert_eq!(embedded.element_count(), db.element_count());
        assert_eq!(
            embedded
                .mu_elam_at("Fe", 10_000.0, xraydb::CrossSectionKind::Total)
                .ok(),
            db.mu_elam_at("Fe", 10_000.0, xraydb::CrossSectionKind::Total)
                .ok(),
        );
    }

    // ── Free functions work off whatever is loaded ─────────────────────────
    let water = xraydb::chemparse("H2O").expect("chemparse uses the loaded database");
    assert_eq!(water.get("H"), Some(2.0));
}
