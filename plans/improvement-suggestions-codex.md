# Broader Perf + Quality Sweep for `xraydb-lib`

## Summary
Improve runtime performance by replacing repeated linear scans with precomputed indexes, and improve code quality by making material/formula calculations fail fast on invalid elements. Add regression tests plus Criterion benchmarks to make gains measurable and prevent regressions.

## Scope
In scope:
- Internal lookup/index refactor for element-driven data access.
- Strict error propagation in material composition math.
- New tests for behavior and edge cases.
- Criterion benchmarks for hot lookup paths.

Out of scope:
- Changing serialized database format.
- Breaking public method signatures.
- Feature-gating strict behavior (strict is default for this pass).

## Public API / Interface Impact
No signature changes to existing public methods.
Behavior change:
- `material_mu` and `xray_delta_beta` paths will now return `Err` for unknown element symbols in formulas instead of silently treating missing molar mass as `0.0`.
- Error shape should be the existing propagated element-resolution error (`XrayDbError::UnknownElement(...)`) unless parsing fails earlier.

## Implementation Plan

1. Add internal indexes in `xraydb-lib/src/db.rs`
- Extend `InitializedDb` with:
  - `element_by_z: HashMap<u16, usize>` (or fixed-size `Vec<Option<usize>>` indexed by Z).
  - `photoabs_by_symbol: HashMap<String, usize>`.
  - `scattering_by_symbol: HashMap<String, usize>`.
  - `chantler_by_symbol: HashMap<String, usize>`.
- Build these once in `db()` during initialization.
- Keep current symbol/name maps, but normalize lookup paths to avoid repeated ad-hoc scanning.

2. Refactor core element resolution/access
- Update `resolve_element` to validate numeric Z via `element_by_z` instead of scanning `elements`.
- Update `element_record` to fetch via indexed position, not `.iter().find(...)`.
- Add small internal helpers for indexed table retrieval (photoabsorption/scattering/chantler) to centralize lookup logic and error mapping.

3. Refactor hot-path callers to use indexed helpers
- `xraydb-lib/src/elam.rs`: replace per-call `.iter().find(...)` with indexed retrieval.
- `xraydb-lib/src/chantler.rs`: replace `.iter().find(...)` with indexed retrieval.
- Keep return values and interpolation behavior unchanged.

4. Enforce strict error propagation in material math
- `xraydb-lib/src/materials.rs`:
  - Replace `self.molar_mass(sym).unwrap_or(0.0)` with `self.molar_mass(sym)?` in:
    - total formula weight computation,
    - mass fraction computation,
    - delta/beta weight computations.
- Preserve existing `zero weight formula` error when appropriate, but unknown symbols should now error early and explicitly.

5. Quality cleanup (while touching code)
- Consolidate duplicated lookup/error branches where possible.
- Keep all refactors internal and documented with short comments only where non-obvious.

6. Add Criterion benchmarks
- Update `xraydb-lib/Cargo.toml` dev-dependencies with `criterion`.
- Add `xraydb-lib/benches/lookup_bench.rs` covering:
  - repeated `atomic_number` / `symbol` / `molar_mass` lookups,
  - `mu_elam` for one element across a moderate energy vector,
  - `f2_chantler` for one element across a moderate energy vector.
- Benchmark both mixed identifier styles (`"Fe"`, `"iron"`, `"26"`) to exercise resolver paths.

## Tests and Scenarios

1. Regression tests for lookup behavior
- File: `xraydb-lib/tests/element_tests.rs`
- Add/extend cases:
  - numeric element resolution valid/invalid,
  - case-insensitive name/symbol resolution still works,
  - unknown element still returns `Err`.

2. Strict material error tests
- File: `xraydb-lib/tests/material_tests.rs`
- Add cases:
  - `material_mu("Xx2O", ...)` returns `Err(UnknownElement)` (or equivalent propagated error).
  - `xray_delta_beta("SiXx", ...)` returns `Err`.
  - Existing valid formulas continue to pass.

3. No-regression coverage for optics/material lookup
- File: `xraydb-lib/tests/optics_tests.rs`
- Ensure `find_material` behavior remains unchanged for known names/formulas and case-insensitive names.

4. Full verification commands
- `cargo test -p xraydb --tests`
- `cargo bench -p xraydb` (or `cargo bench` at workspace level if configured)

## Acceptance Criteria
- All existing tests pass.
- New strict-error tests pass.
- Benchmarks compile and run.
- Indexed lookup code removes repeated `.iter().find(...)` from `db.rs`, `elam.rs`, and `chantler.rs` hot paths.
- No public API signature changes.

## Assumptions and Defaults
- Chosen scope: broader perf sweep.
- Chosen validation mode: strict errors by default (no feature flag).
- Chosen benchmark approach: Criterion-based benches in `xraydb-lib`.
- Internal memory overhead for indexes is acceptable in exchange for faster steady-state queries.
