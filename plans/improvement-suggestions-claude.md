# Performance & Code Quality Improvements for xraydb-rs

## Context

The xraydb-rs codebase is well-structured and fully functional (72 tests passing), but has accumulated code quality issues (29 clippy warnings), performance inefficiencies (O(n) linear scans where O(1) lookups are possible, redundant allocations in hot paths), and gaps in test coverage (no error path tests, no boundary tests). This plan addresses all three areas while keeping the public API unchanged.

---

## Group 1: Fix All Clippy Warnings (29 warnings → 0)

### 1.1 — `xraydb-lib/src/chemparser.rs`: 6 redundant closures
Lines 137, 201, 208, 212, 227, 231: Replace `.map_err(|e| XrayDbError::InvalidFormula(e))` → `.map_err(XrayDbError::InvalidFormula)`

### 1.2 — `xraydb-lib/src/db.rs`: 1 collapsible if
Line 70-80: Collapse nested `if let Ok(z)` + `if .any()` into `if let Ok(z) = ... && ...` (edition 2024 let-chains)

### 1.3 — `xraydb-lib/src/transitions.rs`: 2 collapsible ifs
- Line 82-85: `if let Some(level) = initial_level { if trans.initial_level != level {` → collapse with `&&`
- Line 89-97: `if let Some(level) = self...find(...) { if level.absorption_edge > max_energy {` → collapse with `&&`

### 1.4 — `xraydb-lib/src/optics.rs`: 4 manual assign operations
- Line 289: `ktz = ktz / n` → `ktz /= n`
- Line 295: `r = r * (...)` → `r *= (...)`
- Line 408-409: `r_amp = r_amp * (...)` → `r_amp *= (...)`
- Line 438-439: `r_amp = r_amp * (...)` → `r_amp *= (...)`

### 1.5 — `xraydb-lib/src/optics.rs`, `ionchamber.rs`: 4 too-many-arguments warnings
- `ionchamber_fluxes` (line 59): 8 args
- `darwin_width` (line 91): 9 args
- `multilayer_reflectivity` (line 318): 12 args
- `coated_reflectivity` (line 462): 12 args

Add `#[allow(clippy::too_many_arguments)]` above each — these are physics APIs where all parameters are necessary.

### 1.6 — `xraydb-lib/tests/phase3_tests.rs`: 4 test warnings
- Line 34: `.len() >= 1` → `!is_empty()`
- Line 78: `for (_, line) in &k_lines` → `for line in k_lines.values()`
- Line 88: `for (name, _) in &lines` → `for name in lines.keys()`
- Line 120: `p >= 0.0 && p <= 1.0` → `(0.0..=1.0).contains(&p)`

### 1.7 — `xraydb-generate/src/`: 3 generator warnings
- `chantler.rs:139`: `.map_or(false, ...)` → `.is_some_and(...)`
- `waasmaier.rs:25-26`: `starts_with("#S ") + line[3..]` → `strip_prefix("#S ")`
- `elam.rs:14-20`: Define `type ElamParseResult = (Vec<XrayLevelRecord>, ...)` type alias

---

## Group 2: Performance — O(1) Lookup Indices

Currently, every record lookup (chantler, photoabsorption, scattering, waasmaier) does `iter().find()` — O(n) string comparison. Build HashMap indices at init time in `InitializedDb`.

### 2.1 — Add index maps to `InitializedDb` in `xraydb-lib/src/db.rs`

Add fields:
```rust
z_to_element_idx: HashMap<u16, usize>,         // Z → elements[idx]
symbol_to_chantler_idx: HashMap<String, usize>, // symbol → chantler[idx]
symbol_to_photo_idx: HashMap<String, usize>,    // symbol → photoabsorption[idx]
symbol_to_scatter_idx: HashMap<String, usize>,  // symbol → scattering[idx]
ion_to_waasmaier_idx: HashMap<String, usize>,   // ion → waasmaier[idx]
```

Build during `db()` init with single pass over each table.

### 2.2 — Add private accessor methods to `XrayDb` in `db.rs`

```rust
fn element_by_z(&self, z: u16) -> Option<&ElementRecord>
fn chantler_by_symbol(&self, sym: &str) -> Option<&ChantlerRecord>
fn photo_by_symbol(&self, sym: &str) -> Option<&PhotoabsorptionRecord>
fn scatter_by_symbol(&self, sym: &str) -> Option<&ScatteringRecord>
fn waasmaier_by_ion(&self, ion: &str) -> Option<&WaasmaierRecord>
```

### 2.3 — Update callers to use indices

- `db.rs:92-99` (`element_record`): Use `element_by_z()` instead of `iter().find()`
- `db.rs:70-80` (`resolve_element`): Use `z_to_element_idx.contains_key()` instead of `iter().any()`
- `chantler.rs:17-24` (`chantler_record`): Use `chantler_by_symbol()`
- `elam.rs:54-59, 70-75, 86-91` (`cross_section_elam`): Use `photo_by_symbol()` / `scatter_by_symbol()`
- `waasmaier.rs:31-36` (`f0`): Use `waasmaier_by_ion()`
- `waasmaier.rs:9-17` (`f0_ions`): Use index filter for element-specific ions

---

## Group 3: Performance — Reduce Allocations

### 3.1 — Deduplicate energy clamping in `chantler.rs`

Extract `clamp_energies()` helper used by `f1_chantler`, `f2_chantler`, `mu_chantler` (lines 53-56, 69-72, 90-93).

### 3.2 — Deduplicate safe-value clamping in `chantler.rs`

Extract `safe_for_log()` helper used at lines 75 and 102.

### 3.3 — Reduce `interp_loglog` allocations in `interp.rs`

Currently allocates 3 intermediate Vecs. Reduce to 2 by inlining the `log_x` computation:
```rust
pub fn interp_loglog(x: &[f64], xp: &[f64], fp: &[f64]) -> Vec<f64> {
    let log_xp: Vec<f64> = xp.iter().map(|v| v.ln()).collect();
    let log_fp: Vec<f64> = fp.iter().map(|v| v.ln()).collect();
    x.iter()
        .map(|&xi| interp_one(xi.ln(), &log_xp, &log_fp).exp())
        .collect()
}
```

### 3.4 — Optimize `mu_elam` Total kind in `elam.rs`

Currently `Total` calls `cross_section_elam` 3 times, each doing its own symbol lookup, clamp, and ln. Compute all three in a single pass: 1 symbol lookup, 1 clamp, 1 ln, 3 spline evaluations.

### 3.5 — Eliminate double `molar_mass` lookups in `materials.rs`

In `material_mu` (lines 27-45), `molar_mass()` is called twice per element. Precompute masses once into a local map, then reuse for both total weight and mass fractions.

---

## Group 4: Additional Tests (~25 new tests)

### 4.1 — Error path tests in `element_tests.rs`
- Unknown element by symbol, by out-of-range Z, by invalid name
- Unknown molar_mass, unknown density

### 4.2 — Error path tests in `cross_section_tests.rs`
- `mu_elam` / `f1_chantler` / `f2_chantler` / `mu_chantler` with unknown element

### 4.3 — Error path tests in `phase3_tests.rs`
- Unknown element/edge for `xray_edge`, `xray_lines`, `core_width`
- Unknown gas for `ionization_potential`
- Unknown element for `ck_probability`

### 4.4 — Error path tests in `material_tests.rs`
- Invalid formula, empty formula for `material_mu`
- Invalid formula for `xray_delta_beta`
- Unknown material without density for `material_mu_named`

### 4.5 — Boundary condition tests in `cross_section_tests.rs`
- Single energy point
- Empty energy array → empty result
- Energy at exact edge boundary

### 4.6 — Concurrency test in `element_tests.rs`
- Multi-threaded `XrayDb::new()` + queries to verify `OnceLock` safety

---

## Files to Modify

| File | Groups |
|------|--------|
| `xraydb-lib/src/chemparser.rs` | 1 |
| `xraydb-lib/src/db.rs` | 1, 2 |
| `xraydb-lib/src/transitions.rs` | 1 |
| `xraydb-lib/src/optics.rs` | 1 |
| `xraydb-lib/src/ionchamber.rs` | 1 |
| `xraydb-lib/src/chantler.rs` | 2, 3 |
| `xraydb-lib/src/elam.rs` | 2, 3 |
| `xraydb-lib/src/interp.rs` | 3 |
| `xraydb-lib/src/waasmaier.rs` | 2 |
| `xraydb-lib/src/materials.rs` | 3 |
| `xraydb-lib/tests/phase3_tests.rs` | 1, 4 |
| `xraydb-lib/tests/element_tests.rs` | 4 |
| `xraydb-lib/tests/cross_section_tests.rs` | 4 |
| `xraydb-lib/tests/material_tests.rs` | 4 |
| `xraydb-generate/src/chantler.rs` | 1 |
| `xraydb-generate/src/waasmaier.rs` | 1 |
| `xraydb-generate/src/elam.rs` | 1 |

## Verification

After each group:
```bash
cargo test --all-features          # all existing + new tests pass
cargo clippy --all-targets --all-features  # zero warnings (after Group 1)
```

Expected final state: ~97 tests passing, 0 clippy warnings.
