<img src="https://raw.githubusercontent.com/Ameyanagi/xraydb-rs/main/assets/icon-512.png" alt="" width="88" align="left" hspace="12">

# xraydb-rs

X-ray reference data for the elements in Rust. A pure-Rust port of the [XrayDB](https://github.com/xraypy/XrayDB) project.

The complete database — elements, absorption edges, emission lines, Elam and Chantler cross-sections, Waasmaier–Kirfel form factors — is compiled into your binary as a 3 MB zstd blob and decoded lazily on first use. No runtime data files, no network access.

## Install

```toml
[dependencies]
xraydb = "0.4"
```

## Usage

```rust
use xraydb::{CrossSectionKind, XrayDb};

let db = XrayDb::try_new()?;

// Element facts — symbol, name, or Z all work as identifiers
assert_eq!(db.atomic_number("Fe")?, 26);
assert_eq!(db.symbol("iron")?, "Fe");
assert_eq!(db.atomic_name("26")?, "iron");

// Absorption edges
assert_eq!(db.xray_edge("Fe", "K")?.energy, 7112.0);

// Mass attenuation from the Elam tables (cm²/g)
let mu = db.mu_elam_at("Fe", 10_000.0, CrossSectionKind::Total)?;

// Compound attenuation (1/cm), by formula and density
let mu_water = db.material_mu_at("H2O", 1.0, 10_000.0, CrossSectionKind::Total)?;
assert!((mu_water - 5.33).abs() < 0.05);

// Anomalous scattering factors from the Chantler tables
let f1 = db.f1_chantler_at("Fe", 10_000.0)?;
let f2 = db.f2_chantler_at("Fe", 10_000.0)?;

// Refractive index decrements: n = 1 - delta - i*beta
let n = db.xray_delta_beta("SiO2", 2.2, 10_000.0)?;
println!("delta={:.3e} beta={:.3e} atlen={:.3} cm", n.delta, n.beta, n.attenuation_length_cm);
# Ok::<(), xraydb::XrayDbError>(())
```

### Scalar and batch forms

Every array-valued calculation has an `_at` scalar counterpart that allocates nothing. Use the batch form when you have many energies — it resolves the element and its table rows once for the whole slice.

```rust
# use xraydb::{CrossSectionKind, XrayDb};
# let db = XrayDb::try_new()?;
let energies: Vec<f64> = (0..500).map(|i| 5_000.0 + 20.0 * i as f64).collect();
let mu = db.mu_elam("Fe", &energies, CrossSectionKind::Total)?;
assert_eq!(mu.len(), energies.len());
# Ok::<(), xraydb::XrayDbError>(())
```

### Chemical formulas

One parser serves the whole crate, so anything `material_mu` accepts, `xray_delta_beta` accepts too.

| Notation | Example |
|---|---|
| plain | `H2O` |
| nested groups | `Mn(SO4)2(H2O)7` |
| fractional stoichiometry | `Fe0.7Mg0.3O`, `Fe.7Mg.3O` |
| scientific notation | `Zn1.e-5Fe3O4` |
| fractional group multipliers | `(N2)0.7808(O2)0.2095` |
| weight-percent mixtures | `Ru1wt%SiO2` |
| deuterium alias | `D2O` (counted as hydrogen) |

### Named materials

```rust
# use xraydb::{CrossSectionKind, XrayDb};
# let db = XrayDb::try_new()?;
let kapton = db.find_material("kapton").expect("in the materials table");
assert_eq!(kapton.formula, "C22H10N2O5");
assert_eq!(kapton.density, 1.42);

// Density is taken from the table unless you override it
let mu = db.material_mu_named("kapton", &[10_000.0], CrossSectionKind::Total, None)?;
# Ok::<(), xraydb::XrayDbError>(())
```

### Optics (feature `optics`)

Crystal Darwin widths and mirror/multilayer reflectivity. Parameter structs with `Default` keep call sites readable.

```toml
[dependencies]
xraydb = { version = "0.4", features = ["optics"] }
```

```rust
# #[cfg(feature = "optics")] {
use xraydb::{DarwinParams, MirrorParams, XrayDb};

let db = XrayDb::try_new()?;

let dw = db.darwin_width(DarwinParams {
    energy: 10_000.0,
    crystal: "Si",
    hkl: (1, 1, 1),
    ..Default::default()
})?.expect("Si(111) diffracts at 10 keV");
println!("Darwin width: {:.2} eV FWHM", dw.energy_fwhm);

let angles: Vec<f64> = (1..100).map(|i| i as f64 * 0.1e-3).collect();
let refl = db.mirror_reflectivity(MirrorParams {
    formula: "Pt",
    theta: &angles,
    energy: 10_000.0,
    density: 21.45,
    roughness: 5.0,
    ..Default::default()
})?;
# }
# Ok::<(), xraydb::XrayDbError>(())
```

### Accuracy of `f1_chantler`

f1 is evaluated with an interpolating natural cubic spline through the Chantler grid,
which reproduces upstream XrayDB to within 5e-12.

Upstream fits its spline to a *window* of the grid spanning the requested energies padded
by three points either side, so its answer depends on what else is in the same call —
`f1_chantler('Au', 11919)` alone gives −17.745813, while the same energy inside a wider
batch gives −17.769546. Fitting once, globally, is the limit that windowing converges to
and does not depend on the query.

Validated against 1,727 reference values from upstream spanning 30 elements and every
tabulated absorption edge (`xraydb-lib/tests/data/f1_chantler_reference.csv`).

One practical difference: caesium's grid repeats 11.4 eV, which makes upstream raise
`ValueError: x must be strictly increasing`. Here the spline is fitted to the strictly
increasing subsequence, so every element stays queryable.

### Energy clamping

Following upstream XrayDB, energies outside a table's range are **clamped** to its endpoints rather than rejected. Check first with `XrayDb::elam_energy_range()` and `db.chantler_energy_range(element)`.

## Command-line tool

```sh
cargo install --path xraydb-cli    # installs `xraydb`

xraydb element Fe
xraydb edges Fe
xraydb lines Cu --level K
xraydb mu H2O --density 1.0 --energy 5000:15000/11
xraydb delta-beta SiO2 --density 2.2 --energy 10000
xraydb guess-edge --energy 7100
xraydb materials
```

Energies accept a single value, a comma list, `start:stop:step`, or `start:stop/count` (log-spaced). Add `--json` or `--csv` for machine-readable output.

### For scripts and agents

The CLI is designed to be driven programmatically as well as typed:

```sh
xraydb commands --json      # every subcommand, argument, enum value, and default
```

- `--json` works on **every** subcommand, including `commands` itself, so an agent can
  discover the whole interface in one call instead of parsing `--help` prose.
- With `--json`, **failures are also JSON**, written to stderr as
  `{"error": "...", "context": [...]}`. Without it, errors stay plain text.
- On failure stdout is left empty, so a consumer parsing stdout never sees partial output.
- Exit code is `0` on success and `1` on error. A broken pipe (`| head`) exits `0`.
- Colour is only emitted to a TTY, and `NO_COLOR` disables it.

## npm package

The WASM bindings ship on npm as [`xraydb-wasm`](https://www.npmjs.com/package/xraydb-wasm) —
self-contained (module + database + typed loader):

```sh
npm install xraydb-wasm
```

```js
import initXraydb from 'xraydb-wasm/loader.mjs';
const xraydb = await initXraydb();
xraydb.atomic_number('Fe');   // 26
```

Publish a new version with `./xraydb-wasm/build-pkg.sh && cd web/pkg && npm publish`
(the version follows the workspace's Cargo version).

## Browser demo

A zero-dependency demo page with a periodic-table selector and live cross-section plots:

```sh
./xraydb-wasm/build-pkg.sh     # wasm-pack + bundles the data blob and loader into web/pkg
python3 -m http.server -d web 8080
```

Then open <http://localhost:8080>. See [web/README.md](web/README.md).

## Workspace Structure

| Crate | Description |
|-------|-------------|
| `xraydb-data` | Shared serde data model (`#![no_std]`) |
| `xraydb-generate` | Binary that parses raw data sources into compressed binary format |
| `xraydb-lib` | Main library crate (`xraydb`) with embedded compressed data |
| `xraydb-wasm` | WASM bindings via `wasm-bindgen` |
| `xraydb-cli` | `xraydb` command-line tool |
| `web/` | Browser demo built on `xraydb-wasm` |

## Performance

Typical timings on an M-series Mac, release build (see `cargo bench`):

| Operation | Time |
|---|---|
| `f1_chantler_at` | ~100 ns |
| `f2_chantler_at` | ~120 ns |
| `mu_elam_at` | ~113 ns |
| `xray_edge` | ~78 ns |
| `xray_delta_beta` | ~1.1 µs |
| `mu_elam` batch of 200 | ~7 µs |

## Keeping the data out of your binary

The database is compiled in by default, which costs about 3 MB. Turn the
`embedded-data` feature off and supply the bytes at runtime instead:

```toml
[dependencies]
xraydb = { version = "0.4", default-features = false, features = ["zstd"] }
```

```rust,ignore
// Ship data/xraydb.bin.zst alongside your binary, or fetch it over the network.
let bytes = std::fs::read("xraydb.bin.zst")?;
let db = xraydb::XrayDb::load_compressed(&bytes)?;
```

Measured on a release binary that actually queries the database: **3.72 MB → 0.67 MB,
82% smaller**. This matters most for WebAssembly, where the blob dominates the `.wasm`,
and for embedded targets.

`load_uncompressed` takes an already-decompressed postcard blob, letting you drop
`ruzstd` entirely (`default-features = false` with no `zstd` feature) and decompress
however you like.

The database is global and initialised once — the first successful load wins, and later
calls return that same database rather than replacing it. `XrayDb::current()` returns
whatever is loaded (falling back to the embedded blob when that feature is on), which is
what the crate's free functions use.

### Features

| Feature | Default | Effect |
|---|---|---|
| `embedded-data` | on | Compiles the ~3 MB database in; enables `XrayDb::new`/`try_new` |
| `zstd` | on | zstd decompression; needed by `embedded-data` and `load_compressed` |
| `optics` | off | Darwin widths, mirror and multilayer reflectivity |

## Minimum supported Rust version

**1.87.** The crate's own code compiles on 1.85 (edition 2024's floor); the extra two
versions come from `ruzstd`, which decompresses the embedded database. `xraydb-data`,
which has no such dependency, declares 1.85.

MSRV is checked in CI against the version declared in `Cargo.toml`, so it cannot drift.

## Development

```sh
cargo test --all-features
cargo clippy --all-targets --all-features -- -D warnings
cargo fmt --all --check
```

Install the pre-commit hook (runs fmt, clippy, and tests):

```sh
git config core.hooksPath .githooks
```

The library contains no `unsafe` (`#![forbid(unsafe_code)]`) and no `unwrap`/`expect`
outside test code (lint-enforced).

## Regenerating the data

```sh
git clone https://github.com/xraypy/XrayDB.git XrayDB
cargo run -p xraydb-generate --release
```

## Attribution

The X-ray data used in this project comes from the [XrayDB](https://github.com/xraypy/XrayDB) project by Matt Newville et al., which is placed in the public domain (CC0 1.0). This Rust port is an independent reimplementation.

Data sources include:

- **Elam/Ravel/Sieber tables** — photoabsorption, scattering, and emission line data
- **Chantler tables** — anomalous scattering factors (f', f'') and mass attenuation coefficients
- **Waasmaier-Kirfel coefficients** — elastic (Thomson) scattering factors f0

## License

Dual-licensed under MIT and Apache-2.0. See [LICENSE-MIT](LICENSE-MIT) and [LICENSE-APACHE](LICENSE-APACHE).

The underlying X-ray data is in the public domain (CC0 1.0) from the [XrayDB](https://github.com/xraypy/XrayDB) project.
