# xraydb-rs

X-ray reference data for the elements in Rust. A pure-Rust port of the [XrayDB](https://github.com/xraypy/XrayDB) project.

## Attribution

The X-ray data used in this project comes from the [XrayDB](https://github.com/xraypy/XrayDB) project by Matt Newville et al., which is placed in the public domain (CC0 1.0). This Rust port is an independent reimplementation.

Data sources include:

- **Elam/Ravel/Sieber tables** — photoabsorption, scattering, and emission line data
- **Chantler tables** — anomalous scattering factors (f', f'') and mass attenuation coefficients
- **Waasmaier-Kirfel coefficients** — elastic (Thomson) scattering factors f0

## Workspace Structure

| Crate | Description |
|-------|-------------|
| `xraydb-data` | Shared serde data model (`#![no_std]`) |
| `xraydb-generate` | Binary that parses raw data sources into compressed binary format |
| `xraydb-lib` | Main library crate with embedded compressed data |
| `xraydb-wasm` | WASM bindings via `wasm-bindgen` |

## Usage

```rust
use xraydb::{XrayDb, CrossSectionKind};

let db = XrayDb::new();

// X-ray mass attenuation coefficient for Fe at 10 keV (Elam tables)
let energies = [10000.0]; // eV
let mu = db.mu_elam("Fe", &energies, CrossSectionKind::Total).unwrap();

// Material attenuation (compound formulas supported)
let mu = db.material_mu("Fe2O3", &energies, 5.26, CrossSectionKind::Total).unwrap();

// Anomalous scattering factors (Chantler tables)
let f1 = db.f1_chantler("Fe", &energies).unwrap();
let f2 = db.f2_chantler("Fe", &energies).unwrap();
```

## Features

- **`optics`** — Enables X-ray optics calculations (Darwin width, mirror reflectivity, multilayer reflectivity). Adds a dependency on `num-complex`.

```toml
[dependencies]
xraydb = { version = "0.1", features = ["optics"] }
```

## License

Dual-licensed under MIT and Apache-2.0. See [LICENSE-MIT](LICENSE-MIT) and [LICENSE-APACHE](LICENSE-APACHE).

The underlying X-ray data is in the public domain (CC0 1.0) from the [XrayDB](https://github.com/xraypy/XrayDB) project.
