# xraydb-wasm

X-ray reference data for the elements, in the browser. WebAssembly bindings for the
[`xraydb`](https://crates.io/crates/xraydb) Rust crate — a port of
[XrayDB](https://github.com/xraypy/XrayDB): absorption edges, emission lines,
Elam/Chantler cross-sections, scattering factors, and compound-material calculations
for all 118 elements.

The package is self-contained: it ships a ~350 KB WebAssembly module and the ~3 MB
compressed database as separate files, downloaded in parallel and cached
independently by the browser. No network access after load; all computation is local.

## Usage

```js
import initXraydb from 'xraydb-wasm/loader.mjs';

const xraydb = await initXraydb();

xraydb.atomic_number('Fe');                            // 26
xraydb.xray_edge_energy('Fe', 'K');                    // 7112 (eV)
xraydb.mu_elam('Fe', new Float64Array([10000]), 'total');   // cm²/g
xraydb.material_mu('H2O', 1.0, new Float64Array([10000]), 'total'); // 1/cm
xraydb.xray_delta_beta('SiO2', 2.2, 10000);            // { delta, beta, attenuation_length_cm }
```

Works in plain `<script type="module">` pages and with bundlers (Vite, webpack 5,
Rollup) — the database blob is resolved via `import.meta.url`, so bundlers include it
as an asset automatically. To host the data elsewhere (e.g. a CDN), pass a URL:
`initXraydb('https://cdn.example.com/xraydb.bin.zst')`.

TypeScript definitions are included for the full API and the loader.

## Lower-level use

Skip the loader and drive the module directly:

```js
import init, * as xraydb from 'xraydb-wasm';

await init();
xraydb.init();                        // panic hook
xraydb.load_database(bytes);          // contents of xraydb.bin.zst
```

Every query fails with a "no database loaded" error until `load_database` succeeds.

## Data and licensing

Code is MIT OR Apache-2.0. The underlying X-ray data comes from the
[XrayDB](https://github.com/xraypy/XrayDB) project (public domain, CC0 1.0):
Elam/Ravel/Sieber tables, Chantler tables, and Waasmaier–Kirfel coefficients.

Built from [Ameyanagi/xraydb-rs](https://github.com/Ameyanagi/xraydb-rs); a live demo
lives in that repository's `web/` directory.
