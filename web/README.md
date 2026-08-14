# Browser demo

A dependency-free demo of [`xraydb-wasm`](../xraydb-wasm): periodic-table element
selector, edge and emission-line tables, and live log-log cross-section plots.

## Build and run

From the repository root:

```sh
wasm-pack build --target web --release xraydb-wasm --out-dir ../web/pkg
python3 -m http.server -d web 8080
```

The WASM module is built without the embedded database (~350 KB instead of ~3.5 MB);
the page fetches `data/xraydb.bin.zst` in parallel with it. `web/data/xraydb.bin.zst`
is a committed symlink into `xraydb-lib/data/`, so it needs no copy step and can never
lag a regenerated blob. Code and data are cached separately by the browser: editing
the demo re-downloads 350 KB, not 3.5 MB.

Deploying to a static host that does not follow symlinks? Replace the symlink with the
real file: `cp xraydb-lib/data/xraydb.bin.zst web/data/`.

Then open <http://localhost:8080>.

`pkg/` is generated and git-ignored. If it is missing, the page says so and prints the
command to run rather than failing silently.

![The demo showing Fe with its absorption edges and cross-section plot](screenshot.png)

## What it shows

- **Periodic table** — all 118 elements, coloured by series, keyboard navigable with the
  arrow keys. Elements without tabulated data are dimmed.
- **Element panel** — Z, molar mass, density, every absorption edge (energy, fluorescence
  yield, jump ratio), and the ten strongest emission lines.
- **Cross-section plot** — µ/ρ (total, photoelectric, coherent, incoherent) on a log axis
  plus f₁ and f₂ on a linear right-hand axis, over a user-set energy range. Absorption
  edges are drawn as labelled vertical rules. Hovering reads out the values at that
  energy.
- **Search and deep links** — type `Fe`, `iron`, or `26` in the search box; the selected
  element lives in the URL hash (`#Au`), so views are bookmarkable and shareable.
- **Compounds & materials** — pick from the built-in materials table or type any formula
  (`H2O`, `Fe2O3`, `Ru1wt%SiO2`); shows µ, attenuation length, δ/β, and the critical
  angle at a chosen energy, plus a µ-vs-energy plot. Invalid formulas get an inline
  error rather than a blank panel.

Everything is inline SVG and vanilla ES modules — no charting library, no CDN, no
bundler. The `.wasm` is ~3.5 MB because it embeds the entire compressed database, so the
page works offline once loaded.

## Notes

- Traces are computed with one batched WASM call per trace over a 600-point log-spaced
  grid, not one call per point.
- Energies outside a table's range are clamped by the library; the "Full table" preset
  spans the whole Elam range (100 eV – 800 keV).
