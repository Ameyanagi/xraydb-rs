// xraydb-rs browser demo — no dependencies, no build step beyond wasm-pack.
import init, * as xray from './pkg/xraydb_wasm.js';
import { LAYOUT } from './layout.js';

const SERIES = {
  alkali: 'Alkali metal',
  alkaline: 'Alkaline earth',
  transition: 'Transition metal',
  post: 'Post-transition metal',
  metalloid: 'Metalloid',
  nonmetal: 'Nonmetal',
  halogen: 'Halogen',
  noble: 'Noble gas',
  lanthanide: 'Lanthanide',
  actinide: 'Actinide',
};

const $ = (id) => document.getElementById(id);
const fmt = (v, sig = 4) => {
  if (!isFinite(v)) return '—';
  if (v === 0) return '0';
  const mag = Math.abs(v);
  if (mag >= 1e5 || mag < 1e-3) return v.toExponential(sig - 1);
  return Number(v.toPrecision(sig)).toString();
};

const state = {
  z: 26,
  emin: 1000,
  emax: 40000,
  traces: new Set(['total']),
  showEdges: true,
  ready: false,
};

// ── Boot ────────────────────────────────────────────────────────────────────

async function boot() {
  try {
    await init();
  } catch (err) {
    $('boot-error').hidden = false;
    $('boot-error-detail').textContent = String(err);
    $('app').hidden = true;
    return;
  }
  xray.init();

  const version = xray.data_version();
  $('meta').textContent =
    `${xray.element_count()} elements · XrayDB ${version ?? 'unknown'}`;

  buildPeriodicTable();
  wireControls();
  state.ready = true;
  select(26);
}

// ── Periodic table ──────────────────────────────────────────────────────────

function buildPeriodicTable() {
  const grid = $('ptable');
  const frag = document.createDocumentFragment();

  for (const [z, group, period, series] of LAYOUT) {
    const cell = document.createElement('button');
    cell.type = 'button';
    cell.className = `el el-${series}`;
    cell.style.gridColumn = group;
    cell.style.gridRow = period;
    cell.dataset.z = z;

    let sym = String(z);
    let name = '';
    try {
      sym = xray.symbol(String(z));
      name = xray.atomic_name(String(z));
    } catch {
      cell.classList.add('el-nodata');
    }

    cell.innerHTML = `<span class="el-z">${z}</span><span class="el-sym">${sym}</span>`;
    cell.title = `${name || sym} (Z=${z}) — ${SERIES[series] ?? series}`;
    cell.setAttribute('aria-label', cell.title);
    cell.addEventListener('click', () => select(z));
    frag.appendChild(cell);
  }

  // Spacer labels for the f-block rows.
  for (const [row, text] of [[9, '57–71'], [10, '89–103']]) {
    const label = document.createElement('span');
    label.className = 'fblock-label';
    label.style.gridColumn = 2;
    label.style.gridRow = row;
    label.textContent = text;
    frag.appendChild(label);
  }

  grid.appendChild(frag);

  grid.addEventListener('keydown', (ev) => {
    const step = { ArrowRight: 1, ArrowLeft: -1, ArrowUp: -18, ArrowDown: 18 }[ev.key];
    if (!step) return;
    ev.preventDefault();
    const next = Math.min(118, Math.max(1, state.z + step));
    select(next);
  });
}

function select(z) {
  state.z = z;
  for (const cell of document.querySelectorAll('.el')) {
    cell.classList.toggle('el-active', Number(cell.dataset.z) === z);
  }
  const active = document.querySelector('.el-active');
  if (active) active.focus({ preventScroll: true });
  render();
}

// ── Controls ────────────────────────────────────────────────────────────────

function wireControls() {
  $('emin').addEventListener('change', () => {
    state.emin = Math.max(1, Number($('emin').value) || 1000);
    render();
  });
  $('emax').addEventListener('change', () => {
    state.emax = Math.max(state.emin * 1.01, Number($('emax').value) || 40000);
    render();
  });
  for (const box of document.querySelectorAll('input[data-trace]')) {
    box.addEventListener('change', () => {
      box.checked ? state.traces.add(box.dataset.trace) : state.traces.delete(box.dataset.trace);
      render();
    });
  }
  $('show-edges').addEventListener('change', (ev) => {
    state.showEdges = ev.target.checked;
    render();
  });
  $('preset-full').addEventListener('click', () => {
    $('emin').value = state.emin = 100;
    $('emax').value = state.emax = 800000;
    render();
  });
  $('preset-hard').addEventListener('click', () => {
    $('emin').value = state.emin = 1000;
    $('emax').value = state.emax = 40000;
    render();
  });
}

// ── Rendering ───────────────────────────────────────────────────────────────

function render() {
  if (!state.ready) return;
  const id = String(state.z);

  let info;
  try {
    info = xray.element(id);
  } catch (err) {
    $('facts').innerHTML = `<p class="empty">No data for Z=${state.z}.</p>`;
    $('edges').innerHTML = '';
    $('lines').innerHTML = '';
    $('plot').innerHTML = '';
    return;
  }

  $('el-title').textContent = `${info.symbol} · ${info.name}`;
  $('facts').innerHTML = `
    <div><dt>Atomic number</dt><dd>${info.z}</dd></div>
    <div><dt>Molar mass</dt><dd>${fmt(info.molar_mass, 6)} g/mol</dd></div>
    <div><dt>Density</dt><dd>${fmt(info.density, 5)} g/cm³</dd></div>
  `;

  renderEdges(id);
  renderLines(id);
  renderPlot(id);
}

function renderEdges(id) {
  let edges = [];
  try {
    edges = xray.xray_edges(id);
  } catch {
    /* no edge data */
  }
  if (!edges.length) {
    $('edges').innerHTML = '<p class="empty">No tabulated edges.</p>';
    return;
  }
  const rows = edges
    .map(
      (e) => `<tr><th scope="row">${e.label}</th><td>${fmt(e.energy, 6)}</td>
              <td>${fmt(e.fluorescence_yield, 3)}</td><td>${fmt(e.jump_ratio, 3)}</td></tr>`,
    )
    .join('');
  $('edges').innerHTML = `<table>
    <caption>Absorption edges</caption>
    <thead><tr><th scope="col">Edge</th><th scope="col">Energy (eV)</th>
      <th scope="col">Fluor. yield</th><th scope="col">Jump ratio</th></tr></thead>
    <tbody>${rows}</tbody></table>`;
}

function renderLines(id) {
  let groups = [];
  try {
    groups = xray.xray_lines_by_level(id, undefined);
  } catch {
    /* no line data */
  }
  if (!groups.length) {
    $('lines').innerHTML = '<p class="empty">No tabulated emission lines.</p>';
    return;
  }

  // Elam intensities are normalised within each initial level, so lines are grouped
  // by level rather than ranked globally — a global sort would put weak L lines
  // above Ka1, which is not what a spectrum looks like.
  const body = groups
    .map((g) => {
      const head = `<tr class="group-row"><th scope="rowgroup" colspan="4">${g.level} lines</th></tr>`;
      const rows = g.lines
        .slice(0, 6)
        .map(
          (l) => `<tr><th scope="row">${l.label}</th><td>${fmt(l.energy, 6)}</td>
              <td>${fmt(l.intensity, 3)}</td><td>${l.initial_level}→${l.final_level}</td></tr>`,
        )
        .join('');
      return head + rows;
    })
    .join('');

  $('lines').innerHTML = `<table>
    <caption>Emission lines — intensity is relative within each level</caption>
    <thead><tr><th scope="col">Line</th><th scope="col">Energy (eV)</th>
      <th scope="col">Intensity</th><th scope="col">Transition</th></tr></thead>
    <tbody>${body}</tbody></table>`;
}

// Trace definitions: label, colour, and how to evaluate over an energy grid.
const TRACES = {
  total: { label: 'µ/ρ total', axis: 'mu', color: 'var(--c1)', get: (id, e) => xray.mu_elam(id, e, 'total') },
  photo: { label: 'µ/ρ photo', axis: 'mu', color: 'var(--c2)', get: (id, e) => xray.mu_elam(id, e, 'photo') },
  coherent: { label: 'µ/ρ coherent', axis: 'mu', color: 'var(--c3)', get: (id, e) => xray.mu_elam(id, e, 'coherent') },
  incoherent: { label: 'µ/ρ incoherent', axis: 'mu', color: 'var(--c4)', get: (id, e) => xray.mu_elam(id, e, 'incoherent') },
  f1: { label: "f₁ (f′)", axis: 'f', color: 'var(--c5)', get: (id, e) => xray.f1_chantler(id, e) },
  f2: { label: 'f₂ (f″)', axis: 'f', color: 'var(--c6)', get: (id, e) => xray.f2_chantler(id, e) },
};

const PLOT = { w: 900, h: 460, ml: 74, mr: 62, mt: 18, mb: 48 };

function logspace(a, b, n) {
  const la = Math.log(a);
  const lb = Math.log(b);
  return Array.from({ length: n }, (_, i) => Math.exp(la + ((lb - la) * i) / (n - 1)));
}

function renderPlot(id) {
  const svg = $('plot');
  const energies = logspace(state.emin, state.emax, 600);
  const grid = new Float64Array(energies);

  const series = [];
  for (const key of Object.keys(TRACES)) {
    if (!state.traces.has(key)) continue;
    const spec = TRACES[key];
    try {
      const values = Array.from(spec.get(id, grid));
      if (values.some((v) => isFinite(v))) series.push({ ...spec, key, values });
    } catch {
      /* element lacks data for this trace */
    }
  }

  if (!series.length) {
    svg.innerHTML = '<text x="50%" y="50%" class="empty-text" text-anchor="middle">No data for this element in the selected range.</text>';
    $('legend').innerHTML = '';
    $('readout').textContent = '';
    return;
  }

  const { w, h, ml, mr, mt, mb } = PLOT;
  const px = (e) =>
    ml + ((Math.log(e) - Math.log(state.emin)) / (Math.log(state.emax) - Math.log(state.emin))) * (w - ml - mr);

  // Two independent y-axes: log for µ/ρ (left), linear for f1/f2 (right).
  const muSeries = series.filter((s) => s.axis === 'mu');
  const fSeries = series.filter((s) => s.axis === 'f');

  const muVals = muSeries.flatMap((s) => s.values).filter((v) => isFinite(v) && v > 0);
  const muMin = muVals.length ? Math.min(...muVals) : 1;
  const muMax = muVals.length ? Math.max(...muVals) : 10;
  const muY = (v) =>
    v > 0
      ? h - mb - ((Math.log(v) - Math.log(muMin)) / (Math.log(muMax) - Math.log(muMin) || 1)) * (h - mt - mb)
      : null;

  const fVals = fSeries.flatMap((s) => s.values).filter(isFinite);
  const fMin = fVals.length ? Math.min(...fVals) : 0;
  const fMax = fVals.length ? Math.max(...fVals) : 1;
  const fY = (v) => (h - mb - ((v - fMin) / (fMax - fMin || 1)) * (h - mt - mb));

  const parts = [];

  // Vertical gridlines at decades.
  for (let d = Math.floor(Math.log10(state.emin)); d <= Math.ceil(Math.log10(state.emax)); d++) {
    for (const m of [1, 2, 5]) {
      const e = m * 10 ** d;
      if (e < state.emin || e > state.emax) continue;
      const x = px(e);
      parts.push(`<line class="grid" x1="${x}" y1="${mt}" x2="${x}" y2="${h - mb}"/>`);
      parts.push(`<text class="tick" x="${x}" y="${h - mb + 18}" text-anchor="middle">${e >= 1000 ? `${e / 1000}k` : e}</text>`);
    }
  }

  // Horizontal gridlines for the log µ axis.
  if (muSeries.length) {
    for (let d = Math.floor(Math.log10(muMin)); d <= Math.ceil(Math.log10(muMax)); d++) {
      const v = 10 ** d;
      const y = muY(v);
      if (y === null || y < mt || y > h - mb) continue;
      parts.push(`<line class="grid" x1="${ml}" y1="${y}" x2="${w - mr}" y2="${y}"/>`);
      parts.push(`<text class="tick" x="${ml - 8}" y="${y + 4}" text-anchor="end">1e${d}</text>`);
    }
  }
  if (fSeries.length) {
    for (let i = 0; i <= 4; i++) {
      const v = fMin + ((fMax - fMin) * i) / 4;
      const y = fY(v);
      parts.push(`<text class="tick" x="${w - mr + 8}" y="${y + 4}" text-anchor="start">${fmt(v, 3)}</text>`);
    }
  }

  // Absorption edges as vertical rules.
  if (state.showEdges) {
    let edges = [];
    try {
      edges = xray.xray_edges(id);
    } catch {
      /* none */
    }
    for (const edge of edges) {
      if (edge.energy < state.emin || edge.energy > state.emax) continue;
      const x = px(edge.energy);
      parts.push(`<line class="edge-rule" x1="${x}" y1="${mt}" x2="${x}" y2="${h - mb}"/>`);
      parts.push(`<text class="edge-label" x="${x + 3}" y="${mt + 12}">${edge.label}</text>`);
    }
  }

  // Data paths.
  for (const s of series) {
    const yOf = s.axis === 'mu' ? muY : fY;
    let d = '';
    let pen = false;
    for (let i = 0; i < energies.length; i++) {
      const y = yOf(s.values[i]);
      if (y === null || !isFinite(y)) {
        pen = false;
        continue;
      }
      d += `${pen ? 'L' : 'M'}${px(energies[i]).toFixed(2)},${y.toFixed(2)}`;
      pen = true;
    }
    parts.push(`<path class="trace" d="${d}" style="stroke:${s.color}"/>`);
  }

  // Axes.
  parts.push(`<line class="axis" x1="${ml}" y1="${h - mb}" x2="${w - mr}" y2="${h - mb}"/>`);
  parts.push(`<line class="axis" x1="${ml}" y1="${mt}" x2="${ml}" y2="${h - mb}"/>`);
  parts.push(`<text class="axis-label" x="${(ml + w - mr) / 2}" y="${h - 8}" text-anchor="middle">Energy (eV)</text>`);
  if (muSeries.length) {
    parts.push(`<text class="axis-label" transform="translate(16,${(mt + h - mb) / 2}) rotate(-90)" text-anchor="middle">µ/ρ (cm²/g)</text>`);
  }
  if (fSeries.length) {
    parts.push(`<text class="axis-label" transform="translate(${w - 12},${(mt + h - mb) / 2}) rotate(-90)" text-anchor="middle">f₁ / f₂ (e⁻/atom)</text>`);
  }

  // Crosshair, moved by the pointer handler below.
  parts.push(`<line id="crosshair" class="crosshair" x1="0" y1="${mt}" x2="0" y2="${h - mb}" style="display:none"/>`);
  parts.push(`<rect id="hit" x="${ml}" y="${mt}" width="${w - ml - mr}" height="${h - mt - mb}" fill="transparent"/>`);

  svg.setAttribute('viewBox', `0 0 ${w} ${h}`);
  svg.innerHTML = parts.join('');

  $('legend').innerHTML = series
    .map((s) => `<span class="key"><i style="background:${s.color}"></i>${s.label}</span>`)
    .join('');

  attachCrosshair(svg, energies, series, px);
}

function attachCrosshair(svg, energies, series, px) {
  const hit = svg.querySelector('#hit');
  const hair = svg.querySelector('#crosshair');
  const readout = $('readout');
  if (!hit || !hair) return;

  const move = (ev) => {
    const rect = svg.getBoundingClientRect();
    const sx = ((ev.clientX - rect.left) / rect.width) * PLOT.w;
    hair.setAttribute('x1', sx);
    hair.setAttribute('x2', sx);
    hair.style.display = '';

    // Nearest sample to the pointer.
    let best = 0;
    let bestD = Infinity;
    for (let i = 0; i < energies.length; i++) {
      const d = Math.abs(px(energies[i]) - sx);
      if (d < bestD) {
        bestD = d;
        best = i;
      }
    }
    const bits = series.map((s) => `${s.label} ${fmt(s.values[best], 4)}`);
    readout.textContent = `E = ${fmt(energies[best], 6)} eV — ${bits.join(' · ')}`;
  };

  hit.addEventListener('pointermove', move);
  hit.addEventListener('pointerleave', () => {
    hair.style.display = 'none';
    readout.textContent = '';
  });
}

boot();
