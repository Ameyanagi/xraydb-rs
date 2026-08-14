// One-call setup for the xraydb-wasm package: initialise the module, fetch the
// database blob that ships alongside it, and load it.
//
//     import initXraydb from 'xraydb-wasm/loader.mjs';
//     const xraydb = await initXraydb();
//     xraydb.atomic_number('Fe');   // 26
//
// The blob URL is resolved relative to this module (`import.meta.url`), which both
// plain `<script type="module">` pages and bundlers (Vite, webpack 5, Rollup)
// understand — bundlers will treat the .zst as an asset and include it. Pass an
// explicit URL to load the data from somewhere else.
//
// The .wasm module itself is built without the embedded database (~350 KB instead
// of ~3.5 MB), so code and data download in parallel and cache independently.

import init, * as xraydb from './xraydb_wasm.js';

export * from './xraydb_wasm.js';

/** Initialise the wasm module and load the database. Returns the API namespace. */
export default async function initXraydb(dataUrl) {
  const url = dataUrl ?? new URL('xraydb.bin.zst', import.meta.url);
  const [, resp] = await Promise.all([init(), fetch(url)]);
  if (!resp.ok) {
    throw new Error(`failed to fetch the xraydb database from ${url}: HTTP ${resp.status}`);
  }
  xraydb.init(); // panic hook
  xraydb.load_database(new Uint8Array(await resp.arrayBuffer()));
  return xraydb;
}
