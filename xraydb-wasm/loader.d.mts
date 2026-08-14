// Type definitions for loader.mjs — the one-call entry point.

export * from './xraydb_wasm.js';

import * as xraydb from './xraydb_wasm.js';

/**
 * Initialise the wasm module and load the database that ships with this package.
 *
 * @param dataUrl Where to fetch `xraydb.bin.zst` from. Defaults to the copy
 *                bundled alongside this module, resolved via `import.meta.url`.
 * @returns The full xraydb API namespace, ready to query.
 */
export default function initXraydb(dataUrl?: string | URL): Promise<typeof xraydb>;
