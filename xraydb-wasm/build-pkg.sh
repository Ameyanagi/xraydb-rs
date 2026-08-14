#!/bin/sh
# Build the distributable xraydb-wasm package.
#
# wasm-pack has no mechanism for shipping extra assets, so this wraps it: after the
# build it copies the database blob and the loader into pkg/ and lists them in
# package.json, making the package self-contained — an npm consumer gets the data
# and a one-call `initXraydb()` without sourcing anything else.
#
# Usage: xraydb-wasm/build-pkg.sh [out-dir]     (default: ../web/pkg, for the demo)

set -eu
cd "$(dirname "$0")"
OUT="${1:-../web/pkg}"

wasm-pack build --target web --release --out-dir "$OUT"

cp ../xraydb-lib/data/xraydb.bin.zst "$OUT/xraydb.bin.zst"
cp loader.mjs "$OUT/loader.mjs"
cp loader.d.mts "$OUT/loader.d.mts"

# Register the extras in package.json so `npm publish` from pkg/ would include them.
node - "$OUT/package.json" <<'EOF'
const fs = require('fs');
const path = process.argv[2];
const pkg = JSON.parse(fs.readFileSync(path, 'utf8'));
for (const f of ['xraydb.bin.zst', 'loader.mjs', 'loader.d.mts']) {
    if (!pkg.files.includes(f)) pkg.files.push(f);
}
pkg.homepage = 'https://github.com/Ameyanagi/xraydb-rs';
pkg.keywords = ['xray', 'physics', 'spectroscopy', 'xraydb', 'wasm', 'webassembly',
                'absorption', 'crystallography'];
fs.writeFileSync(path, JSON.stringify(pkg, null, 2) + '\n');
EOF

echo "package ready in $OUT:"
ls -l "$OUT" | awk 'NR>1 {printf "  %-28s %s bytes\n", $NF, $5}'
