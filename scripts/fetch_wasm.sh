#!/usr/bin/env bash
# Fetch the pinned RhizAeon wasm-pack bundle from the HyphAeon monorepo
# release into docs/pkg/.
#
# RHIZAEON_WASM_TAG is the SINGLE source of truth for the engine version
# served by the GitHub Pages app — .github/workflows/deploy-pages.yml
# runs this script on every deploy, and local devs run it once after
# cloning (or after bumping the pin).
#
# Bump by editing the tag below; the next Pages deploy serves the new
# bundle. The app code never changes — it always imports ./pkg/.
set -euo pipefail

RHIZAEON_WASM_TAG="${RHIZAEON_WASM_TAG:-rhizaeon-v0.1.0-beta.2}"
PKG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/docs/pkg"
ASSET="rhizaeon-wasm-${RHIZAEON_WASM_TAG##*-v}-pkg.tar.gz"
URL="https://github.com/veg/HyphAeon/releases/download/${RHIZAEON_WASM_TAG}/${ASSET}"

echo "[*] Fetching ${ASSET} from ${URL}"
# Wipe first: a plain tar overlay would leave stale files behind when the
# previous bundle had entries the new version no longer ships.
rm -rf "$PKG_DIR"
mkdir -p "$PKG_DIR"
curl -fsSL "$URL" | tar -xz -C "$PKG_DIR"
echo "[+] docs/pkg now serves ${RHIZAEON_WASM_TAG}"
