#!/usr/bin/env bash
# MAINTAINER SCRIPT -- build the prebuilt FELIX R package that ships in
# binaries/, so users do not have to compile it.
#
#   pixi run bash tools/build_r_binary.sh
#
# Compiling the R package is the slow, memory-hungry part of a source install
# (~26 C++ files, >=2 GB RAM, several minutes). The tarball produced here is the
# already-compiled package, installed straight into the pixi environment's R
# library in a couple of seconds.
#
# The shared object is linked against libraries in the pixi environment. R's
# build records the build machine's absolute path as the only rpath, which
# would break on any other checkout, so we rewrite it to a path relative to the
# object itself:
#
#   <env>/lib/R/library/SAIGE/libs/SAIGE.so   ->  ../../../..  ==  <env>/lib
#
# The binary is therefore tied to the pixi environment (same R version, same
# pinned dependency builds) but not to any particular directory.
set -euo pipefail

if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "ERROR: run this through pixi (pixi run bash tools/build_r_binary.sh)" >&2
    exit 1
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

VERSION="$(awk '/^Version:/{print $2}' DESCRIPTION)"
case "$(uname -s)/$(uname -m)" in
    Linux/x86_64)              PLATFORM=linux-x86_64 ;;
    Linux/aarch64|Linux/arm64) PLATFORM=linux-aarch64 ;;
    Darwin/arm64)              PLATFORM=macos-arm64 ;;
    Darwin/x86_64)             PLATFORM=macos-x86_64 ;;
    *) echo "ERROR: unsupported platform $(uname -s)/$(uname -m)" >&2; exit 1 ;;
esac

STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT
OUT_DIR="$REPO_ROOT/binaries"
OUT="$OUT_DIR/FELIX_${VERSION}_${PLATFORM}.tgz"
mkdir -p "$OUT_DIR"

echo "==> compiling the FELIX R package (version $VERSION, $PLATFORM)"
rm -f src/*.o src/*.so
R CMD INSTALL --preclean --library="$STAGE" .

SO="$STAGE/SAIGE/libs/SAIGE.so"
[[ -f "$SO" ]] || { echo "ERROR: $SO not produced" >&2; exit 1; }

echo "==> making SAIGE.so relocatable"
if [[ "$(uname -s)" == "Darwin" ]]; then
    while read -r abs; do
        [[ -n "$abs" ]] && install_name_tool -delete_rpath "$abs" "$SO" 2>/dev/null || true
    done < <(otool -l "$SO" | awk '/LC_RPATH/{f=1} f&&/path /{print $2; f=0}' | grep '^/' || true)
    install_name_tool -add_rpath '@loader_path/../../../..' "$SO"
    otool -l "$SO" | awk '/LC_RPATH/{f=1} f&&/path /{print "    rpath: " $2; f=0}'
elif command -v patchelf >/dev/null 2>&1; then
    patchelf --set-rpath '$ORIGIN/../../../..' "$SO"
    readelf -d "$SO" | grep -E 'RPATH|RUNPATH' | sed 's/^/    /'
else
    echo "ERROR: patchelf is required on Linux to make the binary relocatable." >&2
    echo "       Add 'patchelf' to pixi.toml, or install it, then re-run." >&2
    exit 1
fi

echo "==> packing $OUT"
rm -f "$OUT"
tar -czf "$OUT" -C "$STAGE" SAIGE
ls -lh "$OUT" | awk '{print "    " $5, $9}'

# Verify from the *installed* location, not the staging directory: the rpath we
# just wrote is relative to <env>/lib/R/library/SAIGE/libs, so it only resolves
# once the package sits there.
echo "==> verifying the packed package installs and loads"
Rscript -e "install.packages('$OUT', repos = NULL, type = 'source')" >/dev/null 2>&1 \
    || { echo "ERROR: the packed package does not install" >&2; exit 1; }
Rscript -e 'suppressMessages(library(SAIGE)); cat("    loaded SAIGE", as.character(packageVersion("SAIGE")), "\n")' \
    || { echo "ERROR: the packed package installs but does not load" >&2; exit 1; }

echo
echo "Done. Commit binaries/ and users will skip the compile entirely."
