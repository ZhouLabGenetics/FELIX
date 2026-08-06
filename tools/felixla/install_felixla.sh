#!/usr/bin/env bash
# Install the FELIXla command line tool into the active pixi/conda environment.
#
# Run through pixi so that CONDA_PREFIX, the compiler and htslib are set up:
#
#     pixi run install-felixla
#
# The prebuilt binary for the current platform is used when one is shipped in
# tools/felixla/prebuilt/<platform>/ and it actually runs; otherwise felixla is
# compiled from tools/felixla/src (takes under a minute).
set -euo pipefail

if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "ERROR: CONDA_PREFIX is not set. Run this through pixi:" >&2
    echo "         pixi run install-felixla" >&2
    exit 1
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TOOL_DIR="$REPO_ROOT/tools/felixla"
DEST="$CONDA_PREFIX/bin"

# `felixla --help` prints the usage block and exits 0; this is the cheapest
# end-to-end "does this binary load and run here?" test.
runs_here() {
    local out
    out="$("$1" --help 2>&1 || true)"
    [[ "$out" == *"felixla [input flags]"* ]]
}

# --- which prebuilt directory matches this machine? -------------------------
case "$(uname -s)/$(uname -m)" in
    Linux/x86_64)              PLATFORM=linux-64 ;;
    Linux/aarch64|Linux/arm64) PLATFORM=linux-aarch64 ;;
    Darwin/arm64)              PLATFORM=osx-arm64 ;;
    Darwin/x86_64)             PLATFORM=osx-64 ;;
    *)                         PLATFORM=unknown ;;
esac
PREBUILT="$TOOL_DIR/prebuilt/$PLATFORM/felixla"

mkdir -p "$DEST"

use_prebuilt() {
    [[ -f "$PREBUILT" ]] || return 1
    # Copy first: the binary resolves libhts via an rpath relative to its own
    # location ($ORIGIN/../lib, @loader_path/../lib), so it can only be tested
    # once it sits in $CONDA_PREFIX/bin.
    install -m 0755 "$PREBUILT" "$DEST/felixla"
    if ! runs_here "$DEST/felixla"; then
        echo "  the prebuilt binary does not run here; falling back to a source build."
        return 1
    fi
    return 0
}

build_from_source() {
    # CXXFLAGS is set explicitly because the Makefile's default includes
    # -march=native. The rpath is set explicitly rather than relying on the
    # conda compiler wrappers, which differ between Linux and macOS.
    make -C "$TOOL_DIR" clean >/dev/null 2>&1 || true
    make -C "$TOOL_DIR" \
        CXX="${CXX:-c++}" \
        CXXFLAGS="-O3 -std=c++17 -Wall" \
        HTSLIB_CFLAGS="-I$CONDA_PREFIX/include" \
        HTSLIB_LIBS="-L$CONDA_PREFIX/lib -lhts -Wl,-rpath,$CONDA_PREFIX/lib"
    install -m 0755 "$TOOL_DIR/bin/felixla" "$DEST/felixla"
}

echo "Installing FELIXla into $DEST"
if use_prebuilt; then
    echo "  using the prebuilt $PLATFORM binary"
else
    echo "  compiling from source (platform: $PLATFORM)"
    build_from_source
fi

if ! runs_here "$DEST/felixla"; then
    echo "ERROR: felixla was installed but does not run." >&2
    exit 1
fi
echo "  ok: felixla"
