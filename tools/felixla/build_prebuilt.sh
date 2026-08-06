#!/usr/bin/env bash
# MAINTAINER SCRIPT -- regenerate the prebuilt felixla binary that ships in
# tools/felixla/prebuilt/<platform>/.
#
# Users never need to run this: `pixi run setup` uses the prebuilt binary if it
# runs on their machine and otherwise compiles from source.
#
#   # native build for the current platform
#   pixi run bash tools/felixla/build_prebuilt.sh
#
#   # linux-64 build from any machine with Docker
#   see tools/felixla/README.md
#
# The binary is linked with an rpath relative to its own location
# ($ORIGIN/../lib on Linux, @loader_path/../lib on macOS). install_felixla.sh
# copies it into $CONDA_PREFIX/bin, where that resolves to the pixi
# environment's own libhts -- so the binary stays valid wherever the repository
# is checked out.
set -euo pipefail

if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "ERROR: run this through pixi (pixi run bash tools/felixla/build_prebuilt.sh)" >&2
    exit 1
fi

TOOL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# $ORIGIN has to survive two expansions on the way to the linker: make (which
# would read `$O` as an empty make variable) and the recipe's shell (which
# would read $ORIGIN as an empty shell variable). Hence `$$` for make and the
# single quotes for the shell.
LINUX_RPATH="-Wl,-rpath,'\$\$ORIGIN'/../lib"
case "$(uname -s)/$(uname -m)" in
    Linux/x86_64)              PLATFORM=linux-64;      RPATH="$LINUX_RPATH" ;;
    Linux/aarch64|Linux/arm64) PLATFORM=linux-aarch64; RPATH="$LINUX_RPATH" ;;
    Darwin/arm64)              PLATFORM=osx-arm64;     RPATH='-Wl,-rpath,@loader_path/../lib' ;;
    Darwin/x86_64)             PLATFORM=osx-64;        RPATH='-Wl,-rpath,@loader_path/../lib' ;;
    *) echo "ERROR: unsupported platform $(uname -s)/$(uname -m)" >&2; exit 1 ;;
esac

OUT="$TOOL_DIR/prebuilt/$PLATFORM"
mkdir -p "$OUT"

echo "==> building felixla for $PLATFORM"
make -C "$TOOL_DIR" clean >/dev/null 2>&1 || true
make -C "$TOOL_DIR" \
    CXX="${CXX:-c++}" \
    CXXFLAGS="-O3 -std=c++17 -Wall" \
    HTSLIB_CFLAGS="-I$CONDA_PREFIX/include" \
    HTSLIB_LIBS="-L$CONDA_PREFIX/lib -lhts" \
    LDFLAGS="${LDFLAGS:-} $RPATH"

echo "==> running the module's own smoke test"
make -C "$TOOL_DIR" test CXXFLAGS="-O3 -std=c++17 -Wall" >/dev/null

install -m 0755 "$TOOL_DIR/bin/felixla" "$OUT/felixla"
strip -S "$OUT/felixla" 2>/dev/null || true

# Drop the absolute -rpath entries the conda compiler wrappers add, so the
# binary carries only the relocatable one.
if [[ "$(uname -s)" == "Darwin" ]]; then
    while read -r abs; do
        [[ -n "$abs" ]] && install_name_tool -delete_rpath "$abs" "$OUT/felixla" 2>/dev/null || true
    done < <(otool -l "$OUT/felixla" | awk '/LC_RPATH/{f=1} f&&/path /{print $2; f=0}' | grep '^/' || true)
elif command -v patchelf >/dev/null 2>&1; then
    patchelf --set-rpath '$ORIGIN/../lib' "$OUT/felixla"
else
    echo "NOTE: patchelf not found; the build-time -L path stays in RPATH ahead of" >&2
    echo "      \$ORIGIN/../lib. Harmless (it will not exist on other machines)," >&2
    echo "      but install patchelf for a clean binary." >&2
fi

echo
ls -l "$OUT"
echo
echo "Runtime dependencies:"
if [[ "$(uname -s)" == "Darwin" ]]; then
    otool -L "$OUT/felixla"
    otool -l "$OUT/felixla" | awk '/LC_RPATH/{f=1} f&&/path /{print "  rpath: " $2; f=0}'
else
    ldd "$OUT/felixla" || true
    readelf -d "$OUT/felixla" | grep -E 'RPATH|RUNPATH' || true
fi
