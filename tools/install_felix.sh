#!/usr/bin/env bash
# Build and install the FELIX R package (FELIXla reader + FELIXassoc tests)
# and the command-line step wrappers into the active pixi environment.
#
# Run through pixi:
#
#     pixi run install-felix
set -euo pipefail

if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "ERROR: CONDA_PREFIX is not set. Run this through pixi:" >&2
    echo "         pixi run install-felix" >&2
    exit 1
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="$CONDA_PREFIX/bin"
cd "$REPO_ROOT"

WRAPPERS=(step1_fitNULLGLMM.R step2_SPAtests.R step3_LDmat.R createSparseGRM.R)
INSTALLED_SO="$CONDA_PREFIX/lib/R/library/SAIGE/libs/SAIGE.so"

# Skip the (~1 minute) package rebuild when the installed shared object is
# already newer than every package source file. Wrappers are checked separately
# below: a missing or updated wrapper must not force an unnecessary C++ rebuild.
# `pixi run test` / `pixi run check` depend on this task, so they should be
# cheap once the environment is built. FELIX_FORCE=1 rebuilds unconditionally.
package_up_to_date() {
    [[ "${FELIX_FORCE:-0}" == "1" ]] && return 1
    [[ -f "$INSTALLED_SO" ]] || return 1
    [[ -z "$(find src R DESCRIPTION NAMESPACE \
                  -newer "$INSTALLED_SO" -print -quit 2>/dev/null)" ]]
}

wrappers_up_to_date() {
    [[ "${FELIX_FORCE:-0}" == "1" ]] && return 1
    for w in "${WRAPPERS[@]}"; do
        [[ -x "$DEST/$w" && "$DEST/$w" -nt "extdata/$w" ]] || return 1
    done
}

# --- prebuilt R package -----------------------------------------------------
# binaries/FELIX_<version>_<platform>.tgz is the already-compiled package, built
# by tools/build_r_binary.sh against this same pixi environment. Its shared
# object carries a relocatable rpath, so it works from any checkout.
case "$(uname -s)/$(uname -m)" in
    Linux/x86_64)              R_PLATFORM=linux-x86_64 ;;
    Linux/aarch64|Linux/arm64) R_PLATFORM=linux-aarch64 ;;
    Darwin/arm64)              R_PLATFORM=macos-arm64 ;;
    Darwin/x86_64)             R_PLATFORM=macos-x86_64 ;;
    *)                         R_PLATFORM=unknown ;;
esac

install_prebuilt() {
    local tgz
    tgz="$(ls -1 "binaries/FELIX_"*"_${R_PLATFORM}.tgz" 2>/dev/null | head -n1 || true)"
    [[ -n "$tgz" ]] || return 1

    echo "==> Installing the prebuilt FELIX R package"
    echo "    $tgz"
    if ! Rscript -e "install.packages('$tgz', repos = NULL, type = 'source')" >/dev/null 2>&1; then
        echo "    the prebuilt package would not install here; compiling instead."
        return 1
    fi
    # A binary built against a different R or a different dependency build would
    # install but fail to load, so prove it loads before accepting it.
    if ! Rscript -e 'suppressMessages(library(SAIGE))' >/dev/null 2>&1; then
        echo "    the prebuilt package does not load here; compiling instead."
        return 1
    fi
    echo "    ok (skipped the source build)"
    return 0
}

if package_up_to_date; then
    echo "==> FELIX R package is up to date (FELIX_FORCE=1 to rebuild)"
else
    # lintools is used by the null-model fit (lintools::pinv) but is not packaged
    # for conda, so it comes from CRAN into the environment's own R library.
    echo "==> Checking R dependencies"
    Rscript -e 'if (!requireNamespace("lintools", quietly = TRUE)) {
                  install.packages("lintools", repos = "https://cloud.r-project.org")
                }
                if (!requireNamespace("lintools", quietly = TRUE)) {
                  stop("failed to install lintools")
                }
                cat("    lintools", as.character(packageVersion("lintools")), "ok\n")'

    # Prefer the prebuilt package: it installs in seconds, where compiling from
    # source takes several minutes and needs >= 2 GB of RAM (which is what kills
    # the build on cluster login nodes). Set FELIX_FROM_SOURCE=1 to skip it.
    if [[ "${FELIX_FROM_SOURCE:-0}" != "1" ]] && install_prebuilt; then
        :
    else
        # Object files left over from an earlier build (possibly with different
        # include paths) would be silently reused by R CMD INSTALL.
        echo "==> Compiling the FELIX R package from source"
        echo "    (this needs >= 2 GB of RAM and a few minutes; on a cluster run"
        echo "     it on a compute node, not the login node)"
        rm -f src/*.o src/*.so
        R CMD INSTALL --preclean .
    fi
fi

if wrappers_up_to_date; then
    echo "==> FELIX step wrappers are up to date"
else
    echo "==> Installing the step wrappers into $DEST"
    mkdir -p "$DEST"
    for script in "${WRAPPERS[@]}"; do
        # The copies in extdata/ carry the Docker image's shebang; rewrite it so
        # the installed command runs under this environment's Rscript.
        {
            echo '#!/usr/bin/env Rscript'
            tail -n +2 "extdata/$script"
        } > "$DEST/$script"
        chmod 0755 "$DEST/$script"
        echo "    $script"
    done
fi

echo "==> Done. Activate the environment with 'pixi shell', or prefix commands with 'pixi run'."
