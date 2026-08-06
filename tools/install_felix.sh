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

# Skip the (~1 minute) rebuild when the installed artefacts are already newer
# than every source file. `pixi run test` / `pixi run check` depend on this
# task, so they should be cheap once the environment is built.
# FELIX_FORCE=1 rebuilds unconditionally.
up_to_date() {
    [[ "${FELIX_FORCE:-0}" == "1" ]] && return 1
    [[ -f "$INSTALLED_SO" ]] || return 1
    for w in "${WRAPPERS[@]}"; do
        [[ -x "$DEST/$w" ]] || return 1
    done
    # Any source newer than the installed shared object means a rebuild.
    [[ -z "$(find src R extdata DESCRIPTION NAMESPACE tools/install_felix.sh \
                  -newer "$INSTALLED_SO" -print -quit 2>/dev/null)" ]]
}

if up_to_date; then
    echo "==> FELIX R package and step wrappers are up to date (FELIX_FORCE=1 to rebuild)"
    exit 0
fi

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

# Object files left over from an earlier build (possibly with different
# include paths) would be silently reused by R CMD INSTALL.
echo "==> Building the FELIX R package"
rm -f src/*.o src/*.so
R CMD INSTALL --preclean .

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

echo "==> Done. Activate the environment with 'pixi shell', or prefix commands with 'pixi run'."
