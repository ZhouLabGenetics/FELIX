#!/usr/bin/env bash
# Verify that a source (pixi) install of FELIX is complete and working.
#
#     pixi run check
set -euo pipefail

fail=0
ok()   { printf '  \033[32mok\033[0m    %s\n' "$1"; }
bad()  { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; fail=1; }

echo "FELIX install check"
echo "  environment: ${CONDA_PREFIX:-<none>}"
echo

for cmd in bgzip tabix Rscript; do
    if command -v "$cmd" >/dev/null 2>&1; then ok "$cmd  ($(command -v "$cmd"))"; else bad "$cmd not found"; fi
done

if command -v felixla >/dev/null 2>&1 && [[ "$(felixla --help 2>&1 || true)" == *"felixla [input flags]"* ]]; then
    ok "felixla  ($(command -v felixla))"
else
    bad "felixla missing or not runnable  --  run 'pixi run install-felixla'"
fi

for cmd in step1_fitNULLGLMM.R step2_SPAtests.R createSparseGRM.R step3_LDmat.R; do
    if command -v "$cmd" >/dev/null 2>&1; then ok "$cmd  ($(command -v "$cmd"))"; else bad "$cmd not found  --  run 'pixi run install-felix'"; fi
done

if Rscript -e 'suppressMessages(library(SAIGE)); cat("")' >/dev/null 2>&1; then
    ok "R package SAIGE (FELIX core) $(Rscript -e 'cat(as.character(packageVersion("SAIGE")))' 2>/dev/null)"
else
    bad "the FELIX R package does not load  --  run 'pixi run install-felix'"
fi

echo
if [[ $fail -eq 0 ]]; then
    echo "All checks passed. Next: cd tutorial && bash run_example.sh"
else
    echo "Some checks failed (see above)."
    exit 1
fi
