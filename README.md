# FELIX

**FELIX** (Full-cohort Efficient Local ancestry-Integrated miXed-model
framework) is a scalable, local-ancestry-aware GWAS framework for biobank cohorts. 
It has two components:

- **FELIXla** — a compact, streamable storage format that packs phased
  genotypes with per-haplotype local-ancestry calls, giving efficient access
  to ancestry-specific dosages at biobank scale.
- **FELIXassoc** — the association module. From one null mixed-model fit it
  tests, at every variant, both a shared-effect (homogeneous) and an
  ancestry-specific (heterogeneous) model and combines them with a Cauchy
  combination test, while controlling for relatedness and case/control
  imbalance (built on the SAIGE mixed-model core).

> **📘 Start with the [step-by-step tutorial](tutorial/README.md)** — a
> runnable end-to-end example (install → pack → Step 1 → Step 2 → reading the
> output) on a tiny simulated cohort.

## Why FELIX

- Retains admixed and unassigned participants that discrete-ancestry GWAS drops.
- Estimates per-ancestry effects *and* the shared effect from a single run.
- Scales to biobank cohorts via a sparse GRM + variance-ratio approximation.
- Calibrated under severe case/control imbalance (saddlepoint approximation).

## Install

Two supported ways. **Docker is recommended** — it is one command and pins the
whole software stack. The source install exists for machines where Docker is
unavailable (many HPC login nodes) or unwanted, and works on both Linux and
macOS.

### Option A — Docker (recommended)

```bash
docker pull lhu1/felix:latest
```

The image is multi-architecture, so this gives a native build on `linux/amd64`
(servers, HPC, Intel Macs) and on `linux/arm64` (Apple Silicon). Add
`--platform linux/amd64` to pin the architecture.

On HPC, convert it once to a Singularity/Apptainer image:

```bash
singularity pull FELIX_latest.sif docker://lhu1/felix:latest
```

### Option B — source install with pixi (Linux and macOS)

[pixi](https://pixi.sh) is a single-binary package manager that builds a
self-contained environment inside this checkout. You do **not** need conda, R,
a compiler, or admin rights beforehand — pixi fetches all of it.

```bash
# 1. install pixi (once per machine), then reopen your terminal
curl -fsSL https://pixi.sh/install.sh | sh

# 2. get FELIX and build it
git clone https://github.com/ZhouLabGenetics/FELIX.git
cd FELIX
pixi run setup          # ~10 min the first time; installs R, htslib, FELIX
pixi run check          # verifies every command is present and runnable

# 3. run everything through pixi (or `pixi shell` once, then plain commands)
pixi run felixla --help
```

A prebuilt `felixla` for `linux-64` and `osx-arm64` ships in
`tools/felixla/prebuilt/`; on any other platform `pixi run setup` compiles it
from source automatically. See the
[tutorial's installation section](tutorial/README.md#2-installation) for the
full walkthrough, including what to do if you have never used pixi.

## Quick start

```bash
git clone https://github.com/ZhouLabGenetics/FELIX.git # if you have already downloaded it from pixi installation, no need to re-download
cd FELIX/tutorial/
python3 simulate_example.py     # make a tiny simulated admixed cohort
bash run_example.sh             # FELIXla pack -> Step 1 -> Step 2
```

`run_example.sh` uses the Docker image by default, and the source install
automatically when the FELIX commands are already on `PATH` (i.e. inside
`pixi shell`). From the repository root, `pixi run test` runs the same example
against the source install.

This produces `tutorial/output/step2_results.tsv` with per-ancestry and joint
p-values. See [`tutorial/README.md`](tutorial/README.md) for a full, annotated
walkthrough (inputs, each step, and how to read the output).

## Pipeline at a glance

1. **Pack** phased genotypes + local-ancestry calls into FELIXla with `felixla`
   (`--flare-vcf` for FLARE, `--rfmix-msp` for RFMix `.msp.tsv`,
   `--tractor-dosage-vcf` for a TRACTOR dosage VCF).
2. **Step 1** — fit the null GLMM once per cohort/trait (`step1_fitNULLGLMM.R`).
3. **Step 2** — ancestry-aware association test (`step2_SPAtests.R --FELIXlaPrefix ... --is_admixed=TRUE`).

## Repository layout

- `R/`, `src/` — the FELIX R package (FELIXla reader + FELIXassoc tests).
- `tools/felixla/` — the FELIXla module: sources for the single `felixla`
  command (packing, RFMix/TRACTOR conversion, export, query, region subsetting),
  its test fixtures, and prebuilt binaries.
- `tools/` — the pixi source-install scripts (`install_felix.sh`,
  `check_install.sh`).
- `extdata/` — the step wrapper scripts.
- `docker/` — Dockerfile and build/publish instructions.
- `tutorial/` — a runnable end-to-end example on simulated data.
- `pixi.toml` / `pixi.lock` — the pinned source-install environment.

## Citation

FELIX — Hu, L. et al. (in submission).

Questions or bugs: open an issue at <https://github.com/ZhouLabGenetics/FELIX/issues>.
