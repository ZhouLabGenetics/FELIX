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
> runnable end-to-end example (install → pack genotype → Step 1 → Step 2 →
> reading the output) on a tiny simulated cohort.

## Why FELIX

- Retains admixed participants that discrete-global-ancestry GWAS drops.
- Estimates both per-local-ancestry effects and the shared effect from a single run.
- Scales to biobank cohorts.
- Calibrated under severe case/control imbalance.

## Install

Two supported ways. **Docker is recommended** — it is one command and pins the
whole software stack. Pixi installation exists for machines where Docker is
unavailable or unwanted, and works on both Linux and macOS.

### Option A — Docker (recommended)

```bash
docker pull lhu1/felix:latest
```

The image is multi-architecture, so this gives a native build on `linux/amd64`
and on `linux/arm64`. Add `--platform linux/amd64` to pin the architecture.

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
pixi run setup          # downloads the environment; installs prebuilt FELIX
pixi run check          # verifies every command is present and runnable

# 3. run everything through pixi (or `pixi shell` once, then plain commands)
pixi run felixla --help
```

Nothing is compiled on your machine: both the `felixla` command
(`tools/felixla/prebuilt/`) and the FELIX R package (`binaries/`) ship prebuilt
for `linux-x86_64` and `macos-arm64`, and each is verified to load before being
accepted. On any other platform `pixi run setup` falls back to compiling from
source, which needs ≥ 2 GB of RAM.

If your cluster runs an old distro (CentOS/RHEL 7, glibc 2.17), this still
works — the lock file pins that as the supported floor. If `curl` there is too
old for modern TLS, the tutorial has a `wget` fallback.

See the [tutorial's installation section](tutorial/README.md#2-installation) for
the full walkthrough, including what to do if you have never used pixi.

## Quick start

Clone the repository once, then use the block for your installation method.

```bash
git clone https://github.com/ZhouLabGenetics/FELIX.git
cd FELIX
```

### Docker (Linux or macOS)

```bash
docker pull lhu1/felix:latest
cd tutorial
bash run_example.sh             # simulate -> FELIXla pack -> Step 1 -> Step 2
```

### Singularity or Apptainer (Linux / HPC)

```bash
# Run this from the FELIX repository root.
apptainer pull FELIX_latest.sif docker://lhu1/felix:latest
# Or: singularity pull FELIX_latest.sif docker://lhu1/felix:latest

cd tutorial
FELIX_BACKEND=apptainer SIF=../FELIX_latest.sif bash run_example.sh
# Or: FELIX_BACKEND=singularity SIF=../FELIX_latest.sif bash run_example.sh
```

### Pixi source install (Linux or macOS)

```bash
pixi run test                   # from the FELIX repository root
```

The runner uses the Docker image by default, the source install when FELIX
commands are already on `PATH` (for example, inside `pixi shell`), and the
`.sif` image when you select the Singularity/Apptainer backend. It generates
the simulation inside the selected environment, so container users do not
need host Python or NumPy.

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
- `binaries/` — the prebuilt FELIX R package, one tarball per platform, so a
  source install does not have to compile it.
- `tools/` — the pixi source-install scripts (`install_felix.sh`,
  `check_install.sh`).
- `extdata/` — the step wrapper scripts.
- `docker/` — Dockerfile and build/publish instructions.
- `tutorial/` — a runnable end-to-end example on simulated data.
- `pixi.toml` / `pixi.lock` — the pinned source-install environment.

## Citation

FELIX — Hu, L. et al. (in submission).
