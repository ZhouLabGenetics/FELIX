#!/usr/bin/env bash
# End-to-end FELIX example on the simulated data: FELIXla pack -> Step 1 -> Step 2.
# Produces ./output/step2_results.tsv. Run from this tutorial/ directory.
#
#   bash run_example.sh          # uses the Docker image (recommended)
#   pixi run test                # uses a source install, from the repo root
#
# The backend is auto-detected: if the FELIX commands are already on PATH (as
# they are inside `pixi shell` / `pixi run`), they are used directly; otherwise
# every command is run inside the Docker image. Force one with
# FELIX_BACKEND=native or FELIX_BACKEND=docker.
#
# Set LA_SOURCE=rfmix to pack from the RFMix .msp.tsv instead of the FLARE VCF.
# Both produce the same association results.
#
# The image is multi-arch, so Docker picks the native build. Set
# DOCKER_PLATFORM=linux/amd64 to pin the architecture (that is the build the
# committed example_outputs/ reference was generated with).
set -euo pipefail

IMAGE="${IMAGE:-lhu1/felix:latest}"
LA_SOURCE="${LA_SOURCE:-flare}"
DOCKER_PLATFORM="${DOCKER_PLATFORM:-}"

if [[ -z "${FELIX_BACKEND:-}" ]]; then
    if command -v felixla >/dev/null 2>&1 && command -v step1_fitNULLGLMM.R >/dev/null 2>&1; then
        FELIX_BACKEND=native
    else
        FELIX_BACKEND=docker
    fi
fi

case "$FELIX_BACKEND" in
    native)
        echo "== FELIX example (native / pixi install, local ancestry from: $LA_SOURCE) =="
        # Empty prefix: commands run directly. Expanded as ${RUN[@]+...} below
        # so it also works under bash 3.2 (macOS system bash) with `set -u`.
        RUN=()
        PYTHON="${PYTHON:-python3}"
        ;;
    docker)
        echo "== FELIX example (Docker image $IMAGE, local ancestry from: $LA_SOURCE) =="
        RUN=(docker run --rm)
        [ -n "$DOCKER_PLATFORM" ] && RUN+=(--platform "$DOCKER_PLATFORM")
        RUN+=(-v "${PWD}":/work -w /work "${IMAGE}")
        PYTHON="${PYTHON:-python3}"
        ;;
    *)
        echo "ERROR: FELIX_BACKEND must be 'native' or 'docker' (got '$FELIX_BACKEND')" >&2
        exit 1
        ;;
esac

mkdir -p packed output
[ -f example_inputs/genotypes.phased.vcf ] || "$PYTHON" simulate_example.py

# Step 0 — bgzip + index the phased genotype VCF, then pack it together with
# the local-ancestry calls into the FELIXla format.
case "$LA_SOURCE" in
    flare)
        ${RUN[@]+"${RUN[@]}"} bash -c '
          bgzip -c example_inputs/genotypes.phased.vcf > example_inputs/genotypes.phased.vcf.gz
          bgzip -c example_inputs/localanc.flare.vcf  > example_inputs/localanc.flare.vcf.gz
          tabix -f -p vcf example_inputs/genotypes.phased.vcf.gz
          tabix -f -p vcf example_inputs/localanc.flare.vcf.gz
          felixla \
            --phase-vcf example_inputs/genotypes.phased.vcf.gz \
            --flare-vcf example_inputs/localanc.flare.vcf.gz \
            --n-ancestries 3 \
            --mac-threshold 20 \
            --make-felixla \
            --out packed/example
        '
        ;;
    rfmix)
        ${RUN[@]+"${RUN[@]}"} bash -c '
          bgzip -c example_inputs/genotypes.phased.vcf > example_inputs/genotypes.phased.vcf.gz
          tabix -f -p vcf example_inputs/genotypes.phased.vcf.gz
          felixla \
            --phase-vcf example_inputs/genotypes.phased.vcf.gz \
            --rfmix-msp example_inputs/localanc.rfmix.msp.tsv \
            --n-ancestries 3 \
            --mac-threshold 20 \
            --make-felixla \
            --out packed/example
        '
        ;;
    *)
        echo "ERROR: LA_SOURCE must be 'flare' or 'rfmix' (got '$LA_SOURCE')" >&2
        exit 1
        ;;
esac

# Step 1 — fit the null GLMM (does not use local ancestry).
${RUN[@]+"${RUN[@]}"} step1_fitNULLGLMM.R \
  --plinkFile=example_inputs/example_grm \
  --phenoFile=example_inputs/phenotype.tsv \
  --phenoCol=Y \
  --traitType=quantitative \
  --invNormalize=FALSE \
  --covarColList=age,sex,PC1,PC2,PC3 \
  --qCovarColList=sex \
  --sampleIDColinphenoFile=IID \
  --outputPrefix=output/step1_example \
  --nThreads=1 \
  --LOCO=FALSE \
  --IsOverwriteVarianceRatioFile=TRUE

# Step 2 — ancestry-aware association test on the FELIXla input.
${RUN[@]+"${RUN[@]}"} step2_SPAtests.R \
  --FELIXlaPrefix=packed/example \
  --chrom=chr22 \
  --is_admixed=TRUE \
  --number_of_ancestry=3 \
  --pvalcutoff_of_haplotype=0.000005 \
  --GMMATmodelFile=output/step1_example.rda \
  --varianceRatioFile=output/step1_example.varianceRatio.txt \
  --SAIGEOutputFile=output/step2_results.tsv \
  --LOCO=FALSE \
  --minMAF=0 \
  --minMAC=1

echo "Done. Results: output/step2_results.tsv"
