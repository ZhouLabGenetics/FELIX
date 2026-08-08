#!/usr/bin/env bash
# End-to-end FELIX example on the simulated data: FELIXla pack -> Step 1 -> Step 2.
# Produces ./output/step2_results.tsv. Run from this tutorial/ directory.
#
#   bash run_example.sh                         # uses the Docker image
#   FELIX_BACKEND=singularity bash run_example.sh # uses FELIX_latest.sif
#   pixi run test                               # uses a source install, from the repo root
#
# The backend is auto-detected: if the FELIX commands are already on PATH (as
# they are inside `pixi shell` / `pixi run`), they are used directly; otherwise
# every command is run inside the Docker image. Force one with
# FELIX_BACKEND=native, docker, singularity, or apptainer. For the latter two,
# SIF points to the image file (default: FELIX_latest.sif) and CONTAINER_BIN
# optionally selects the runtime explicitly.
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
SIF="${SIF:-${SINGULARITY_IMAGE:-FELIX_latest.sif}}"
CONTAINER_BIN="${CONTAINER_BIN:-}"

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
    singularity|apptainer)
        if [[ -z "$CONTAINER_BIN" ]]; then
            if [[ "$FELIX_BACKEND" == "apptainer" ]]; then
                CONTAINER_BIN=apptainer
            else
                CONTAINER_BIN=singularity
            fi
        fi
        if ! command -v "$CONTAINER_BIN" >/dev/null 2>&1; then
            echo "ERROR: '$CONTAINER_BIN' is not on PATH. Set CONTAINER_BIN to your Singularity/Apptainer command." >&2
            exit 1
        fi
        if [[ ! -f "$SIF" ]]; then
            echo "ERROR: SIF image not found: $SIF" >&2
            echo "       Create it with: $CONTAINER_BIN pull $SIF docker://lhu1/felix:latest" >&2
            exit 1
        fi
        echo "== FELIX example ($CONTAINER_BIN image $SIF, local ancestry from: $LA_SOURCE) =="
        RUN=("$CONTAINER_BIN" exec --bind "${PWD}:/work" --pwd /work "$SIF")
        PYTHON="${PYTHON:-python3}"
        ;;
    *)
        echo "ERROR: FELIX_BACKEND must be 'native', 'docker', 'singularity', or 'apptainer' (got '$FELIX_BACKEND')" >&2
        exit 1
        ;;
esac

run() {
    if [[ ${#RUN[@]} -eq 0 ]]; then
        "$@"
    else
        "${RUN[@]}" "$@"
    fi
}

mkdir -p packed output
# Run the simulator through the selected backend too. Docker and SIF users
# therefore need no host Python packages; the image provides Python + NumPy.
[ -f example_inputs/genotypes.phased.vcf ] || run "$PYTHON" simulate_example.py

# Step 0 — bgzip + index the phased genotype VCF, then pack it together with
# the local-ancestry calls into the FELIXla format.
case "$LA_SOURCE" in
    flare)
        run bash -c '
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
        run bash -c '
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
run step1_fitNULLGLMM.R \
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
run step2_SPAtests.R \
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

if [[ ! -s output/step2_results.tsv ]]; then
    echo "ERROR: Step 2 finished without producing output/step2_results.tsv" >&2
    exit 1
fi

echo "Done. Results: output/step2_results.tsv"
