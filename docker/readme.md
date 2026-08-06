# Building and publishing the FELIX image

The image bundles the FELIX R package (FELIXassoc), the `felixla` command
(packing, RFMix/TRACTOR conversion, export, query), and the step wrappers
(`step1_fitNULLGLMM.R`, `step2_SPAtests.R`, `step3_LDmat.R`,
`createSparseGRM.R`).

> This page is for **maintainers building the image**. Users should either
> `docker pull lhu1/felix:latest` or do the pixi source install — see
> [`../README.md`](../README.md#install) and
> [the tutorial](../tutorial/README.md#2-installation).
>
> The image and the source install share the same `pixi.toml`/`pixi.lock`, so a
> change to the dependency list affects both. After editing `pixi.toml`, run
> `pixi run check` locally **and** rebuild the image before pushing.
>
> `felixla` is linked against the **pixi environment's** htslib inside the
> image, not the distro `libhts-dev`: the base image is Ubuntu 20.04, whose
> htslib is 1.10, and FELIXla needs 1.11+ for `hts_idx_nseq()`. That is why the
> build step runs under `pixi run` with explicit `HTSLIB_*` flags and an rpath
> to `/app/.pixi/envs/default/lib`.

## Recommended: linux/amd64

`linux/amd64` runs natively on Linux and, under Docker Desktop, on both Intel
and Apple Silicon macOS (Apple Silicon runs it via built-in emulation). Always
pass `--platform linux/amd64` when running on a Mac.

```bash
docker login                      # as lhu1

docker buildx build \
  --platform linux/amd64 \
  -t lhu1/felix:latest \
  -f docker/Dockerfile \
  --push .
```

Build a local test image before pushing:

```bash
docker build --platform linux/amd64 -t lhu1/felix:latest -f docker/Dockerfile .
```

## Native multi-arch (amd64 + arm64)

To run natively (no emulation) on both Linux/amd64 and Apple Silicon, publish a
multi-arch manifest. Both arches have been built from this Dockerfile and each
runs the tutorial end to end; within each arch, packing from a FLARE VCF and
from an RFMix `.msp.tsv` gives bit-identical Step 2 output. The arches differ
from each other only in the last digits (see the reproducibility note in the
tutorial).

```bash
docker buildx build \
  --platform linux/amd64,linux/arm64 \
  -t lhu1/felix:latest \
  -f docker/Dockerfile \
  --push .
```

Docker automatically serves the right architecture on `docker pull`, so the
tutorial's `--platform linux/amd64` flag becomes optional once this is pushed
(it still works, and still forces the amd64 variant under emulation).

`buildx --push` rebuilds both arches inside the buildx builder, which does not
share the default builder's cache. To publish images you have **already** built
and tested locally, tag and push them per-arch and join them with a manifest
instead:

```bash
docker tag felix:amd64 lhu1/felix:v0.0-amd64 && docker push lhu1/felix:v0.0-amd64
docker tag felix:arm64 lhu1/felix:v0.0-arm64 && docker push lhu1/felix:v0.0-arm64
docker manifest create lhu1/felix:v0.0 \
    --amend lhu1/felix:v0.0-amd64 --amend lhu1/felix:v0.0-arm64
docker manifest push lhu1/felix:v0.0
```

## Singularity / Apptainer

On HPC, convert the pushed image to the `.sif` referenced in the tutorial:

```bash
singularity pull FELIX_latest.sif docker://lhu1/felix:latest
# or
apptainer pull FELIX_latest.sif docker://lhu1/felix:latest
```
