# Prebuilt `felixla`

Ready-to-run builds of the FELIXla command, so a source install does not have
to compile it.

| directory | target |
|---|---|
| `linux-64/` | Linux, x86-64 — servers and essentially all HPC clusters |
| `osx-arm64/` | macOS, Apple Silicon (M1/M2/M3/M4) |

`pixi run setup` (via `../install_felixla.sh`) picks the directory matching
`uname -s`/`uname -m`, copies the binary into `.pixi/envs/default/bin`, and runs
it. **If there is no directory for your platform — Intel macOS (`osx-64`), ARM
Linux (`linux-aarch64`) — or if the binary fails to run, setup falls back to
compiling from `../src/` automatically.** No action is needed either way.

## Why these are portable

They are built inside a pixi environment and linked with a runtime search path
relative to the executable itself:

- Linux: `RUNPATH = $ORIGIN/../lib`
- macOS: `LC_RPATH = @loader_path/../lib`

Installed into `<env>/bin`, that resolves to `<env>/lib`, where the pixi
environment's `libhts` lives. So the binary works regardless of where the
repository is cloned, and never depends on a system-wide htslib.

They link `libhts.so.3` / `libhts.3.dylib` from the pixi environment
(htslib 1.24; FELIXla needs 1.11+). They are also compiled **without** the Makefile's default `-march=native`, so
they are not tied to the CPU that built them. The Linux build uses the
conda-forge sysroot (glibc 2.17) and runs on any distribution from CentOS 7
onwards.

## Rebuilding

See [`../README.md`](../README.md#regenerating-the-shipped-prebuilt-binary-maintainers).
