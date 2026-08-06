# Legacy conda recipes — not the FELIX install path

The files in this directory are inherited from upstream SAIGE. They build a
conda environment for **SAIGE**, not for FELIX, and
`install_SAIGE_conda_linux.sh` clones <https://github.com/saigegit/SAIGE>.
They are kept only for reference.

To install FELIX, use one of the two supported paths:

- **Docker** (recommended): `docker pull --platform linux/amd64 lhu1/felix:latest`
- **Source, via pixi** (Linux and macOS): `pixi run setup` from the repository root

Both are documented in [`../README.md`](../README.md#install) and, step by
step, in [the tutorial](../tutorial/README.md#2-installation). The pixi
environment is defined by `../pixi.toml` and pinned by `../pixi.lock`; it
supersedes `environment-RSAIGE.yml`.
