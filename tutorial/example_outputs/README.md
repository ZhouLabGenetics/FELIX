# Example outputs

`step2_results.tsv` here is the reference output, produced with the Docker
image (`lhu1/felix:latest`, linux/amd64) by packing from the FLARE VCF and
fitting Step 1 with `--nThreads=1`. Regenerate it from the `tutorial/`
directory with:

```bash
bash run_example.sh
```

Running the tutorial's commands verbatim with the Docker image reproduces it
bit-for-bit, from either the FLARE VCF or the RFMix `.msp.tsv`. A pixi source
install (`linux-64` or `osx-arm64`) may land on Step 1's other solution and
report `P_cct_admixed` = 6.34e-12 instead of 5.89e-12 — same ranking, same
conclusions. See the reproducibility note in
[section 2B of the tutorial](../README.md#2b-source-install-with-pixi-linux-and-macos).

This writes, into `../output/`:

- `step1_example.rda`, `step1_example.varianceRatio.txt` — the Step 1 null model.
- `step2_results.tsv` — the Step 2 association results (one row per variant;
  per-ancestry `BETA_anc{j}` / `p.value_anc{j}`, plus `P_het_admixed`,
  `P_hom_admixed`, `P_cct_admixed`).

For the causal variant **v31**, expect a strong AFR-specific signal
(`BETA_anc1` ≈ 1.4, small `p.value_anc1`), null EUR/AMR arms, and a
`P_het_admixed` far more significant than the collapsed `P_hom_admixed`. See
section 9 of `../README.md` for the exact command to print that row.
