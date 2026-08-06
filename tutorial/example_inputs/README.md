# Example inputs

Simulated, shareable inputs produced by `../simulate_example.py` (seeded, so they
reproduce exactly). No real human data. 200 individuals, 3 ancestries, 60 variants
on `chr22`; variant **v31** is the simulated AFR-haplotype-specific causal variant.

| file | description |
|---|---|
| `samples.txt` | 200 sample IDs, one per line |
| `phenotype.tsv` | quantitative trait `Y` + covariates `age, sex, PC1, PC2, PC3` (`IID` matches the VCFs) |
| `genotypes.phased.vcf` | phased genotype VCF (`0|1`), input to the FELIXla converter |
| `localanc.flare.vcf` | FLARE local-ancestry VCF with per-haplotype `AN1`/`AN2` ancestry labels |
| `localanc.rfmix.msp.tsv` | the same ancestry calls as an RFMix `.msp.tsv` (interval form); an alternative input to the FLARE VCF |
| `example_grm.{bed,bim,fam}` | PLINK fileset used by Step 1 for the GRM and variance-ratio markers |

Step 0 of the tutorial bgzips and tabix-indexes the VCFs before packing.
`felixla` reads either `localanc.flare.vcf` (`--flare-vcf`) or
`localanc.rfmix.msp.tsv` (`--rfmix-msp`) and produces the same packed output.
