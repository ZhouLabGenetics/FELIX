# SHAPEIT5 b38 4 cM Chunks

This directory vendors GRCh38 4 cM autosome chunk files for driving
`flare_subset_to_tractor_hybrid` region conversion.

Files are in `b38_4cM/chunks_chr1.txt` through `chunks_chr22.txt`. They use
chromosome names without the `chr` prefix, matching the original SHAPEIT5 b38
resource naming convention.

Each line has the SHAPEIT5/GLIMPSE chunk layout:

```text
chunk_index  chrom  buffered_region  core_region  ...
```

For tractor_hybrid conversion, use `core_region` (column 4). The buffered
region in column 3 is intended for phasing/ligation and overlaps neighboring
chunks; using it for conversion would duplicate variants at chunk boundaries.

Source and attribution:

- These files were copied from
  `Atkinson-Lab/Tractor/resources/genomic_chunks/chunks_b38_4cM`.
- The Tractor resource README states that its b38 4 cM chunks were sourced
  directly from SHAPEIT5 GitHub resources at
  `odelaneau/shapeit5/resources/chunks/b38/4cM`.
- Atkinson-Lab/Tractor is MIT licensed.
