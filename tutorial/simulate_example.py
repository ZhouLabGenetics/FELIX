"""
Simulate a tiny admixed-cohort dataset for the FELIX tutorial.

Produces a self-contained, fully simulated, publicly shareable example with:
  - 200 individuals from a 3-way admixed population (AFR / EUR / AMR)
  - 60 variants on a small chr22 region
  - one "causal" variant with an AFR-haplotype-specific effect
  - a quantitative phenotype + covariates
  - phased genotype VCF (input to the FELIXla converter)
  - FLARE-format local-ancestry VCF (per-haplotype AN1/AN2 calls)
  - PLINK .bed/.bim/.fam (input to SAIGE Step 1 for the GRM)
  - sample list

Nothing here is real human data; all draws are pseudo-random.

Run:
    python simulate_example.py            # writes everything into ./example_inputs/

No external dependencies beyond numpy.
"""

from __future__ import annotations
import os
import struct
import numpy as np

# --------------------------------------------------------- configuration
SEED        = 7
N           = 200          # samples
M           = 60           # variants
K           = 3            # ancestries (0=AFR, 1=EUR, 2=AMR)
ANC_LABELS  = ["AFR", "EUR", "AMR"]
CHR         = "chr22"
POS_START   = 16000000     # arbitrary chr22 position
POS_STEP    = 5000
CAUSAL_IDX  = 30           # index of the causal variant
OUT_DIR     = os.path.join(os.path.dirname(__file__), "example_inputs")

# per-individual ancestry-proportion archetypes (sums to 1)
ARCHETYPES = np.array([
    [0.80, 0.15, 0.05],   # AFR-dominant
    [0.10, 0.85, 0.05],   # EUR-dominant
    [0.05, 0.20, 0.75],   # AMR-dominant
    [0.45, 0.40, 0.15],   # AFR/EUR admixed
    [0.30, 0.35, 0.35],   # 3-way admixed
])

# per-ancestry allele frequencies for the *causal* variant
# (mimics a variant common in AFR, rare in EUR/AMR)
CAUSAL_AF_BY_ANC = [0.30, 0.04, 0.06]

# AFR-haplotype-specific effect on Y for the causal variant
CAUSAL_BETA_AFR = 1.4

# noise level for Y
NOISE_SD = 1.0


# --------------------------------------------------------- helpers
def make_dir():
    os.makedirs(OUT_DIR, exist_ok=True)


def sample_ancestry_proportions(rng: np.random.Generator) -> np.ndarray:
    """For each individual, pick an archetype and jitter."""
    arch_idx = rng.integers(0, len(ARCHETYPES), size=N)
    base = ARCHETYPES[arch_idx]
    jitter = rng.dirichlet([8.0, 8.0, 8.0], size=N)
    p = 0.7 * base + 0.3 * jitter
    p = p / p.sum(axis=1, keepdims=True)
    return p


def paint_haplotypes(rng: np.random.Generator,
                     ancprop: np.ndarray) -> np.ndarray:
    """
    Paint each variant of each haplotype with an ancestry label.

    Returns ancestries with shape (N, 2, M), entries in {0,1,2}.

    Uses a simple block-tract model: for each haplotype we sample 1-3
    breakpoints and assign ancestries to the resulting tracts in
    proportion to that individual's global ancestry mix.
    """
    anc = np.zeros((N, 2, M), dtype=np.int8)
    for i in range(N):
        for h in range(2):
            n_tracts = rng.integers(1, 4)  # 1, 2, or 3 tracts
            if n_tracts == 1:
                breaks = []
            else:
                breaks = sorted(rng.choice(np.arange(5, M - 5),
                                           size=n_tracts - 1, replace=False).tolist())
            segs = [0] + breaks + [M]
            tract_ancs = rng.choice(K, size=n_tracts, p=ancprop[i])
            for t, (a, b) in enumerate(zip(segs[:-1], segs[1:])):
                anc[i, h, a:b] = tract_ancs[t]
    return anc


def simulate_haplotypes(rng: np.random.Generator,
                        anc: np.ndarray) -> np.ndarray:
    """
    Draw 0/1 haplotype alleles given the ancestry label at each variant.

    Returns hap with shape (N, 2, M).

    For non-causal variants: per-ancestry AF drawn from Beta(0.4, 4) to give
    realistic MAF spread.  For the causal variant: use CAUSAL_AF_BY_ANC.
    """
    # per-ancestry AF for each variant
    af = rng.beta(0.4, 4.0, size=(K, M))
    af = np.clip(af, 0.01, 0.5)
    af[:, CAUSAL_IDX] = CAUSAL_AF_BY_ANC

    hap = np.zeros((N, 2, M), dtype=np.int8)
    for k in range(K):
        mask = (anc == k)                                # (N, 2, M)
        u = rng.random(size=mask.shape)
        af_full = np.broadcast_to(af[k][np.newaxis, np.newaxis, :],
                                  mask.shape)
        hap[mask] = (u[mask] < af_full[mask]).astype(np.int8)
    return hap


def simulate_phenotype(rng: np.random.Generator,
                       anc: np.ndarray,
                       hap: np.ndarray,
                       ancprop: np.ndarray) -> dict:
    """
    Simulate covariates and a quantitative phenotype.

    Y = 0.04 * age  +  0.6 * sex  +  1.2 * PC1  -  0.5 * PC2
        + CAUSAL_BETA_AFR * G_AFR_at_causal   + N(0, NOISE_SD^2)

    where G_AFR_at_causal is the number of *AFR-painted* alleles carrying
    the alt at the causal variant.  Only AFR haplotypes contribute to the
    genetic component, so the signal is ancestry-specific.
    """
    age = rng.integers(35, 80, size=N).astype(float)
    sex = rng.integers(0, 2, size=N).astype(float)

    # PCs informed by ancestry proportion (so they capture global ancestry)
    pcs = ancprop @ rng.normal(0, 1, size=(K, 3)) + rng.normal(0, 0.1, size=(N, 3))

    afr_mask  = (anc[:, :, CAUSAL_IDX] == 0)             # (N, 2)
    g_afr_caus = (hap[:, :, CAUSAL_IDX] * afr_mask).sum(axis=1).astype(float)

    y = (0.04 * (age - 55.0)
         + 0.6  * sex
         + 1.2  * pcs[:, 0]
         - 0.5  * pcs[:, 1]
         + CAUSAL_BETA_AFR * g_afr_caus
         + rng.normal(0, NOISE_SD, size=N))

    return dict(age=age, sex=sex, pcs=pcs, y=y, g_afr_caus=g_afr_caus)


# --------------------------------------------------------- writers
def write_samples(samples):
    path = os.path.join(OUT_DIR, "samples.txt")
    with open(path, "w") as f:
        for s in samples:
            f.write(s + "\n")
    print(f"wrote {path}  ({N} samples)")


def write_phenotype(samples, pheno):
    path = os.path.join(OUT_DIR, "phenotype.tsv")
    cols = ["IID", "age", "sex", "PC1", "PC2", "PC3", "Y"]
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for i, s in enumerate(samples):
            row = [s,
                   f"{pheno['age'][i]:.0f}",
                   f"{int(pheno['sex'][i])}",
                   f"{pheno['pcs'][i,0]:.4f}",
                   f"{pheno['pcs'][i,1]:.4f}",
                   f"{pheno['pcs'][i,2]:.4f}",
                   f"{pheno['y'][i]:.4f}"]
            f.write("\t".join(row) + "\n")
    print(f"wrote {path}")


def write_phased_vcf(samples, hap, positions, refs, alts):
    path = os.path.join(OUT_DIR, "genotypes.phased.vcf")
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##contig=<ID={CHR}>\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotype\">\n")
        f.write("\t".join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"] + samples) + "\n")
        for v in range(M):
            row = [CHR, str(positions[v]), f"v{v+1}",
                   refs[v], alts[v], ".", "PASS", ".", "GT"]
            for i in range(N):
                row.append(f"{hap[i,0,v]}|{hap[i,1,v]}")
            f.write("\t".join(row) + "\n")
    print(f"wrote {path}  ({M} variants, {N} samples)")


def write_flare_vcf(samples, hap, anc, positions, refs, alts):
    path = os.path.join(OUT_DIR, "localanc.flare.vcf")
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##contig=<ID={CHR}>\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotype\">\n")
        f.write("##FORMAT=<ID=AN1,Number=1,Type=Integer,Description=\"Ancestry label of first haplotype\">\n")
        f.write("##FORMAT=<ID=AN2,Number=1,Type=Integer,Description=\"Ancestry label of second haplotype\">\n")
        for k, name in enumerate(ANC_LABELS):
            f.write(f"##ANCESTRY=<ID={k},Name={name}>\n")
        f.write("\t".join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"] + samples) + "\n")
        for v in range(M):
            row = [CHR, str(positions[v]), f"v{v+1}",
                   refs[v], alts[v], ".", "PASS", ".", "GT:AN1:AN2"]
            for i in range(N):
                row.append(f"{hap[i,0,v]}|{hap[i,1,v]}:{anc[i,0,v]}:{anc[i,1,v]}")
            f.write("\t".join(row) + "\n")
    print(f"wrote {path}")


def write_rfmix_msp(samples, anc, positions):
    """
    Write the same local-ancestry calls in RFMix `.msp.tsv` form.

    RFMix reports ancestry per *interval* rather than per variant, so the
    per-variant labels are collapsed into runs where every haplotype keeps the
    same label. Two columns per sample (`ID.0`, `ID.1`), in the same haplotype
    order as the `|` separator in the phased VCF.

    This is an alternative to localanc.flare.vcf: FELIXla reads either one and
    produces the same packed output.
    """
    path = os.path.join(OUT_DIR, "localanc.rfmix.msp.tsv")

    rows = []
    for v in range(M):
        labels = []
        for i in range(N):
            labels += [int(anc[i, 0, v]), int(anc[i, 1, v])]
        rows.append((positions[v], labels))

    blocks = []
    for pos, labels in rows:
        if blocks and blocks[-1]["labels"] == labels:
            blocks[-1]["epos"] = pos
            blocks[-1]["n"] += 1
        else:
            blocks.append({"spos": pos, "epos": pos, "labels": labels, "n": 1})

    with open(path, "w") as f:
        f.write("#Subpopulation order/codes: "
                + "\t".join(f"{name}={k}" for k, name in enumerate(ANC_LABELS)) + "\n")
        hap_cols = []
        for s in samples:
            hap_cols += [f"{s}.0", f"{s}.1"]
        f.write("\t".join(["#chm", "spos", "epos", "sgpos", "egpos", "n snps"] + hap_cols) + "\n")
        for b in blocks:
            f.write("\t".join(
                [CHR, str(b["spos"]), str(b["epos"]),
                 f"{b['spos'] / 1e6:.6f}", f"{b['epos'] / 1e6:.6f}", str(b["n"])]
                + [str(x) for x in b["labels"]]) + "\n")
    print(f"wrote {path}  ({len(blocks)} ancestry intervals)")


def write_plink(samples, hap, pheno, positions, refs, alts):
    """
    PLINK 1 binary format writer (SNP-major, no external deps).
    Encoding per 2 bits:  00 = hom A1 (ref/ref),
                          01 = missing,
                          10 = het,
                          11 = hom A2 (alt/alt).
    """
    # collapse the two haplotypes to genotype counts {0,1,2} of ALT
    geno = hap.sum(axis=1).astype(np.int8)        # (N, M)

    # .fam
    fam = os.path.join(OUT_DIR, "example_grm.fam")
    with open(fam, "w") as f:
        for i, s in enumerate(samples):
            sex_code = int(pheno["sex"][i]) + 1   # 1=male, 2=female
            f.write(f"{s}\t{s}\t0\t0\t{sex_code}\t-9\n")

    # .bim
    # SAIGE step 1 requires a numeric chromosome in the PLINK bim (the GRM
    # markers are independent of the VCF), so strip the "chr" prefix here.
    chrom_num = CHR.replace("chr", "")
    bim = os.path.join(OUT_DIR, "example_grm.bim")
    with open(bim, "w") as f:
        for v in range(M):
            # PLINK convention: A1 = minor (alt), A2 = major (ref) in .bim is
            # NOT enforced; here we just write REF as A2 and ALT as A1.
            f.write(f"{chrom_num}\tv{v+1}\t0\t{positions[v]}\t{alts[v]}\t{refs[v]}\n")

    # .bed
    bed = os.path.join(OUT_DIR, "example_grm.bed")
    bytes_per_var = (N + 3) // 4
    with open(bed, "wb") as f:
        f.write(struct.pack("BBB", 0x6c, 0x1b, 0x01))    # magic + SNP-major
        for v in range(M):
            buf = bytearray(bytes_per_var)
            for i in range(N):
                g = geno[i, v]
                # encode 0->00 (hom A2? PLINK convention: A1 is allele1)
                # We write so that:
                #   0 alt copies -> 11 (hom A2 = REF/REF)
                #   1 alt copies -> 10 (het)
                #   2 alt copies -> 00 (hom A1 = ALT/ALT)
                if   g == 0: bits = 0b11
                elif g == 1: bits = 0b10
                elif g == 2: bits = 0b00
                else:        bits = 0b01    # missing
                byte_i = i // 4
                bit_off = (i % 4) * 2
                buf[byte_i] |= (bits & 0b11) << bit_off
            f.write(bytes(buf))

    print(f"wrote {fam}, {bim}, {bed}")


# --------------------------------------------------------- main
def main():
    rng = np.random.default_rng(SEED)
    make_dir()

    samples   = [f"S{i:04d}" for i in range(N)]
    positions = [POS_START + v * POS_STEP for v in range(M)]
    refs      = list(rng.choice(list("ACGT"), size=M))
    alts      = []
    for r in refs:
        other = [x for x in "ACGT" if x != r]
        alts.append(rng.choice(other))

    ancprop = sample_ancestry_proportions(rng)
    anc     = paint_haplotypes(rng, ancprop)
    hap     = simulate_haplotypes(rng, anc)
    pheno   = simulate_phenotype(rng, anc, hap, ancprop)

    write_samples(samples)
    write_phenotype(samples, pheno)
    write_phased_vcf(samples, hap, positions, refs, alts)
    write_flare_vcf(samples, hap, anc, positions, refs, alts)
    write_rfmix_msp(samples, anc, positions)
    write_plink(samples, hap, pheno, positions, refs, alts)

    print("\nSanity check:")
    print(f"  causal variant idx        = {CAUSAL_IDX}  (chr22:{positions[CAUSAL_IDX]})")
    print(f"  AFR-only carriers of alt  = {int((pheno['g_afr_caus'] > 0).sum())}/{N}")
    print(f"  simulated AFR-effect      = {CAUSAL_BETA_AFR}")
    print("\nDone.  All outputs in:", OUT_DIR)


if __name__ == "__main__":
    main()
