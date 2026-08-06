#!/usr/bin/env bash
# designed by Kai, implemented by codex
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN_DIR="${BIN_DIR:-$ROOT_DIR/bin}"
OUT_DIR="$(mktemp -d "${TMPDIR:-/tmp}/tractor-hybrid-tiny.XXXXXX")"
trap 'rm -rf "$OUT_DIR"' EXIT

python3 - "$BIN_DIR" <<'PY'
import pathlib
import sys

bin_dir = pathlib.Path(sys.argv[1])
unexpected = sorted(
    path.name for path in bin_dir.iterdir()
    if path.name != "felixla"
)
assert not unexpected, unexpected
assert (bin_dir / "felixla").is_file(), bin_dir
PY

"$BIN_DIR/felixla" \
  from-flare \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  "$ROOT_DIR/testdata/tiny.flare.vcf" \
  2 \
  1 \
  "$OUT_DIR/tiny" >/dev/null

"$BIN_DIR/felixla" \
  to-vcf \
  "$OUT_DIR/tiny" \
  "$OUT_DIR/tiny.roundtrip.vcf.gz" >/dev/null

"$BIN_DIR/felixla" \
  --felixla \
  "$OUT_DIR/tiny" \
  --query chr1:160 \
  --ref A \
  --alt T \
  >"$OUT_DIR/tiny.common.query.tsv"

"$BIN_DIR/felixla" \
  --felixla \
  "$OUT_DIR/tiny" \
  --query chr1:100 \
  --ref A \
  --alt G \
  --nonzero-only \
  >"$OUT_DIR/tiny.rare.query.tsv"

python3 - "$OUT_DIR/tiny.common.query.tsv" "$OUT_DIR/tiny.rare.query.tsv" <<'PY'
import pathlib
import sys

common = pathlib.Path(sys.argv[1]).read_text().strip().splitlines()
rare = pathlib.Path(sys.argv[2]).read_text().strip().splitlines()
header = "global_variant_index\tchr\tpos\tid\tref\talt\tsample\tDSALL\tDS1\tDS2"
assert common == [
    header,
    "4\tchr1\t160\tcommon_A_T\tA\tT\ts1\t2\t1\t1",
    "4\tchr1\t160\tcommon_A_T\tA\tT\ts2\t0\t0\t0",
], common
assert rare == [
    header,
    "2\tchr1\t100\tmulti_A_G\tA\tG\ts2\t1\t0\t1",
], rare
PY

"$BIN_DIR/felixla" --version >"$OUT_DIR/felixla.version.txt"
grep -q "FELIXla CLI v0" "$OUT_DIR/felixla.version.txt"

"$BIN_DIR/felixla" --recommend-mac-threshold --n-samples 2 >"$OUT_DIR/mac_threshold.tsv"
grep -q $'recommended_mac_threshold\t1' "$OUT_DIR/mac_threshold.tsv"

"$BIN_DIR/felixla" \
  --phase-vcf \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  --flare-vcf \
  "$ROOT_DIR/testdata/tiny.flare.vcf" \
  --n-ancestries 2 \
  --make-felixla \
  --out \
  "$OUT_DIR/tiny.felixla_cli" >/dev/null

"$BIN_DIR/felixla" \
  --felixla \
  "$OUT_DIR/tiny.felixla_cli" \
  --export vcf \
  --out \
  "$OUT_DIR/tiny.felixla_cli.roundtrip.vcf.gz" >/dev/null

"$BIN_DIR/felixla" \
  --compare-vcfs \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  "$OUT_DIR/tiny.felixla_cli.roundtrip.vcf.gz" \
  --split-multiallelic \
  --out "$OUT_DIR/felixla_cli.compare.txt"

grep -q "Differences:             0" "$OUT_DIR/felixla_cli.compare.txt"

cat >"$OUT_DIR/keep.samples" <<'EOF'
s2
EOF

cat >"$OUT_DIR/extract.sites.vcf" <<'EOF'
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ignored_sample
chr1	100	extract_id_is_ignored_1	A	C,G	.	PASS	.	GT	0|1
chr1	160	extract_id_is_ignored_2	A	T	.	PASS	.	GT	0|0
EOF

python3 - "$OUT_DIR/extract.sites.vcf" "$OUT_DIR/extract.sites.vcf.gz" <<'PY'
import gzip
import pathlib
import sys

gzip.open(sys.argv[2], "wb").write(pathlib.Path(sys.argv[1]).read_bytes())
PY

"$BIN_DIR/felixla" \
  --phase-vcf \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  --flare-vcf \
  "$ROOT_DIR/testdata/tiny.flare.vcf" \
  --n-ancestries 2 \
  --keep "$OUT_DIR/keep.samples" \
  --extract "$OUT_DIR/extract.sites.vcf.gz" \
  --make-felixla \
  --out \
  "$OUT_DIR/tiny.keep_extract" >/dev/null

"$BIN_DIR/felixla" \
  to-vcf \
  "$OUT_DIR/tiny.keep_extract" \
  "$OUT_DIR/tiny.keep_extract.roundtrip.vcf.gz" >/dev/null

python3 - "$OUT_DIR/tiny.keep_extract" <<'PY'
import gzip
import pathlib
import sys

prefix = pathlib.Path(sys.argv[1])
meta = pathlib.Path(str(prefix) + ".meta").read_text()
samples = pathlib.Path(str(prefix) + ".samples").read_text().strip().splitlines()
assert "n_samples\t1" in meta, meta
assert "keep_samples\t" in meta, meta
assert "extract_sites\t" in meta, meta
assert samples == ["s2"], samples

with gzip.open(str(prefix) + ".roundtrip.vcf.gz", "rt") as fh:
    lines = [line.rstrip() for line in fh]
header = [line for line in lines if line.startswith("#CHROM")][0].split("\t")
assert header[9:] == ["s2"], header
rows = [line.split("\t") for line in lines if line and not line.startswith("#")]
assert [(r[0], r[1], r[2], r[3], r[4], r[9]) for r in rows] == [
    ("chr1", "100", "multi_A_C", "A", "C", "0|0"),
    ("chr1", "100", "multi_A_G", "A", "G", "1|0"),
    ("chr1", "160", "common_A_T", "A", "T", "0|0"),
], rows
PY

cat >"$OUT_DIR/extract.bad_ref.vcf" <<'EOF'
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	multi	T	G	.	PASS	.
EOF

if "$BIN_DIR/felixla" \
  --phase-vcf \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  --flare-vcf \
  "$ROOT_DIR/testdata/tiny.flare.vcf" \
  --n-ancestries 2 \
  --extract "$OUT_DIR/extract.bad_ref.vcf" \
  --make-felixla \
  --out \
  "$OUT_DIR/tiny.bad_extract_ref" >/dev/null 2>"$OUT_DIR/tiny.bad_extract_ref.err"; then
  echo "expected bad REF in --extract fixture to fail" >&2
  exit 1
fi

grep -q "REF mismatch" "$OUT_DIR/tiny.bad_extract_ref.err"

"$BIN_DIR/felixla" \
  --felixla \
  "$OUT_DIR/tiny.felixla_cli" \
  --region chr1:100-160 \
  --make-felixla \
  --out \
  "$OUT_DIR/tiny.felixla_cli.chr1_100_160" >/dev/null

"$BIN_DIR/felixla" \
  --vcf \
  "$ROOT_DIR/testdata/tiny.flare.vcf" \
  --admixture \
  --progress-every 0 \
  >"$OUT_DIR/felixla_admixture.tsv" \
  2>"$OUT_DIR/felixla_admixture.err"

grep -q $'s1\t6\t0.3333333333\t0.6666666667\t2\t4' "$OUT_DIR/felixla_admixture.tsv"

cat >"$OUT_DIR/shapeit.chunks.txt" <<'EOF'
0	chr1	chr1:1-180	chr1:1-150	4.0	150	10	3
1	chr1	chr1:120-260	chr1:151-250	4.0	100	10	3
EOF

"$BIN_DIR/felixla" \
  --shapeit-args \
  --chunks "$OUT_DIR/shapeit.chunks.txt" \
  --phase-template "$OUT_DIR/source/{chrom}.phased.vcf.gz" \
  --flare-template "$OUT_DIR/source/{chrom}.flare.vcf.gz" \
  --n-ancestries 2 \
  --n-samples 2 \
  --out-prefix-template "$OUT_DIR/{chrom}.shapeit4cM.chunk{chunk0}" \
  >"$OUT_DIR/shapeit.args.tsv"

python3 - "$OUT_DIR/shapeit.args.tsv" "$OUT_DIR" <<'PY'
import pathlib
import sys

args_file = pathlib.Path(sys.argv[1])
out_dir = pathlib.Path(sys.argv[2])
rows = [line.rstrip().split("\t") for line in args_file.read_text().splitlines()]
assert [[pathlib.Path(row[0]), pathlib.Path(row[1]), row[2], row[3], pathlib.Path(row[4]), row[5]] for row in rows] == [
    [
        out_dir / "source" / "chr1.phased.vcf.gz",
        out_dir / "source" / "chr1.flare.vcf.gz",
        "2",
        "1",
        out_dir / "chr1.shapeit4cM.chunk0001",
        "chr1:1-150",
    ],
    [
        out_dir / "source" / "chr1.phased.vcf.gz",
        out_dir / "source" / "chr1.flare.vcf.gz",
        "2",
        "1",
        out_dir / "chr1.shapeit4cM.chunk0002",
        "chr1:151-250",
    ],
], rows
PY

"$BIN_DIR/felixla" \
  compare-vcfs \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  "$OUT_DIR/tiny.roundtrip.vcf.gz" \
  --split-multiallelic \
  >"$OUT_DIR/compare.txt"

grep -q "Differences:             0" "$OUT_DIR/compare.txt"

"$BIN_DIR/felixla" \
  extract \
  "$OUT_DIR/tiny" \
  chr1:100-160 \
  "$OUT_DIR/tiny.chr1_100_160" >/dev/null

"$BIN_DIR/felixla" \
  to-vcf \
  "$OUT_DIR/tiny.chr1_100_160" \
  "$OUT_DIR/tiny.chr1_100_160.roundtrip.vcf.gz" >/dev/null

"$BIN_DIR/felixla" \
  from-flare \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  "$ROOT_DIR/testdata/tiny.flare.vcf" \
  2 \
  1 \
  "$OUT_DIR/tiny.direct_region" \
  chr1:100-160 >/dev/null

grep -q $'selected_region\tchr1:100-160' "$OUT_DIR/tiny.direct_region.meta"

"$BIN_DIR/felixla" \
  to-vcf \
  "$OUT_DIR/tiny.direct_region" \
  "$OUT_DIR/tiny.direct_region.roundtrip.vcf.gz" >/dev/null

"$BIN_DIR/felixla" \
  compare-vcfs \
  "$OUT_DIR/tiny.chr1_100_160.roundtrip.vcf.gz" \
  "$OUT_DIR/tiny.direct_region.roundtrip.vcf.gz" \
  --split-multiallelic \
  >"$OUT_DIR/direct_region.compare.txt"

grep -q "Differences:             0" "$OUT_DIR/direct_region.compare.txt"

grep -v '^##contig=' "$ROOT_DIR/testdata/tiny.genotypes.vcf" >"$OUT_DIR/tiny.genotypes.no_contig.vcf"
grep -v '^##contig=' "$ROOT_DIR/testdata/tiny.flare.vcf" >"$OUT_DIR/tiny.flare.no_contig.vcf"

"$BIN_DIR/felixla" \
  from-flare \
  "$OUT_DIR/tiny.genotypes.no_contig.vcf" \
  "$OUT_DIR/tiny.flare.no_contig.vcf" \
  2 \
  1 \
  "$OUT_DIR/tiny.no_contig" >/dev/null

"$BIN_DIR/felixla" \
  to-vcf \
  "$OUT_DIR/tiny.no_contig" \
  "$OUT_DIR/tiny.no_contig.roundtrip.vcf.gz" >/dev/null

"$BIN_DIR/felixla" \
  compare-vcfs \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  "$OUT_DIR/tiny.no_contig.roundtrip.vcf.gz" \
  --split-multiallelic \
  >"$OUT_DIR/no_contig.compare.txt"

grep -q "Differences:             0" "$OUT_DIR/no_contig.compare.txt"

"$BIN_DIR/felixla" \
  from-flare \
  "$OUT_DIR/tiny.genotypes.no_contig.vcf" \
  "$OUT_DIR/tiny.flare.no_contig.vcf" \
  2 \
  1 \
  "$OUT_DIR/tiny.no_contig_region" \
  chr1:100-160 >/dev/null

grep -q $'selected_region\tchr1:100-160' "$OUT_DIR/tiny.no_contig_region.meta"

"$BIN_DIR/felixla" \
  to-vcf \
  "$OUT_DIR/tiny.no_contig_region" \
  "$OUT_DIR/tiny.no_contig_region.roundtrip.vcf.gz" >/dev/null

"$BIN_DIR/felixla" \
  compare-vcfs \
  "$OUT_DIR/tiny.direct_region.roundtrip.vcf.gz" \
  "$OUT_DIR/tiny.no_contig_region.roundtrip.vcf.gz" \
  --split-multiallelic \
  >"$OUT_DIR/no_contig_region.compare.txt"

grep -q "Differences:             0" "$OUT_DIR/no_contig_region.compare.txt"

"$BIN_DIR/felixla" \
  --phase-vcf \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  --rfmix-msp \
  "$ROOT_DIR/testdata/tiny.rfmix.msp.tsv" \
  --n-ancestries 2 \
  --make-felixla \
  --out \
  "$OUT_DIR/rfmix" >/dev/null

"$BIN_DIR/felixla" \
  to-vcf \
  "$OUT_DIR/rfmix" \
  "$OUT_DIR/rfmix.roundtrip.vcf.gz" >/dev/null

"$BIN_DIR/felixla" \
  compare-vcfs \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  "$OUT_DIR/rfmix.roundtrip.vcf.gz" \
  --split-multiallelic \
  >"$OUT_DIR/rfmix.compare.txt"

grep -q "Differences:             0" "$OUT_DIR/rfmix.compare.txt"

"$BIN_DIR/felixla" \
  --tractor-dosage-vcf \
  "$ROOT_DIR/testdata/tiny.tractor_dosage.vcf" \
  --n-ancestries 2 \
  --make-felixla \
  --out \
  "$OUT_DIR/dosage" >/dev/null

"$BIN_DIR/felixla" \
  to-vcf \
  "$OUT_DIR/dosage" \
  "$OUT_DIR/dosage.roundtrip.vcf.gz" >/dev/null

python3 - "$OUT_DIR/tiny" <<'PY'
import gzip
import pathlib
import struct
import sys

prefix = pathlib.Path(sys.argv[1])

def read_string(data, off):
    (n,) = struct.unpack_from("<I", data, off)
    off += 4
    value = data[off:off + n].decode()
    off += n
    return value, off

def read_common_idx(path):
    data = pathlib.Path(path).read_bytes()
    assert data[:8] == b"TRCMMKS1", data[:8]
    off = 8
    records = []
    while off < len(data):
        mks_offset = off
        common_index, global_index = struct.unpack_from("<QI", data, off)
        off += 12
        chrom, off = read_string(data, off)
        (pos,) = struct.unpack_from("<q", data, off)
        off += 8
        vid, off = read_string(data, off)
        ref, off = read_string(data, off)
        alt, off = read_string(data, off)
        alt_index, block_id = struct.unpack_from("<II", data, off)
        off += 8
        geno_offset, = struct.unpack_from("<Q", data, off)
        off += 8
        mac, = struct.unpack_from("<I", data, off)
        off += 4
        records.append({
            "mks_offset": mks_offset,
            "common_index": common_index,
            "global_index": global_index,
            "chr": chrom,
            "pos": pos,
            "id": vid,
            "ref": ref,
            "alt": alt,
            "alt_index": alt_index,
            "block_id": block_id,
            "geno_offset": geno_offset,
            "mac": mac,
        })
    return records

def read_rare_idx(path):
    data = pathlib.Path(path).read_bytes()
    assert data[:8] == b"TRRAMKS1", data[:8]
    off = 8
    records = []
    while off < len(data):
        mks_offset = off
        rare_index, global_index = struct.unpack_from("<QI", data, off)
        off += 12
        chrom, off = read_string(data, off)
        (pos,) = struct.unpack_from("<q", data, off)
        off += 8
        vid, off = read_string(data, off)
        ref, off = read_string(data, off)
        alt, off = read_string(data, off)
        (alt_index,) = struct.unpack_from("<I", data, off)
        off += 4
        carrier_offset, = struct.unpack_from("<Q", data, off)
        off += 8
        n_carriers, mac = struct.unpack_from("<II", data, off)
        off += 8
        records.append({
            "mks_offset": mks_offset,
            "rare_index": rare_index,
            "global_index": global_index,
            "chr": chrom,
            "pos": pos,
            "id": vid,
            "ref": ref,
            "alt": alt,
            "alt_index": alt_index,
            "carrier_offset": carrier_offset,
            "n_carriers": n_carriers,
            "mac": mac,
        })
    return records

def read_anc_idx(path):
    data = pathlib.Path(path).read_bytes()
    assert data[:8] == b"TRANMKS1", data[:8]
    off = 8
    records = []
    while off < len(data):
        mks_offset = off
        (block_id,) = struct.unpack_from("<I", data, off)
        off += 4
        chrom, off = read_string(data, off)
        start, end = struct.unpack_from("<qq", data, off)
        off += 16
        anc_offset, = struct.unpack_from("<Q", data, off)
        off += 8
        records.append({
            "mks_offset": mks_offset,
            "block_id": block_id,
            "chr": chrom,
            "start": start,
            "end": end,
            "anc_offset": anc_offset,
        })
    return records

def read_common_offset_idx(path):
    data = pathlib.Path(path).read_bytes()
    assert data[:8] == b"TRCMIDX2", data[:8]
    return [
        {
            "common_index": struct.unpack_from("<Q", data, off)[0],
            "global_index": struct.unpack_from("<I", data, off + 8)[0],
            "mks_offset": struct.unpack_from("<Q", data, off + 12)[0],
            "geno_offset": struct.unpack_from("<Q", data, off + 20)[0],
        }
        for off in range(8, len(data), 28)
    ]

def read_rare_offset_idx(path):
    data = pathlib.Path(path).read_bytes()
    assert data[:8] == b"TRRAIDX2", data[:8]
    return [
        {
            "rare_index": struct.unpack_from("<Q", data, off)[0],
            "global_index": struct.unpack_from("<I", data, off + 8)[0],
            "mks_offset": struct.unpack_from("<Q", data, off + 12)[0],
            "carrier_offset": struct.unpack_from("<Q", data, off + 20)[0],
            "n_carriers": struct.unpack_from("<I", data, off + 28)[0],
        }
        for off in range(8, len(data), 32)
    ]

def read_anc_offset_idx(path):
    data = pathlib.Path(path).read_bytes()
    assert data[:8] == b"TRANIDX2", data[:8]
    return [
        {
            "block_id": struct.unpack_from("<I", data, off)[0],
            "mks_offset": struct.unpack_from("<Q", data, off + 4)[0],
            "anc_offset": struct.unpack_from("<Q", data, off + 12)[0],
        }
        for off in range(8, len(data), 20)
    ]

common_idx = read_common_idx(str(prefix) + ".common.variant.mks")
rare_idx = read_rare_idx(str(prefix) + ".rare.variant.mks")
anc_idx = read_anc_idx(str(prefix) + ".ancblock.mks")
common_offsets = read_common_offset_idx(str(prefix) + ".common.variant.idx")
rare_offsets = read_rare_offset_idx(str(prefix) + ".rare.variant.idx")
anc_offsets = read_anc_offset_idx(str(prefix) + ".ancblock.idx")
meta = pathlib.Path(str(prefix) + ".meta").read_text().strip().splitlines()
samples = pathlib.Path(str(prefix) + ".samples").read_text().strip().splitlines()
roundtrip_path = pathlib.Path(str(prefix) + ".roundtrip.vcf.gz")
assert roundtrip_path.exists(), roundtrip_path
assert pathlib.Path(str(roundtrip_path) + ".tbi").exists(), str(roundtrip_path) + ".tbi"
with gzip.open(roundtrip_path, "rt") as fh:
    roundtrip_vcf = fh.read().strip().splitlines()

assert len(common_idx) == 1, common_idx
assert len(rare_idx) == 8, rare_idx
assert len(anc_idx) == 2, anc_idx
assert common_offsets == [
    {
        "common_index": r["common_index"],
        "global_index": r["global_index"],
        "mks_offset": r["mks_offset"],
        "geno_offset": r["geno_offset"],
    }
    for r in common_idx
], common_offsets
assert rare_offsets == [
    {
        "rare_index": r["rare_index"],
        "global_index": r["global_index"],
        "mks_offset": r["mks_offset"],
        "carrier_offset": r["carrier_offset"],
        "n_carriers": r["n_carriers"],
    }
    for r in rare_idx
], rare_offsets
assert anc_offsets == [
    {
        "block_id": r["block_id"],
        "mks_offset": r["mks_offset"],
        "anc_offset": r["anc_offset"],
    }
    for r in anc_idx
], anc_offsets
assert "n_samples\t2" in meta, meta
assert "n_words\t1" in meta, meta
assert samples == ["s1", "s2"], samples

assert {k: common_idx[0][k] for k in ("global_index", "chr", "pos", "id", "ref", "alt", "alt_index", "block_id")} == {
    "global_index": 4,
    "chr": "chr1",
    "pos": 160,
    "id": "common_A_T",
    "ref": "A",
    "alt": "T",
    "alt_index": 1,
    "block_id": 0,
}, common_idx[0]

rare_globals = [record["global_index"] for record in rare_idx]
assert rare_globals == [0, 1, 2, 3, 5, 6, 7, 8], rare_globals

assert {k: anc_idx[0][k] for k in ("block_id", "chr", "start", "end")} == {
    "block_id": 0, "chr": "chr1", "start": 1, "end": 250
}, anc_idx[0]
assert {k: anc_idx[1][k] for k in ("block_id", "chr", "start", "end")} == {
    "block_id": 1, "chr": "chr2", "start": 1, "end": 150
}, anc_idx[1]

common_bin = pathlib.Path(str(prefix) + ".common.geno.bin").read_bytes()
assert struct.unpack("<Q", common_bin) == (0b0011,), common_bin

rare_bin = pathlib.Path(str(prefix) + ".rare.carrier.bin").read_bytes()
records = [
    struct.unpack("<II", rare_bin[i:i + 8])
    for i in range(0, len(rare_bin), 8)
]

expected = [
    (0, (1 << 27) | 1),
    (1, (1 << 27) | 1),
    (2, (1 << 27) | 2),
    (3, 0),
    (5, 3),
    (6, (1 << 27) | 0),
    (7, 2),
    (8, (1 << 27) | 1),
]
assert records == expected, records

vcf_records = [
    line.split("\t")
    for line in roundtrip_vcf
    if line and not line.startswith("#")
]
assert len(vcf_records) == 9, vcf_records
assert [(r[0], r[1], r[2], r[3], r[4], r[9], r[10]) for r in vcf_records] == [
    ("chr1", "50", "pre_A_T", "A", "T", "0|1", "0|0"),
    ("chr1", "100", "multi_A_C", "A", "C", "0|1", "0|0"),
    ("chr1", "100", "multi_A_G", "A", "G", "0|0", "1|0"),
    ("chr1", "150", "rare_A_T", "A", "T", "1|0", "0|0"),
    ("chr1", "160", "common_A_T", "A", "T", "1|1", "0|0"),
    ("chr1", "250", "tail1_A_T", "A", "T", "0|0", "0|1"),
    ("chr2", "50", "should_skip_A_T", "A", "T", "1|0", "0|0"),
    ("chr2", "100", "chr2var_A_T", "A", "T", "0|0", "1|0"),
    ("chr2", "150", "tail2_A_T", "A", "T", "0|1", "0|0"),
], vcf_records
PY

python3 - "$OUT_DIR/tiny.chr1_100_160" <<'PY'
import gzip
import pathlib
import struct
import sys

prefix = pathlib.Path(sys.argv[1])
meta = pathlib.Path(str(prefix) + ".meta").read_text()
assert "extracted_region\tchr1:100-160" in meta, meta
assert "source_hybrid_prefix" in meta, meta

with gzip.open(str(prefix) + ".roundtrip.vcf.gz", "rt") as fh:
    rows = [line.rstrip().split("\t") for line in fh if not line.startswith("#")]
assert [(r[0], r[1], r[2], r[3], r[4], r[9], r[10]) for r in rows] == [
    ("chr1", "100", "multi_A_C", "A", "C", "0|1", "0|0"),
    ("chr1", "100", "multi_A_G", "A", "G", "0|0", "1|0"),
    ("chr1", "150", "rare_A_T", "A", "T", "1|0", "0|0"),
    ("chr1", "160", "common_A_T", "A", "T", "1|1", "0|0"),
], rows

rare_bin = pathlib.Path(str(prefix) + ".rare.carrier.bin").read_bytes()
records = [
    struct.unpack("<II", rare_bin[i:i + 8])
    for i in range(0, len(rare_bin), 8)
]
assert [r[0] for r in records] == [0, 1, 2], records
PY

python3 - "$OUT_DIR/dosage" <<'PY'
import gzip
import pathlib
import sys

prefix = pathlib.Path(sys.argv[1])
meta = pathlib.Path(str(prefix) + ".meta").read_text()
samples = pathlib.Path(str(prefix) + ".samples").read_text().strip().splitlines()
assert "source_tractor_dosage_vcf" in meta, meta
assert samples == ["s1", "s2"], samples
with gzip.open(str(prefix) + ".roundtrip.vcf.gz", "rt") as fh:
    rows = [line.rstrip().split("\t") for line in fh if not line.startswith("#")]
assert [(r[0], r[1], r[2], r[3], r[4], r[9], r[10]) for r in rows] == [
    ("chr1", "100", "v1", "A", "G", "1|0", "0|1"),
    ("chr1", "160", "v2", "C", "T", "1|1", "1|0"),
    ("chr1", "220", "v3", "G", "A", "1|0", "1|0"),
], rows
PY

python3 - "$OUT_DIR/rfmix" <<'PY'
import pathlib
import sys

prefix = pathlib.Path(sys.argv[1])
meta = pathlib.Path(str(prefix) + ".meta").read_text()
samples = pathlib.Path(str(prefix) + ".samples").read_text().strip().splitlines()
assert "source_rfmix_msp" in meta, meta
assert "rfmix_msp_interval_note" in meta, meta
assert "rfmix_subpopulation_order_codes" in meta, meta
assert samples == ["s1", "s2"], samples
PY

if "$BIN_DIR/felixla" \
  from-flare \
  "$ROOT_DIR/testdata/tiny.missing_gt.vcf" \
  "$ROOT_DIR/testdata/tiny.flare.vcf" \
  2 \
  1 \
  "$OUT_DIR/missing_gt" >/dev/null 2>"$OUT_DIR/missing_gt.err"; then
  echo "expected missing GT fixture to fail" >&2
  exit 1
fi

grep -q "missing genotype" "$OUT_DIR/missing_gt.err"

if "$BIN_DIR/felixla" \
  from-flare \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  "$ROOT_DIR/testdata/tiny.missing_flare.vcf" \
  2 \
  1 \
  "$OUT_DIR/missing_flare" >/dev/null 2>"$OUT_DIR/missing_flare.err"; then
  echo "expected missing FLARE fixture to fail" >&2
  exit 1
fi

grep -q "missing FORMAT/AN1" "$OUT_DIR/missing_flare.err"

"$BIN_DIR/felixla" \
  from-flare \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  "$ROOT_DIR/testdata/tiny.flare_swapped_samples.vcf" \
  2 \
  1 \
  "$OUT_DIR/swapped_samples" >/dev/null 2>"$OUT_DIR/swapped_samples.err"

grep -q "Using genotype/FLARE sample intersection" "$OUT_DIR/swapped_samples.err"

"$BIN_DIR/felixla" \
  --felixla "$OUT_DIR/swapped_samples" \
  --query chr1:100 \
  --ref A \
  --alt C \
  > "$OUT_DIR/swapped_samples.query.tsv"

python3 - "$OUT_DIR/swapped_samples" "$OUT_DIR/swapped_samples.query.tsv" <<'PY'
import pathlib
import sys

prefix = pathlib.Path(sys.argv[1])
query = pathlib.Path(sys.argv[2]).read_text().strip().splitlines()
samples = pathlib.Path(str(prefix) + ".samples").read_text().strip().splitlines()

assert samples == ["s1", "s2"], samples
assert query == [
    "global_variant_index\tchr\tpos\tid\tref\talt\tsample\tDSALL\tDS1\tDS2",
    "1\tchr1\t100\tmulti_A_C\tA\tC\ts1\t1\t0\t1",
    "1\tchr1\t100\tmulti_A_C\tA\tC\ts2\t0\t0\t0",
], query
PY

"$BIN_DIR/felixla" \
  from-flare \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  "$ROOT_DIR/testdata/tiny.duplicate_flare.vcf" \
  2 \
  1 \
  "$OUT_DIR/duplicate_flare" >/dev/null

python3 - "$OUT_DIR/duplicate_flare" <<'PY'
import pathlib
import struct
import sys

prefix = pathlib.Path(sys.argv[1])
rare_bin = pathlib.Path(str(prefix) + ".rare.carrier.bin").read_bytes()

def read_string(data, off):
    (n,) = struct.unpack_from("<I", data, off)
    off += 4
    value = data[off:off + n].decode()
    off += n
    return value, off

def read_rare_idx(path):
    data = pathlib.Path(path).read_bytes()
    assert data[:8] == b"TRRAMKS1", data[:8]
    off = 8
    records = []
    while off < len(data):
        rare_index, global_index = struct.unpack_from("<QI", data, off)
        off += 12
        chrom, off = read_string(data, off)
        (pos,) = struct.unpack_from("<q", data, off)
        off += 8
        vid, off = read_string(data, off)
        ref, off = read_string(data, off)
        alt, off = read_string(data, off)
        (alt_index,) = struct.unpack_from("<I", data, off)
        off += 4
        carrier_offset, = struct.unpack_from("<Q", data, off)
        off += 8
        n_carriers, mac = struct.unpack_from("<II", data, off)
        off += 8
        records.append({"global_index": global_index, "n_carriers": n_carriers, "mac": mac})
    return records

def read_anc_idx(path):
    data = pathlib.Path(path).read_bytes()
    assert data[:8] == b"TRANMKS1", data[:8]
    off = 8
    records = []
    while off < len(data):
        (block_id,) = struct.unpack_from("<I", data, off)
        off += 4
        chrom, off = read_string(data, off)
        start, end = struct.unpack_from("<qq", data, off)
        off += 16
        anc_offset, = struct.unpack_from("<Q", data, off)
        off += 8
        records.append({"block_id": block_id, "chr": chrom, "start": start, "end": end})
    return records

anc_idx = read_anc_idx(str(prefix) + ".ancblock.mks")
rare_idx = read_rare_idx(str(prefix) + ".rare.variant.mks")

assert len(anc_idx) == 1, anc_idx
assert {k: anc_idx[0][k] for k in ("block_id", "chr", "start", "end")} == {
    "block_id": 0, "chr": "chr1", "start": 1, "end": 250
}, anc_idx[0]
assert [record["global_index"] for record in rare_idx] == [0, 1, 2, 3, 5], rare_idx

records = [
    struct.unpack("<II", rare_bin[i:i + 8])
    for i in range(0, len(rare_bin), 8)
]
assert records == [
    (0, 1),
    (1, 1),
    (2, 2),
    (3, (1 << 27) | 0),
    (5, (1 << 27) | 3),
], records
PY

echo "tiny smoke test passed"
