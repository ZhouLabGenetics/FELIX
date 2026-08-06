#!/usr/bin/env python3
"""Intense keep/extract regression tests for FELIXla FLARE conversion."""

from __future__ import annotations

import argparse
import gzip
import math
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile


N_ANCESTRIES = 3
ALLELES = ["A", "C", "G", "T"]


def fail(message: str) -> None:
    raise AssertionError(message)


def run(cmd: list[str], *, expect_fail: bool = False, contains: str | None = None) -> subprocess.CompletedProcess:
    result = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if expect_fail:
        if result.returncode == 0:
            fail(f"expected failure but command succeeded: {' '.join(cmd)}")
        if contains and contains not in (result.stdout + result.stderr):
            fail(
                f"expected failed command output to contain {contains!r}\n"
                f"command: {' '.join(cmd)}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
            )
        return result

    if result.returncode != 0:
        fail(
            f"command failed: {' '.join(cmd)}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result


def opener(path: pathlib.Path, mode: str = "rt"):
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return open(path, mode)


def alt_list(ref: str, count: int, variant_index: int) -> list[str]:
    choices = [a for a in ALLELES if a != ref]
    offset = variant_index % len(choices)
    rotated = choices[offset:] + choices[:offset]
    return rotated[:count]


def build_inputs(work: pathlib.Path):
    samples = [f"s{i:03d}" for i in range(1, 38)]
    endpoints = {
        "chr1": [250, 500, 750, 1000],
        "chr2": [230, 520, 770, 1000],
    }
    chrom_offset = {"chr1": 0, "chr2": 1}

    def ancestry(chrom: str, pos: int, sample_index: int, hap: int) -> int:
        endpoint_index = next(i for i, end in enumerate(endpoints[chrom]) if pos <= end)
        if hap == 0:
            return (sample_index + endpoint_index + chrom_offset[chrom]) % N_ANCESTRIES
        return (2 * sample_index + endpoint_index + chrom_offset[chrom]) % N_ANCESTRIES

    records = []
    positions = {
        "chr1": [40, 85, 130, 175, 220, 265, 310, 355, 410, 465, 530, 585, 640, 695, 760, 825, 890, 955],
        "chr2": [35, 95, 155, 215, 275, 335, 395, 455, 535, 595, 655, 715, 775, 835, 895, 955],
    }
    variant_index = 0
    for chrom in ["chr1", "chr2"]:
        for pos in positions[chrom]:
            ref = ALLELES[(variant_index + (1 if chrom == "chr2" else 0)) % len(ALLELES)]
            alt_count = [1, 2, 1, 3, 1, 2, 1][variant_index % 7]
            alts = alt_list(ref, alt_count, variant_index)
            raw_id = "." if variant_index % 10 == 4 else f"v{variant_index:03d}"
            n_alleles = 1 + len(alts)
            gts: list[tuple[int, int]] = []
            for sample_index, _sample in enumerate(samples):
                if variant_index % 13 == 0:
                    a0 = 1 if sample_index == 0 else 0
                    a1 = 0
                elif variant_index % 17 == 0 and n_alleles > 2:
                    a0 = 2 if sample_index == len(samples) - 1 else 0
                    a1 = 0
                else:
                    a0 = (sample_index + 2 * variant_index) % n_alleles
                    a1 = (2 * sample_index + variant_index + 1) % n_alleles
                gts.append((a0, a1))
            records.append(
                {
                    "chrom": chrom,
                    "pos": pos,
                    "id": raw_id,
                    "ref": ref,
                    "alts": alts,
                    "gts": gts,
                    "variant_index": variant_index,
                }
            )
            variant_index += 1

    genotype_path = work / "synthetic.phased.vcf"
    with genotype_path.open("w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##contig=<ID=chr1>\n")
        out.write("##contig=<ID=chr2>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        out.write("\t".join(samples))
        out.write("\n")
        for record in records:
            sample_gts = [f"{a0}|{a1}" for a0, a1 in record["gts"]]
            out.write(
                f"{record['chrom']}\t{record['pos']}\t{record['id']}\t"
                f"{record['ref']}\t{','.join(record['alts'])}\t.\tPASS\t.\tGT\t"
                + "\t".join(sample_gts)
                + "\n"
            )

    flare_path = work / "synthetic.flare.vcf"
    with flare_path.open("w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##contig=<ID=chr1>\n")
        out.write("##contig=<ID=chr2>\n")
        out.write('##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="First haplotype local ancestry">\n')
        out.write('##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="Second haplotype local ancestry">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        out.write("\t".join(samples))
        out.write("\n")
        for chrom in ["chr1", "chr2"]:
            for pos in endpoints[chrom]:
                sample_ans = [
                    f"{ancestry(chrom, pos, sample_index, 0)}:{ancestry(chrom, pos, sample_index, 1)}"
                    for sample_index, _sample in enumerate(samples)
                ]
                out.write(f"{chrom}\t{pos}\t.\tA\tC\t.\tPASS\t.\tAN1:AN2\t")
                out.write("\t".join(sample_ans))
                out.write("\n")

    return samples, records, ancestry, genotype_path, flare_path


def split_id(record: dict, alt: str) -> str:
    if record["id"] == ".":
        return f"{record['chrom']}:{record['pos']}:{record['ref']}:{alt}"
    return f"{record['id']}_{record['ref']}_{alt}"


def split_records(
    records: list[dict],
    sample_indices: list[int],
    *,
    selected: set[tuple[str, int, str, str]] | None = None,
    region: tuple[str, int, int] | None = None,
) -> list[dict]:
    result = []
    for record in records:
        if region:
            r_chrom, r_start, r_end = region
            if record["chrom"] != r_chrom or not (r_start <= record["pos"] <= r_end):
                continue
        for alt_index, alt in enumerate(record["alts"], start=1):
            key = (record["chrom"], record["pos"], record["ref"], alt)
            if selected is not None and key not in selected:
                continue
            gts = []
            for sample_index in sample_indices:
                a0, a1 = record["gts"][sample_index]
                gts.append(f"{1 if a0 == alt_index else 0}|{1 if a1 == alt_index else 0}")
            result.append(
                {
                    "global_index": len(result),
                    "chrom": record["chrom"],
                    "pos": record["pos"],
                    "id": split_id(record, alt),
                    "ref": record["ref"],
                    "alt": alt,
                    "alt_index": alt_index,
                    "gts": gts,
                    "record": record,
                }
            )
    return result


def read_roundtrip_vcf(path: pathlib.Path):
    with opener(path) as fh:
        lines = [line.rstrip("\n") for line in fh]
    header = next(line.split("\t") for line in lines if line.startswith("#CHROM"))
    rows = []
    for line in lines:
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        rows.append(
            {
                "chrom": fields[0],
                "pos": int(fields[1]),
                "id": fields[2],
                "ref": fields[3],
                "alt": fields[4],
                "gts": fields[9:],
            }
        )
    return header[9:], rows


def assert_vcf_matches(path: pathlib.Path, expected_samples: list[str], expected: list[dict]) -> None:
    actual_samples, actual_rows = read_roundtrip_vcf(path)
    if actual_samples != expected_samples:
        fail(f"samples mismatch for {path}: {actual_samples} != {expected_samples}")
    compact_actual = [
        (r["chrom"], r["pos"], r["id"], r["ref"], r["alt"], tuple(r["gts"]))
        for r in actual_rows
    ]
    compact_expected = [
        (r["chrom"], r["pos"], r["id"], r["ref"], r["alt"], tuple(r["gts"]))
        for r in expected
    ]
    if compact_actual != compact_expected:
        fail(f"VCF rows mismatch for {path}\nactual={compact_actual[:8]}\nexpected={compact_expected[:8]}")


def read_meta(prefix: pathlib.Path) -> dict[str, str]:
    meta = {}
    with pathlib.Path(str(prefix) + ".meta").open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t", 1)
            if len(parts) == 2:
                meta[parts[0]] = parts[1]
    return meta


def read_samples(prefix: pathlib.Path) -> list[str]:
    return pathlib.Path(str(prefix) + ".samples").read_text().splitlines()


def export_prefix(bin_dir: pathlib.Path, prefix: pathlib.Path) -> pathlib.Path:
    out_vcf = pathlib.Path(str(prefix) + ".roundtrip.vcf.gz")
    run([str(bin_dir / "felixla"), "to-vcf", str(prefix), str(out_vcf)])
    return out_vcf


def make_selected_alleles(records: list[dict]) -> set[tuple[str, int, str, str]]:
    selected: set[tuple[str, int, str, str]] = set()
    multi_records = [record for record in records if len(record["alts"]) > 1]
    first_multi = multi_records[0]
    selected.add((first_multi["chrom"], first_multi["pos"], first_multi["ref"], first_multi["alts"][1]))

    second_multi = multi_records[3]
    for alt in second_multi["alts"][:2]:
        selected.add((second_multi["chrom"], second_multi["pos"], second_multi["ref"], alt))

    for record in records[5::6][:8]:
        selected.add((record["chrom"], record["pos"], record["ref"], record["alts"][0]))

    zero_after_keep = next(record for record in records if record["variant_index"] % 17 == 0 and len(record["alts"]) > 1)
    selected.add((zero_after_keep["chrom"], zero_after_keep["pos"], zero_after_keep["ref"], zero_after_keep["alts"][1]))
    return selected


def write_extract_files(
    work: pathlib.Path,
    selected: set[tuple[str, int, str, str]],
    stem: str = "extract.sites",
) -> tuple[pathlib.Path, pathlib.Path]:
    by_site: dict[tuple[str, int, str], list[str]] = {}
    for chrom, pos, ref, alt in sorted(selected):
        by_site.setdefault((chrom, pos, ref), []).append(alt)

    pvar = work / f"{stem}.pvar"
    with pvar.open("w") as out:
        out.write("#CHROM\tPOS\tID\tREF\tALT\n")
        for (chrom, pos, ref), alts in sorted(by_site.items()):
            out.write(
                f"{chrom}\t{pos}\tpvar_id_must_be_ignored_{chrom}_{pos}"
                f"\t{ref}\t{','.join(alts)}\n"
            )

    vcf = work / f"{stem}.vcf"
    with vcf.open("w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tignored\n")
        for (chrom, pos, ref), alts in sorted(by_site.items()):
            out.write(
                f"{chrom}\t{pos}\tvcf_id_must_be_ignored_{chrom}_{pos}"
                f"\t{ref}\t{','.join(alts)}\t.\tPASS\t.\tGT\t0|0\n"
            )

    vcfgz = work / f"{stem}.vcf.gz"
    with vcf.open("rb") as src, gzip.open(vcfgz, "wb") as dst:
        shutil.copyfileobj(src, dst)

    return pvar, vcfgz


def write_keep_file(work: pathlib.Path, samples: list[str], keep_indices: list[int]) -> pathlib.Path:
    keep = work / "keep.samples"
    scrambled = [keep_indices[i] for i in [5, 0, 16, 2, 9, 20, 1, 12, 26, 3, 7, 10, 30, 4, 6, 8, 11, 13, 14, 15, 17, 18, 19, 21, 22, 23, 24, 25, 27, 28, 29, 31, 32]]
    keep.write_text("\n".join(samples[i] for i in scrambled) + "\n")
    return keep


def write_reordered_subset_flare(
    work: pathlib.Path,
    flare_path: pathlib.Path,
    sample_order_indices: list[int],
) -> pathlib.Path:
    out_path = work / "synthetic.flare.reordered_subset.vcf"
    with flare_path.open() as src, out_path.open("w") as out:
        for line in src:
            line = line.rstrip("\n")
            if line.startswith("#CHROM"):
                fields = line.split("\t")
                header_samples = fields[9:]
                reordered_samples = [header_samples[i] for i in sample_order_indices]
                out.write("\t".join(fields[:9] + reordered_samples) + "\n")
            elif line.startswith("#"):
                out.write(line + "\n")
            else:
                fields = line.split("\t")
                sample_values = fields[9:]
                reordered_values = [sample_values[i] for i in sample_order_indices]
                out.write("\t".join(fields[:9] + reordered_values) + "\n")
    return out_path


def query_rows(
    bin_dir: pathlib.Path,
    prefix: pathlib.Path,
    split: dict,
    samples: list[str],
    sample_indices: list[int],
    ancestry,
    *,
    nonzero_only: bool = False,
) -> list[list[str]]:
    cmd = [
        str(bin_dir / "felixla"),
        "--felixla",
        str(prefix),
        "--query",
        f"{split['chrom']}:{split['pos']}",
        "--ref",
        split["ref"],
        "--alt",
        split["alt"],
    ]
    if nonzero_only:
        cmd.append("--nonzero-only")
    actual = run(cmd).stdout.strip().splitlines()
    expected_header = "global_variant_index\tchr\tpos\tid\tref\talt\tsample\tDSALL\tDS1\tDS2\tDS3"
    if not actual or actual[0] != expected_header:
        fail(f"unexpected query header for {prefix}: {actual[:1]}")
    actual_rows = [line.split("\t") for line in actual[1:]]

    expected_rows = []
    for out_i, sample_index in enumerate(sample_indices):
        gt = split["gts"][out_i]
        h0 = 1 if gt[0] == "1" else 0
        h1 = 1 if gt[2] == "1" else 0
        ds = [0, 0, 0]
        ds[ancestry(split["chrom"], split["pos"], sample_index, 0)] += h0
        ds[ancestry(split["chrom"], split["pos"], sample_index, 1)] += h1
        dsall = h0 + h1
        if nonzero_only and dsall == 0:
            continue
        expected_rows.append(
            [
                str(split["global_index"]),
                split["chrom"],
                str(split["pos"]),
                split["id"],
                split["ref"],
                split["alt"],
                samples[sample_index],
                str(dsall),
                *(str(x) for x in ds),
            ]
        )
    if actual_rows != expected_rows:
        fail(f"query mismatch for {prefix} {split['chrom']}:{split['pos']} {split['ref']}>{split['alt']}")
    return actual_rows


def check_queries(
    bin_dir: pathlib.Path,
    prefix: pathlib.Path,
    expected: list[dict],
    samples: list[str],
    sample_indices: list[int],
    ancestry,
) -> None:
    targets = expected[:3] + expected[len(expected) // 2: len(expected) // 2 + 2] + expected[-3:]
    seen = set()
    unique_targets = []
    for target in targets:
        key = (target["chrom"], target["pos"], target["ref"], target["alt"])
        if key not in seen:
            seen.add(key)
            unique_targets.append(target)
    for target in unique_targets:
        query_rows(bin_dir, prefix, target, samples, sample_indices, ancestry)

    nonzero_target = next(target for target in expected if any(gt != "0|0" for gt in target["gts"]))
    query_rows(bin_dir, prefix, nonzero_target, samples, sample_indices, ancestry, nonzero_only=True)


def build_with_cli(
    bin_dir: pathlib.Path,
    genotype_path: pathlib.Path,
    flare_path: pathlib.Path,
    prefix: pathlib.Path,
    *extra: str,
) -> pathlib.Path:
    run(
        [
            str(bin_dir / "felixla"),
            "--phase-vcf",
            str(genotype_path),
            "--flare-vcf",
            str(flare_path),
            "--n-ancestries",
            str(N_ANCESTRIES),
            *extra,
            "--make-felixla",
            "--out",
            str(prefix),
        ]
    )
    return export_prefix(bin_dir, prefix)


def bounded_extract_expected_from_split_records(
    split_records_all: list[dict],
    selected: set[tuple[str, int, str, str]],
) -> list[dict]:
    expected = []
    for split in split_records_all:
        key = (split["chrom"], split["pos"], split["ref"], split["alt"])
        if key not in selected:
            continue
        copied = dict(split)
        copied["global_index"] = len(expected)
        copied["alt_index"] = 1
        copied["id"] = f"{split['id']}_{split['ref']}_{split['alt']}"
        expected.append(copied)
    return expected


def bounded_extract_scanned_records(
    split_records_all: list[dict],
    selected: set[tuple[str, int, str, str]],
) -> int:
    bounds: dict[str, list[int]] = {}
    for chrom, pos, _ref, _alt in selected:
        if chrom not in bounds:
            bounds[chrom] = [pos, pos]
        else:
            bounds[chrom][0] = min(bounds[chrom][0], pos)
            bounds[chrom][1] = max(bounds[chrom][1], pos)

    return sum(
        1 for split in split_records_all
        if split["chrom"] in bounds and
        bounds[split["chrom"]][0] <= split["pos"] <= bounds[split["chrom"]][1]
    )


def last_progress_record_count(stderr: str) -> int:
    matches = re.findall(r"Progress: records (\d+)", stderr)
    if not matches:
        fail(f"missing progress records line in stderr:\n{stderr}")
    return int(matches[-1])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin-dir", type=pathlib.Path, default=pathlib.Path("bin"))
    parser.add_argument("--work-dir", type=pathlib.Path)
    parser.add_argument("--keep-work", action="store_true")
    args = parser.parse_args()

    owned_tmp = None
    if args.work_dir:
        work = args.work_dir
        work.mkdir(parents=True, exist_ok=True)
    else:
        owned_tmp = tempfile.TemporaryDirectory(prefix="felixla-keep-extract-intense.")
        work = pathlib.Path(owned_tmp.name)

    try:
        bin_dir = args.bin_dir.resolve()
        if not (bin_dir / "felixla").exists():
            fail(f"missing tool: {bin_dir / 'felixla'}")
        unexpected = sorted(
            path.name for path in bin_dir.iterdir()
            if path.name != "felixla"
        )
        if unexpected:
            fail(f"unexpected files in {bin_dir}; only felixla should be present: {unexpected}")

        samples, records, ancestry, genotype_path, flare_path = build_inputs(work)
        all_indices = list(range(len(samples)))
        keep_indices = [i for i in all_indices if i not in {5, 10, 26, 36}]
        keep_samples = [samples[i] for i in keep_indices]
        keep_path = write_keep_file(work, samples, keep_indices)
        selected = make_selected_alleles(records)
        pvar_path, vcfgz_path = write_extract_files(work, selected)

        legacy_prefix = work / "full.legacy"
        run(
            [
                str(bin_dir / "felixla"),
                "from-flare",
                str(genotype_path),
                str(flare_path),
                str(N_ANCESTRIES),
                "auto",
                str(legacy_prefix),
            ]
        )
        legacy_vcf = export_prefix(bin_dir, legacy_prefix)
        full_expected = split_records(records, all_indices)
        assert_vcf_matches(legacy_vcf, samples, full_expected)
        check_queries(bin_dir, legacy_prefix, full_expected, samples, all_indices, ancestry)
        legacy_meta = read_meta(legacy_prefix)
        assert legacy_meta["n_samples"] == str(len(samples)), legacy_meta
        assert legacy_meta["n_words"] == "2", legacy_meta
        assert legacy_meta["rare_threshold"] == str(math.ceil(len(samples) / 32)), legacy_meta
        assert "keep_samples" not in legacy_meta and "extract_sites" not in legacy_meta, legacy_meta

        cli_prefix = work / "full.cli"
        cli_vcf = build_with_cli(bin_dir, genotype_path, flare_path, cli_prefix)
        assert_vcf_matches(cli_vcf, samples, full_expected)
        if read_roundtrip_vcf(cli_vcf) != read_roundtrip_vcf(legacy_vcf):
            fail("full CLI output differs from legacy positional output")

        subset_expected = split_records(records, keep_indices, selected=selected)
        subset_prefix = work / "subset.vcfgz"
        subset_vcf = build_with_cli(
            bin_dir,
            genotype_path,
            flare_path,
            subset_prefix,
            "--keep",
            str(keep_path),
            "--extract",
            str(vcfgz_path),
        )
        assert_vcf_matches(subset_vcf, keep_samples, subset_expected)
        check_queries(bin_dir, subset_prefix, subset_expected, samples, keep_indices, ancestry)
        subset_meta = read_meta(subset_prefix)
        assert subset_meta["n_samples"] == str(len(keep_indices)), subset_meta
        assert subset_meta["n_words"] == "2", subset_meta
        assert subset_meta["rare_threshold"] == str(math.ceil(len(keep_indices) / 32)), subset_meta
        assert subset_meta["keep_samples"] == str(keep_path), subset_meta
        assert subset_meta["extract_sites"] == str(vcfgz_path), subset_meta
        assert read_samples(subset_prefix) == keep_samples

        pvar_prefix = work / "subset.pvar"
        pvar_vcf = build_with_cli(
            bin_dir,
            genotype_path,
            flare_path,
            pvar_prefix,
            "--keep",
            str(keep_path),
            "--extract",
            str(pvar_path),
        )
        assert_vcf_matches(pvar_vcf, keep_samples, subset_expected)
        if read_roundtrip_vcf(pvar_vcf) != read_roundtrip_vcf(subset_vcf):
            fail("PVAR extract output differs from gzipped VCF extract output")

        missing_from_flare = {keep_indices[1], keep_indices[8], keep_indices[-2]}
        flare_subset_indices = [i for i in all_indices if i not in missing_from_flare]
        reordered_flare_indices = (
            flare_subset_indices[3::5]
            + flare_subset_indices[0::5]
            + flare_subset_indices[4::5]
            + flare_subset_indices[1::5]
            + flare_subset_indices[2::5]
        )
        reordered_flare_path = write_reordered_subset_flare(work, flare_path, reordered_flare_indices)
        intersection_indices = [i for i in keep_indices if i not in missing_from_flare]
        intersection_samples = [samples[i] for i in intersection_indices]
        intersection_expected = split_records(records, intersection_indices, selected=selected)
        intersection_prefix = work / "subset.reordered_flare_intersection"
        intersection_vcf = build_with_cli(
            bin_dir,
            genotype_path,
            reordered_flare_path,
            intersection_prefix,
            "--keep",
            str(keep_path),
            "--extract",
            str(pvar_path),
        )
        assert_vcf_matches(intersection_vcf, intersection_samples, intersection_expected)
        check_queries(bin_dir, intersection_prefix, intersection_expected, samples, intersection_indices, ancestry)
        intersection_meta = read_meta(intersection_prefix)
        assert intersection_meta["n_samples"] == str(len(intersection_indices)), intersection_meta
        assert intersection_meta["n_words"] == str(math.ceil(2 * len(intersection_indices) / 64)), intersection_meta
        assert intersection_meta["rare_threshold"] == str(math.ceil(len(intersection_indices) / 32)), intersection_meta
        assert read_samples(intersection_prefix) == intersection_samples

        indexed_selected = {
            (split["chrom"], split["pos"], split["ref"], split["alt"])
            for split in [full_expected[1], full_expected[7], full_expected[19], full_expected[-3]]
        }
        indexed_extract_pvar, _indexed_extract_vcfgz = write_extract_files(
            work,
            indexed_selected,
            "indexed.extract.sites",
        )
        indexed_expected = bounded_extract_expected_from_split_records(full_expected, indexed_selected)
        indexed_prefix = work / "subset.indexed_extract"
        indexed_result = run(
            [
                str(bin_dir / "felixla"),
                "--phase-vcf",
                str(legacy_vcf),
                "--flare-vcf",
                str(flare_path),
                "--n-ancestries",
                str(N_ANCESTRIES),
                "--extract",
                str(indexed_extract_pvar),
                "--make-felixla",
                "--out",
                str(indexed_prefix),
            ]
        )
        if "Using indexed --extract bounded genotype reader" not in indexed_result.stderr:
            fail(f"indexed bounded --extract reader was not used:\n{indexed_result.stderr}")
        scanned_records = last_progress_record_count(indexed_result.stderr)
        expected_scanned_records = bounded_extract_scanned_records(full_expected, indexed_selected)
        if scanned_records != expected_scanned_records:
            fail(
                f"indexed --extract scanned {scanned_records} records; "
                f"expected exactly {expected_scanned_records} records inside bounded extract intervals"
            )
        if scanned_records >= len(full_expected):
            fail(f"bounded indexed --extract scanned the full VCF: {scanned_records} of {len(full_expected)} records")
        indexed_vcf = export_prefix(bin_dir, indexed_prefix)
        assert_vcf_matches(indexed_vcf, samples, indexed_expected)

        region = ("chr1", 200, 650)
        region_expected = split_records(records, keep_indices, selected=selected, region=region)
        region_prefix = work / "subset.region"
        region_vcf = build_with_cli(
            bin_dir,
            genotype_path,
            flare_path,
            region_prefix,
            "--region",
            f"{region[0]}:{region[1]}-{region[2]}",
            "--keep",
            str(keep_path),
            "--extract",
            str(vcfgz_path),
        )
        assert_vcf_matches(region_vcf, keep_samples, region_expected)
        check_queries(bin_dir, region_prefix, region_expected, samples, keep_indices, ancestry)

        bad_keep_dup = work / "bad.keep.dup"
        bad_keep_dup.write_text(f"{samples[0]}\n{samples[0]}\n")
        build_bad_base = [
            str(bin_dir / "felixla"),
            "--phase-vcf",
            str(genotype_path),
            "--flare-vcf",
            str(flare_path),
            "--n-ancestries",
            str(N_ANCESTRIES),
            "--make-felixla",
            "--out",
        ]
        run([*build_bad_base, str(work / "bad.keep.dup.out"), "--keep", str(bad_keep_dup)], expect_fail=True, contains="duplicate sample ID")

        bad_keep_extra = work / "bad.keep.extra"
        bad_keep_extra.write_text(f"{samples[0]} extra\n")
        run([*build_bad_base, str(work / "bad.keep.extra.out"), "--keep", str(bad_keep_extra)], expect_fail=True, contains="one sample ID per line")

        bad_keep_missing = work / "bad.keep.missing"
        bad_keep_missing.write_text("not_a_sample\n")
        run([*build_bad_base, str(work / "bad.keep.missing.out"), "--keep", str(bad_keep_missing)], expect_fail=True, contains="absent from genotype VCF")

        bad_ref_unknown = work / "bad.ref.unknown.pvar"
        first_key = next(iter(selected))
        bad_ref_unknown.write_text(f"#CHROM\tPOS\tID\tREF\tALT\n{first_key[0]}\t{first_key[1]}\t.\t.\t{first_key[3]}\n")
        run([*build_bad_base, str(work / "bad.ref.unknown.out"), "--extract", str(bad_ref_unknown)], expect_fail=True, contains="requires known REF")

        bad_ref_mismatch = work / "bad.ref.mismatch.pvar"
        wrong_ref = next(a for a in ALLELES if a != first_key[2])
        bad_ref_mismatch.write_text(f"#CHROM\tPOS\tID\tREF\tALT\n{first_key[0]}\t{first_key[1]}\t.\t{wrong_ref}\t{first_key[3]}\n")
        run([*build_bad_base, str(work / "bad.ref.mismatch.out"), "--extract", str(bad_ref_mismatch)], expect_fail=True, contains="REF mismatch")

        bad_ref_conflict = work / "bad.ref.conflict.pvar"
        bad_ref_conflict.write_text(
            f"#CHROM\tPOS\tID\tREF\tALT\n"
            f"{first_key[0]}\t{first_key[1]}\t.\t{first_key[2]}\t{first_key[3]}\n"
            f"{first_key[0]}\t{first_key[1]}\t.\t{wrong_ref}\t{first_key[3]}\n"
        )
        run([*build_bad_base, str(work / "bad.ref.conflict.out"), "--extract", str(bad_ref_conflict)], expect_fail=True, contains="conflicting REF")

        bad_dup_allele = work / "bad.dup.allele.pvar"
        bad_dup_allele.write_text(
            f"#CHROM\tPOS\tID\tREF\tALT\n"
            f"{first_key[0]}\t{first_key[1]}\tfirst_id\t{first_key[2]}\t{first_key[3]}\n"
            f"{first_key[0]}\t{first_key[1]}\tsecond_id\t{first_key[2]}\t{first_key[3]}\n"
        )
        run([*build_bad_base, str(work / "bad.dup.allele.out"), "--extract", str(bad_dup_allele)], expect_fail=True, contains="duplicate allele")

        run(
            [
                str(bin_dir / "felixla"),
                "--felixla",
                str(cli_prefix),
                "--export",
                "vcf",
                "--keep",
                str(keep_path),
                "--out",
                str(work / "bad.unsupported.vcf.gz"),
            ],
            expect_fail=True,
            contains="supported only for --phase-vcf + --flare-vcf --make-felixla",
        )

        print(f"intense keep/extract regression passed in {work}")
    finally:
        if owned_tmp and not args.keep_work:
            owned_tmp.cleanup()

    return 0


if __name__ == "__main__":
    sys.exit(main())
