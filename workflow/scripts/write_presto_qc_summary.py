#!/usr/bin/env python3
"""Write one-sample pRESTO QC metrics for the tcr-seek MultiQC plugin."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def count_fastq_records(path: Path) -> int | None:
    if not path.exists():
        return None
    with open_text(path) as handle:
        return sum(1 for _ in handle) // 4


def find_fastq(sample_dir: Path, *stems: str) -> Path | None:
    for stem in stems:
        for suffix in (".fastq", ".fastq.gz", ".fq", ".fq.gz"):
            path = sample_dir / f"{stem}{suffix}"
            if path.exists():
                return path
    return None


def count_named_fastq(sample_dir: Path, *stems: str) -> int | None:
    path = find_fastq(sample_dir, *stems)
    if path is None:
        return None
    return count_fastq_records(path)


def count_table_rows(path: Path) -> int | None:
    if not path.exists():
        return None
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            next(reader)
        except StopIteration:
            return 0
        return sum(1 for _ in reader)


def read_final_tsv(path: Path) -> tuple[int | None, int | None, int | None]:
    if not path.exists():
        return None, None, None
    rows = 0
    total_conscount = 0
    c_regions: set[str] = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows += 1
            value = row.get("CONSCOUNT", "")
            try:
                total_conscount += int(float(value))
            except ValueError:
                pass
            c_region = row.get("CREGION", "")
            if c_region:
                c_regions.add(c_region)
    return rows, total_conscount, len(c_regions)


def pct(numerator: int | None, denominator: int | None) -> float | None:
    if numerator is None or denominator in (None, 0):
        return None
    return numerator / denominator


def value_or_blank(value: int | float | str | None) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--sample-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    sample = args.sample
    sample_dir = Path(args.sample_dir)
    parsed_logs = sample_dir / "parsed_logs"

    metrics: dict[str, int | float | str | None] = {
        "Sample": sample,
        "R1InputReads": count_named_fastq(sample_dir, f"{sample}_R1.input"),
        "R2InputReads": count_named_fastq(sample_dir, f"{sample}_R2.input"),
        "R1Q20Reads": count_named_fastq(sample_dir, f"{sample}_R1.q20"),
        "R2Q20Reads": count_named_fastq(sample_dir, f"{sample}_R2.q20"),
        "R1PrimerMaskedReads": count_named_fastq(sample_dir, f"{sample}_R1.mask-pass", f"{sample}_R1.mask"),
        "R2PrimerCutReads": count_named_fastq(sample_dir, f"{sample}_R2.cut-pass", f"{sample}_R2.cut"),
        "InitialPairPassReads": count_named_fastq(sample_dir, f"{sample}_R1.mask_pair-pass"),
        "R1ConsensusReads": count_named_fastq(sample_dir, f"{sample}_R1.consensus-pass", f"{sample}_R1.consensus"),
        "R2ConsensusReads": count_named_fastq(sample_dir, f"{sample}_R2.consensus-pass", f"{sample}_R2.consensus"),
        "ConsensusPairPassReads": count_named_fastq(sample_dir, f"{sample}_R1.consensus_pair-pass"),
        "AssembledReads": count_named_fastq(sample_dir, f"{sample}.assembled-pass", f"{sample}.assembled"),
        "CRegionReads": count_named_fastq(sample_dir, f"{sample}.Cprimer-pass", f"{sample}.Cprimer"),
        "QualityLogR1Records": count_table_rows(parsed_logs / "quality-1_table.tab"),
        "QualityLogR2Records": count_table_rows(parsed_logs / "quality-2_table.tab"),
        "PrimerLogR1Records": count_table_rows(parsed_logs / "consensus-1_table.tab"),
        "PrimerLogR2Records": count_table_rows(parsed_logs / "consensus-2_table.tab"),
        "AssemblyLogRecords": count_table_rows(parsed_logs / "assemble_table.tab"),
        "CRegionLogRecords": count_table_rows(parsed_logs / "cregion_table.tab"),
    }
    final_rows, final_abundance, c_region_count = read_final_tsv(sample_dir / f"{sample}_uniqC_atleast-2.tsv")
    metrics["FinalConsCountGe2Sequences"] = final_rows
    metrics["FinalConsCountGe2Abundance"] = final_abundance
    metrics["CRegionCount"] = c_region_count
    metrics["R1Q20PassFraction"] = pct(metrics["R1Q20Reads"], metrics["R1InputReads"])  # type: ignore[arg-type]
    metrics["R2Q20PassFraction"] = pct(metrics["R2Q20Reads"], metrics["R2InputReads"])  # type: ignore[arg-type]
    metrics["InitialPairPassFraction"] = pct(metrics["InitialPairPassReads"], metrics["R1Q20Reads"])  # type: ignore[arg-type]
    metrics["AssemblyFraction"] = pct(metrics["AssembledReads"], metrics["ConsensusPairPassReads"])  # type: ignore[arg-type]
    metrics["FinalRecoveryFraction"] = pct(metrics["FinalConsCountGe2Sequences"], metrics["R1InputReads"])  # type: ignore[arg-type]

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(metrics))
        writer.writeheader()
        writer.writerow({key: value_or_blank(value) for key, value in metrics.items()})


if __name__ == "__main__":
    main()
