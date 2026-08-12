#!/usr/bin/env python3
"""Aggregate tcr-seek QC metrics for the MultiQC plugin."""

from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def as_number(value: str) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def fraction(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return ""
    return f"{numerator / denominator:.6g}"


def summarize_changeo(path: Path) -> dict[str, str]:
    if not path.exists():
        return {
            "ChangeoDbPassSequences": "",
            "ProductiveFraction": "",
            "InFrameFraction": "",
            "StopCodonFraction": "",
            "MeanConsensusCount": "",
            "MedianConsensusCount": "",
            "TopVCall": "",
            "TopJCall": "",
            "TopLocus": "",
        }

    total = productive = in_frame = stop_codon = 0
    consensus_counts: list[float] = []
    v_counts: dict[str, int] = {}
    j_counts: dict[str, int] = {}
    locus_counts: dict[str, int] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            total += 1
            if row.get("productive", "").upper() == "T":
                productive += 1
            if row.get("vj_in_frame", "").upper() == "T":
                in_frame += 1
            if row.get("stop_codon", "").upper() == "T":
                stop_codon += 1
            count = as_number(row.get("consensus_count", ""))
            if count is not None:
                consensus_counts.append(count)
            for column, counts in (("v_call", v_counts), ("j_call", j_counts), ("locus", locus_counts)):
                value = row.get(column, "")
                if value:
                    value = value.split(",")[0]
                    counts[value] = counts.get(value, 0) + 1

    def top_value(counts: dict[str, int]) -> str:
        if not counts:
            return ""
        return sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0]

    return {
        "ChangeoDbPassSequences": str(total),
        "ProductiveFraction": fraction(productive, total),
        "InFrameFraction": fraction(in_frame, total),
        "StopCodonFraction": fraction(stop_codon, total),
        "MeanConsensusCount": f"{statistics.mean(consensus_counts):.6g}" if consensus_counts else "",
        "MedianConsensusCount": f"{statistics.median(consensus_counts):.6g}" if consensus_counts else "",
        "TopVCall": top_value(v_counts),
        "TopJCall": top_value(j_counts),
        "TopLocus": top_value(locus_counts),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--samples", required=True, help="Comma-separated sample names")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    samples = [sample for sample in args.samples.split(",") if sample]
    rows: list[dict[str, str]] = []
    for sample in samples:
        presto_qc = run_dir / "01_run" / sample / f"{sample}.presto_qc.tsv"
        row = read_tsv(presto_qc)[0] if presto_qc.exists() else {"Sample": sample}
        row.update(summarize_changeo(run_dir / "02_changeo" / sample / f"{sample}.al2_db-pass.tsv"))
        final_sequences = as_number(row.get("FinalConsCountGe2Sequences", ""))
        db_pass = as_number(row.get("ChangeoDbPassSequences", ""))
        if final_sequences is not None and db_pass is not None and final_sequences > 0:
            row["ChangeoPassFraction"] = f"{db_pass / final_sequences:.6g}"
        else:
            row["ChangeoPassFraction"] = ""
        rows.append(row)

    fieldnames = [
        "Sample",
        "R1InputReads",
        "R2InputReads",
        "R1Q20Reads",
        "R2Q20Reads",
        "R1Q20PassFraction",
        "R2Q20PassFraction",
        "InitialPairPassReads",
        "InitialPairPassFraction",
        "R1ConsensusReads",
        "R2ConsensusReads",
        "ConsensusPairPassReads",
        "AssembledReads",
        "AssemblyFraction",
        "CRegionReads",
        "FinalConsCountGe2Sequences",
        "FinalConsCountGe2Abundance",
        "FinalRecoveryFraction",
        "ChangeoDbPassSequences",
        "ChangeoPassFraction",
        "ProductiveFraction",
        "InFrameFraction",
        "StopCodonFraction",
        "MeanConsensusCount",
        "MedianConsensusCount",
        "CRegionCount",
        "TopVCall",
        "TopJCall",
        "TopLocus",
    ]
    write_tsv(Path(args.output), rows, fieldnames)


if __name__ == "__main__":
    main()
