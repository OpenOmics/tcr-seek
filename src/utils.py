from __future__ import annotations

import csv
import os
import re
from pathlib import Path


def sanitize_sample(sample: str) -> str:
    sample = re.sub(r"[^A-Za-z0-9_.-]+", "_", sample)
    sample = sample.strip("._-")
    if not sample:
        raise ValueError("Empty sample name after FASTQ parsing")
    return sample


def normalize_fastq_name(filename: str) -> str:
    """Normalize FASTQ names to the cell-seek style: sample.R1.fastq.gz."""
    name = Path(filename).name
    if name.endswith((".R1.fastq.gz", ".R2.fastq.gz")):
        return name

    # Mirrors cell-seek's supported conventions while keeping the match anchored.
    patterns = [
        (r"[._-]R1(?:[._-]001)?\.f(?:ast)?q\.gz$", ".R1.fastq.gz"),
        (r"[._-]R2(?:[._-]001)?\.f(?:ast)?q\.gz$", ".R2.fastq.gz"),
        (r"[._-]R1(?:[._-][A-Za-z0-9]+)*\.f(?:ast)?q\.gz$", ".R1.fastq.gz"),
        (r"[._-]R2(?:[._-][A-Za-z0-9]+)*\.f(?:ast)?q\.gz$", ".R2.fastq.gz"),
        (r"_1\.f(?:ast)?q\.gz$", ".R1.fastq.gz"),
        (r"_2\.f(?:ast)?q\.gz$", ".R2.fastq.gz"),
    ]
    for pattern, replacement in patterns:
        if re.search(pattern, name):
            return re.sub(pattern, replacement, name)

    raise ValueError(
        "Could not normalize FASTQ name. Expected examples include "
        "sample.R1.fastq.gz, sample_R1_001.fastq.gz, or sample_1.fastq.gz: "
        f"{filename}"
    )


def sample_from_normalized_fastq(filename: str, preserve_samples: set[str] | None = None) -> tuple[str, str]:
    name = Path(filename).name
    base_match = re.match(r"(.+?)\.R([12])\.fastq\.gz$", name)
    if not base_match:
        raise ValueError(f"FASTQ is not in normalized sample.R1/R2.fastq.gz form: {filename}")
    prefix, read_number = base_match.group(1), base_match.group(2)
    preserve_samples = preserve_samples or set()
    if prefix in preserve_samples:
        return sanitize_sample(prefix), f"R{read_number}"

    match = re.match(r"(.+?)(?:_S\d+)?(?:_L\d{3})?$", prefix)
    if not match:
        raise ValueError(f"FASTQ is not in normalized sample.R1/R2.fastq.gz form: {filename}")
    return sanitize_sample(match.group(1)), f"R{read_number}"


def read_sample_map(path: str | Path | None) -> dict[str, str]:
    """Read an optional TSV with fastq and sample columns."""
    if not path:
        return {}
    map_path = Path(path).expanduser().absolute()
    if not map_path.exists():
        raise FileNotFoundError(f"Sample map TSV does not exist: {path}")

    sample_map: dict[str, str] = {}
    with map_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Sample map TSV is empty: {path}")
        columns = {name.lower(): name for name in reader.fieldnames}
        missing = [name for name in ["fastq", "sample"] if name not in columns]
        if missing:
            raise ValueError(
                f"Sample map TSV must contain columns fastq and sample. Missing: {', '.join(missing)}"
            )
        fastq_col = columns["fastq"]
        sample_col = columns["sample"]
        for line_number, row in enumerate(reader, start=2):
            fastq = (row.get(fastq_col) or "").strip()
            sample = sanitize_sample((row.get(sample_col) or "").strip())
            if not fastq:
                raise ValueError(f"Sample map TSV has an empty fastq value on line {line_number}")
            for key in sample_map_keys(fastq):
                existing = sample_map.get(key)
                if existing and existing != sample:
                    raise ValueError(
                        f"Sample map TSV assigns conflicting sample names to {fastq}: {existing}, {sample}"
                    )
                sample_map[key] = sample
    return sample_map


def read_sample_fastq_inputs(path: str | Path) -> list[str]:
    """Read FASTQ paths from a sample_fastq-style TSV with a fastq column."""
    table_path = Path(path).expanduser().absolute()
    if not table_path.exists():
        raise FileNotFoundError(f"Sample FASTQ TSV does not exist: {path}")

    inputs: list[str] = []
    with table_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Sample FASTQ TSV is empty: {path}")
        columns = {name.lower(): name for name in reader.fieldnames}
        if "fastq" not in columns:
            raise ValueError("Sample FASTQ TSV must contain a fastq column")
        fastq_col = columns["fastq"]
        for line_number, row in enumerate(reader, start=2):
            fastq = (row.get(fastq_col) or "").strip()
            if not fastq:
                raise ValueError(f"Sample FASTQ TSV has an empty fastq value on line {line_number}")
            inputs.append(fastq)

    if not inputs:
        raise ValueError(f"Sample FASTQ TSV contains no FASTQ rows: {path}")
    return inputs


def sample_map_keys(path: str | Path) -> list[str]:
    item = Path(path).expanduser()
    keys = [item.name, str(item)]
    if item.is_absolute():
        keys.append(str(item.absolute()))
        try:
            keys.append(str(item.resolve()))
        except OSError:
            pass
    return list(dict.fromkeys(keys))


def sample_override_for(source: Path, sample_map: dict[str, str]) -> str | None:
    for key in sample_map_keys(source):
        if key in sample_map:
            return sample_map[key]
    return None


def write_concatenated_fastq(sources: list[Path], target: Path) -> None:
    """Concatenate lane-split gzip FASTQs without decompressing them."""
    manifest = target.with_name(f"{target.name}.sources.txt")
    source_list = [str(source.resolve()) for source in sources]
    manifest_text = "\n".join(source_list) + "\n"

    if target.exists() and manifest.exists() and manifest.read_text() == manifest_text:
        return

    tmp = target.with_name(f"{target.name}.tmp")
    if tmp.exists() or tmp.is_symlink():
        tmp.unlink()
    with tmp.open("wb") as output:
        for source in sources:
            with source.open("rb") as handle:
                while True:
                    chunk = handle.read(1024 * 1024 * 8)
                    if not chunk:
                        break
                    output.write(chunk)
    os.replace(tmp, target)
    manifest.write_text(manifest_text)


def stage_fastqs(inputs: list[str], target: str | Path, sample_map: dict[str, str] | None = None) -> list[str]:
    """Create normalized FASTQs, merging lane-split inputs when needed."""
    target = Path(target)
    target.mkdir(parents=True, exist_ok=True)
    sample_map = sample_map or {}
    grouped: dict[tuple[str, str], list[tuple[Path, str]]] = {}

    for item in inputs:
        source = Path(item).expanduser().absolute()
        if not source.exists():
            raise FileNotFoundError(f"Input FASTQ does not exist: {item}")
        normalized_name = normalize_fastq_name(source.name)
        override = sample_override_for(source, sample_map)
        if override:
            _, read = sample_from_normalized_fastq(normalized_name)
            sample = sanitize_sample(override)
            normalized_name = f"{sample}.{read}.fastq.gz"
        else:
            sample, read = sample_from_normalized_fastq(normalized_name)
        grouped.setdefault((sample, read), []).append((source, normalized_name))

    staged: list[str] = []
    for (sample, read), records in sorted(grouped.items()):
        unique_records = list(dict.fromkeys(records))
        if len(unique_records) == 1:
            source, normalized_name = unique_records[0]
            normalized = target / normalized_name
            if normalized.exists() or normalized.is_symlink():
                if normalized.resolve() != source.resolve():
                    raise FileExistsError(
                        f"Staged FASTQ already exists and points elsewhere: {normalized}"
                    )
            else:
                os.symlink(source, normalized)
            staged.append(str(normalized))
            continue

        normalized = target / f"{sample}.{read}.fastq.gz"
        if normalized.is_symlink():
            normalized.unlink()
        elif normalized.exists() and not normalized.is_file():
            raise FileExistsError(f"Staged FASTQ path exists but is not a file: {normalized}")
        sources = [source for source, _ in unique_records]
        write_concatenated_fastq(sources, normalized)
        staged.append(str(normalized))

    return sorted(staged)


def build_sample_table(inputs: list[str], preserve_samples: set[str] | None = None) -> dict[str, dict[str, str]]:
    samples: dict[str, dict[str, str]] = {}
    preserve_samples = preserve_samples or set()
    for item in inputs:
        sample, read = sample_from_normalized_fastq(item, preserve_samples=preserve_samples)
        if read in samples.setdefault(sample, {}):
            raise ValueError(f"Duplicate {read} FASTQ for sample {sample}: {item}")
        samples[sample][read] = str(Path(item).absolute())

    incomplete = {s: reads for s, reads in samples.items() if set(reads) != {"R1", "R2"}}
    if incomplete:
        details = ", ".join(f"{s}: {sorted(reads)}" for s, reads in incomplete.items())
        raise ValueError(f"Each sample must have exactly one R1 and one R2 FASTQ. Incomplete: {details}")
    return dict(sorted(samples.items()))
