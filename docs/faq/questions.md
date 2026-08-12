# General Questions

## Should I use local or Slurm mode?

Use local mode for dry-runs, testing, and tiny analyses. Use Slurm mode for real
projects so per-sample jobs can run in parallel on the HPC system.

## How should I provide the container?

For Slurm mode, either pass an exact SIF with `--container-image` or pass a
shareable cache directory with `--sif-cache`. `--container-image` is the most
explicit option and overrides the cache. `--sif-cache` uses the default SIF
recorded in `config/containers.json` and pulls it if the file is missing.

Pre-populate the cache with:

```bash
tcr-seek cache --sif-cache /path/to/SIFs
```

## Can I skip rarefaction and run it later?

Yes. Run the full pipeline with `--run-rarefaction false` to skip rarefied
richness and rarefaction curves. After the run completes, rarefaction can be run
from the existing `03_immunarch/Immdata.rds` by enabling rarefaction in the run
configuration and forcing the `repertoire_diversity` and `qc_report` Snakemake
rules.

## Does `tcr-seek` merge lane-split FASTQs?

Yes. If multiple FASTQs map to the same sample and read through `--sample-fastq`
or `--sample-map`, `tcr-seek` writes one concatenated staged FASTQ and a
`.sources.txt` manifest. Single FASTQs are staged as symlinks.

## Where are the reports?

The immunarch/Quarto QC report is written to:

```text
<output>/03_immunarch/qc_report.html
```

The aggregate MultiQC report is written to:

```text
<output>/04_multiqc/multiqc_report.html
```

## Does `tcr-seek` run RNA-seq or HLA analysis?

Not currently. The present workflow is focused on bulk TCR processing, diversity
metrics, immunarch QC reporting, and aggregate MultiQC reporting. RNA-seq and
HLA stages can be added later as separate workflow modules.

## Where do worker logs go in Slurm mode?

Slurm worker logs are written to `<output>/logfiles/slurm/`. The controller log
is written under `<output>/logfiles/controller-<jobid>.out`.
