# Slurm Mode

## 1. About

Slurm mode lets Snakemake submit independent rule jobs to the scheduler. For
bulk TCR-seq, each sample can run pRESTO independently, then each sample can run
Change-O/IgBLAST independently, followed by shared immunarch, diversity, QC
report, and MultiQC jobs.

`tcr-seek run --mode slurm` prepares the run directory, writes `config.json`,
and submits a lightweight controller job with `sbatch`. The controller job runs
Snakemake on the host. Snakemake submits each rule job with `sbatch` and uses
native `--use-singularity` support so the rule commands run inside the resolved
SIF image.

Use `--no-submit-controller` only when you want to keep the Snakemake controller
in the foreground for debugging. Dry-runs print the Snakemake plan to the
terminal and also write it to `<output>/logfiles/dry-run.log`; add `--quiet` to
write only the log file.

## 2. Container Options

Use one of these approaches:

- `--container-image /path/to/image.sif`: pin an exact SIF path.
- `--sif-cache /path/to/SIFs`: use the default SIF from a shareable cache
  directory, pulling it from the default container URI if it is missing.

The default container is defined in `config/containers.json`. The current default
is `tcrseq-rnaseq_1.3.sif` from `docker://pauls85/tcrseq-rnaseq:1.3`.

Pre-cache the image when desired:

```bash
module load singularity

./tcr-seek/tcr-seek cache \
  --sif-cache /data/NIAMS_SPA/processed/BulkTCR/containers
```

## 3. Example

```bash
module purge
module load snakemake singularity

./tcr-seek/tcr-seek run \
  --sample-fastq /path/to/sample_fastq.tsv \
  --output /path/to/tcr_output \
  --refdir /path/to/01_ref \
  --sample-info /path/to/sample_info.tsv \
  --threads 2 \
  --mode slurm \
  --sif-cache /data/NIAMS_SPA/processed/BulkTCR/containers \
  --run-rarefaction false \
  --cleanup-intermediates false \
  --jobs 100
```

## 4. Resource Model

Default resources are declared in the Snakemake rules:

| Rule | Parallelism | Memory | Walltime | Partition |
|---|---:|---:|---:|---|
| `raw_fastq_qc` | one job per sample | 4 GB | 1 hour | `norm` |
| `presto_process` | one job per sample | 32 GB | 1 day | `norm` |
| `changeo_igblast` | one job per sample | 16 GB | 8 hours | `norm` |
| `immunarch_object` | one final job | 100 GB | 3 days | `norm` |
| `repertoire_diversity` | one final job | 100 GB | 4 hours | `norm` |
| `qc_report` | one final job | 16 GB | 4 hours | `norm` |
| `multiqc_input` | one final job | 4 GB | 1 hour | `norm` |
| `multiqc_report` | one final job | 8 GB | 1 hour | `norm` |

`--threads` controls CPUs requested for per-sample pRESTO and Change-O jobs.
`--jobs` controls how many Slurm jobs Snakemake may have active at once.

## 5. Logs

Controller logs are written under:

```text
<output>/logfiles/controller-<jobid>.out
```

Worker job logs are written under:

```text
<output>/logfiles/slurm/
```

The copied workflow and resources are stored inside the run directory, so each
output folder keeps the exact workflow logic used for that run.
