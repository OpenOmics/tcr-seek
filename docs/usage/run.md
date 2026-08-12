# <code>tcr-seek <b>run</b></code>

## 1. About

The `run` command stages paired FASTQ files, writes a reproducible `config.json`,
copies the workflow into the output directory, and starts Snakemake. It accepts
FASTQs directly with `--input` or a cell-seek-style table through
`--sample-fastq`. When multiple FASTQs map to the same sample and read, the
runner concatenates the gzip streams into one staged FASTQ and writes a
`.sources.txt` manifest beside it.

## 2. Synopsis

```text
$ tcr-seek run [--help] \
    [--input INPUT [INPUT ...] | --sample-fastq SAMPLE_FASTQ] \
    --output OUTPUT \
    --refdir REFDIR \
    [--metadata METADATA] [--sample-info SAMPLE_INFO] [--sample-map SAMPLE_MAP] \
    [--mixcr-license-file MIXCR_LICENSE_FILE] \
    [--species SPECIES] [--receptor RECEPTOR] \
    [--clonotype-definition CLONOTYPE_DEFINITION] \
    [--run-rarefaction {true,false}] [--cleanup-intermediates {true,false}] \
    [--qc-variables QC_VARIABLES] \
    [--threads THREADS] [--cores CORES] \
    [--mode {local,slurm}] \
    [--container-image CONTAINER_IMAGE] [--sif-cache SIF_CACHE] \
    [--jobs JOBS] [--local-cores LOCAL_CORES] [--job-name JOB_NAME] \
    [--controller-partition CONTROLLER_PARTITION] \
    [--controller-time CONTROLLER_TIME] \
    [--controller-mem-mb CONTROLLER_MEM_MB] \
    [--no-submit-controller] [--quiet] \
    [--latency-wait LATENCY_WAIT] [--restart-times RESTART_TIMES] \
    [--singularity-bind SINGULARITY_BIND] [--dry-run]
```

## 3. Required Arguments

`--input INPUT [INPUT ...]` or `--sample-fastq SAMPLE_FASTQ`  
> Provide paired gzipped FASTQ files directly, or provide a TSV with a `fastq`
> column and a `sample` column. With `--sample-fastq`, the runner reads FASTQ
> paths from the `fastq` column and uses the same file as the sample map unless
> `--sample-map` is also supplied.

`--output OUTPUT`  
> Output working directory. `tcr-seek` creates this directory if it does not
> exist and writes `input_fastqs/`, `config.json`, copied workflow files, logs,
> and analysis outputs under it.

`--refdir REFDIR`  
> Reference directory containing `IS_Human_R1_Primers.fasta`,
> `IS_Human_R2_Primers.fasta`, `IS_Human_C-Region.fasta`, and
> `Immune_Ref.fasta`.

## 4. Common Options

`--sample-map SAMPLE_MAP`  
> Optional TSV with `fastq` and `sample` columns. Use this when FASTQ basenames
> should map to custom sample names.

`--sample-info SAMPLE_INFO`  
> Optional tab-delimited sample information file. It must contain `Sample`,
> `SampleID`, `sample`, or `sample_id`. Columns are merged into immunarch
> metadata and can be used by the QC report.

`--qc-variables QC_VARIABLES`  
> Comma-separated sample metadata variables to use as QC report grouping tabs,
> for example `Condition,Sex,B27,B38,Batch,Diagnosis`. If omitted, all
> non-sample metadata columns are considered.

`--clonotype-definition CLONOTYPE_DEFINITION`  
> Clonotype definition preset or comma-separated immunarch columns for diversity
> metrics. Presets are `gene`, `nt`, `aa`, and `strict`. Default: `strict`.

`--run-rarefaction {true,false}`  
> Whether to calculate rarefied richness and rarefaction curves for the QC
> report. Default: `true`. Use `false` for large runs when rarefaction should be
> skipped and run later.

`--cleanup-intermediates {true,false}`  
> After successful final outputs, gzip retained pRESTO handoff FASTQs and remove
> reproducible intermediates. Default: `false`.

`--threads THREADS`  
> Threads per sample-level pRESTO and Change-O job. Default: `2`.

`--mode {local,slurm}`  
> Execution mode. Use `local` for testing and `slurm` for full HPC runs. In
> Slurm mode, real runs submit a lightweight controller job and return after
> printing the Slurm job id.

`--container-image CONTAINER_IMAGE`  
> Explicit SIF path for Slurm worker jobs. This overrides `--sif-cache`.

`--sif-cache SIF_CACHE`  
> Directory containing or receiving the default tcr-seek SIF. When no explicit
> `--container-image` is supplied, the runner looks for the default SIF recorded
> in `config/containers.json` and pulls it if missing.

`--dry-run`  
> Build the run directory and show the Snakemake DAG without executing jobs.
> Slurm dry-runs write the full plan to `<output>/logfiles/dry-run.log`.

## 5. Local Example

Local mode is usually run from inside the SIF.

```bash
singularity exec \
  --bind /data:/data \
  --bind /vf:/vf \
  --bind /tmp:/tmp \
  /data/NIAMS_SPA/processed/BulkTCR/containers/tcrseq-rnaseq_1.3.sif \
  ./tcr-seek/tcr-seek run \
    --sample-fastq /path/to/sample_fastq.tsv \
    --output /path/to/tcr_output \
    --refdir /path/to/01_ref \
    --sample-info /path/to/sample_info.tsv \
    --qc-variables Condition,Sex,B27,B38,Batch,Diagnosis \
    --threads 2 \
    --cores 8
```

## 6. Slurm Example With SIF Cache

```bash
module load snakemake singularity

./tcr-seek/tcr-seek run \
  --sample-fastq /path/to/sample_fastq.tsv \
  --output /path/to/tcr_output \
  --refdir /path/to/01_ref \
  --sample-info /path/to/sample_info.tsv \
  --qc-variables Condition,Sex,B27,B38,Batch,Diagnosis \
  --threads 2 \
  --mode slurm \
  --sif-cache /data/NIAMS_SPA/processed/BulkTCR/containers \
  --run-rarefaction false \
  --cleanup-intermediates false \
  --jobs 100 \
  --dry-run
```

Remove `--dry-run` to submit the controller job. Use
`tcr-seek cache --sif-cache /path/to/SIFs` to pre-populate the SIF cache before
running on systems with limited external network access.

## 7. Main Outputs

Outputs are written under the selected output directory:

- `input_fastqs/`: normalized staged FASTQs and `.sources.txt` manifests for merged inputs.
- `00_fastqc/`: upstream raw FASTQ QC from fastp and, when available, FastQC.
- `01_run/`: per-sample pRESTO outputs.
- `02_changeo/`: per-sample Change-O/IgBLAST AIRR tables.
- `03_immunarch/`: `Immdata.rds`, diversity outputs, and `qc_report.html`.
- `04_multiqc/`: `multiqc_report.html` and `tcr_seek_multiqc.tsv`.
