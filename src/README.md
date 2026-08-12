# src

This directory implements the `tcr-seek` command-line interface. The runner
parses user inputs, identifies paired FASTQs, creates an isolated output working
directory, writes `config.json`, copies workflow resources into the output
folder, and launches Snakemake.

The FASTQ parser accepts common paired-end names containing `R1`/`R2` or `_1`/`_2`. Users can pass FASTQs directly with `--input`, or pass a cell-seek-style `--sample-fastq` TSV. The `--sample-fastq` table must contain a `fastq` column and normally also contains a `sample` column; when no separate `--sample-map` is provided, the runner uses the same table for sample-name overrides. If multiple input FASTQs resolve to the same sample and read, the runner concatenates the gzip streams into one staged FASTQ and writes a `.sources.txt` manifest beside it so the merge is reproducible.

## Execution modes

The command-line runner supports `--mode local` and `--mode slurm`. Local mode creates the run directory and executes Snakemake directly with `--cores`. Slurm mode creates the same run directory, writes `config.json`, copies `workflow/` and `resources/`, and starts Snakemake with a cluster submission command.

In Slurm mode, Snakemake submits jobs with `sbatch` and uses its native `--use-singularity` support. The host Snakemake jobscript stays in the host environment, while each rule shell command runs inside the resolved SIF image declared by the workflow `container:` directive. The image can be supplied explicitly with `--container-image`, resolved from a shareable SIF directory with `--sif-cache`, or pre-pulled with the `tcr-seek cache` subcommand.

## Container cache resolution

`src/run.py` reads `config/containers.json` to identify the default runtime image.
`--container-image` remains the highest-priority override. When `--sif-cache` is
provided and no explicit image is supplied, the runner looks for the default SIF
in that cache directory and pulls it with Singularity or Apptainer if it is
missing. The resolved SIF path is written to the run-specific `config.json` so
Snakemake sees a concrete local image path.
