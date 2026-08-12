# resources

This directory is reserved for small reusable resources. Project-specific primer
and reference FASTA files are not stored here; pass them with `--refdir` so the
same pipeline can be reused across projects.

## Slurm execution resources

Slurm mode follows the `cell-seek` style. The CLI prepares the output directory and submits a lightweight controller job with `sbatch`. That controller runs Snakemake from the copied workflow and uses Snakemake's native `--cluster` and `--use-singularity` options to submit one worker job per runnable rule instance.

The generated controller script is written into each output directory as `run_tcr_seek_controller.sh`. Controller logs are written under `logfiles/`, and worker logs are written under `logfiles/slurm/`.

`slurm-job-wrapper.sh` is retained as an experimental helper from the first validation attempt, but it is not used by `tcr-seek run --mode slurm`.
