# workflow

This directory contains the Snakemake workflow. `Snakefile` loads the generated
`config.json`, defines the sample list, includes the rule files, and sets the
final targets: raw FastQC reports, processed TCR FASTQs, Change-O AIRR tables, the immunarch
object, vegan diversity summaries, optional rarefaction summaries, the
immunarch/Quarto QC report, and the separate MultiQC report.

The workflow is copied into each output directory before execution so every run
keeps a record of the pipeline logic used for that run.

## Method implemented in this directory

The workflow keeps the existing immunarch Quarto QC report as an independent
`03_immunarch/qc_report.html` target. It also runs raw R1/R2 QC into `00_fastqc/<sample>/`, always using paired-end fastp and also using FastQC when available, then builds
`04_multiqc/multiqc_report.html`
from the raw-QC outputs and `04_multiqc/tcr_seek_multiqc.tsv` using the
run-local tcr-seek MultiQC plugin under `workflow/multiqc_plugins/tcr_seek`.
The MultiQC rule writes a runtime config so reports show FastQC, tcr-seek, and
fastp when FastQC exists, or fastp and tcr-seek when FastQC is absent. When
`params.cleanup_intermediates` is true in `config.json`, `Snakefile` adds the
`cleanup_intermediates` rule as a final target after the main reports complete;
cleanup gzips retained pRESTO handoff FASTQs and removes reproducible large
intermediates.
