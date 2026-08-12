# workflow scripts

This directory contains helper scripts called by Snakemake rules.

`create_immunarch_object.R` gathers Change-O `*_db-pass.tsv` files, writes the
`ImmunarchInput` directory expected by immunarch, optionally filters user
metadata to sequenced samples, runs `immunarch::repLoad()`, and saves
`Immdata.rds`.

`calculate_diversity_per_sample.R` reads `Immdata.rds`, calculates diversity
metrics independently per sample using the `strict` preset by default, and
writes observed richness, Shannon diversity, vegan Simpson diversity `(1 - D)`,
Pielou evenness, library size, effective clonotypes, clone dominance metrics,
and, when `run_rarefaction` is true, scalar rarefied richness. It writes
`DiversityMethod.txt`, `diversity_by_clonotype_definition.rds`, and the Quarto
HTML QC report. The QC variables section creates one tab for each requested
metadata variable and nested tabs for each available clonotype definition. QC
variable names are resolved case-insensitively, and the report derives a
`fastq.delivery` grouping variable from `sample_fastq.tsv` when available,
falling back to `input_fastqs/*.sources.txt` sidecars. Delivery labels are the
parent directory names of the original FASTQ paths; the current NIAMS delivery
`20260710_SH01446_0005_ASC2252557-SC3_From412` is reported as `batch5`. The
primary `diversity_summary.rds` uses the selected clonotype definition, while
`diversity_by_clonotype_definition.rds` recalculates the same metrics for
`gene`, `nt`, `aa`, and `strict` so the report can compare definitions side by
side. When `run_rarefaction` is false, rarefaction columns are retained with
`NA` values and the no-rarefaction QC report template is rendered. The
no-rarefaction template also summarizes FASTQ delivery batches by matching
`input_fastqs/*.R1.fastq.gz.sources.txt` and
`input_fastqs/*.R2.fastq.gz.sources.txt` sidecar files with escaped R regular
expressions.

`calculate_diversity_vegan.R` is retained as a backup implementation that builds the full sample-by-clonotype matrix and can produce rarefaction curves, but it is memory intensive for large cohorts.

`write_methods_markdown.R` records a short publication-ready methods paragraph
for a completed or in-progress run. It embeds software versions directly in the
paragraph and can include both the Singularity/SIF path and the Docker or OCI
source image used to build it. If `--container-source` is omitted, the script
tries to read the source image from `singularity inspect` metadata. Run it
inside the same container used for the analysis so command-line tools and R
package versions match the workflow environment:

```bash
Rscript workflow/scripts/write_methods_markdown.R \
  --output 06_immunarch/tcr_seek_methods.md \
  --pipeline-version "$(cat /data/NIAMS_IDSS/projects/NIAMS-47/tcr-seek/VERSION)" \
  --container-image /data/NIAMS_IDSS/projects/NIAMS-47/containers/tcrseq-rnaseq_latest.sif \
  --container-source docker://pauls85/tcrseq-rnaseq:latest \
  --clonotype-definition strict
```

When `--run-rarefaction true` is used, `calculate_diversity_per_sample.R` also writes per-sample rarefaction curves for the primary clonotype definition so `qc_report.qmd` can render the rarefaction curve tabs without building a global sample-by-clonotype matrix. The rarefaction section includes QC-variable tabs; for each selected QC variable, sample curves are colored by group and overlaid with a thicker loess-smoothed curve fit to the group mean at each observed clone/read abundance depth.

`write_presto_qc_summary.py` runs at the end of each pRESTO sample job and writes
`01_run/<sample>/<sample>.presto_qc.tsv`. The table captures read counts across
input, Q20 filtering, pairing, consensus, assembly, C-region detection, final
`CONSCOUNT >= 2` sequences, final abundance, and derived pass fractions.

`write_multiqc_input.py` aggregates the per-sample pRESTO QC rows with
Change-O/IgBLAST `*_db-pass.tsv` summaries and writes
`04_multiqc/tcr_seek_multiqc.tsv` for the tcr-seek MultiQC plugin. The aggregate
includes db-pass counts, Change-O pass fraction, productive fraction, in-frame
fraction, stop-codon fraction, consensus-count summaries, and top V/J/locus
calls.

`multiqc_config.yaml` provides the base MultiQC search pattern for
`tcr_seek_multiqc.tsv`. The `multiqc_report` rule writes a runtime config that
orders modules as FastQC, tcr-seek, fastp when FastQC reports are present, or
fastp, tcr-seek when only fastp reports are present. The tcr-seek plugin adds targeted bulk TCR-seq interpretation
notes for raw FASTQ duplication and GC deviation so those upstream QC values
are reviewed as cohort-comparison metrics rather than generic pass/fail calls.

