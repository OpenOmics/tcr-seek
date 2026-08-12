# workflow rules

This directory contains modular Snakemake rules.

- `tcr.smk` runs upstream raw FASTQ QC on each sample, always using paired-end fastp and additionally using FastQC when available. It also implements the pRESTO processing and Change-O/IgBLAST steps, and writes one per-sample pRESTO QC TSV for MultiQC.
- `immunarch.smk` creates immunarch input files, saves the imported immunarch
  object, calculates vegan diversity metrics, optionally calculates rarefaction metrics, renders the matching Quarto QC report, aggregates tcr-seek metrics for MultiQC, and runs the tcr-seek MultiQC plugin.

The rules can run locally through Snakemake or through `tcr-seek run --mode slurm`. In Slurm mode, the per-sample `presto_process` and `changeo_igblast` jobs are independently submitted to the scheduler, while `immunarch_object` runs after all sample-level Change-O outputs are present, followed by `repertoire_diversity`. The `run_rarefaction` config flag selects the report template with or without rarefaction curves.

Default Slurm resources are declared directly in the rules:

- `raw_fastq_qc`: 2 threads, 4 GB memory, 1 hour walltime, `norm` partition; writes fastp reports for every sample and FastQC reports when `fastqc` is available.
- `presto_process`: 2 configurable threads by default, 32 GB memory, 1 day walltime, `norm` partition.
- `changeo_igblast`: 2 configurable threads by default, 16 GB memory, 8 hours walltime, `norm` partition.
- `immunarch_object`: 1 thread, 100 GB memory, 3 days walltime, `norm` partition.
- `repertoire_diversity`: 1 thread, 100 GB memory, 4 hours walltime, `norm` partition.
- `qc_report`: 1 thread, 16 GB memory, 4 hours walltime, `norm` partition.
- `multiqc_input`: 1 thread, 4 GB memory, 1 hour walltime, `norm` partition.
- `multiqc_report`: 1 thread, 8 GB memory, 1 hour walltime, `norm` partition.
