# FASTQ Naming

By default, `tcr-seek` follows the `cell-seek` style of staging FASTQs as
normalized files named `sample.R1.fastq.gz` and `sample.R2.fastq.gz` under the
run directory `input_fastqs/`.

Automatic parsing strips common Illumina read suffixes and sample-index tokens.
For manual control, provide a TSV with `fastq` and `sample` columns:

```text
fastq	sample
294_S29_R1_001.fastq.gz	294_S29
294_S29_R2_001.fastq.gz	294_S29
```

Mapped FASTQs use the supplied sample name. Unmapped FASTQs continue to use
automatic parsing.

## Multiple FASTQs Per Sample

When more than one input FASTQ maps to the same sample and read, the runner
concatenates the gzip streams without decompressing them. The merged output is
written as one staged FASTQ, and a `.sources.txt` file beside it records the
original source files in order.

For example, if four lane files map to `SampleA` read 1 and four lane files map
to `SampleA` read 2, the staged inputs become:

```text
input_fastqs/SampleA.R1.fastq.gz
input_fastqs/SampleA.R1.fastq.gz.sources.txt
input_fastqs/SampleA.R2.fastq.gz
input_fastqs/SampleA.R2.fastq.gz.sources.txt
```

When only one FASTQ maps to a sample/read, the staged FASTQ is a symlink to the
original input. Source FASTQs are never edited in place.
