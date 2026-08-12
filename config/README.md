# config

This directory contains reusable configuration files and local placeholders used
by `tcr-seek`. Run-specific `config.json` files are written into each selected
output directory from command-line arguments.

## Method implemented in this directory

`mi.license` is the project-local MiXCR license file placeholder. Paste the
MiXCR license text directly into `mi.license` when MiXCR-based rules are added
or when you want `tcr-seek` to pass the license location to MiXCR through the
`MI_LICENSE_FILE` environment variable. Leave the file empty if MiXCR is not
being used.

Keep `mi.license` private because it contains license text. You can also pass a
different license path at runtime with `tcr-seek run --mixcr-license-file`.

## Sample map template

`sample_map.tsv.example` shows the optional TSV format accepted by
`tcr-seek run --sample-map`. The file must contain `fastq` and `sample` columns.
The `fastq` value can be a full path or a FASTQ basename. FASTQs listed in the
map use the provided sample name; unlisted FASTQs use automatic parsing.

## Container manifest

`containers.json` records the default tcr-seek runtime container. The manifest
contains the logical container name, version, expected SIF filename, and remote
container URI used by `tcr-seek cache` and by `tcr-seek run --sif-cache` when the
SIF is not already present. Keeping this small manifest in the pipeline source
lets run directories record both the resolved SIF path and the default manifest
that selected it.

The current default SIF is `tcrseq-rnaseq_1.3.sif`. Users can bypass the default
resolver with `tcr-seek run --container-image /path/to/image.sif` when they need
to pin a different image explicitly.
