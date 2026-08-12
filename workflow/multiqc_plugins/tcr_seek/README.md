# tcr-seek MultiQC plugin

This directory contains the Python package for the `tcr_seek` MultiQC module.

## Method implemented in this directory

The plugin searches for `tcr_seek_multiqc.tsv`, parses one row per sample, adds
selected metrics to the MultiQC general statistics table, and renders tcr-seek
summary sections for read processing, Change-O/IgBLAST pass rates, productivity,
and a downloadable sample-level QC table.
