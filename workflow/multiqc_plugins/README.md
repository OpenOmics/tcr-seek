# MultiQC plugins

This directory contains run-local MultiQC plugin packages used by the Snakemake
workflow.

## Method implemented in this directory

The `tcr_seek` plugin is packaged as a standard Python distribution with a
MultiQC module entry point. The workflow installs it into `04_multiqc/plugin_env`
with `pip --target` immediately before running `multiqc`, so the plugin does not
need to be installed globally in the container or user environment.
