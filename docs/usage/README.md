# usage

This directory contains command-oriented documentation for `tcr-seek`. Pages
describe the user-facing CLI, required inputs, container cache management,
execution modes, Slurm operation, and operational recovery commands.

## Method implemented in this directory

The usage pages are source Markdown files for the MkDocs site. They document the
current command behavior implemented in `src/run.py`, including `tcr-seek run`,
`tcr-seek cache`, `tcr-seek unlock`, Slurm controller submission, SIF cache
resolution, rarefaction toggling, cleanup toggling, and expected output
locations.
