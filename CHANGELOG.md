# Changelog

All notable changes to this project will be documented in this file.

This changelog is automatically updated by [release-please](https://github.com/googleapis/release-please) when contributors follow [conventional commit](https://www.conventionalcommits.org/) git messages. If you are using conventional commit messages, you should never need to edit this file manually. A github action will automatically update this file when a new release is created.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-20

### Added

- Initial release of `tcr-seek`, an end-to-end bulk TCR-seq analysis pipeline.
- Paired-end FASTQ quality control and pRESTO-based read processing, including quality filtering, primer processing, consensus generation, read assembly, and duplicate collapsing.
- Change-O/IgBLAST TCR gene assignment with AIRR-formatted output.
- Immunarch repertoire object generation and integration of sample metadata.
- Repertoire diversity analysis, including richness, Shannon diversity, Simpson diversity, evenness, and optional rarefaction.
- Quarto and MultiQC reports for sequencing, processing, repertoire, and diversity quality control.
- Support for local and Slurm-based workflow execution.
