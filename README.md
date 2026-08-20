<div align="center">
   
  <h1>tcr-seek 🔬</h1>
  
  **_A reproducible bulk TCR-sequencing pipeline_**

  [![Docs](https://github.com/OpenOmics/tcr-seek/actions/workflows/docs.yml/badge.svg)](https://github.com/OpenOmics/tcr-seek/actions/workflows/docs.yml)
  [![Tests](https://github.com/OpenOmics/tcr-seek/actions/workflows/tests.yml/badge.svg)](https://github.com/OpenOmics/tcr-seek/actions/workflows/tests.yml)
  [![Issues](https://img.shields.io/github/issues/OpenOmics/tcr-seek)](https://github.com/OpenOmics/tcr-seek/issues)
  [![License](https://img.shields.io/github/license/OpenOmics/tcr-seek)](https://github.com/OpenOmics/tcr-seek/blob/main/LICENSE)</br>
  [![Docker Pulls](https://img.shields.io/docker/pulls/pauls85/tcrseq-rnaseq)](https://hub.docker.com/r/pauls85/tcrseq-rnaseq)
  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21908454.svg)](https://doi.org/10.5281/zenodo.21908454)
  
  <i>
    This is the home of the pipeline, tcr-seek. Its long-term goals: to perform bulk TCR-seq processing, repertoire analysis, and integrated quality control!
  </i>
</div>

## Overview

Welcome to tcr-seek! Before getting started, we highly recommend reading through [tcr-seek's documentation](https://openomics.github.io/tcr-seek/).

The **`./tcr-seek`** pipeline is composed several inter-related sub commands to setup and run the pipeline across different systems. Each of the available sub commands perform different functions: 

 * [<code>tcr-seek <b>run</b></code>](https://openomics.github.io/tcr-seek/usage/run/): Run the tcr-seek pipeline with your input files.
 * [<code>tcr-seek <b>unlock</b></code>](https://openomics.github.io/tcr-seek/usage/unlock/): Unlocks a previous runs output directory.
 * [<code>tcr-seek <b>cache</b></code>](https://openomics.github.io/tcr-seek/usage/cache/): Cache software containers locally.

**tcr-seek** is a reproducible bulk T-cell receptor sequencing (TCR-seq) pipeline for processing paired-end FASTQ files, generating immune repertoire assignments, and performing repertoire-level analysis and quality control. The workflow stages and normalizes input FASTQ files, evaluates raw sequencing quality with fastp and optional FastQC, and processes reads with pRESTO using quality filtering, primer masking, read pairing, consensus generation, pair assembly, C-region masking, duplicate collapsing, and consensus-count filtering.

Processed sequences are then analyzed with Change-O/IgBLAST to generate AIRR-formatted repertoire tables. These tables are imported into immunarch to create a reusable repertoire object for downstream diversity analysis. The pipeline computes repertoire diversity and rarefaction metrics and produces both repertoire-focused Quarto reports and integrated MultiQC reports that summarize sequencing and pipeline quality-control results.

The pipeline supports both local execution and Slurm-based cluster execution, with containerized worker jobs for reproducible analysis at project scale. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

As input, it accepts a set of paired-end FastQ files and can be run locally on a compute instance or on-premise using a cluster. A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM (more coming soon!). A hybrid approach ensures the pipeline is accessible to all users. Before getting started, we highly recommend reading through the [usage](https://openomics.github.io/tcr-seek/usage/run/) section of each available sub command.

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](https://openomics.github.io/tcr-seek/faq/questions/) prior to [opening an issue on Github](https://github.com/OpenOmics/tcr-seek/issues).

## Dependencies

**Requires:** `singularity>=3.5`  `snakemake<=7.32.4`

At the current moment, the pipeline will only uses docker images. With that being said, [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) and [singularity](https://singularity.lbl.gov/all-releases) must be installed on the target system. Snakemake orchestrates the execution of each step in the pipeline. To guarantee the highest level of reproducibility, each step of the pipeline will rely on versioned images from [DockerHub](https://hub.docker.com). Snakemake uses singularity to pull these images onto the local filesystem prior to job execution, and as so, snakemake and singularity will be the only two dependencies in the future.

## Installation

Please clone this repository to your local filesystem using the following command:
```bash
# Clone Repository from Github
git clone https://github.com/OpenOmics/tcr-seek.git
# Change your working directory
cd tcr-seek/
# Add dependencies to $PATH
# Biowulf users should run
module load snakemake singularity
# Get usage information
./tcr-seek -h
```

## Contribute

This site is a living document, created for and by members like you. tcr-seek is maintained by the members of [OpenOmics](https://openomics.github.io) and is improved by continous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [GitHub repository](https://github.com/OpenOmics/tcr-seek).


## Cite

If you use this software, please cite it as below:  

<details>
  <summary><b><i>@BibText</i></b></summary>
 
```text
@software{subrata_paul_2026_21908454,
  author       = {Paul, Subrata and Schaughency, Paul},
  title        = {OpenOmics/tcr-seek},
  month        = aug,
  year         = 2026,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.21908453},
  url          = {https://doi.org/10.5281/zenodo.21908453},
}
```

</details>

<details>
  <summary><b><i>@APA</i></b></summary>

```text
Paul, S., & Schaughency, P. (2026). OpenOmics/tcr-seek: tcr-seek [Computer software]. Zenodo. 10.5281/zenodo.21908453
```

</details>

## References

<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  