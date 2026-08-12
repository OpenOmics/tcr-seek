# <code>tcr-seek <b>cache</b></code>

## 1. About

The `cache` command pulls the default `tcr-seek` runtime container into a local
SIF cache directory. This mirrors the OpenOmics pattern used by pipelines such
as `cell-seek`: the cache directory can be shared on a filesystem and reused by
future `tcr-seek run --sif-cache` commands.

The default container is recorded in `config/containers.json` in the pipeline
source. The current default is:

```text
tcrseq-rnaseq_1.3.sif
docker://pauls85/tcrseq-rnaseq:1.3
```

## 2. Synopsis

```text
$ tcr-seek cache --sif-cache SIF_CACHE [--dry-run]
```

## 3. Required Arguments

`--sif-cache SIF_CACHE`  
> Directory where the default tcr-seek SIF is stored. If the SIF is already
> present, the command reports the cached path and does not pull again. If it is
> missing, the command pulls the default container URI into this directory.

## 4. Options

`--dry-run`  
> Show the pull that would be performed without downloading the image.

## 5. Example

```bash
module load singularity

tcr-seek cache \
  --sif-cache /data/NIAMS_SPA/processed/BulkTCR/containers
```

After the cache exists, pass the same directory to `run`:

```bash
tcr-seek run ... \
  --mode slurm \
  --sif-cache /data/NIAMS_SPA/processed/BulkTCR/containers
```

Use `--container-image /path/to/image.sif` with `tcr-seek run` when an exact SIF
path should override the default cache resolver.
