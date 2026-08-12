# <code>tcr-seek <b>unlock</b></code>

## 1. About

Snakemake locks a working directory while a workflow is running. If a run is interrupted, the directory may need to be unlocked before resuming.

Only unlock a run after confirming that no controller or Slurm worker jobs are still active for that output directory.

## 2. Synopsis

```text
$ tcr-seek unlock --output OUTPUT
```

## 3. Example

```bash
./tcr-seek/tcr-seek unlock --output /path/to/tcr_output
```
