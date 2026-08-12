#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
slurm-job-wrapper.sh submits one Snakemake rule job to Slurm and runs the job body inside a Singularity/SIF container.

Required arguments:
  --container PATH     SIF image used for the worker job
  --bind PATHS         Singularity bind paths, comma-separated
  --logdir PATH        Directory for Slurm stdout/stderr logs
  --partition NAME     Slurm partition
  --time TIME          Slurm runtime, for example 0-08:00:00
  --mem MB             Memory in MB
  --cpus N             CPUs per task
  --job-name NAME      Slurm job name
  JOBSCRIPT            Snakemake-generated jobscript, supplied by Snakemake
EOF
}

container=""
bind_paths=""
logdir=""
partition="norm"
runtime="0-08:00:00"
mem_mb="16000"
cpus="1"
job_name="tcr-seek"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --container) container="$2"; shift 2 ;;
    --bind) bind_paths="$2"; shift 2 ;;
    --logdir) logdir="$2"; shift 2 ;;
    --partition) partition="$2"; shift 2 ;;
    --time) runtime="$2"; shift 2 ;;
    --mem) mem_mb="$2"; shift 2 ;;
    --cpus) cpus="$2"; shift 2 ;;
    --job-name) job_name="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    --*) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
    *) jobscript="$1"; shift ;;
  esac
done

if [[ -z "${container}" || -z "${bind_paths}" || -z "${logdir}" || -z "${jobscript:-}" ]]; then
  usage >&2
  exit 2
fi

mkdir -p "${logdir}"
if [[ ! -s "${container}" ]]; then
  echo "Container image does not exist or is empty: ${container}" >&2
  exit 2
fi

# Slurm receives this wrapper command from Snakemake. The actual Snakemake
# jobscript is executed inside the container so all pRESTO, Change-O, IgBLAST,
# and R package versions come from the validated SIF.
sbatch --parsable   --partition="${partition}"   --time="${runtime}"   --mem="${mem_mb}"   --cpus-per-task="${cpus}"   --job-name="${job_name}"   --output="${logdir}/slurm-%j_${job_name}.out"   --error="${logdir}/slurm-%j_${job_name}.out"   --wrap="singularity exec --bind ${bind_paths} ${container} bash ${jobscript}"
