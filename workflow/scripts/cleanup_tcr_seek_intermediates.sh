#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage:
  cleanup_tcr_seek_intermediates.sh \
    --output-dir DIR \
    --samples "S1 S2 ..." \
    --run-dir 01_run \
    --changeo-dir 02_changeo \
    --immunarch-dir 03_immunarch \
    --marker cleanup/tcr_seek_cleanup.done

The script is intended as the final Snakemake target after Change-O and
immunarch outputs have completed. It keeps final AIRR tables, immunarch outputs,
small per-sample HTML reports, and gzipped Presto handoff FASTQs:

  <run-dir>/<sample>/<sample>.uniqC_atleast-2.fastq.gz

It removes reproducible large intermediates from <run-dir> and Change-O debug
outputs that can be regenerated from raw FASTQs and the recorded workflow.
USAGE
}

output_dir=""
samples=""
run_dir="01_run"
changeo_dir="02_changeo"
immunarch_dir="03_immunarch"
marker=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-dir)
            output_dir="$2"
            shift 2
            ;;
        --samples)
            samples="$2"
            shift 2
            ;;
        --run-dir)
            run_dir="$2"
            shift 2
            ;;
        --changeo-dir)
            changeo_dir="$2"
            shift 2
            ;;
        --immunarch-dir)
            immunarch_dir="$2"
            shift 2
            ;;
        --marker)
            marker="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ -z "$output_dir" || -z "$samples" || -z "$marker" ]]; then
    usage >&2
    exit 2
fi

cd "$output_dir"
mkdir -p "$(dirname "$marker")"
summary="$(dirname "$marker")/tcr_seek_cleanup_summary.tsv"
tmp_summary="${summary}.tmp"

require_file() {
    local file="$1"
    if [[ ! -s "$file" ]]; then
        echo "Required final output is missing or empty: $file" >&2
        exit 1
    fi
}

size_or_zero() {
    local file="$1"
    if [[ -e "$file" ]]; then
        stat -c '%s' "$file"
    else
        printf '0'
    fi
}

record_remove() {
    local file="$1"
    if [[ -f "$file" ]]; then
        printf 'remove\t%s\t%s\n' "$(size_or_zero "$file")" "$file" >> "$tmp_summary"
        rm -f -- "$file"
    fi
}

record_gzip() {
    local file="$1"
    if [[ -f "$file" ]]; then
        printf 'gzip\t%s\t%s\n' "$(size_or_zero "$file")" "$file" >> "$tmp_summary"
        gzip -f -- "$file"
    fi
}

require_file "${immunarch_dir}/Immdata.rds"
require_file "${immunarch_dir}/ImmunarchInput/metadata.txt"
require_file "${immunarch_dir}/diversity_summary.rds"
require_file "${immunarch_dir}/rarefaction.rds"
require_file "${immunarch_dir}/rarefaction_curve.rds"
require_file "${immunarch_dir}/qc_report.html"

for sample in $samples; do
    require_file "${changeo_dir}/${sample}/${sample}.al2_db-pass.tsv"

    handoff="${run_dir}/${sample}/${sample}.uniqC_atleast-2.fastq"
    handoff_gz="${handoff}.gz"
    if [[ ! -f "$handoff" && ! -s "$handoff_gz" ]]; then
        echo "Expected Presto handoff FASTQ or gzip is missing: $handoff" >&2
        exit 1
    fi
done

: > "$tmp_summary"
printf 'action\tbytes_before\tpath\n' >> "$tmp_summary"

for sample in $samples; do
    handoff="${run_dir}/${sample}/${sample}.uniqC_atleast-2.fastq"
    handoff_gz="${handoff}.gz"
    if [[ -f "$handoff" ]]; then
        record_gzip "$handoff"
    fi

    if [[ -d "${run_dir}/${sample}" ]]; then
        while IFS= read -r -d '' file; do
            case "$file" in
                "$handoff"|"$handoff_gz"|"${run_dir}/${sample}/parsed_logs/${sample}_report.html")
                    ;;
                *)
                    record_remove "$file"
                    ;;
            esac
        done < <(find "${run_dir}/${sample}" -type f -print0)
    fi

    record_remove "${changeo_dir}/${sample}/${sample}.al2_igblast.fmt7"
    record_remove "${changeo_dir}/${sample}/${sample}.al2.fasta"
    record_remove "${changeo_dir}/${sample}/${sample}.al2_db-fail.tsv"
done

mv -f "$tmp_summary" "$summary"
{
    echo "Cleanup completed: $(date -Iseconds)"
    echo "Summary: $summary"
} > "$marker"
