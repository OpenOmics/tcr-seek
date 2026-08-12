rule immunarch_object:
    container: WORKER_CONTAINER
    input:
        dbs=expand("02_changeo/{sample}/{sample}.al2_db-pass.tsv", sample=SAMPLES),
        script="workflow/scripts/create_immunarch_object.R"
    output:
        rds="03_immunarch/Immdata.rds",
        metadata="03_immunarch/ImmunarchInput/metadata.txt",
        missing="03_immunarch/SampleNames_notSequenced.txt",
        session="03_immunarch/SessionInfo.txt"
    params:
        metadata=lambda wildcards: config["paths"].get("metadata", ""),
        samples=",".join(SAMPLES),
        sample_info=lambda wildcards: config["paths"].get("sample_info", "")
    threads: 1
    resources:
        mem_mb=100000,
        walltime="3-00:00:00",
        partition="norm"
    shell:
        r'''
        set -euo pipefail
        mkdir -p 03_immunarch/ImmunarchInput
        Rscript {input.script} \
            --changeo-dir 02_changeo \
            --output-dir 03_immunarch \
            --metadata "{params.metadata}" \
            --sample-info "{params.sample_info}" \
            --samples "{params.samples}"
        '''


rule repertoire_diversity:
    container: WORKER_CONTAINER
    input:
        rds="03_immunarch/Immdata.rds",
        script="workflow/scripts/calculate_diversity_per_sample.R"
    output:
        summary="03_immunarch/diversity_summary.rds",
        comparison="03_immunarch/diversity_by_clonotype_definition.rds",
        rarefaction="03_immunarch/rarefaction.rds",
        curve="03_immunarch/rarefaction_curve.rds",
        method="03_immunarch/DiversityMethod.txt"
    params:
        outdir="03_immunarch",
        clonotype=lambda wildcards: config["params"].get("clonotype_definition", "strict"),
        run_rarefaction=lambda wildcards: str(config["params"].get("run_rarefaction", True)).lower()
    threads: 1
    resources:
        mem_mb=100000,
        walltime="0-04:00:00",
        partition="norm"
    shell:
        r"""
        set -euo pipefail
        Rscript {input.script} \
            --immdata {input.rds} \
            --output-dir {params.outdir} \
            --clonotype "{params.clonotype}" \
            --run-rarefaction "{params.run_rarefaction}"
        """


rule qc_report:
    container: WORKER_CONTAINER
    input:
        immdata="03_immunarch/Immdata.rds",
        diversity="03_immunarch/diversity_summary.rds",
        comparison="03_immunarch/diversity_by_clonotype_definition.rds",
        rarefaction_curve="03_immunarch/rarefaction_curve.rds",
        qmd=lambda wildcards: "workflow/scripts/qc_report.qmd" if RUN_RAREFACTION else "workflow/scripts/qc_report_no_rarefaction.qmd",
        script="workflow/scripts/render_qc_report.R",
        logo="resources/tcr-seek.svg"
    output:
        html="03_immunarch/qc_report.html"
    params:
        outdir="03_immunarch",
        qc_variables=lambda wildcards: config["params"].get("qc_variables", ""),
        run_rarefaction=lambda wildcards: str(config["params"].get("run_rarefaction", True)).lower()
    threads: 1
    resources:
        mem_mb=16000,
        walltime="0-04:00:00",
        partition="norm"
    shell:
        r"""
        set -euo pipefail
        Rscript {input.script} \
            --qmd {input.qmd} \
            --output-dir {params.outdir} \
            --immdata {input.immdata} \
            --diversity {input.diversity} \
            --rarefaction-curve {input.rarefaction_curve} \
            --run-rarefaction "{params.run_rarefaction}" \
            --qc-variables "{params.qc_variables}" \
            --logo {input.logo}
        """

rule multiqc_input:
    input:
        presto_qc=expand("01_run/{sample}/{sample}.presto_qc.tsv", sample=SAMPLES),
        dbs=expand("02_changeo/{sample}/{sample}.al2_db-pass.tsv", sample=SAMPLES),
        script="workflow/scripts/write_multiqc_input.py"
    output:
        summary="04_multiqc/tcr_seek_multiqc.tsv"
    params:
        samples=",".join(SAMPLES)
    threads: 1
    resources:
        mem_mb=4000,
        walltime="0-01:00:00",
        partition="norm"
    shell:
        r"""
        set -euo pipefail
        python {input.script} \
            --run-dir . \
            --samples "{params.samples}" \
            --output {output.summary}
        """


rule multiqc_report:
    input:
        summary="04_multiqc/tcr_seek_multiqc.tsv",
        raw_qc=expand("00_fastqc/{sample}/{sample}.raw_fastq_qc.done", sample=SAMPLES),
        config="workflow/scripts/multiqc_config.yaml",
        plugin_py="workflow/multiqc_plugins/tcr_seek/multiqc_tcr_seek/modules/tcr_seek.py",
        plugin_project="workflow/multiqc_plugins/tcr_seek/pyproject.toml"
    output:
        html="04_multiqc/multiqc_report.html"
    params:
        outdir="04_multiqc",
        plugin_dir="workflow/multiqc_plugins/tcr_seek",
        plugin_target="04_multiqc/plugin_env"
    threads: 1
    resources:
        mem_mb=8000,
        walltime="0-01:00:00",
        partition="norm"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} {params.plugin_target}
        if ! command -v multiqc >/dev/null 2>&1 && command -v module >/dev/null 2>&1; then
            module load multiqc || true
        fi
        if ! command -v multiqc >/dev/null 2>&1; then
            echo "MultiQC executable was not found in PATH, and module load multiqc did not provide it." >&2
            exit 1
        fi
        python -m pip install --quiet --no-deps --upgrade --target {params.plugin_target} {params.plugin_dir}
        rm -rf {params.plugin_dir}/build {params.plugin_dir}/*.egg-info
        runtime_config="{params.outdir}/multiqc_runtime_config.yaml"
        if find 00_fastqc -name '*_fastqc.zip' -print -quit | grep -q .; then
            cat > "$runtime_config" <<'YAML'
title: "tcr-seek MultiQC report"
subtitle: "Bulk TCR-seq processing and repertoire QC"
module_order:
  - fastqc
  - tcr_seek
  - fastp
sp:
  tcr_seek/tcr_seek_multiqc:
    fn: "tcr_seek_multiqc.tsv"
YAML
        else
            cat > "$runtime_config" <<'YAML'
title: "tcr-seek MultiQC report"
subtitle: "Bulk TCR-seq processing and repertoire QC"
module_order:
  - fastp
  - tcr_seek
sp:
  tcr_seek/tcr_seek_multiqc:
    fn: "tcr_seek_multiqc.tsv"
YAML
        fi
        PYTHONPATH="{params.plugin_target}:${{PYTHONPATH:-}}" \
            multiqc {input.summary} 00_fastqc \
                --config "$runtime_config" \
                --outdir {params.outdir} \
                --filename multiqc_report.html \
                --force
        test -s {output.html}
        """

