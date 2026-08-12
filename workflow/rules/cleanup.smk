rule cleanup_intermediates:
    input:
        dbs=expand(f"{CHANGEO_DIR}/{{sample}}/{{sample}}.al2_db-pass.tsv", sample=SAMPLES),
        immdata=f"{IMMUNARCH_DIR}/Immdata.rds",
        metadata=f"{IMMUNARCH_DIR}/ImmunarchInput/metadata.txt",
        diversity=f"{IMMUNARCH_DIR}/diversity_summary.rds",
        rarefaction=f"{IMMUNARCH_DIR}/rarefaction.rds",
        rarefaction_curve=f"{IMMUNARCH_DIR}/rarefaction_curve.rds",
        qc_report=f"{IMMUNARCH_DIR}/qc_report.html",
        script="workflow/scripts/cleanup_tcr_seek_intermediates.sh"
    output:
        marker="cleanup/tcr_seek_cleanup.done"
    params:
        samples=" ".join(SAMPLES),
        run_dir=RUN_DIR,
        changeo_dir=CHANGEO_DIR,
        immunarch_dir=IMMUNARCH_DIR
    threads: 1
    resources:
        mem_mb=1000,
        walltime="0-08:00:00",
        partition="norm"
    shell:
        r'''
        set -euo pipefail
        bash {input.script} \
          --output-dir . \
          --samples "{params.samples}" \
          --run-dir {params.run_dir} \
          --changeo-dir {params.changeo_dir} \
          --immunarch-dir {params.immunarch_dir} \
          --marker {output.marker}
        '''
