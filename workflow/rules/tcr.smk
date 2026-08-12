rule raw_fastq_qc:
    container: WORKER_CONTAINER
    input:
        r1=lambda wildcards: config["samples"][wildcards.sample]["R1"],
        r2=lambda wildcards: config["samples"][wildcards.sample]["R2"]
    output:
        done="00_fastqc/{sample}/{sample}.raw_fastq_qc.done"
    params:
        outdir="00_fastqc/{sample}",
        staged_r1="00_fastqc/{sample}/{sample}_R1.fastq.gz",
        staged_r2="00_fastqc/{sample}/{sample}_R2.fastq.gz",
        fastp_json="00_fastqc/{sample}/{sample}_fastp.json",
        fastp_html="00_fastqc/{sample}/{sample}_fastp.html"
    threads: 2
    resources:
        mem_mb=4000,
        walltime="0-01:00:00",
        partition="norm"
    shell:
        r'''
        set -euo pipefail
        mkdir -p {params.outdir}
        ln -sf {input.r1} {params.staged_r1}
        ln -sf {input.r2} {params.staged_r2}

        fastp             --thread {threads}             --in1 {params.staged_r1}             --in2 {params.staged_r2}             --stdout             --disable_adapter_trimming             --disable_quality_filtering             --disable_length_filtering             --disable_trim_poly_g             --overrepresentation_analysis             --json {params.fastp_json}             --html {params.fastp_html}             --report_title "{wildcards.sample} raw FASTQ QC"             > /dev/null
        test -s {params.fastp_json}
        test -s {params.fastp_html}

        if command -v fastqc >/dev/null 2>&1; then
            fastqc --quiet --threads {threads} --outdir {params.outdir} {params.staged_r1} {params.staged_r2}
            test -s {params.outdir}/{wildcards.sample}_R1_fastqc.zip
            test -s {params.outdir}/{wildcards.sample}_R2_fastqc.zip
            printf 'tools\tfastp,fastqc\n' > {output.done}
        else
            printf 'tools\tfastp\n' > {output.done}
        fi
        '''


rule presto_process:
    container: WORKER_CONTAINER
    input:
        r1=lambda wildcards: config["samples"][wildcards.sample]["R1"],
        r2=lambda wildcards: config["samples"][wildcards.sample]["R2"],
        r1_primers=lambda wildcards: os.path.join(config["paths"]["refdir"], "IS_Human_R1_Primers.fasta"),
        r2_primers=lambda wildcards: os.path.join(config["paths"]["refdir"], "IS_Human_R2_Primers.fasta"),
        c_region=lambda wildcards: os.path.join(config["paths"]["refdir"], "IS_Human_C-Region.fasta"),
        immune_ref=lambda wildcards: os.path.join(config["paths"]["refdir"], "Immune_Ref.fasta"),
        qc_script="workflow/scripts/write_presto_qc_summary.py"
    output:
        uniq="01_run/{sample}/{sample}.uniqC_atleast-2.fastq",
        cregions="01_run/{sample}/{sample}_uniqC_atleast-2.tsv",
        qc="01_run/{sample}/{sample}.presto_qc.tsv"
    params:
        outdir="01_run/{sample}",
        parsed="01_run/{sample}/parsed_logs",
        threads=lambda wildcards: int(config["params"].get("threads", 2))
    threads: lambda wildcards: int(config["params"].get("threads", 2))
    resources:
        mem_mb=32000,
        walltime="1-00:00:00",
        partition="norm"
    shell:
        r'''
        set -euo pipefail
        mkdir -p {params.outdir} {params.parsed}
        cd {params.outdir}

        gunzip -c {input.r1} > {wildcards.sample}_R1.input.fastq
        gunzip -c {input.r2} > {wildcards.sample}_R2.input.fastq

        FilterSeq.py quality --nproc {params.threads} -s {wildcards.sample}_R1.input.fastq -o {wildcards.sample}_R1.q20.fastq -q 20 --log qual_R1.log
        FilterSeq.py quality --nproc {params.threads} -s {wildcards.sample}_R2.input.fastq -o {wildcards.sample}_R2.q20.fastq -q 20 --log qual_R2.log
        ParseLog.py -l qual_R1.log -f ID QUALITY -o parsed_logs/quality-1_table.tab
        ParseLog.py -l qual_R2.log -f ID QUALITY -o parsed_logs/quality-2_table.tab

        MaskPrimers.py score --nproc {params.threads} -s {wildcards.sample}_R1.q20.fastq -o {wildcards.sample}_R1.mask.fastq -p {input.r1_primers} --maxerror 0.2 --start 0 --log mask_R1.log
        MaskPrimers.py score --mode cut --barcode --nproc {params.threads} -s {wildcards.sample}_R2.q20.fastq -o {wildcards.sample}_R2.cut.fastq -p {input.r2_primers} --maxerror 0.5 --start 17 --log cut_R2.log
        ParseLog.py -l mask_R1.log -f ID BARCODE PRIMER ERROR -o parsed_logs/consensus-1_table.tab
        ParseLog.py -l cut_R2.log -f ID BARCODE PRIMER ERROR -o parsed_logs/consensus-2_table.tab

        PairSeq.py -1 {wildcards.sample}_R1.mask.fastq -2 {wildcards.sample}_R2.cut.fastq --2f BARCODE --coord illumina

        BuildConsensus.py --nproc {params.threads} -s {wildcards.sample}_R1.mask_pair-pass.fastq -o {wildcards.sample}_R1.consensus.fastq -n 1 --bf BARCODE -q 0 --freq 0.6 --maxgap 0.5 --pf PRIMER --prcons 0.6 --maxerror 0.1 --log consensus_R1.log
        BuildConsensus.py --nproc {params.threads} -s {wildcards.sample}_R2.cut_pair-pass.fastq -o {wildcards.sample}_R2.consensus.fastq -n 1 --bf BARCODE -q 0 --freq 0.6 --maxgap 0.5 --pf PRIMER --prcons 0.6 --maxerror 0.1 --log consensus_R2.log
        ParseLog.py -l consensus_R1.log -f BARCODE SEQCOUNT CONSCOUNT PRIMER PRCONS PRCOUNT PRFREQ ERROR -o parsed_logs/consensus-1_table.tab
        ParseLog.py -l consensus_R2.log -f BARCODE SEQCOUNT CONSCOUNT PRIMER PRCONS PRCOUNT PRFREQ ERROR -o parsed_logs/consensus-2_table.tab

        PairSeq.py -1 {wildcards.sample}_R1.consensus.fastq -2 {wildcards.sample}_R2.consensus.fastq --coord presto

        AssemblePairs.py sequential --nproc {params.threads} -1 {wildcards.sample}_R1.consensus_pair-pass.fastq -2 {wildcards.sample}_R2.consensus_pair-pass.fastq -r {input.immune_ref} -o {wildcards.sample}.assembled.fastq --rc tail --1f CONSCOUNT --2f PRCONS CONSCOUNT --alpha 1e-05 --maxerror 0.3 --minlen 8 --maxlen 1000 --scanrev --minident 0.5 --evalue 1e-05 --maxhits 100 --aligner blastn --coord presto --log assemble.log
        ParseLog.py -l assemble.log -f ID REFID LENGTH OVERLAP GAP ERROR PVALUE EVALUE1 EVALUE2 IDENTITY FIELDS1 FIELDS2 -o parsed_logs/assemble_table.tab

        FilterSeq.py maskqual --nproc {params.threads} -s {wildcards.sample}.assembled.fastq -o {wildcards.sample}.q0.fastq -q 0 --l mask.log

        MaskPrimers.py align --nproc {params.threads} -s {wildcards.sample}.q0.fastq -o {wildcards.sample}.Cprimer.fastq -p {input.c_region} --maxlen 50 --gap 1 1 --maxerror 0.2 --l CPrimer.log
        ParseLog.py -l CPrimer.log -f ID PRIMER ERROR -o parsed_logs/cregion_table.tab

        Rscript -e 'library(prestor); report_abseq3("parsed_logs", sample="{wildcards.sample}", output_dir="parsed_logs", output_file="{wildcards.sample}_report.html", format="html", quiet=FALSE)'

        ParseHeaders.py rename -s {wildcards.sample}.Cprimer.fastq -o {wildcards.sample}.rename.fastq -f PRIMER -k CREGION --act first
        ParseHeaders.py collapse -s {wildcards.sample}.rename.fastq -o {wildcards.sample}.collapse.fastq -f CONSCOUNT --act min
        ParseHeaders.py table -s {wildcards.sample}.collapse.fastq -o {wildcards.sample}_cregions.tsv -f ID PRCONS CREGION CONSCOUNT
        CollapseSeq.py -s {wildcards.sample}.collapse.fastq -o {wildcards.sample}.uniqC.fastq -n 0 --uf PRCONS CREGION --cf CONSCOUNT --act sum --inner --keepmiss
        ParseHeaders.py table -s {wildcards.sample}.uniqC.fastq -o {wildcards.sample}_uniqcregions.tsv -f ID PRCONS CREGION CONSCOUNT
        SplitSeq.py group -s {wildcards.sample}.uniqC.fastq -f CONSCOUNT --num 2
        ParseHeaders.py table -s {wildcards.sample}.uniqC_atleast-2.fastq -o {wildcards.sample}_uniqC_atleast-2.tsv -f ID PRCONS CREGION CONSCOUNT
        python ../../{input.qc_script} --sample {wildcards.sample} --sample-dir . --output {wildcards.sample}.presto_qc.tsv

        rm -f {wildcards.sample}_R1.input.fastq {wildcards.sample}_R2.input.fastq
        '''

rule changeo_igblast:
    container: WORKER_CONTAINER
    input:
        fastq="01_run/{sample}/{sample}.uniqC_atleast-2.fastq"
    output:
        db="02_changeo/{sample}/{sample}.al2_db-pass.tsv"
    params:
        outdir="02_changeo/{sample}",
        threads=lambda wildcards: int(config["params"].get("threads", 2)),
        species=lambda wildcards: config["params"].get("species", "human"),
        receptor=lambda wildcards: config["params"].get("receptor", "tr"),
        igdata="/usr/local/share/igblast",
        germlines=lambda wildcards: "/usr/local/share/germlines/imgt/{species}/vdj".format(
            species=config["params"].get("species", "human")
        )
    threads: lambda wildcards: int(config["params"].get("threads", 2))
    resources:
        mem_mb=16000,
        walltime="0-08:00:00",
        partition="norm"
    shell:
        r'''
        set -euo pipefail
        mkdir -p {params.outdir}
        seqkit fq2fa --quiet --threads {params.threads} {input.fastq} -o {params.outdir}/{wildcards.sample}.al2.fasta
        AssignGenes.py igblast \
            -s {params.outdir}/{wildcards.sample}.al2.fasta \
            -b {params.igdata} \
            --organism {params.species} \
            --loci {params.receptor} \
            --format blast \
            --outdir {params.outdir} \
            --outname {wildcards.sample}.al2 \
            --nproc {params.threads}
        MakeDb.py igblast \
            -i {params.outdir}/{wildcards.sample}.al2_igblast.fmt7 \
            -s {params.outdir}/{wildcards.sample}.al2.fasta \
            -r {params.germlines} \
            --outdir {params.outdir} \
            --outname {wildcards.sample}.al2 \
            --format airr \
            --extended \
            --failed \
            --nproc {params.threads}
        test -s {output.db}
        '''
