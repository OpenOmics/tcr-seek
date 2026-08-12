# Workflow Method

`tcr-seek` currently implements a TCR-only bulk repertoire workflow.

## pRESTO

For each sample, paired FASTQs are quality-filtered, primer-masked, paired, consensus-collapsed, assembled, C-region annotated, deduplicated, and split to retain `CONSCOUNT >= 2` reads.

## Change-O / IgBLAST

Each `*.uniqC_atleast-2.fastq` file is converted to FASTA, assigned with `AssignGenes.py igblast`, and converted to AIRR-style Change-O tables with `MakeDb.py igblast --format airr --extended --failed`.

## immunarch

All per-sample `*_db-pass.tsv` files are copied into `03_immunarch/ImmunarchInput`, loaded with `immunarch::repLoad()`, and saved as `Immdata.rds`. The downstream `repertoire_diversity` rule then calculates sample-level diversity metrics from the imported object.

## Repertoire Diversity

The diversity rule converts `Immdata.rds` to a vegan community matrix where rows are samples, columns are clonotypes, and values are `Clones` abundance. By default, the rule uses the `strict` preset: CDR3 nucleotide plus available V/D/J/C gene calls. The definition can be changed with `--clonotype-definition`; supported presets are `gene`, `nt`, `aa`, and `strict`, and a custom comma-separated immunarch column list can also be supplied.

The workflow writes `diversity_summary.rds`, `diversity_by_clonotype_definition.rds`, `rarefaction.rds`, `rarefaction_curve.rds`, `DiversityMethod.txt`, and the Quarto HTML report `qc_report.html` under `03_immunarch`. `diversity_summary.rds` contains the selected clonotype definition used for the primary QC metrics and rarefaction curves. `diversity_by_clonotype_definition.rds` contains the same diversity metrics recalculated for the `gene`, `nt`, `aa`, and `strict` presets, and the report displays those rows in a clonotype-definition comparison section. Observed richness is the number of detected clonotypes before depth adjustment. Shannon diversity summarizes richness and evenness and decreases when a repertoire is dominated by expanded clones. vegan Simpson diversity is reported as `1 - D`, so larger values indicate higher diversity and lower dominance. Pielou evenness is Shannon divided by `log(richness)` and measures how evenly abundance is distributed across observed clonotypes. Library size is the total `Clones` abundance used by vegan. Effective clonotypes, top clone fraction, and top-10 clone fraction are also reported as practical QC checks for repertoire dominance and expansion.

Rarefaction estimates expected richness after subsampling each sample to a common abundance depth. The pipeline uses the minimum sample library size, calculated as `rowSums` of the clonotype abundance matrix. This is the theoretically appropriate input for `vegan::rarefy()`; using the minimum observed clonotype count would rarefy to richness rather than to sequencing or clone abundance depth. Rarefaction curves show how expected richness accumulates with increasing depth and are useful for identifying samples that remain undersampled.

`RichnessPerMillion` is retained only as a descriptive rough normalization. Richness is not linear with sequencing depth: newly sampled clonotypes accumulate quickly at first and then asymptotically. Dividing observed richness by reads-per-million can therefore exaggerate shallow samples, compress deep samples, and create artifacts when library sizes differ greatly. For cross-sample inference, prefer rarefied richness, coverage-aware estimators, or models that include library depth. Optional `--sample-info` annotations are merged into `immdata`; optional `--qc-variables` selects which merged variables receive their own report tabs.
