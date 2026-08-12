from __future__ import annotations

import csv

try:
    from multiqc.modules.base_module import BaseMultiqcModule
except ImportError:  # MultiQC < 1.13 compatibility
    from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, table


def parse_value(value: str):
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return value


class MultiqcModule(BaseMultiqcModule):
    """MultiQC module for tcr-seek aggregate QC tables."""

    def __init__(self):
        super().__init__(
            name="tcr-seek",
            anchor="tcr_seek",
            href="https://multiqc.info/",
            info="summarizes bulk TCR-seq pRESTO, Change-O/IgBLAST, and repertoire QC metrics.",
        )

        self.tcr_seek_data = {}
        for log_file in self.find_log_files("tcr_seek/tcr_seek_multiqc", filehandles=True):
            reader = csv.DictReader(log_file["f"], delimiter="\t")
            for row in reader:
                sample = row.get("Sample")
                if not sample:
                    continue
                self.tcr_seek_data[sample] = {
                    key: parse_value(value)
                    for key, value in row.items()
                    if key != "Sample" and value != ""
                }

        if not self.tcr_seek_data:
            raise UserWarning

        self.write_data_file(self.tcr_seek_data, "multiqc_tcr_seek")
        self.add_general_stats()
        self.add_fastq_interpretation_sections()
        self.add_processing_section()
        self.add_changeo_section()
        self.add_table_section()

    def add_general_stats(self):
        headers = {
            "R1InputReads": {
                "title": "R1 reads",
                "description": "Input R1 reads",
                "format": "{:,.0f}",
                "scale": "YlGnBu",
            },
            "FinalConsCountGe2Sequences": {
                "title": "Final seqs",
                "description": "Collapsed sequences with CONSCOUNT >= 2",
                "format": "{:,.0f}",
                "scale": "PuBu",
            },
            "ChangeoPassFraction": {
                "title": "Change-O pass",
                "description": "Fraction of final sequences in Change-O db-pass output",
                "format": "{:.1%}",
                "min": 0,
                "max": 1,
                "scale": "RdYlGn",
            },
            "ProductiveFraction": {
                "title": "Productive",
                "description": "Fraction of Change-O db-pass sequences marked productive",
                "format": "{:.1%}",
                "min": 0,
                "max": 1,
                "scale": "RdYlGn",
            },
        }
        self.general_stats_addcols(self.tcr_seek_data, headers)

    def add_fastq_interpretation_sections(self):
        self.add_section(
            name="Raw FASTQ Duplication Note",
            anchor="tcr_seek_fastq_duplication_note",
            content=(
                '<div class="alert alert-info">'
                "<strong>FASTQ duplication interpretation for bulk TCR-seq:</strong> "
                "High sequence duplication is expected in targeted repertoire libraries because "
                "PCR amplification, UMI/consensus generation, and true clonal expansion all increase "
                "repeated sequence observations. Use FastQC or fastp duplication primarily to compare libraries "
                "within the same run or protocol, and interpret it together with downstream consensus "
                "counts, recovered clonotypes, and productivity metrics."
                "</div>"
            ),
            autoformat=False,
        )
        self.add_section(
            name="Raw FASTQ GC Note",
            anchor="tcr_seek_fastq_gc_note",
            content=(
                '<div class="alert alert-info">'
                "<strong>FASTQ GC interpretation for bulk TCR-seq:</strong> "
                "GC deviation is useful for detecting unusual libraries, contamination, severe primer "
                "or adapter carryover, or sample swaps. Because this assay is targeted rather than "
                "whole-transcriptome or whole-genome sequencing, GC profiles may not follow generic "
                "FastQC expectations; compare samples against the cohort and sequencing batch."
                "</div>"
            ),
            autoformat=False,
        )

    def add_processing_section(self):
        categories = {
            "R1InputReads": {"name": "R1 input"},
            "R1Q20Reads": {"name": "R1 Q20"},
            "InitialPairPassReads": {"name": "Initial pair pass"},
            "ConsensusPairPassReads": {"name": "Consensus pair pass"},
            "AssembledReads": {"name": "Assembled"},
            "CRegionReads": {"name": "C-region"},
            "FinalConsCountGe2Sequences": {"name": "Final CONSCOUNT >= 2"},
        }
        self.add_section(
            name="Read Processing",
            anchor="tcr_seek_read_processing",
            description="Read and sequence counts across the pRESTO processing path.",
            plot=bargraph.plot(
                self.tcr_seek_data,
                categories,
                pconfig={
                    "id": "tcr_seek_read_processing",
                    "title": "tcr-seek: read processing counts",
                    "ylab": "Sequences / reads",
                },
            ),
        )

    def add_changeo_section(self):
        categories = {
            "ChangeoPassFraction": {"name": "Change-O pass"},
            "ProductiveFraction": {"name": "Productive"},
            "InFrameFraction": {"name": "In frame"},
            "StopCodonFraction": {"name": "Stop codon"},
        }
        self.add_section(
            name="Change-O / IgBLAST",
            anchor="tcr_seek_changeo",
            description="Assignment and productivity fractions from Change-O db-pass tables.",
            plot=bargraph.plot(
                self.tcr_seek_data,
                categories,
                pconfig={
                    "id": "tcr_seek_changeo_fractions",
                    "title": "tcr-seek: Change-O / IgBLAST fractions",
                    "ylab": "Fraction",
                    "ymin": 0,
                    "ymax": 1,
                },
            ),
        )

    def add_table_section(self):
        headers = {
            "R1InputReads": {"title": "R1 input", "format": "{:,.0f}"},
            "R1Q20PassFraction": {"title": "R1 Q20", "format": "{:.1%}"},
            "InitialPairPassFraction": {"title": "Pair pass", "format": "{:.1%}"},
            "AssemblyFraction": {"title": "Assembly", "format": "{:.1%}"},
            "FinalConsCountGe2Sequences": {"title": "Final seqs", "format": "{:,.0f}"},
            "FinalConsCountGe2Abundance": {"title": "Final abundance", "format": "{:,.0f}"},
            "ChangeoDbPassSequences": {"title": "Change-O db-pass", "format": "{:,.0f}"},
            "ChangeoPassFraction": {"title": "Change-O pass", "format": "{:.1%}"},
            "ProductiveFraction": {"title": "Productive", "format": "{:.1%}"},
            "MeanConsensusCount": {"title": "Mean consensus", "format": "{:.2f}"},
            "TopVCall": {"title": "Top V"},
            "TopJCall": {"title": "Top J"},
            "TopLocus": {"title": "Top locus"},
        }
        self.add_section(
            name="Sample QC Table",
            anchor="tcr_seek_sample_qc",
            description="Sample-level tcr-seek metrics collected for MultiQC.",
            plot=table.plot(
                self.tcr_seek_data,
                headers,
                pconfig={
                    "id": "tcr_seek_sample_qc_table",
                    "title": "tcr-seek: sample QC metrics",
                },
            ),
        )
