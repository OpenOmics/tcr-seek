suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(immunarch)
})

parser <- ArgumentParser(description = "Create immunarch input files and an immunarch object from Change-O AIRR outputs.")
parser$add_argument("--changeo-dir", required = TRUE, help = "Directory containing per-sample Change-O outputs")
parser$add_argument("--output-dir", required = TRUE, help = "Output directory for ImmunarchInput and Immdata.rds")
parser$add_argument("--metadata", default = "", help = "Optional metadata TSV with a Sample column")
parser$add_argument("--sample-info", default = "", help = "Optional sample information TSV with Sample, SampleID, sample, or sample_id column")
parser$add_argument("--samples", required = TRUE, help = "Comma-separated expected sample names")
args <- parser$parse_args()

read_sample_info <- function(path) {
  info <- read.table(path, sep = "\t", header = TRUE, check.names = FALSE, quote = "", comment.char = "")
  sample_col <- intersect(c("Sample", "SampleID", "sample", "sample_id"), colnames(info))
  if (!length(sample_col)) {
    stop("Sample info file must contain one of these columns: Sample, SampleID, sample, sample_id")
  }
  info <- info %>% rename(Sample = all_of(sample_col[[1]]))
  info %>% distinct(.data$Sample, .keep_all = TRUE)
}

merge_sample_info <- function(metadata, sample_info) {
  overlap <- setdiff(intersect(colnames(metadata), colnames(sample_info)), "Sample")
  merged <- metadata %>% left_join(sample_info, by = "Sample", suffix = c(".metadata", ".sample_info"))
  for (column in overlap) {
    metadata_col <- paste0(column, ".metadata")
    sample_info_col <- paste0(column, ".sample_info")
    merged[[column]] <- dplyr::coalesce(merged[[sample_info_col]], merged[[metadata_col]])
    merged[[metadata_col]] <- NULL
    merged[[sample_info_col]] <- NULL
  }
  merged
}

expected_samples <- strsplit(args$samples, ",", fixed = TRUE)[[1]]
expected_samples <- expected_samples[nzchar(expected_samples)]

input_dir <- file.path(args$output_dir, "ImmunarchInput")
dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)

airr_files <- Sys.glob(file.path(args$changeo_dir, "*", "*_db-pass.tsv"))
if (!length(airr_files)) {
  stop("No Change-O db-pass TSV files found under ", args$changeo_dir)
}

sample_from_file <- function(path) {
  name <- basename(path)
  name <- sub("_db-pass\\.tsv$", "", name)
  name <- sub("\\.al2$", "", name)
  name
}

observed_samples <- sample_from_file(airr_files)
for (i in seq_along(airr_files)) {
  file.copy(airr_files[[i]], file.path(input_dir, paste0(observed_samples[[i]], ".tsv")), overwrite = TRUE)
}

if (nzchar(args$metadata)) {
  metadata <- read.table(args$metadata, sep = "\t", header = TRUE, check.names = FALSE, quote = "", comment.char = "")
  if (!"Sample" %in% colnames(metadata)) {
    stop("Metadata file must contain a Sample column")
  }
} else {
  metadata <- tibble(Sample = expected_samples)
}

missing_samples <- setdiff(metadata$Sample, observed_samples)
write.table(missing_samples, file.path(args$output_dir, "SampleNames_notSequenced.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

metadata <- metadata %>% filter(.data$Sample %in% observed_samples)
if (nzchar(args$sample_info)) {
  sample_info <- read_sample_info(args$sample_info)
  metadata <- merge_sample_info(metadata, sample_info)
}
if (!nrow(metadata)) {
  stop("No metadata samples matched Change-O db-pass files")
}
write.table(metadata, file.path(input_dir, "metadata.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

immdata <- repLoad(input_dir)
immdata$meta <- metadata %>%
  filter(.data$Sample %in% names(immdata$data)) %>%
  arrange(match(.data$Sample, names(immdata$data)))
saveRDS(immdata, file.path(args$output_dir, "Immdata.rds"))
capture.output(sessionInfo(), file = file.path(args$output_dir, "SessionInfo.txt"))
