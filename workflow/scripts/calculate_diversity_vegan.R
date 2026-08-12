suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(vegan)
})

parser <- ArgumentParser(description = "Calculate repertoire diversity and rarefaction metrics from an immunarch Immdata.rds object.")
parser$add_argument("--immdata", required = TRUE, help = "Path to Immdata.rds produced by immunarch::repLoad().")
parser$add_argument("--output-dir", required = TRUE, help = "Directory where diversity_summary.rds, rarefaction.rds, and rarefaction_curve.rds will be written.")
parser$add_argument("--clonotype", default = "strict", help = "Clonotype definition preset or comma-separated immunarch columns. Presets: gene, nt, aa, strict. Default: strict")
parser$add_argument("--run-rarefaction", default = "true", help = "Whether to calculate rarefaction metrics and curves. Accepted values: true/false. Default: true")
args <- parser$parse_args()

parse_bool <- function(value) {
  normalized <- tolower(trimws(as.character(value)))
  if (normalized %in% c("true", "t", "yes", "y", "1")) {
    return(TRUE)
  }
  if (normalized %in% c("false", "f", "no", "n", "0")) {
    return(FALSE)
  }
  stop("--run-rarefaction must be true or false")
}

run_rarefaction <- parse_bool(args$run_rarefaction)

first_repertoire_columns <- function(immdata) {
  if (!length(immdata$data)) {
    stop("Immdata object has no repertoire tables")
  }
  names(immdata$data[[1]])
}

resolve_clonotype_definition <- function(value, available_columns) {
  raw_value <- trimws(value)
  lower_value <- tolower(raw_value)
  gene_candidates <- c("V.name", "D.name", "J.name", "C.name")
  preset_columns <- switch(
    lower_value,
    gene = intersect(gene_candidates, available_columns),
    nt = intersect(c("CDR3.nt", "junction", "cdr3"), available_columns)[1],
    aa = intersect(c("CDR3.aa", "junction_aa"), available_columns)[1],
    strict = c(intersect(c("CDR3.nt", "junction", "cdr3"), available_columns)[1], intersect(gene_candidates, available_columns)),
    NULL
  )
  if (!is.null(preset_columns)) {
    columns <- preset_columns[!is.na(preset_columns) & nzchar(preset_columns)]
    if (!length(columns)) {
      stop("Could not resolve clonotype preset '", raw_value, "' from available columns: ", paste(available_columns, collapse = ", "))
    }
    list(mode = lower_value, columns = unique(columns))
  } else {
    columns <- trimws(strsplit(raw_value, ",", fixed = TRUE)[[1]])
    columns <- columns[nzchar(columns)]
    if (!length(columns)) {
      stop("At least one clonotype column or preset must be supplied with --clonotype")
    }
    list(mode = "custom", columns = columns)
  }
}

immdata_to_vegan <- function(immdata, clonotype_columns, clonotype_mode){
  samples <- names(immdata$data)
  names(samples) <- samples
  immdata2 <- lapply(samples,
                     function(sample){
                       missing_cols <- setdiff(c(clonotype_columns, "Clones"), names(immdata$data[[sample]]))
                       if (length(missing_cols) > 0) {
                         stop("Missing expected columns in sample ", sample, ": ", paste(missing_cols, collapse = ", "))
                       }
                       immdata$data[[sample]] %>%
                         dplyr::select(all_of(c(clonotype_columns, "Clones"))) %>%
                         dplyr::filter(dplyr::if_all(all_of(clonotype_columns), ~ !is.na(.x) & .x != "")) %>%
                         tidyr::unite("Clonotype", all_of(clonotype_columns), remove = FALSE) %>%
                         dplyr::group_by(Clonotype) %>%
                         dplyr::summarise(count = sum(Clones), .groups = "drop") %>%
                         dplyr::rename(!!sample := count)
                     }
  )

  immdata3 <- immdata2[[1]]
  for(sample in samples[-1]){
    immdata3 <- dplyr::full_join(immdata3, immdata2[[sample]], by = "Clonotype")
  }
  immdata4 <- t(immdata3[,-1])
  colnames(immdata4) <- immdata3$Clonotype
  immdata4[is.na(immdata4)] <- 0
  storage.mode(immdata4) <- "numeric"
  sample.depths <- rowSums(immdata4)
  observed.richness <- vegan::specnumber(immdata4)
  list(
    immdata = immdata4,
    depth = sample.depths,
    observed_richness = observed.richness,
    clonotype_mode = clonotype_mode,
    clonotype_definition = paste(clonotype_columns, collapse = " + ")
  )
}

immdata <- readRDS(args$immdata)
clonotype <- resolve_clonotype_definition(args$clonotype, first_repertoire_columns(immdata))
vegan.data <- immdata_to_vegan(immdata = immdata, clonotype_columns = clonotype$columns, clonotype_mode = clonotype$mode)
vd_matrix <- vegan.data$immdata

if (any(vegan.data$depth <= 0)) {
  zero_samples <- names(vegan.data$depth)[vegan.data$depth <= 0]
  stop("Cannot calculate diversity for samples with zero clonotype abundance: ", paste(zero_samples, collapse = ", "))
}

samples <- rownames(vd_matrix)
names(samples) <- samples
rarefaction_depth <- if (run_rarefaction) min(vegan.data$depth) else NA_real_

if (run_rarefaction) {
  rarefaction <- vegan::rarefy(x = vd_matrix, sample = rarefaction_depth, se = FALSE, MARGIN = 1)
  rarefaction.curve <- vegan::rarecurve(vd_matrix, step = 100, sample = rarefaction_depth)
  names(rarefaction.curve) <- samples
} else {
  rarefaction <- rep(NA_real_, length(samples))
  names(rarefaction) <- samples
  rarefaction.curve <- list()
}

richness <- vegan::specnumber(vd_matrix)
shannon <- vegan::diversity(vd_matrix, index = "shannon")
simpson <- vegan::diversity(vd_matrix, index = "simpson")
evenness <- shannon / log(richness)
library_size <- rowSums(vd_matrix)
richness_per_million <- richness / (library_size / 1e6)
top_clone_fraction <- apply(vd_matrix, 1, max) / library_size
top10_clone_fraction <- apply(vd_matrix, 1, function(x) sum(sort(x, decreasing = TRUE)[seq_len(min(10, length(x)))]) ) / library_size
effective_clonotypes <- exp(shannon)

diversity_summary <- data.frame(
  Sample = rownames(vd_matrix),
  Richness = richness,
  Shannon = shannon,
  Simpson = simpson,
  Evenness = evenness,
  LibrarySize = library_size,
  RichnessPerMillion = richness_per_million,
  EffectiveClonotypes = effective_clonotypes,
  TopCloneFraction = top_clone_fraction,
  Top10CloneFraction = top10_clone_fraction,
  Rarefaction = as.numeric(rarefaction),
  RarefactionDepth = rarefaction_depth,
  ClonotypeMode = vegan.data$clonotype_mode,
  ClonotypeDefinition = vegan.data$clonotype_definition,
  stringsAsFactors = FALSE
)

definition_modes <- c("gene", "nt", "aa", "strict")
calculate_summary_for_definition <- function(clonotype_value) {
  resolved <- resolve_clonotype_definition(clonotype_value, first_repertoire_columns(immdata))
  mode_data <- immdata_to_vegan(immdata = immdata, clonotype_columns = resolved$columns, clonotype_mode = resolved$mode)
  mode_matrix <- mode_data$immdata

  if (any(mode_data$depth <= 0)) {
    zero_samples <- names(mode_data$depth)[mode_data$depth <= 0]
    stop("Cannot calculate diversity for samples with zero clonotype abundance: ", paste(zero_samples, collapse = ", "))
  }

  mode_richness <- vegan::specnumber(mode_matrix)
  mode_shannon <- vegan::diversity(mode_matrix, index = "shannon")
  mode_simpson <- vegan::diversity(mode_matrix, index = "simpson")
  mode_evenness <- mode_shannon / log(mode_richness)
  mode_library_size <- rowSums(mode_matrix)
  mode_richness_per_million <- mode_richness / (mode_library_size / 1e6)
  mode_top_clone_fraction <- apply(mode_matrix, 1, max) / mode_library_size
  mode_top10_clone_fraction <- apply(mode_matrix, 1, function(x) sum(sort(x, decreasing = TRUE)[seq_len(min(10, length(x)))]) ) / mode_library_size
  mode_effective_clonotypes <- exp(mode_shannon)
  mode_rarefaction_depth <- if (run_rarefaction) min(mode_data$depth) else NA_real_
  if (run_rarefaction) {
    mode_rarefaction <- vegan::rarefy(x = mode_matrix, sample = mode_rarefaction_depth, se = FALSE, MARGIN = 1)
  } else {
    mode_rarefaction <- rep(NA_real_, nrow(mode_matrix))
  }

  data.frame(
    Sample = rownames(mode_matrix),
    Richness = mode_richness,
    Shannon = mode_shannon,
    Simpson = mode_simpson,
    Evenness = mode_evenness,
    LibrarySize = mode_library_size,
    RichnessPerMillion = mode_richness_per_million,
    EffectiveClonotypes = mode_effective_clonotypes,
    TopCloneFraction = mode_top_clone_fraction,
    Top10CloneFraction = mode_top10_clone_fraction,
    Rarefaction = as.numeric(mode_rarefaction),
    RarefactionDepth = mode_rarefaction_depth,
    ClonotypeMode = mode_data$clonotype_mode,
    ClonotypeDefinition = mode_data$clonotype_definition,
    stringsAsFactors = FALSE
  )
}

diversity_by_definition <- purrr::map_dfr(definition_modes, calculate_summary_for_definition)
if (!all(diversity_summary$ClonotypeMode %in% diversity_by_definition$ClonotypeMode)) {
  diversity_by_definition <- bind_rows(diversity_by_definition, diversity_summary)
}

dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(diversity_summary, file.path(args$output_dir, "diversity_summary.rds"))
saveRDS(diversity_by_definition, file.path(args$output_dir, "diversity_by_clonotype_definition.rds"))
saveRDS(rarefaction, file.path(args$output_dir, "rarefaction.rds"))
saveRDS(rarefaction.curve, file.path(args$output_dir, "rarefaction_curve.rds"))

method_lines <- c(
  paste0("Clonotype mode: ", vegan.data$clonotype_mode),
  paste0("Clonotype definition: ", vegan.data$clonotype_definition),
  paste0("Rarefaction run: ", run_rarefaction),
  paste0("Rarefaction depth: ", ifelse(is.na(rarefaction_depth), "not calculated", rarefaction_depth)),
  "Preset definitions: gene = V/D/J/C genes when present; nt = CDR3 nucleotide; aa = CDR3 amino acid; strict = CDR3 nucleotide plus V/D/J/C genes when present.",
  "Rarefaction depth is the minimum sample library size, calculated as rowSums of the clonotype abundance matrix.",
  "RichnessPerMillion is retained only as a rough descriptive normalization because observed richness is nonlinear with sequencing depth. Prefer rarefied richness or models that account for depth for cross-sample comparisons."
)
writeLines(method_lines, file.path(args$output_dir, "DiversityMethod.txt"))
