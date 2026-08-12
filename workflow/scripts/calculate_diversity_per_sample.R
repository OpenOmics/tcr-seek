suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(vegan)
})

parser <- ArgumentParser(description = "Calculate repertoire diversity metrics per sample without building a global sample-by-clonotype matrix.")
parser$add_argument("--immdata", required = TRUE, help = "Path to Immdata.rds produced by immunarch::repLoad().")
parser$add_argument("--output-dir", required = TRUE, help = "Directory where diversity_summary.rds and diversity_by_clonotype_definition.rds will be written.")
parser$add_argument("--clonotype", default = "strict", help = "Clonotype definition preset or comma-separated immunarch columns. Presets: gene, nt, aa, strict. Default: strict")
parser$add_argument("--run-rarefaction", default = "true", help = "Whether to calculate rarefied richness. Accepted values: true/false. Default: true")
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

collapse_sample_counts <- function(sample_name, repertoire, clonotype_columns) {
  missing_cols <- setdiff(c(clonotype_columns, "Clones"), names(repertoire))
  if (length(missing_cols) > 0) {
    stop("Missing expected columns in sample ", sample_name, ": ", paste(missing_cols, collapse = ", "))
  }

  collapsed <- repertoire %>%
    dplyr::select(all_of(c(clonotype_columns, "Clones"))) %>%
    dplyr::filter(dplyr::if_all(all_of(clonotype_columns), ~ !is.na(.x) & .x != "")) %>%
    tidyr::unite("Clonotype", all_of(clonotype_columns), remove = FALSE) %>%
    dplyr::group_by(Clonotype) %>%
    dplyr::summarise(count = sum(Clones), .groups = "drop")

  if (!nrow(collapsed)) {
    stop("Cannot calculate diversity for sample with zero clonotype abundance: ", sample_name)
  }
  collapsed$count
}

summarise_counts <- function(sample_name, counts, clonotype_mode, clonotype_definition, rarefaction_depth, run_rarefaction) {
  library_size <- sum(counts)
  richness <- vegan::specnumber(counts)
  shannon <- vegan::diversity(counts, index = "shannon")
  simpson <- vegan::diversity(counts, index = "simpson")
  evenness <- shannon / log(richness)
  top_counts <- sort(counts, decreasing = TRUE)
  rarefaction <- if (run_rarefaction) {
    as.numeric(vegan::rarefy(counts, sample = rarefaction_depth, se = FALSE))
  } else {
    NA_real_
  }

  data.frame(
    Sample = sample_name,
    Richness = richness,
    Shannon = shannon,
    Simpson = simpson,
    Evenness = evenness,
    LibrarySize = library_size,
    RichnessPerMillion = richness / (library_size / 1e6),
    EffectiveClonotypes = exp(shannon),
    TopCloneFraction = top_counts[1] / library_size,
    Top10CloneFraction = sum(top_counts[seq_len(min(10, length(top_counts)))]) / library_size,
    Rarefaction = rarefaction,
    RarefactionDepth = rarefaction_depth,
    ClonotypeMode = clonotype_mode,
    ClonotypeDefinition = clonotype_definition,
    stringsAsFactors = FALSE
  )
}

make_rarefaction_curve <- function(counts, rarefaction_depth, step = 100) {
  library_size <- sum(counts)
  depths <- unique(sort(c(seq.int(1, library_size, by = step), rarefaction_depth, library_size)))
  depths <- depths[!is.na(depths) & depths >= 1 & depths <= library_size]
  values <- vapply(depths, function(depth) {
    as.numeric(vegan::rarefy(counts, sample = depth, se = FALSE))
  }, numeric(1))
  stats::setNames(values, paste0("N", depths))
}

calculate_counts_for_definition <- function(immdata, resolved) {
  samples <- names(immdata$data)
  names(samples) <- samples
  purrr::map(samples, function(sample_name) {
    collapse_sample_counts(sample_name, immdata$data[[sample_name]], resolved$columns)
  })
}

calculate_summary_for_definition <- function(immdata, clonotype_value, run_rarefaction) {
  resolved <- resolve_clonotype_definition(clonotype_value, first_repertoire_columns(immdata))
  samples <- names(immdata$data)
  names(samples) <- samples
  sample_counts <- lapply(samples, function(sample_name) {
    collapse_sample_counts(sample_name, immdata$data[[sample_name]], resolved$columns)
  })
  depths <- vapply(sample_counts, sum, numeric(1))
  if (any(depths <= 0)) {
    stop("Cannot calculate diversity for samples with zero clonotype abundance: ", paste(names(depths)[depths <= 0], collapse = ", "))
  }
  rarefaction_depth <- if (run_rarefaction) min(depths) else NA_real_
  clonotype_definition <- paste(resolved$columns, collapse = " + ")

  purrr::map2_dfr(
    names(sample_counts),
    sample_counts,
    ~ summarise_counts(.x, .y, resolved$mode, clonotype_definition, rarefaction_depth, run_rarefaction)
  )
}

immdata <- readRDS(args$immdata)
diversity_summary <- calculate_summary_for_definition(immdata, args$clonotype, run_rarefaction)
definition_modes <- c("gene", "nt", "aa", "strict")
diversity_by_definition <- purrr::map_dfr(definition_modes, ~ calculate_summary_for_definition(immdata, .x, run_rarefaction))
if (!all(diversity_summary$ClonotypeMode %in% diversity_by_definition$ClonotypeMode)) {
  diversity_by_definition <- bind_rows(diversity_by_definition, diversity_summary)
}

rarefaction <- stats::setNames(diversity_summary$Rarefaction, diversity_summary$Sample)
if (run_rarefaction) {
  primary_resolved <- resolve_clonotype_definition(args$clonotype, first_repertoire_columns(immdata))
  primary_counts <- calculate_counts_for_definition(immdata, primary_resolved)
  rarefaction_depth <- diversity_summary$RarefactionDepth[1]
  rarefaction.curve <- purrr::map(primary_counts, ~ make_rarefaction_curve(.x, rarefaction_depth))
} else {
  rarefaction.curve <- list()
}

dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(diversity_summary, file.path(args$output_dir, "diversity_summary.rds"))
saveRDS(diversity_by_definition, file.path(args$output_dir, "diversity_by_clonotype_definition.rds"))
saveRDS(rarefaction, file.path(args$output_dir, "rarefaction.rds"))
saveRDS(rarefaction.curve, file.path(args$output_dir, "rarefaction_curve.rds"))

method_lines <- c(
  paste0("Clonotype mode: ", diversity_summary$ClonotypeMode[1]),
  paste0("Clonotype definition: ", diversity_summary$ClonotypeDefinition[1]),
  paste0("Rarefaction run: ", run_rarefaction),
  paste0("Rarefaction depth: ", ifelse(is.na(diversity_summary$RarefactionDepth[1]), "not calculated", diversity_summary$RarefactionDepth[1])),
  "Diversity metrics were calculated independently per sample without constructing a global sample-by-clonotype matrix."
)
writeLines(method_lines, file.path(args$output_dir, "DiversityMethod.txt"))
