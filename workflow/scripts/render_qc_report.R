suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description = "Render the tcr-seek diversity QC Quarto report.")
parser$add_argument("--qmd", required = TRUE, help = "Path to qc_report.qmd")
parser$add_argument("--output-dir", required = TRUE, help = "Directory for rendered HTML")
parser$add_argument("--immdata", required = TRUE, help = "Path to Immdata.rds")
parser$add_argument("--diversity", required = TRUE, help = "Path to diversity_summary.rds")
parser$add_argument("--rarefaction-curve", default = "", help = "Path to rarefaction_curve.rds")
parser$add_argument("--run-rarefaction", default = "true", help = "Whether the report should include rarefaction curves")
parser$add_argument("--qc-variables", default = "", help = "Comma-separated metadata variables for report tabs")
parser$add_argument("--logo", default = "", help = "Optional logo path")
args <- parser$parse_args()

parse_bool <- function(value) {
  normalized <- tolower(trimws(as.character(value)))
  normalized %in% c("true", "t", "yes", "y", "1")
}
run_rarefaction <- parse_bool(args$run_rarefaction)

quarto <- Sys.which("quarto")
if (!nzchar(quarto)) {
  stop("The quarto command-line executable is required to render the QC report")
}

output_dir <- args$output_dir
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
render_qmd <- file.path(output_dir, "qc_report.qmd")
file.copy(args$qmd, render_qmd, overwrite = TRUE)
if (nzchar(args$logo) && file.exists(args$logo)) {
  file.copy(args$logo, file.path(output_dir, "tcr-seek.svg"), overwrite = TRUE)
}

render_args <- c(
  "render", render_qmd,
  "--output-dir", output_dir,
  "-P", paste0("immdata_rds:", normalizePath(args$immdata, mustWork = TRUE)),
  "-P", paste0("diversity_rds:", normalizePath(args$diversity, mustWork = TRUE)),
  "-P", paste0("qc_variables:", args$qc_variables),
  "-P", "logo:tcr-seek.svg"
)
if (run_rarefaction) {
  render_args <- c(render_args, "-P", paste0("rarefaction_curve_rds:", normalizePath(args$rarefaction_curve, mustWork = TRUE)))
}

status <- system2(quarto, args = render_args)
if (!identical(status, 0L)) {
  stop("quarto render failed with exit status ", status)
}
