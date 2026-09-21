# No network or Shiny server required. Test the exported artifact after relocation.
root <- normalizePath(getwd(), winslash = "/")
rscript <- file.path(R.home("bin"), "Rscript")
if (.Platform$OS.type == "windows") rscript <- paste0(rscript, ".exe")
scratch <- tempfile("clinreport cli ")
dir.create(scratch)
source("R/config/parameters.R")
invoke <- function(script, args, expected = 0L) {
  output <- suppressWarnings(system2(rscript, c("--vanilla", shQuote(script), shQuote(args)), stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  if (status != expected) stop(paste(output, collapse = "\n"))
  invisible(output)
}
dir.create(file.path(scratch, "inputs"))
writeLines("fixture", file.path(scratch, "inputs", "same name.tsv"))
dir.create(file.path(scratch, "inputs", "drugs"))
dir.create(file.path(scratch, "inputs", "drugs", "nested"))
writeLines("{}", file.path(scratch, "inputs", "drugs", "nested", "drug.json"))
input_file <- normalizePath(file.path(scratch, "inputs", "same name.tsv"), winslash = "/")
input_drugs <- normalizePath(file.path(scratch, "inputs", "drugs"), winslash = "/")
input_args <- c(
  "--drug-network-file", input_file,
  "--clinreport-dir", input_drugs,
  "--gene-file", input_file,
  "--wes-vcf-file", input_file,
  "--interaction-network-file", input_file,
  "--tx2gene-file", input_file,
  "--gsea-carcinoma-report-file", input_file,
  "--gsea-normal-report-file", input_file,
  "--pathway-gmt-file", input_file,
  "--tf-file=", "--dorothea-file=",
  "--patient-id", "CLI patient",
  "--padj-cutoff", "0.00000123456789",
  "--network-display-name", "Drug–gene сеть"
)
bundle <- file.path(scratch, "bundle")
invoke(file.path(root, "run.R"), c("--prepare", bundle, input_args))
invoke(file.path(root, "run.R"), c("--prepare", bundle, input_args), 1L)
invoke(file.path(root, "run.R"), c("--check", input_args))
moved <- file.path(scratch, "moved bundle")
stopifnot(file.rename(bundle, moved))
# Move original inputs away: the export must not depend on them.
stopifnot(file.rename(file.path(scratch, "inputs"), file.path(scratch, "unavailable")))
source(file.path(moved, "R", "config", "prepared_parameters.R"))
resolved <- clinreport_prepared_parameters
stopifnot(identical(resolved$patient_id, "CLI patient"), resolved$padj_cutoff == 0.00000123456789,
  identical(resolved$network_display_name, "Drug–gene сеть"),
  identical(resolved$wes_display_name, "WES · PCGR"),
  file.exists(file.path(moved, resolved$clinreport_dir, "nested", "drug.json")),
  file.exists(file.path(moved, "www", "css", "app.css")))
invoke(file.path(root, "run.R"), c("--check", input_args, "--gene-file", "missing.tsv"), 1L)
invoke(file.path(root, "run.R"), c("--check", input_args, "--padj-cutoff", "NaN"), 1L)
invoke(file.path(root, "run.R"), c("--check", input_args, "--interaction-network-max-neighbors", "1.5"), 1L)
invoke(file.path(root, "run.R"), c("--check", input_args, "--unknown", "value"), 1L)
invoke(file.path(root, "run.R"), "--check", 1L)
cat("PASS: explicit CLI, relative paths, copied nested data, relocation, invalid inputs, overwrite protection\n")
