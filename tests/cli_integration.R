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
writeLines('{"genes":{"ENSG00000000003":{"name":"tetraspanin 6","summary":"Example summary"}}}',
  file.path(scratch, "inputs", "gene_annotations.json"))
for (name in c("case.first.gsea_report_for_carcinoma.tsv",
               "case.second.gsea_report_for_normal.tsv", "first.gmt", "second.gmt")) {
  writeLines("fixture", file.path(scratch, "inputs", name))
}
writeLines(c("name,url", "First,https://example.org/first.gmt", "Second,https://example.org/second.gmt"),
  file.path(scratch, "inputs", "gmt_list.csv"))
dir.create(file.path(scratch, "inputs", "drugs"))
dir.create(file.path(scratch, "inputs", "drugs", "nested"))
writeLines("{}", file.path(scratch, "inputs", "drugs", "nested", "drug.json"))
input_file <- normalizePath(file.path(scratch, "inputs", "same name.tsv"), winslash = "/")
input_gsea <- normalizePath(file.path(scratch, "inputs", c(
  "case.first.gsea_report_for_carcinoma.tsv",
  "case.second.gsea_report_for_normal.tsv")), winslash = "/")
input_gmts <- normalizePath(file.path(scratch, "inputs", c("first.gmt", "second.gmt")), winslash = "/")
input_gmt_list <- normalizePath(file.path(scratch, "inputs", "gmt_list.csv"), winslash = "/")
input_drugs <- normalizePath(file.path(scratch, "inputs", "drugs"), winslash = "/")
input_annotations <- normalizePath(file.path(scratch, "inputs", "gene_annotations.json"), winslash = "/")
required_args <- c(
  "--drug-network-file", input_file,
  "--gene-file", input_file,
  "--interaction-network-file", input_file,
  "--tx2gene-file", input_file,
  "--gsea-files", input_gsea,
  "--gmt-files", input_gmts,
  "--gmt-list-file", input_gmt_list
)
input_args <- c(
  required_args,
  "--clinreport-dir", input_drugs,
  "--gene-annotations-file", input_annotations,
  "--wes-vcf-file", input_file,
  "--tf-file=", "--dorothea-file=",
  "--patient-id", "CLI patient",
  "--padj-cutoff", "0.00000123456789",
  "--network-display-name", "Drug–gene сеть"
)
bundle <- file.path(scratch, "bundle")
invoke(file.path(root, "run.R"), c("--prepare", bundle, input_args))
invoke(file.path(root, "run.R"), c("--prepare", bundle, input_args), 1L)
invoke(file.path(root, "run.R"), c("--check", input_args))
core_bundle <- file.path(scratch, "core bundle")
invoke(file.path(root, "run.R"), c("--prepare", core_bundle, required_args))
source(file.path(core_bundle, "R", "config", "prepared_parameters.R"))
core_parameters <- clinreport_prepared_parameters
stopifnot(
  identical(core_parameters$clinreport_dir, ""),
  identical(core_parameters$gene_annotations_file, ""),
  identical(core_parameters$wes_vcf_file, ""),
  identical(core_parameters$tf_file, ""),
  length(core_parameters$gsea_files) == 2L,
  length(core_parameters$gmt_files) == 2L,
  all(file.exists(file.path(core_bundle, core_parameters$gsea_files))),
  all(file.exists(file.path(core_bundle, core_parameters$gmt_files)))
)
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
  file.exists(file.path(moved, resolved$gene_annotations_file)),
  file.exists(file.path(moved, "www", "css", "app.css")))
invoke(file.path(root, "run.R"), c("--check", input_args, "--gene-file", "missing.tsv"), 1L)
invoke(file.path(root, "run.R"), c("--check", input_args, "--padj-cutoff", "NaN"), 1L)
invoke(file.path(root, "run.R"), c("--check", input_args, "--interaction-network-max-neighbors", "1.5"), 1L)
invoke(file.path(root, "run.R"), c("--check", input_args, "--unknown", "value"), 1L)
writeLines(c("name,url", "Wrong,https://example.org/missing.gmt"),
  file.path(scratch, "unavailable", "bad_gmt_list.csv"))
invoke(file.path(root, "run.R"), c("--check", input_args,
  "--gmt-list-file", file.path(scratch, "unavailable", "bad_gmt_list.csv")), 1L)
invoke(file.path(root, "run.R"), "--check", 1L)
cat("PASS: explicit CLI, optional-input omission, relative paths, copied nested data, relocation, invalid inputs, overwrite protection\n")
