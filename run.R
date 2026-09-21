#!/usr/bin/env Rscript
# Launch or export ClinReport. CLI paths are relative to the calling directory.
main <- function() {
  if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
    Sys.setlocale("LC_CTYPE", "English_United States.utf8")
  }
  script <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
  app_dir <- dirname(normalizePath(sub("^--file=", "", script), winslash = "/", mustWork = TRUE))
  caller_dir <- getwd()
  source(file.path(app_dir, "R/config/parameters.R"), local = TRUE)
  cli <- clinreport_parse_args(commandArgs(TRUE))
  if (isTRUE(cli$help)) {
    cat("Usage: clinreport --prepare OUTPUT_DIR [parameters]\n",
        "       clinreport [--host 127.0.0.1] [--port 3838]\n",
        "       clinreport --check [parameters]\n\n",
        "Required with --prepare and --check:\n",
        "  --drug-network-file PATH\n",
        "  --clinreport-dir PATH\n",
        "  --gene-file PATH\n",
        "  --wes-vcf-file PATH\n",
        "  --interaction-network-file PATH\n",
        "  --tx2gene-file PATH\n",
        "  --gsea-carcinoma-report-file PATH\n",
        "  --gsea-normal-report-file PATH\n",
        "  --pathway-gmt-file PATH\n\n",
        "Optional data:\n",
        "  --tf-file PATH                 Default: raw/carcinoma_vs_normal_significant_tfs.tsv\n",
        "  --dorothea-file PATH           Default: download DoRothEA at application startup\n\n",
        "Optional metadata and thresholds:\n",
        "  --patient-id TEXT              Default: PATIENT-001\n",
        "  --wes-display-name TEXT        Default: WES · PCGR\n",
        "  --network-display-name TEXT    Default: Drug–gene interaction network\n",
        "  --interaction-network-name TEXT Default: STRING\n",
        "  --pathway-collection-name TEXT Default: Hallmark\n",
        "  --pathway-id-prefix TEXT       Default: HALLMARK_\n",
        "  --padj-cutoff NUMBER           Default: 0.05\n",
        "  --log2fc-cutoff NUMBER         Default: 1.0\n",
        "  --interaction-network-max-neighbors INTEGER  Default: 40\n\n",
        "Launch-only options:\n",
        "  --host ADDRESS                 Default: 127.0.0.1\n",
        "  --port NUMBER                  Default: 3838\n\n",
        "Paths are relative to the calling working directory.\n",
        "--prepare copies code and data and embeds the supplied settings; --check validates without starting Shiny.\n")
    return(invisible(NULL))
  }
  controls <- c("prepare", "check", "host", "port")
  unknown <- setdiff(names(cli), c(controls, clinreport_keys))
  if (length(unknown)) stop("Unknown option(s): ", paste(unknown, collapse = ", "))
  if (!is.null(cli$prepare) && isTRUE(cli$check)) stop("Choose --prepare or --check")
  values <- list()
  for (key in intersect(names(cli), clinreport_keys)) {
    value <- cli[[key]]
    if (key %in% clinreport_numbers) value <- suppressWarnings(as.numeric(value))
    if (key %in% clinreport_paths) value <- clinreport_resolve(value, caller_dir)
    values[[key]] <- value
  }
  setwd(app_dir)
  on.exit(setwd(caller_dir), add = TRUE)
  env <- new.env(parent = globalenv())
  source("R/config/config.R", local = env, encoding = "UTF-8")
  for (key in names(values)) assign(key, values[[key]], envir = env)
  if (!is.null(cli$prepare) || isTRUE(cli$check)) {
    required <- setdiff(clinreport_paths, c("tf_file", "dorothea_file"))
    missing <- setdiff(required, names(cli))
    if (length(missing)) stop("Explicit input required: ", paste(paste0("--", gsub("_", "-", missing)), collapse = ", "))
  }
  source("R/config/validation.R", local = env, encoding = "UTF-8")
  if (isTRUE(cli$check)) {
    cat("PASS: input paths and parameters\n")
    return(invisible(NULL))
  }
  if (!is.null(cli$prepare)) {
    output <- clinreport_resolve(cli$prepare, caller_dir)
    if (file.exists(output) || dir.exists(output)) stop("Output already exists: ", output)
    for (key in clinreport_paths) {
      path <- get(key, env)
      if (nzchar(path) && dir.exists(path) && startsWith(tolower(paste0(output, "/")),
          tolower(paste0(clinreport_resolve(path, app_dir), "/")))) stop("Output is inside an input directory")
    }
    if (!dir.create(output, recursive = TRUE)) stop("Cannot create output: ", output)
    copy_checked <- function(from, to, recursive = FALSE) {
      if (!all(file.copy(from, to, recursive = recursive, copy.mode = TRUE))) stop("Failed to copy: ", from)
    }
    copy_checked(file.path(app_dir, c("app.R", "run.R", "R", "www")), output, TRUE)
    bundled <- mget(clinreport_keys, envir = env)
    for (key in clinreport_paths) {
      path <- bundled[[key]]
      if (!nzchar(path) || !file.exists(path)) {
        bundled[[key]] <- ""
        next
      }
      relative_dir <- file.path("data", key)
      target_dir <- file.path(output, relative_dir)
      dir.create(target_dir, recursive = TRUE)
      copy_checked(path, target_dir, dir.exists(path))
      bundled[[key]] <- gsub("\\\\", "/", file.path(relative_dir, basename(path)))
    }
    prepared_file <- file.path(output, "R", "config", "prepared_parameters.R")
    value_lines <- vapply(names(bundled), function(key) {
      paste0("  ", key, " = ", paste(deparse(bundled[[key]], control = "all"), collapse = ""))
    }, character(1))
    writeLines(c(
      "# Generated by clinreport --prepare. Do not edit manually.",
      "clinreport_prepared_parameters <- structure(list(",
      paste(value_lines, collapse = ",\n"),
      "))"
    ), prepared_file, useBytes = TRUE)
    writeLines(c("Launch locally in an R environment containing ClinReport's dependencies:",
      "Rscript run.R", "", "Apptainer, from this directory:",
      'apptainer exec --bind "$PWD:/report" IMAGE.sif clinreport --host 127.0.0.1',
      "", "Open http://127.0.0.1:3838 on the same machine.",
      "DoRothEA needs internet unless dorothea_file was supplied during preparation."), file.path(output, "README.txt"))
    cat("Prepared ClinReport: ", output, "\n", sep = "")
    return(invisible(NULL))
  }
  host <- if (is.null(cli$host)) "127.0.0.1" else cli$host
  port <- if (is.null(cli$port)) 3838 else suppressWarnings(as.numeric(cli$port))
  if (length(port) != 1 || !is.finite(port) || port != floor(port) || port < 1 || port > 65535) stop("Invalid port")
  shiny::runApp(app_dir, host = host, port = as.integer(port), launch.browser = FALSE)
}
tryCatch(main(), error = function(e) {
  message("ClinReport: ", conditionMessage(e))
  quit(status = 1)
})
