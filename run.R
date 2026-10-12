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
        "  --gene-file PATH\n",
        "  --interaction-network-file PATH\n",
        "  --tx2gene-file PATH\n",
        "  --gsea-files PATH...         One or more GSEA reports\n",
        "  --gmt-files PATH...          One or more GMT files\n",
        "  --gmt-list-file PATH         GMT metadata CSV\n\n",
        "Optional data:\n",
        "  --clinreport-dir PATH         Drug annotations; omitted means unavailable\n",
        "  --gene-annotations-file PATH  Gene descriptions (MyGene.info); optional\n",
        "  --wes-vcf-file PATH           Variant evidence; omitted means unavailable\n",
        "  --tf-file PATH                 TF evidence; omitted means unavailable\n",
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
        "  --report-top-n INTEGER         Default: 15\n",
        "  --report-module-top-n INTEGER Default: 5\n\n",
        "  --report-language en|ru        Default: en\n\n",
        "  --use-tf-activity TRUE|FALSE  Include TF activity in Drug module prioritisation (default: FALSE)\n",
        "  --use-wes TRUE|FALSE          Include WES evidence in Drug module prioritisation (default: FALSE)\n",
        "  --skip-network-processing TRUE|FALSE  Disable network expansion (default: FALSE)\n\n",
        "Launch-only options:\n",
        "  --host ADDRESS                 Default: 127.0.0.1\n",
        "  --port NUMBER                  Default: 3838\n\n",
        "Paths are relative to the calling working directory.\n",
        "--prepare copies code and data and embeds the supplied settings; --check validates without starting Shiny.\n")
    return(invisible(NULL))
  }
  controls <- c("prepare", "check", "generate_default_pdf", "host", "port")
  unknown <- setdiff(names(cli), c(controls, clinreport_keys))
  if (length(unknown)) stop("Unknown option(s): ", paste(unknown, collapse = ", "))
  if (!is.null(cli$prepare) && isTRUE(cli$check)) stop("Choose --prepare or --check")
  values <- list()
  for (key in intersect(names(cli), clinreport_keys)) {
    value <- cli[[key]]
    if (key %in% clinreport_numbers) value <- suppressWarnings(as.numeric(value))
    if (key %in% clinreport_paths || key == "gmt_list_file") value <- clinreport_resolve(value, caller_dir)
    if (key %in% clinreport_list_paths) value <- vapply(value, clinreport_resolve, character(1), base = caller_dir)
    values[[key]] <- value
  }
  setwd(app_dir)
  on.exit(setwd(caller_dir), add = TRUE)

  if (isTRUE(cli$generate_default_pdf)) {
    source("R/config/config.R", local = .GlobalEnv, encoding = "UTF-8")
    source("R/report/report_pdf.R", local = .GlobalEnv, encoding = "UTF-8")
    required_report_packages <- c("shiny", "DT", "jsonlite", "dplyr", "purrr", "stringr", "htmltools", "plotly", "visNetwork")
    missing_report_packages <- required_report_packages[!vapply(required_report_packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing_report_packages)) {
      generate_clinreport_fallback_pdf(
        file.path(app_dir, "default_report.pdf"),
        patient = patient_id,
        reason = paste("missing R packages:", paste(missing_report_packages, collapse = ", ")),
        top_n = report_top_n,
        module_top_n = report_module_top_n,
        language = report_language,
        state_json_path = file.path(app_dir, "default_selection.json")
      )
    } else {
      source("app.R", encoding = "UTF-8")
      default_ids <- head(report_default_order(drug_network), as.integer(report_top_n))
      generate_clinreport_pdf(
        file.path(app_dir, "default_report.pdf"),
        final_ids = default_ids,
        baseline_ids = default_ids,
        top_n = length(default_ids),
        module_top_n = min(as.integer(report_module_top_n), length(default_ids)),
        language = report_language,
        state_json_path = file.path(app_dir, "default_selection.json")
      )
    }
    cat("Generated default PDF: ", file.path(app_dir, "default_report.pdf"), "\n", sep = "")
    return(invisible(NULL))
  }

  env <- new.env(parent = globalenv())
  source("R/config/config.R", local = env, encoding = "UTF-8")
  if (!is.null(cli$prepare) || isTRUE(cli$check)) {
    # In a pipeline, an omitted optional analysis must mean "not supplied",
    # not "use the local development fixture from raw/".
    for (key in c("clinreport_dir", "gene_annotations_file", "wes_vcf_file", "tf_file")) {
      if (!key %in% names(values)) values[[key]] <- ""
    }
  }
  for (key in names(values)) assign(key, values[[key]], envir = env)
  if (!is.null(cli$prepare) || isTRUE(cli$check)) {
    required <- c("drug_network_file", "gene_file", "interaction_network_file", "tx2gene_file",
      "gsea_files", "gmt_files", "gmt_list_file")
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
    for (key in c(clinreport_paths, clinreport_list_paths, "gmt_list_file")) {
      path <- get(key, env)
      if (length(path) && any(vapply(path, function(item) {
        nzchar(item) && dir.exists(item) && startsWith(tolower(paste0(output, "/")),
          tolower(paste0(clinreport_resolve(item, app_dir), "/")))
      }, logical(1)))) stop("Output is inside an input directory")
    }
    if (!dir.create(output, recursive = TRUE)) stop("Cannot create output: ", output)
    copy_checked <- function(from, to, recursive = FALSE) {
      if (!all(file.copy(from, to, recursive = recursive, copy.mode = TRUE))) stop("Failed to copy: ", from)
    }
    copy_checked(file.path(app_dir, c("app.R", "run.R", "R", "www")), output, TRUE)
    bundled <- mget(clinreport_keys, envir = env)
    effective_command_args <- c(
      commandArgs(TRUE),
      paste0("--use-tf-activity=", get("use_tf_activity", env)),
      paste0("--use-wes=", get("use_wes", env)),
      paste0("--skip-network-processing=", get("skip_network_processing", env)),
      paste0("--report-top-n=", get("report_top_n", env)),
      paste0("--report-module-top-n=", get("report_module_top_n", env)),
      paste0("--report-language=", get("report_language", env))
    )
    bundled$analysis_command <- paste(
      c("clinreport", effective_command_args),
      collapse = " "
    )
    for (key in c(clinreport_paths, clinreport_list_paths, "gmt_list_file")) {
      path <- bundled[[key]]
      if (is.null(path) || !length(path) || (length(path) == 1 && !nzchar(path))) {
        bundled[[key]] <- ""
        next
      }
      if (length(path) > 1) {
        target_dir <- file.path(output, "data", key)
        dir.create(target_dir, recursive = TRUE)
        for (item in path) copy_checked(item, file.path(target_dir, basename(item)), dir.exists(item))
        bundled[[key]] <- gsub("\\\\", "/", file.path("data", key, basename(path)))
        next
      }
      if (!file.exists(path)) { bundled[[key]] <- ""; next }
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
      ",\n  analysis_command = ", paste(deparse(bundled$analysis_command, control = "all"), collapse = ""),
      "))"
    ), prepared_file, useBytes = TRUE)
    writeLines(c("Launch locally in an R environment containing ClinReport's dependencies:",
      "Rscript run.R", "", "Apptainer, from this directory:",
      'apptainer exec --bind "$PWD:/report" IMAGE.sif clinreport --host 127.0.0.1',
      "", "Open http://127.0.0.1:3838 on the same machine.",
      "DoRothEA needs internet unless dorothea_file was supplied during preparation."), file.path(output, "README.txt"))

    rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
    generation_output <- system2(
      rscript,
      c(shQuote(file.path(output, "run.R")), "--generate-default-pdf"),
      stdout = TRUE,
      stderr = TRUE
    )
    if (!is.null(attr(generation_output, "status")) && attr(generation_output, "status") != 0) {
      stop("Default PDF generation failed: ", paste(generation_output, collapse = "\n"))
    }
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
