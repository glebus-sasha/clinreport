# Shared input contract for launch and preparation.
clinreport_paths <- c("drug_network_file", "clinreport_dir", "gene_file", "gene_annotations_file",
  "wes_vcf_file", "interaction_network_file", "tx2gene_file",
  "gsea_carcinoma_report_file", "gsea_normal_report_file", "pathway_gmt_file",
  "tf_file", "dorothea_file")
clinreport_list_paths <- c("gsea_files", "gmt_files")
clinreport_numbers <- c("padj_cutoff", "log2fc_cutoff", "interaction_network_max_neighbors",
  "report_top_n", "report_module_top_n")
clinreport_strings <- c("patient_id", "wes_display_name", "network_display_name",
  "interaction_network_name", "pathway_collection_name", "pathway_id_prefix",
  "use_tf_activity", "use_wes", "skip_network_processing", "report_language")
clinreport_keys <- c(clinreport_paths, clinreport_list_paths, clinreport_numbers, clinreport_strings,
  "gmt_list_file")

clinreport_resolve <- function(path, base) {
  if (!nzchar(path)) return(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)) path <- file.path(base, path)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

clinreport_validate_values <- function(values) {
  if (!is.list(values) || (length(values) &&
      (is.null(names(values)) || anyDuplicated(names(values))))) stop("Parameters must have unique names")
  unknown <- setdiff(names(values), clinreport_keys)
  if (length(unknown)) stop("Unknown parameter(s): ", paste(unknown, collapse = ", "))
  for (key in names(values)) {
    x <- values[[key]]
    if (key %in% clinreport_list_paths) {
      if (!is.character(x) || !length(x) || anyNA(x) || any(!nzchar(x))) {
        stop("Expected one or more file paths for ", key)
      }
      next
    }
    if (length(x) != 1 || is.na(x)) stop("Expected a scalar for ", key)
    if (key %in% clinreport_numbers) {
      if (!is.numeric(x) || !is.finite(x)) stop("Expected a finite number for ", key)
    } else if (!is.character(x)) stop("Expected a string for ", key)
  }
  values
}

clinreport_parse_args <- function(args) {
  values <- list()
  while (length(args)) {
    if (!startsWith(args[1], "--")) stop("Expected --option, got: ", args[1])
    option <- substring(args[1], 3)
    has_equals <- grepl("=", option, fixed = TRUE)
    key <- gsub("-", "_", sub("=.*$", "", option), fixed = TRUE)
    if (key %in% c("help", "check", "generate_default_pdf")) {
      if (has_equals) stop("Flag does not take a value: --", key)
      values[[key]] <- TRUE
      args <- args[-1]
    } else if (has_equals) {
      values[[key]] <- substring(option, regexpr("=", option, fixed = TRUE) + 1)
      args <- args[-1]
    } else if (key %in% c("gsea_files", "gmt_files")) {
      next_option <- which(startsWith(args[-1], "--"))[1]
      if (!length(next_option) || is.na(next_option)) next_option <- length(args)
      value_count <- next_option - 1L
      values[[key]] <- args[seq.int(2L, value_count + 1L)]
      if (!length(values[[key]])) stop("Missing value for --", key)
      args <- args[-seq_len(value_count + 1L)]
    } else {
      if (length(args) < 2 || startsWith(args[2], "--")) stop("Missing value for --", key)
      values[[key]] <- args[2]
      args <- args[-c(1, 2)]
    }
  }
  values
}
