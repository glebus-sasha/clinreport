# Configuration and input validation

# 2. VALIDATION
# ============================================================

required_file_inputs <- setdiff(clinreport_paths, c("clinreport_dir", "wes_vcf_file", "tf_file", "dorothea_file", "gene_annotations_file"))
for (input_key in required_file_inputs) {
  input_path <- get(input_key)
  if (nzchar(input_path) && dir.exists(input_path)) stop("Expected a file for ", input_key, ": ", input_path)
}

if (!file.exists(drug_network_file)) {
  stop("Drug network file does not exist: ", drug_network_file)
}

gsea_input_files <- if (exists("gsea_files") && length(gsea_files)) gsea_files else
  c(gsea_carcinoma_report_file, gsea_normal_report_file)
gmt_input_files <- if (exists("gmt_files") && length(gmt_files)) gmt_files else pathway_gmt_file
gsea_input_files <- c(gsea_input_files, gmt_input_files, if (exists("gmt_list_file")) gmt_list_file else character())

missing_gsea_files <- gsea_input_files[!file.exists(gsea_input_files)]
if (length(missing_gsea_files) > 0) {
  stop("Missing GSEA input files: ", paste(missing_gsea_files, collapse = ", "))
}
if (exists("gmt_list_file") && nzchar(gmt_list_file)) {
  collection_list <- read.csv(gmt_list_file, stringsAsFactors = FALSE)
  if (!all(c("name", "url") %in% names(collection_list)) || !nrow(collection_list)) {
    stop("GMT collection list must contain name and url columns")
  }
  listed_gmts <- basename(sub("[?#].*$", "", collection_list$url))
  if (anyNA(collection_list$name) || any(!nzchar(collection_list$name)) ||
      anyDuplicated(collection_list$name) || anyDuplicated(listed_gmts) ||
      !setequal(listed_gmts, basename(gmt_input_files)) ||
      length(listed_gmts) != length(gmt_input_files)) {
    stop("GMT files do not match the collection list")
  }
  report_names <- basename(if (exists("gsea_files") && length(gsea_files)) gsea_files else
    c(gsea_carcinoma_report_file, gsea_normal_report_file))
  if (anyDuplicated(report_names)) stop("GSEA report filenames must be unique")
  stems <- tools::file_path_sans_ext(listed_gmts)
  report_matches <- vapply(report_names, function(report) sum(vapply(stems, function(stem) {
    grepl(paste0(stem, ".gsea_report_for_"), report, fixed = TRUE)
  }, logical(1))), integer(1))
  if (any(report_matches != 1L) ||
      any(!grepl("\\.gsea_report_for_[[:alnum:]_-]+\\.tsv$", report_names, ignore.case = TRUE))) {
    stop("Each GSEA report must match exactly one GMT collection and name its group")
  }
  for (stem in stems) {
    if (!any(grepl(paste0(stem, ".gsea_report_for_"), report_names, fixed = TRUE))) {
      stop("No GSEA report for GMT collection: ", stem)
    }
  }
}

has_drug_details <- nzchar(clinreport_dir) && dir.exists(clinreport_dir)
if (nzchar(clinreport_dir) && !has_drug_details) {
  stop("ClinReport directory does not exist: ", clinreport_dir)
}

if (!file.exists(gene_file)) {
  stop("Gene expression file does not exist: ", gene_file)
}
if (nzchar(gene_annotations_file) && !file.exists(gene_annotations_file)) {
  stop("Gene annotations file does not exist: ", gene_annotations_file)
}
if (nzchar(gene_annotations_file) && dir.exists(gene_annotations_file)) {
  stop("Expected a JSON file for gene_annotations_file: ", gene_annotations_file)
}

has_wes <- nzchar(wes_vcf_file) && file.exists(wes_vcf_file)
if (nzchar(wes_vcf_file) && !has_wes) {
  stop("WES VCF file does not exist: ", wes_vcf_file)
}

if (
  !is.numeric(padj_cutoff) ||
  length(padj_cutoff) != 1 ||
  is.na(padj_cutoff) ||
  padj_cutoff <= 0 ||
  padj_cutoff >= 1
) {
  stop("padj_cutoff must be a single number between 0 and 1.")
}

if (
  !is.numeric(log2fc_cutoff) ||
  length(log2fc_cutoff) != 1 ||
  is.na(log2fc_cutoff) ||
  !is.finite(log2fc_cutoff) ||
  log2fc_cutoff < 0
) {
  stop("log2fc_cutoff must be a single non-negative number.")
}



# Network visualization inputs
if (nzchar(dorothea_file) && !file.exists(dorothea_file)) {
  stop("DoRothEA resource not found: ", dorothea_file)
}
if (nzchar(tf_file) && !file.exists(tf_file)) {
  stop("TF input not found: ", tf_file, "; set tf_file to an empty string to disable it")
}
if (!file.exists(interaction_network_file)) {
  stop(interaction_network_name, " edge file not found: ", interaction_network_file)
}

if (!file.exists(tx2gene_file)) {
  stop("tx2gene file not found: ", tx2gene_file)
}

if (!is.numeric(interaction_network_max_neighbors) || length(interaction_network_max_neighbors) != 1 ||
    is.na(interaction_network_max_neighbors) || !is.finite(interaction_network_max_neighbors) ||
    interaction_network_max_neighbors != floor(interaction_network_max_neighbors) || interaction_network_max_neighbors < 1) {
  stop("interaction_network_max_neighbors must be a positive number")
}

if (!is.numeric(report_top_n) || length(report_top_n) != 1 || is.na(report_top_n) ||
    !is.finite(report_top_n) || report_top_n != floor(report_top_n) || report_top_n < 1) {
  stop("report_top_n must be a positive integer")
}
if (!is.numeric(report_graph_top_n) || length(report_graph_top_n) != 1 || is.na(report_graph_top_n) ||
    !is.finite(report_graph_top_n) || report_graph_top_n != floor(report_graph_top_n) || report_graph_top_n < 0) {
  stop("report_graph_top_n must be a non-negative integer")
}
if (!is.character(report_language) || length(report_language) != 1 ||
    !report_language %in% c("en", "ru")) {
  stop("report_language must be either 'en' or 'ru'")
}
