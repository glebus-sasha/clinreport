# Configuration and input validation

# 2. VALIDATION
# ============================================================

for (input_key in setdiff(clinreport_paths, "clinreport_dir")) {
  input_path <- get(input_key)
  if (nzchar(input_path) && dir.exists(input_path)) stop("Expected a file for ", input_key, ": ", input_path)
}

if (!file.exists(drug_network_file)) {
  stop("Drug network file does not exist: ", drug_network_file)
}

gsea_input_files <- c(
  gsea_carcinoma_report_file,
  gsea_normal_report_file,
  pathway_gmt_file
)

missing_gsea_files <- gsea_input_files[!file.exists(gsea_input_files)]
if (length(missing_gsea_files) > 0) {
  stop("Missing GSEA input files: ", paste(missing_gsea_files, collapse = ", "))
}

if (!dir.exists(clinreport_dir)) {
  stop("ClinReport directory does not exist: ", clinreport_dir)
}

if (!file.exists(gene_file)) {
  stop("Gene expression file does not exist: ", gene_file)
}

if (!file.exists(wes_vcf_file)) {
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
