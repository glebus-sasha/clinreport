# Configuration and input validation

# 2. VALIDATION
# ============================================================

if (!file.exists(drug_network_file)) {
  stop("Drug network file does not exist: ", drug_network_file)
}

if (!dir.exists(clinreport_dir)) {
  stop("ClinReport directory does not exist: ", clinreport_dir)
}

if (!file.exists(gene_file)) {
  stop("Gene expression file does not exist: ", gene_file)
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
  log2fc_cutoff < 0
) {
  stop("log2fc_cutoff must be a single non-negative number.")
}



# Network visualization inputs
if (!file.exists(string_network_file)) {
  stop("STRING edge file not found: ", string_network_file)
}

if (!file.exists(tx2gene_file)) {
  stop("tx2gene file not found: ", tx2gene_file)
}

if (!is.numeric(network_max_neighbors) || length(network_max_neighbors) != 1 ||
    is.na(network_max_neighbors) || network_max_neighbors < 1) {
  stop("network_max_neighbors must be a positive number")
}
