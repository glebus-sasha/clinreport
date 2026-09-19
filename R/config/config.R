# Application configuration

# Keep machine-specific paths here. For a shared deployment, move these
# values to environment variables or a config file.

# ============================================================
# 1. INPUT PARAMETERS
# ============================================================

drug_network_file <- file.path(
  "C:/projects/clinreport/raw",
  "all_samples_string_human_links_v12_0_min900_Ensembl_diamond_trustrank.csv"
)

clinreport_dir <- file.path(
  "C:/projects/clinreport/raw",
  "clinreport"
)

gene_file <- file.path(
  "C:/projects/clinreport/raw",
  "carcinoma_vs_normal_gene_names_added.tsv"
)

padj_cutoff <- 0.05
log2fc_cutoff <- 1.0

