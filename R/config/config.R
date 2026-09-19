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

# Protein/gene interaction network used for the drug subgraph visualization.
string_network_file <- file.path(
  "C:/projects/clinreport/raw/network",
  "string.human_links_v12_0_min900.Ensembl.edges.tsv"
)

# Transcript -> gene mapping. Used to enrich/normalize gene identifiers and names.
tx2gene_file <- file.path(
  "C:/projects/clinreport/raw",
  "tx2gene.tsv"
)

# Hallmark GSEA inputs. These remain explicit configuration values so a later
# parameterized run can replace them without changing UI or server logic.
gsea_carcinoma_report_file <- file.path(
  "C:/projects/clinreport/raw/gsea",
  "carcinoma_vs_normal_h_all_v2026_1_Hs_symbols_gsea_report_for_carcinoma.tsv"
)
gsea_normal_report_file <- file.path(
  "C:/projects/clinreport/raw/gsea",
  "carcinoma_vs_normal_h_all_v2026_1_Hs_symbols_gsea_report_for_normal.tsv"
)
hallmark_gmt_file <- file.path(
  "C:/projects/clinreport/raw/gsea",
  "carcinoma_vs_normal_h_all_v2026_1_Hs_symbols_h_all_v2026_1_Hs_symbols.gmt"
)

# Maximum number of STRING neighbors displayed per drug target.
network_max_neighbors <- 40

# Display name for the selected interaction network. Keep this as a distinct
# input so a future dataset metadata loader can provide it at startup.
network_display_name <- "Drug–gene interaction network"

padj_cutoff <- 0.05
log2fc_cutoff <- 1.0
