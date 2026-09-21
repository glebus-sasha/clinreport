# Application configuration

# Local demo defaults. Prepared reports override these with their supplied inputs.

# ============================================================
# 1. INPUT PARAMETERS
# ============================================================

drug_network_file <- file.path(
  "raw",
  "all_samples_string_human_links_v12_0_min900_Ensembl_diamond_trustrank.csv"
)

clinreport_dir <- file.path(
  "raw",
  "clinreport"
)

gene_file <- file.path(
  "raw",
  "carcinoma_vs_normal_gene_names_added.tsv"
)

# Whole-exome sequencing variants annotated by PCGR for the current patient.
wes_vcf_file <- file.path(
  "raw",
  "R_PTA_22.pcgr.grch38.pass.vcf.gz"
)
wes_display_name <- "WES · PCGR"

# Protein/gene interaction network used for the drug subgraph visualization.
interaction_network_file <- file.path(
  "raw/network",
  "string.human_links_v12_0_min900.Ensembl.edges.tsv"
)

# Transcript -> gene mapping. Used to enrich/normalize gene identifiers and names.
tx2gene_file <- file.path(
  "raw",
  "tx2gene.tsv"
)

# Hallmark GSEA inputs. These remain explicit configuration values so a later
# parameterized run can replace them without changing UI or server logic.
gsea_carcinoma_report_file <- file.path(
  "raw/gsea",
  "carcinoma_vs_normal_h_all_v2026_1_Hs_symbols_gsea_report_for_carcinoma.tsv"
)
gsea_normal_report_file <- file.path(
  "raw/gsea",
  "carcinoma_vs_normal_h_all_v2026_1_Hs_symbols_gsea_report_for_normal.tsv"
)
pathway_gmt_file <- file.path(
  "raw/gsea",
  "carcinoma_vs_normal_h_all_v2026_1_Hs_symbols_h_all_v2026_1_Hs_symbols.gmt"
)

# Maximum number of interaction-network neighbours displayed per drug target.
interaction_network_max_neighbors <- 40

# Display metadata. These are explicit inputs so a later dataset metadata
# loader can replace the patient, interaction network or GMT collection
# without changing presentation code.
patient_id <- "PATIENT-001"
network_display_name <- "Drug–gene interaction network"
interaction_network_name <- "STRING"
pathway_collection_name <- "Hallmark"
pathway_id_prefix <- "HALLMARK_"

display_pathway_name <- function(pathways) {
  pathways <- as.character(pathways)
  has_prefix <- startsWith(pathways, pathway_id_prefix)
  ifelse(
    has_prefix,
    substr(pathways, nchar(pathway_id_prefix) + 1, nchar(pathways)),
    pathways
  )
}

padj_cutoff <- 0.05
log2fc_cutoff <- 1.0

tf_file <- file.path("raw", "carcinoma_vs_normal_significant_tfs.tsv")
dorothea_file <- ""

# A prepared application embeds the exact paths and metadata supplied to
# `clinreport --prepare`; no external configuration file is required.
source("R/config/parameters.R", local = TRUE, encoding = "UTF-8")
prepared_parameters_file <- "R/config/prepared_parameters.R"
if (file.exists(prepared_parameters_file)) {
  source(prepared_parameters_file, local = TRUE, encoding = "UTF-8")
  list2env(clinreport_prepared_parameters, envir = environment())
}
