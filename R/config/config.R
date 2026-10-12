# Application configuration

# Local demo defaults. Prepared reports override these with their supplied inputs.

# ============================================================
# 1. INPUT PARAMETERS
# ============================================================

drug_network_file <- file.path(
  "raw", "pta_27", "drug_network.csv"
)

clinreport_dir <- file.path(
  "raw", "pta_27", "drugs"
)

gene_file <- file.path(
  "raw", "pta_27", "gene_data.tsv"
)
gene_annotations_file <- file.path("raw", "pta_27", "annotations", "gene_annotations.json")

# Whole-exome sequencing variants annotated by PCGR for the current patient.
wes_vcf_file <- file.path(
  "raw", "pta_27", "wes.pcgr.grch38.pass.vcf.gz"
)
wes_display_name <- "WES · PCGR"

# Protein/gene interaction network used for the drug subgraph visualization.
interaction_network_file <- file.path(
  "raw", "pta_27", "network", "string.human_links_v12_0_min900.Ensembl_edges.tsv"
)

# Transcript -> gene mapping. Used to enrich/normalize gene identifiers and names.
tx2gene_file <- file.path(
  "raw", "pta_27", "tx2gene.tsv"
)

# The CLI supplies these paths in prepared reports. Local development discovers
# the current fixture without assuming a fixed number of GMT collections.
local_gsea_dir <- file.path("raw", "pta_27", "gsea")
gsea_files <- list.files(local_gsea_dir, pattern = "gsea_report_for_.*\\.tsv$", full.names = TRUE)
gmt_files <- list.files(local_gsea_dir, pattern = "\\.gmt$", full.names = TRUE)
gmt_list_file <- file.path(local_gsea_dir, "gmt_list.csv")
gsea_carcinoma_report_file <- ""
gsea_normal_report_file <- ""
pathway_gmt_file <- ""

# Maximum number of interaction-network neighbours displayed per drug target.
interaction_network_max_neighbors <- 40
report_top_n <- 15
report_module_top_n <- 5
report_language <- "ru"
use_tf_activity <- FALSE
use_wes <- FALSE
skip_network_processing <- FALSE

# Display metadata. These are explicit inputs so a later dataset metadata
# loader can replace the patient, interaction network or GMT collection
# without changing presentation code.
patient_id <- "PTA_27"
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

tf_file <- file.path("raw", "pta_27", "tf.tsv")
dorothea_file <- ""
analysis_command <- ""

# A prepared application embeds the exact paths and metadata supplied to
# `clinreport --prepare`; no external configuration file is required.
source("R/config/parameters.R", local = TRUE, encoding = "UTF-8")
prepared_parameters_file <- "R/config/prepared_parameters.R"
if (file.exists(prepared_parameters_file)) {
  source(prepared_parameters_file, local = TRUE, encoding = "UTF-8")
  list2env(clinreport_prepared_parameters, envir = environment())
}
