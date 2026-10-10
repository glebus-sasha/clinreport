library(jsonlite)
library(purrr)
library(stringr)
library(htmltools)
library(shiny)
library(dplyr)
library(plotly)

source("R/utils/nested.R")
source("R/utils/formatting.R")
source("R/domain/identifiers.R")
source("R/domain/drug_description.R")
source("R/utils/urls.R")
source("R/ui/components/source_buttons.R")
source("R/ui/components/molecular_metric.R")
source("R/ui/components/structure.R")
source("R/ui/tabs/general_tab.R")

drug <- list(
  chembl = list(chembl_id = list(chembl_id = "CHEMBL597365"),
    chembl_molecule = list(pref_name = "SUNITINIB")),
  pubchem = list(
    pubchem_cids = list(IdentifierList = list(CID = 5329102)),
    pubchem_description = jsonlite::fromJSON(
      "raw/pta_27/drugs/R_PTA_27_DB01268_pubchem/pubchem_description.json",
      simplifyVector = FALSE)))
description <- get_drug_description(drug)
stopifnot(grepl("Sunitinib", description$text, fixed = TRUE),
  identical(description$source, "ChEBI"))
stopifnot(is.null(get_drug_description(list(pubchem = list()))))
html <- as.character(render_general("DB01268", drug))
stopifnot(grepl("https://go.drugbank.com/drugs/DB01268", html, fixed = TRUE),
  grepl(description$text, html, fixed = TRUE),
  grepl("Source: ", html, fixed = TRUE),
  grepl("ChEBI", html, fixed = TRUE))

fixture <- tempfile("gene-annotations-", fileext = ".json")
writeLines('{"metadata":{"source":"MyGene.info"},"genes":{"ENSG00000000003":{"name":"tetraspanin 6","summary":"A gene summary"}}}', fixture)
gene_annotations_file <- fixture
clinreport_dir <- ""
source("R/data/load_gene_annotations.R")
stopifnot(identical(gene_annotations[["ENSG00000000003"]]$summary, "A gene summary"))
stopifnot(identical(gene_annotation_source, "MyGene.info"))
source("R/ui/gene_page.R")
wes_variants <- data.frame(gene_symbol = character())
drug_network <- data.frame(drugId = character())
has_wes <- FALSE
padj_cutoff <- 0.05
log2fc_cutoff <- 1
gene <- data.frame(gene_name = "TSPAN6", gene_id_clean = "ENSG00000000003",
  log2FoldChange = 1.5, padj = 0.01, baseMean = 100)
gene_html <- as.character(render_gene_page(gene))
stopifnot(grepl("A gene summary", gene_html, fixed = TRUE),
  grepl("MyGene.info", gene_html, fixed = TRUE))
staged_dir <- tempfile("staged-annotations-")
dir.create(staged_dir)
stopifnot(file.copy(fixture, file.path(staged_dir, "gene_annotations.json")))
gene_annotations_file <- ""
clinreport_dir <- staged_dir
source("R/data/load_gene_annotations.R")
stopifnot(identical(gene_annotations[["ENSG00000000003"]]$summary, "A gene summary"))
cat("PASS: drug description, DrugBank link, source labels, gene annotation loading\n")
