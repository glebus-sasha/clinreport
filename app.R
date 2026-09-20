# ClinReport application entry point

# Dependencies
library(shiny)
library(DT)
library(jsonlite)
library(dplyr)
library(purrr)
library(stringr)
library(htmltools)
library(plotly)
library(visNetwork)

# Configuration and validation
source("R/config/config.R")
source("R/config/validation.R")

# Data layer
source("R/data/load_drug_network.R")
source("R/data/load_gene_data.R")
source("R/data/load_wes.R")
source("R/data/load_string_network.R")
source("R/data/load_tx2gene.R")
source("R/data/load_gsea.R")
source("R/data/drug_directories.R")
source("R/data/json_reader.R")

# Domain layer
source("R/utils/nested.R")
source("R/utils/formatting.R")
source("R/domain/identifiers.R")
source("R/domain/drug_names.R")
source("R/utils/urls.R")
source("R/domain/genes.R")
source("R/domain/drug_network.R")

# UI components
source("R/ui/components/source_buttons.R")
source("R/ui/components/molecular_metric.R")
source("R/ui/components/structure.R")
source("R/ui/tabs/general_tab.R")
source("R/ui/tabs/pharmacology_tab.R")
source("R/ui/tabs/safety_tab.R")
source("R/ui/components/target_cards.R")
source("R/ui/components/wes_variants.R")
source("R/ui/drug_page.R")
source("R/ui/gene_page.R")
source("R/ui/layout.R")

# Server
source("R/server/volcano_plot.R")
source("R/server/drug_network_plot.R")
source("R/server/gsea_panel.R")
source("R/server/server.R")

shinyApp(
  ui = ui,
  server = server
)
