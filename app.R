# ClinReport application entry point

# Dependencies
if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) {
  Sys.setlocale("LC_CTYPE", "English_United States.utf8")
}
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
source("R/config/config.R", encoding = "UTF-8")
source("R/config/validation.R", encoding = "UTF-8")

# Data layer
source("R/data/load_drug_network.R", encoding = "UTF-8")
source("R/data/load_gene_data.R", encoding = "UTF-8")
source("R/data/load_wes.R", encoding = "UTF-8")
source("R/data/load_string_network.R", encoding = "UTF-8")
source("R/data/load_tx2gene.R", encoding = "UTF-8")
source("R/data/load_gsea.R", encoding = "UTF-8")
source("R/data/load_tf.R", encoding = "UTF-8")
source("R/data/drug_directories.R", encoding = "UTF-8")
source("R/data/json_reader.R", encoding = "UTF-8")

# Domain layer
source("R/utils/nested.R", encoding = "UTF-8")
source("R/utils/formatting.R", encoding = "UTF-8")
source("R/domain/identifiers.R", encoding = "UTF-8")
source("R/utils/urls.R", encoding = "UTF-8")
source("R/domain/genes.R", encoding = "UTF-8")
source("R/domain/drug_network.R", encoding = "UTF-8")
source("R/domain/tf_network.R", encoding = "UTF-8")

# UI components
source("R/ui/components/source_buttons.R", encoding = "UTF-8")
source("R/ui/components/molecular_metric.R", encoding = "UTF-8")
source("R/ui/components/structure.R", encoding = "UTF-8")
source("R/ui/tabs/general_tab.R", encoding = "UTF-8")
source("R/ui/tabs/pharmacology_tab.R", encoding = "UTF-8")
source("R/ui/tabs/safety_tab.R", encoding = "UTF-8")
source("R/ui/components/target_cards.R", encoding = "UTF-8")
source("R/ui/components/wes_variants.R", encoding = "UTF-8")
source("R/ui/components/network_overlay.R", encoding = "UTF-8")
source("R/ui/drug_page.R", encoding = "UTF-8")
source("R/ui/gene_page.R", encoding = "UTF-8")
source("R/ui/layout.R", encoding = "UTF-8")

# Server
source("R/server/volcano_plot.R", encoding = "UTF-8")
source("R/server/drug_network_plot.R", encoding = "UTF-8")
source("R/server/gsea_panel.R", encoding = "UTF-8")
source("R/server/tf_panel.R", encoding = "UTF-8")
source("R/server/server.R", encoding = "UTF-8")

shinyApp(
  ui = ui,
  server = server
)
