# Drug page composition

# 19. DRUG PAGE
# ============================================================

render_drug_page <- function(
    drug_id,
    drug_name,
    data,
    genes,
    network_name
) {
  
  chembl_id <- get_chembl_id(data)
  pubchem_cid <- get_pubchem_cid(data)
  
  
  tagList(
    
    div(
      class = "drug-header",
      
      div(
        class = "drug-title-row",
        
        div(
          
          div(
            class = "drug-title",
            drug_name
          ),
          
          div(
            class = "drug-id",
            drug_id
          )
        ),
        
        source_buttons(
          drug_id,
          chembl_id,
          pubchem_cid
        )
      )
    ),

    div(
      class = "drug-network-panel combined-network-panel",
      div(
        class = "network-content-grid",
        div(
          class = "network-graph-column",
        div(
          class = "network-visualization-header",
          div(
            class = "network-visualization-title",
            network_name
          ),
          div(
            class = "network-visualization-subtitle",
            "Selected drug, direct targets, and one-hop STRING neighbors"
          )
        ),
        uiOutput("drug_network_summary"),
        visNetworkOutput("drug_network_graph", height = "500px")
        ),
        div(
          class = "gsea-panel embedded-gsea-panel",
        div(
          class = "gsea-panel-heading",
          div(class = "gsea-panel-title", "Hallmark pathways"),
          actionButton("clear_gsea_pathways", "Clear", class = "gsea-clear-button")
        ),
        div(class = "gsea-panel-subtitle", "All tested gene sets · click to highlight DE genes"),
        div(
          class = "gsea-sort",
          radioButtons(
            "gsea_sort",
            NULL,
            choices = c("Enrichment" = "nes", "FDR" = "fdr"),
            selected = "nes",
            inline = TRUE
          )
        ),
          uiOutput("gsea_matching_pathways"),
          div(class = "gsea-all-pathways-label", "All Hallmark pathways"),
          div(class = "gsea-pathway-list", uiOutput("gsea_pathway_tiles"))
        )
      ),
      div(class = "gsea-details-strip", uiOutput("gsea_pathway_details"))
    ),
    
    
    tabsetPanel(
      
      id = "drug_tabs",
      
      
      tabPanel(
        "General",
        
        render_general(
          drug_id,
          data
        )
      ),
      
      
      tabPanel(
        "Indications & Mechanism",
        
        render_pharmacology(
          data
        )
      ),
      
      
      tabPanel(
        "Safety",
        
        render_safety(
          data
        )
      ),
      
      
      tabPanel(
        "Targets",
        
        div(
          class = "targets-header",
          
          div(
            
            div(
              class = "targets-title",
              "Target genes"
            ),
            
            div(
              class = "targets-criteria",
              
              tags$span(
                class = "criteria-chip",
                paste0(
                  "padj < ",
                  padj_cutoff
                )
              ),
              
              tags$span(
                class = "criteria-chip",
                paste0(
                  "|log2FC| ≥ ",
                  log2fc_cutoff
                )
              )
            )
          ),
          
          div(
            class = "targets-count",
            
            strong(
              nrow(genes)
            ),
            
            " ",
            
            ifelse(
              nrow(genes) == 1,
              "gene",
              "genes"
            )
          )
        ),
        
        div(
          class = "targets-list",
          
          render_target_cards(
            genes
          )
        )
      )
    )
  )
}
