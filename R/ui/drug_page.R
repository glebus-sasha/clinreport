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
            class = "network-header-copy",
            div(
              class = "network-visualization-title",
              network_name
            ),
            div(
              class = "network-visualization-subtitle",
              paste("Drug targets, differential expression, and", pathway_collection_name, "pathway coverage")
            )
          )
        ),
        uiOutput("drug_network_summary"),
        div(
          class = "gsea-sort network-color-toggle",
          tags$label(
            class = "radio-inline",
            tags$input(
              id = "gsea_color_by_direction",
              type = "checkbox",
              class = "shiny-input-checkbox"
            ),
            "Direction"
          )
        ),
        visNetworkOutput("drug_network_graph", height = "500px")
        ),
        div(
          class = "gsea-panel embedded-gsea-panel",
        render_network_overlay_panel()
        )
      ),
      conditionalPanel("input.overlay_mode !== 'tf'", div(class = "gsea-details-strip", uiOutput("gsea_pathway_details")))
    ),
    
    
    tabsetPanel(
      
      id = "drug_tabs",
      selected = "General",
      
      
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
      ),

      tabPanel(
        "WES variants",
        render_wes_variants()
      ),
      render_tf_tab()
    )
  )
}
