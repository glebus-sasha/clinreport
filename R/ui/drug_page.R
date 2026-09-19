# Drug page composition

# 19. DRUG PAGE
# ============================================================

render_drug_page <- function(
    drug_id,
    drug_name,
    data,
    genes
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


