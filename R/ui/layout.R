# UI layout

ui <- fluidPage(
  tags$head(
    includeCSS("www/css/app.css")
  ),
  div(
    class = "app-title",
    "ClinReport — Drug Reference"
  ),
  
  div(
    class = "app-subtitle",
    "Drug network · pharmacology · safety · molecular properties · target genes"
  ),
  
  
  div(
    class = "main-layout",
    
    
    # ========================================================
    # LEFT ZONE
    # ========================================================
    
    div(
      class = "drug-list",
      
      div(
        class = "network-header",
        
        div(
          class = "network-title",
          "Drug network"
        ),
        
        div(
          class = "network-subtitle",
          paste(
            nrow(drug_network),
            "records · select a drug to inspect details"
          )
        )
      ),
      
      DTOutput(
        "drug_table"
      )
    ),
    
    
    # ========================================================
    # RIGHT ZONE
    # ========================================================
    
    div(
      class = "drug-details",
      
      uiOutput(
        "right_panel"
      )
    )
  )
)


# ============================================================
