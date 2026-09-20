# UI layout

render_input_overview <- function() {
  input_card <- function(label, value, note) {
    div(
      class = "input-overview-card",
      div(class = "input-overview-label", label),
      div(class = "input-overview-value", title = value, value),
      div(class = "input-overview-note", note)
    )
  }

  div(
    class = "input-overview",
    div(
      class = "input-overview-heading",
      h2(class = "input-overview-title", "Analysis inputs"),
      div(
        class = "input-overview-patient",
        div(class = "input-overview-patient-label", "Patient ID"),
        div(class = "input-overview-patient-value", patient_id)
      )
    ),
    div(
      class = "input-overview-subtitle",
      "Select a drug from the list to inspect its pharmacology, targets, expression profile and interaction network."
    ),
    div(
      class = "input-overview-section",
      div(class = "input-overview-section-title", "Data sources"),
      div(
        class = "input-overview-grid",
        input_card(
          "Drug–gene network",
          basename(drug_network_file),
          paste(nrow(drug_network), "drug records")
        ),
        input_card(
          "Differential expression",
          basename(gene_file),
          paste(nrow(gene_data), "genes")
        ),
        input_card(
          paste(interaction_network_name, "interactions"),
          basename(interaction_network_file),
          paste(nrow(interaction_edges), "high-confidence edges")
        ),
        input_card(
          paste(pathway_collection_name, "GSEA"),
          basename(pathway_gmt_file),
          paste(nrow(pathway_sets), "gene sets")
        )
      )
    ),
    div(
      class = "input-overview-section",
      div(class = "input-overview-section-title", "Analysis settings"),
      div(
        class = "input-overview-settings",
        span(class = "input-overview-setting", paste0("Network: ", network_display_name)),
        span(class = "input-overview-setting", paste0(interaction_network_name, " neighbours: ", interaction_network_max_neighbors)),
        span(class = "input-overview-setting", paste0("padj < ", padj_cutoff)),
        span(class = "input-overview-setting", paste0("|log2FC| ≥ ", log2fc_cutoff))
      )
    )
  )
}

ui <- fluidPage(
  tags$head(
    tags$link(
      rel = "stylesheet",
      type = "text/css",
      href = "css/app.css"
    ),
    tags$script(src = "js/app.js")
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

      uiOutput("focused_drugs"),
      
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
