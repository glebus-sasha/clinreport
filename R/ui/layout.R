# UI layout

render_input_overview <- function() {
  input_card <- function(label, value, note, status = "available") {
    div(
      class = paste("input-overview-card", status),
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
      "Select a drug to inspect the evidence available for its targets, expression profile and interaction network."
    ),
    div(
      class = "input-overview-section",
      div(class = "input-overview-section-title", "Evidence available for drug assessment"),
      div(
        class = "input-overview-grid evidence-overview-grid",
        input_card(
          "Differential expression",
          "Available",
          paste(nrow(gene_data), "genes; targets are assessed against RNA differential expression")
        ),
        input_card(
          paste(interaction_network_name, "network expansion"),
          "Available",
          paste("One-hop", interaction_network_name, "neighbours extend direct drug targets")
        ),
        input_card(
          "TF regulatory context",
          if (has_tf_analysis) "Available" else "Not available",
          if (has_tf_analysis) {
            paste(nrow(tf_results), "TF activities with DoRothEA target links")
          } else {
            "No usable TF activity input and DoRothEA reference pair was supplied"
          },
          if (has_tf_analysis) "available" else "unavailable"
        ),
        input_card(
          "WES variants",
          if (has_wes) "Available" else "Not supplied",
          if (has_wes) {
            paste(nrow(wes_variants), "PASS variants across", nrow(wes_gene_summary), "genes")
          } else {
            "Variant evidence is excluded from this report"
          },
          if (has_wes) "available" else "unavailable"
        ),
        input_card(
          paste(pathway_collection_name, "GSEA analysis"),
          "Available",
          paste(nrow(pathway_sets), "tested gene sets provide pathway context for targets")
        )
      )
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
        ),
        input_card(
          "Drug annotations",
          if (has_drug_details) basename(clinreport_dir) else "Not supplied",
          if (has_drug_details) "Pharmacology, safety and structure details are available" else "Drug annotation tabs are hidden",
          if (has_drug_details) "available" else "unavailable"
        ),
        input_card(
          wes_display_name,
          if (has_wes) basename(wes_vcf_file) else "Not supplied",
          if (has_wes) paste(nrow(wes_variants), "PASS variants ·", nrow(wes_gene_summary), "genes") else "WES-specific views are hidden",
          if (has_wes) "available" else "unavailable"
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
