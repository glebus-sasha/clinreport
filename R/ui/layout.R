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
  n_degs <- sum(!is.na(gene_data$padj) & gene_data$padj < padj_cutoff &
    !is.na(gene_data$log2FoldChange) &
    abs(gene_data$log2FoldChange) >= log2fc_cutoff)

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
      div(class = "input-overview-section-title", "Interpretation of drug selection"),
      p(class = "input-overview-subtitle",
        paste0(
          "The recommended medicines were prioritised using differential expression",
          if (isTRUE(as.logical(use_tf_activity))) ", transcription-factor activity" else "",
          if (!isTRUE(as.logical(skip_network_processing))) paste0(", and expansion of the ", interaction_network_name, " interaction network") else "",
          ". ",
          if (isTRUE(as.logical(use_wes))) "WES evidence was also included in prioritisation. " else "WES findings are presented as an additional analysis and were not used to prioritise medicines. ",
          "Pathway enrichment (GSEA) provides metabolic and functional context, while variant-level exome results are reported as supplementary evidence."
        )
      )
    ),
    div(
      class = "input-overview-section",
      div(class = "input-overview-section-title", "Analyses used to build the Drug module"),
      div(
        class = "input-overview-grid evidence-overview-grid",
        input_card(
          "Differential expression",
          "Used",
          paste(n_degs, "DEGs passed the configured padj and log2FC thresholds")
        ),
        if (!isTRUE(as.logical(skip_network_processing)))
          input_card(paste(interaction_network_name, "network expansion"), "Used",
            paste("Expanded drug-linked targets using", interaction_network_name)),
        if (isTRUE(as.logical(use_tf_activity)) && has_tf_analysis)
          input_card("TF regulatory context", "Used",
            paste(nrow(tf_results), "TF activities with DoRothEA target links")),
        if (isTRUE(as.logical(use_wes)) && has_wes)
          input_card("WES variants", "Used",
            paste(nrow(wes_variants), "PASS variants across", nrow(wes_gene_summary), "genes"))
      )
    ),
    div(
      class = "input-overview-section",
      div(class = "input-overview-section-title", "Additional analyses and context"),
      div(
        class = "input-overview-grid evidence-overview-grid",
        if (has_tf_analysis && !isTRUE(as.logical(use_tf_activity)))
          input_card("TF regulatory context", "Available",
            "Shown as supplementary regulatory context", "available"),
        if (has_wes && !isTRUE(as.logical(use_wes)))
          input_card("WES variants", "Available",
            "Shown as supplementary variant evidence", "available"),
        if (isTRUE(as.logical(skip_network_processing)))
          input_card(paste(interaction_network_name, "network expansion"), "Not used",
            "Network expansion was disabled for this report", "unavailable"),
        input_card(
          "GSEA analysis",
          "Additional",
          paste(nrow(pathway_sets), "tested gene sets provide pathway context for targets")
        )
      )
    ),
    tags$details(
      class = "input-overview-metadata",
      tags$summary("Analysis metadata"),
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
          "GSEA collections",
          paste(length(gmt_files), "GMT files"),
          paste(nrow(pathway_sets), "reported gene sets")
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
    ),
    div(
      class = "input-overview-section",
      div(class = "input-overview-section-title", "Preparation command"),
      div(
        class = "input-overview-command",
        title = if (nzchar(analysis_command)) analysis_command else "Not recorded",
        if (nzchar(analysis_command)) analysis_command else "Not recorded"
      )
    )
    )
  )
}

render_report_builder <- function() {
  max_drugs <- max(1L, nrow(drug_network))
  default_top_n <- min(as.integer(report_top_n), max_drugs)
  default_module_n <- min(as.integer(report_module_top_n), default_top_n)
  div(
    class = "report-builder",
    div(class = "report-builder-heading",
      div(class = "report-builder-title", "PDF Report"),
      div(class = "report-builder-subtitle",
        "Build the drug list and adjust it manually when needed.")),
    div(class = "report-title-field",
      textInput("report_title", "Report title",
        value = paste(patient_id, format(Sys.time(), "%Y-%m-%d %H-%M-%S")))),
    div(class = "report-builder-controls",
      numericInput("report_top_n", "Number of drugs (Top N)", default_top_n,
        min = 1, max = max_drugs, step = 1),
      numericInput("report_module_top_n", "Drug modules shown (Top K)", default_module_n,
        min = 0, max = max_drugs, step = 1),
      textInput("report_expert_name", "Expert (optional)", value = ""),
      textAreaInput("report_expert_comment", "Expert comment (optional)", value = "",
        rows = 2, placeholder = "For example: excluded because of a contraindication.")),
    uiOutput("report_status"),
    uiOutput("report_editor"),
    div(class = "report-builder-actions",
      actionButton("report_reset", "Reset to automatic list", class = "btn-default"),
      downloadButton("download_report_pdf", "Generate PDF", class = "btn-primary"))
  )
}

ui <- fluidPage(
  tags$head(
    tags$link(
      rel = "stylesheet",
      type = "text/css",
      href = "css/app.css"
    ),
    tags$script(src = "js/app.js"),
    tags$script(src = "js/i18n.js?v=2")
  ),
  div(
    class = "app-heading",
    div(class = "app-heading-copy",
      div(class = "app-title", "O2T ClinReport - Drug Prioritization Explorer"),
      div(class = "app-subtitle",
        "Drug network - pharmacology - safety - molecular properties - target genes")),
    div(class = "language-switcher",
      selectInput("language", "Language",
        choices = c("English" = "en", "Русский" = "ru"),
        selected = report_language, width = "150px"))
  ),

  render_report_builder(),
  
  
  div(
    class = "main-layout",
    
    
    # ========================================================
    # LEFT ZONE
    # ========================================================
    
    div(
      class = "drug-list",
      
      div(
        class = "network-header",
        tags$button(
          id = "drug_list_toggle",
          type = "button",
          class = "drug-list-collapse",
          title = "Collapse drug list",
          "‹"
        ),
        
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
