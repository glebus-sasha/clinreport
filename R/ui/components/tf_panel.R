render_tf_selector <- function() {
  tagList(
    div(class = "gsea-panel-title", "Transcription factors"),
    div(class = "gsea-panel-subtitle", "TFs from the imported contrast · select regulons"),
    div(class = "gsea-sort",
        actionButton("tf_all", "All", class = "gsea-clear-button"),
        actionButton("tf_clear", "Clear", class = "gsea-clear-button")),
    selectizeInput("tf_selection", NULL, choices = tf_results$TF, multiple = TRUE,
                   options = list(placeholder = "Select TFs")),
    checkboxGroupInput("tf_confidence", "DoRothEA confidence", c("A", "B", "C"), c("A", "B", "C"), inline = TRUE),
    checkboxInput("tf_expression_only", "Targets present in RNA table", TRUE),
    checkboxInput("tf_deg_only", "Only DE targets", FALSE),
    uiOutput("tf_reference_status"),
    actionButton("tf_retry", "Download / retry reference", class = "gsea-clear-button"),
    div(class = "gsea-panel-subtitle", "Diamond: selected TF (activity colour). Circle: target (RNA colour). Green arrow: activation; dashed purple: repression. Pink centre: WES variant."),
    div(class = "gsea-panel-subtitle", "Reference edges describe prior knowledge. Activity is not recalculated with display filters.")
  )
}

render_tf_table <- function() {
  tagList(div(class = "targets-header", div(class = "targets-title", "TF activity")),
          div(class = "gsea-panel-subtitle", "Imported results; logFC is an activity difference, not RNA log2FC. File naming does not define significance."),
          DTOutput("tf_table"))
}

render_regulation_panel <- function() {
  div(class = "drug-network-panel combined-network-panel",
    div(class = "network-content-grid",
      div(class = "network-graph-column",
        conditionalPanel("!input.regulation_mode || input.regulation_mode === 'pathways'",
          div(class = "network-visualization-title", network_display_name),
          uiOutput("drug_network_summary"),
          checkboxInput("gsea_color_by_direction", "Direction", FALSE),
          visNetworkOutput("drug_network_graph", height = "500px")),
        conditionalPanel("input.regulation_mode === 'tf'",
          div(class = "network-visualization-title", "TF regulatory network · DoRothEA"),
          uiOutput("tf_network_summary"),
          visNetworkOutput("tf_network_graph", height = "500px"))
      ),
      div(class = "gsea-panel embedded-gsea-panel",
        radioButtons("regulation_mode", NULL, c("Pathways" = "pathways", "TF" = "tf"), inline = TRUE),
        conditionalPanel("!input.regulation_mode || input.regulation_mode === 'pathways'",
          div(class = "gsea-panel-heading",
            div(class = "gsea-panel-title", paste(pathway_collection_name, "pathways")),
            actionButton("clear_gsea_pathways", "Clear", class = "gsea-clear-button")),
          div(class = "gsea-panel-subtitle", paste("All tested", pathway_collection_name, "gene sets · click to highlight DE genes")),
          radioButtons("gsea_sort", NULL, c("Enrichment" = "nes", "FDR" = "fdr"), inline = TRUE),
          uiOutput("gsea_matching_pathways"),
          div(class = "gsea-pathway-section-title", paste("All", pathway_collection_name, "pathways")),
          div(class = "gsea-pathway-list gsea-all-pathway-list", uiOutput("gsea_pathway_tiles"))),
        conditionalPanel("input.regulation_mode === 'tf'", render_tf_selector())
      )
    ),
    conditionalPanel("!input.regulation_mode || input.regulation_mode === 'pathways'",
      div(class = "gsea-details-strip", uiOutput("gsea_pathway_details")))
  )
}
