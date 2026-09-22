render_network_overlay_panel <- function() {
  overlay_choices <- c("Pathways" = "pathways")
  if (has_tf_analysis) overlay_choices <- c(overlay_choices, "TF" = "tf")
  tagList(
    div(class = "overlay-mode-switch",
      radioButtons("overlay_mode", NULL, overlay_choices, inline = TRUE)),
    conditionalPanel("input.overlay_mode !== 'tf'",
      div(class = "gsea-panel-heading",
        div(class = "gsea-panel-title", paste(pathway_collection_name, "pathways")),
        actionButton("clear_gsea_pathways", "Clear", class = "gsea-clear-button")),
      div(class = "tf-actions",
        div(class = "gsea-sort",
          radioButtons("gsea_sort", NULL, c("Enrichment" = "nes", "FDR" = "fdr"), inline = TRUE)),
        actionButton("select_all_pathways", "All target-linked pathways", class = "gsea-clear-button")),
      uiOutput("gsea_matching_pathways"),
      div(class = "gsea-pathway-section-title", "All imported pathways"),
      div(class = "gsea-pathway-list gsea-all-pathway-list", uiOutput("gsea_pathway_tiles"))),
    conditionalPanel("input.overlay_mode === 'tf'",
      div(class = "gsea-panel-heading",
        div(class = "gsea-panel-title", "Transcription factors"),
        actionButton("clear_tfs", "Clear", class = "gsea-clear-button")),
      div(class = "gsea-panel-subtitle", "Imported TFs · click to toggle · DoRothEA A/B/C · targets present in RNA table"),
      uiOutput("tf_drug_summary"),
      div(class = "tf-actions",
        div(class = "gsea-sort", title = "Enrichment: absolute activity difference; FDR: adjusted p-value",
          radioButtons("tf_sort", NULL, c("Enrichment" = "activity", "FDR" = "fdr"), inline = TRUE)),
        actionButton("select_all_tfs", "All target-linked TFs", class = "gsea-clear-button")),
      uiOutput("tf_status"), uiOutput("tf_overlay_summary"),
      uiOutput("tf_matching_factors"),
      div(class = "tf-legend", "▲ TF · blue: drug–target · grey dashed: STRING",
        tags$br(), "Green → activation · pink dashed ⊣ repression",
        tags$br(), "Direction: TF colour = activity difference; other genes = RNA expression."),
      div(class = "gsea-pathway-section-title", "All imported TFs"),
      div(class = "gsea-pathway-list gsea-all-pathway-list", uiOutput("tf_tiles")))
  )
}
render_tf_tab <- function() {
  if (!has_tf_analysis) return(NULL)
  tabPanel("TF", div(class = "wes-variants-panel",
    div(class = "targets-header",
      div(div(class = "targets-title", "Transcription factor activity"),
        div(class = "targets-criteria",
          span(class = "criteria-chip", "DoRothEA A/B/C"),
          span(class = "criteria-chip", paste(sum(tf_results$adj.P.Val < 0.05, na.rm = TRUE), "FDR < 0.05")))),
      div(class = "targets-count", strong(nrow(tf_results)), " transcription factors")),
    div(class = "wes-variants-note",
      "Activity difference is not RNA log2FC. Select a row to toggle its regulon. Drug targets refer to the selected drug. ",
      tags$a(href = "https://saezlab.github.io/dorothea/", target = "_blank", rel = "noopener noreferrer", "DoRothEA source ↗")),
    div(class = "wes-variants-table", DTOutput("tf_table"))))
}
