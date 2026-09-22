register_tf_outputs <- function(input, output, session, selected_tfs,
                                drug_target_symbols, selected_drug, selected_network_gene,
                                significant_ids, target_connected_only) {
  toggle_tf <- function(tf) {
    if (!tf %in% tf_results$TF) return(invisible())
    selected_tfs(if (tf %in% selected_tfs()) setdiff(selected_tfs(), tf) else c(selected_tfs(), tf))
  }
  observeEvent(input$selected_tf_index, {
    i <- suppressWarnings(as.integer(input$selected_tf_index))
    if (length(i) == 1 && !is.na(i) && i >= 1 && i <= nrow(tf_results)) toggle_tf(tf_results$TF[i])
  })
  observeEvent(input$select_all_tfs, selected_tfs(tf_results$TF))
  observeEvent(input$clear_tfs, selected_tfs(character()))
  output$tf_status <- renderUI({
    if (!is.null(tf_input$error)) return(div(class = "network-empty-message", tf_input$error))
    tagList(
      if (!is.null(dorothea_resource$error)) div(class = "network-empty-message",
        paste("DoRothEA download failed. Activity table is available; restart to retry.", dorothea_resource$error)),
      div(class = "gsea-panel-subtitle", paste(sum(tf_results$adj.P.Val < 0.05, na.rm = TRUE),
        "of", nrow(tf_results), "imported TFs pass FDR < 0.05.")))
  })
  make_tf_tiles <- function(results, links = NULL) {
    ordered <- results |> mutate(input_index = match(TF, tf_results$TF))
    ordered <- if (identical(input$tf_sort, "fdr")) {
      ordered |> arrange(adj.P.Val, desc(abs(logFC)), TF)
    } else {
      ordered |> arrange(desc(abs(logFC)), adj.P.Val, TF)
    }
    lapply(seq_len(nrow(ordered)), function(i) {
      x <- ordered[i, ]
      div(class = paste0("gsea-pathway-tile", if (x$TF %in% selected_tfs()) " active" else ""),
        style = "border-left-color:#8b5cf6",
        onclick = sprintf("Shiny.setInputValue('selected_tf_index', %d, {priority:'event'})", x$input_index),
        div(class = "gsea-pathway-name", x$TF),
        div(class = "gsea-pathway-metrics", span(sprintf("Δ activity %+.2f", x$logFC)),
          span(sprintf("FDR %.4g", x$adj.P.Val))),
        if (!is.null(links)) div(class = "gsea-pathway-metrics",
          paste("Drug targets:", paste(sort(unique(links$target[links$tf == x$TF])), collapse = ", "))))
    })
  }
  output$tf_tiles <- renderUI({ make_tf_tiles(tf_results) })
  output$tf_matching_factors <- renderUI({
    drug <- selected_drug()
    req(drug)
    links <- get_target_tf_links(drug_target_symbols())
    related <- tf_results |> filter(TF %in% links$tf)
    section <- function(title, results, edges = NULL) div(
      class = "gsea-pathway-section gsea-related-section",
      div(class = "gsea-pathway-section-title", title),
      div(class = "gsea-pathway-list gsea-related-pathway-list", make_tf_tiles(results, edges)))
    focused <- selected_network_gene()
    focused_links <- if (!is.null(focused)) get_target_tf_links(focused$symbol) else tf_regulons[0, ]
    tagList(
      if (nrow(related)) section(paste0("Regulators of ", drug$label, " targets"), related, links)
      else div(class = "gsea-panel-subtitle", if (!is.null(dorothea_resource$error))
        "Related TFs unavailable: DoRothEA download failed." else "No imported TF regulates these drug targets in DoRothEA A/B/C."),
      if (nrow(focused_links)) section(paste0("Regulators of selected gene: ", focused$symbol),
        tf_results |> filter(TF %in% focused_links$tf)))
  })
  output$tf_drug_summary <- renderUI({
    drug <- selected_drug()
    req(drug)
    targets <- unique(drug_target_symbols())
    target_rows <- gene_data |> filter(gene_name %in% targets)
    wes_symbols <- if (has_wes) unique(wes_gene_summary$gene_symbol) else character()
    related <- get_target_tf_links(targets)
    div(class = "network-summary",
      span(class = "network-summary-drug", drug$label),
      span(paste0(" · ", paste(c(
        paste0(length(targets), " direct targets"),
        paste0(sum(target_rows$gene_id_clean %in% significant_ids), " DE-significant targets"),
        if (has_wes) paste0(sum(targets %in% wes_symbols), " WES-mutated targets"),
        paste0(sum(vapply(pathway_sets$genes, function(gs) any(targets %in% gs), logical(1))),
          " ", pathway_collection_name, " pathways contain drug targets")
      ), collapse = " · "))),
      span(class = "network-summary-context", paste0(
        "Graph context: ", nrow(related), " DoRothEA TF-target links connect imported TFs to direct targets.")))
  })
  output$tf_overlay_summary <- renderUI({
    edges <- tf_regulons |> filter(tf %in% selected_tfs(), target %in% gene_data$gene_name)
    if (target_connected_only()) {
      graph <- get_drug_interaction_subgraph(selected_drug()$drugId)
      linked_tfs <- get_target_tf_links(drug_target_symbols())$tf
      edges <- edges |> filter(tf %in% c(graph$nodes$gene_name, linked_tfs),
                               target %in% graph$nodes$gene_name)
    }
    missing <- setdiff(selected_tfs(), edges$tf)
    div(class = "gsea-panel-subtitle", paste(length(selected_tfs()), "selected TFs ·",
      n_distinct(edges$target), "regulon targets ·", nrow(edges), "regulatory links"),
      if (length(missing)) div(paste("No eligible links:", paste(missing, collapse = ", "))))
  })
  output$tf_table <- renderDT({
    links <- get_target_tf_links(drug_target_symbols())
    source_link <- function(tf) paste0(
      "<a href='https://www.ncbi.nlm.nih.gov/gene/?term=",
      utils::URLencode(paste0(tf, "[Symbol] AND Homo sapiens[Organism]"), reserved = TRUE),
      "' target='_blank' rel='noopener noreferrer'>NCBI Gene ↗</a>")
    tf_results |>
      transmute(TF, `Activity difference` = logFC, `Mean activity` = AveExpr,
                t, `P-value` = P.Value, FDR = adj.P.Val,
                `FDR < 0.05` = if_else(!is.na(adj.P.Val) & adj.P.Val < 0.05, "Yes", "No"),
                `Drug targets` = vapply(TF, function(tf) paste(sort(unique(links$target[links$tf == tf])), collapse = ", "), character(1)),
                Links = vapply(TF, source_link, character(1))) |>
      datatable(rownames = FALSE, selection = "single", escape = 1:8,
        filter = "top", class = "compact stripe hover",
        options = list(pageLength = 15, lengthMenu = c(15, 25, 50), scrollX = TRUE,
          autoWidth = FALSE, columnDefs = list(list(targets = 8, orderable = FALSE)))) |>
      formatSignif(columns = 2:6, digits = 4)
  })
  observeEvent(input$tf_table_rows_selected, {
    i <- input$tf_table_rows_selected
    if (length(i) == 1 && i >= 1 && i <= nrow(tf_results)) {
      toggle_tf(tf_results$TF[i])
      updateRadioButtons(session, "overlay_mode", selected = "tf")
      selectRows(dataTableProxy("tf_table", session = session), NULL)
    }
  })
}
