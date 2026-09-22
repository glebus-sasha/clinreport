# Overlay signed regulatory edges without replacing drug or STRING edges.
get_target_tf_links <- function(target_symbols) {
  tf_regulons |>
    filter(tf %in% tf_results$TF, target %in% target_symbols,
           target %in% gene_data$gene_name) |>
    distinct(tf, target, mor, confidence)
}

build_tf_network <- function(graph, drug, selected, significant_ids, color_by_direction,
                             target_connected_only = FALSE) {
  selected <- intersect(selected, tf_results$TF)
  links <- tf_regulons |>
    filter(tf %in% selected, target %in% gene_data$gene_name)
  if (target_connected_only) {
    target_symbols <- graph$nodes$gene_name[graph$nodes$gene_id %in% graph$targets]
    linked_tfs <- get_target_tf_links(target_symbols)$tf
    selected <- intersect(selected, unique(c(graph$nodes$gene_name, linked_tfs)))
    links <- links |> filter(tf %in% selected, target %in% graph$nodes$gene_name)
  }
  symbols <- unique(c(selected, links$target))
  lookup <- gene_data |>
    transmute(gene_name, gene_id = gene_id_clean) |>
    filter(!is.na(gene_name), nzchar(gene_name)) |>
    distinct(gene_name, .keep_all = TRUE)
  extras <- tibble(gene_name = setdiff(symbols, graph$nodes$gene_name)) |>
    left_join(lookup, by = "gene_name") |>
    mutate(gene_id = coalesce(gene_id, paste0("__TF_GENE__", gene_name)),
           type = if_else(gene_name %in% selected, "TF regulator", "TF regulon target"), degree = 0L)
  expanded <- graph
  expanded$nodes <- bind_rows(graph$nodes, extras) |> distinct(gene_id, .keep_all = TRUE)
  data <- prepare_vis_network(expanded, drug, setNames(character(), character()),
                              significant_ids, color_by_direction)
  source_rows <- which(data$nodes$label %in% selected & data$nodes$id != "__DRUG__")
  activity <- tf_results[match(data$nodes$label[source_rows], tf_results$TF), ]
  data$nodes$shape[source_rows] <- "triangle"
  data$nodes$size[source_rows] <- pmax(data$nodes$size[source_rows], 23)
  data$nodes$color[source_rows] <- if (color_by_direction) {
    ifelse(is.na(activity$logFC), pathway_direction_colors[["unknown"]],
           ifelse(activity$logFC >= 0, pathway_direction_colors[["up"]], pathway_direction_colors[["down"]]))
  } else rep("#8b5cf6", length(source_rows))
  if (length(source_rows)) data$nodes$title[source_rows] <- paste0(data$nodes$title[source_rows],
    "<br><b>Transcription factor</b><br>Activity difference: ", signif(activity$logFC, 4),
    "<br>Activity FDR: ", signif(activity$adj.P.Val, 4),
    "<br>Activity differs from this gene's RNA expression")
  regulatory_edges <- links |>
    mutate(from = data$nodes$id[match(tf, data$nodes$label)],
           to = data$nodes$id[match(target, data$nodes$label)]) |>
    transmute(id = paste("__TF_EDGE__", tf, target, mor, confidence, sep = "::"), from, to,
      title = paste0("DoRothEA: ", htmlEscape(tf), " → ", htmlEscape(target),
                     if_else(mor > 0, " · activation (+)", " · repression (−)"),
                     "<br>Confidence: ", confidence, "<br>Prior regulatory relationship"),
      color = if_else(mor > 0, "#059669", "#db2777"), width = 2,
      dashes = mor < 0, arrows.to.enabled = TRUE,
      arrows.to.type = if_else(mor > 0, "arrow", "bar"),
      smooth.enabled = TRUE, smooth.type = "curvedCW", smooth.roundness = 0.12)
  data$edges <- bind_rows(data$edges, regulatory_edges)
  data$extra_ids <- extras$gene_id
  data$tf_edges <- regulatory_edges
  data
}
