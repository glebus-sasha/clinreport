# Drug -> interaction-network subgraph business logic

count_visible_wes_genes <- function(nodes) {
  hidden <- if (is.null(nodes$hidden)) rep(FALSE, nrow(nodes)) else nodes$hidden %in% TRUE
  length(unique(nodes$label[nodes$mutated %in% TRUE & !hidden]))
}

target_evidence_summary <- function(target_ids, target_symbols, significant_ids,
                                    wes_symbols = NULL, imported_tfs = NULL) {
  deg <- target_ids %in% significant_ids
  wes <- target_symbols %in% wes_symbols
  tf <- target_symbols %in% imported_tfs
  paste0(length(target_ids), " genes linked to this drug in the input network: ",
    sum(deg), " with significant expression changes (DE)",
    if (!is.null(wes_symbols)) paste0("; ", sum(wes), " with WES variants") else "; WES not supplied",
    if (!is.null(imported_tfs)) paste0("; ", sum(tf), " TFs in the imported analysis") else "; TF analysis not supplied",
    "; ", sum(!(deg | wes | tf)), " without these flags in the supplied data. Groups may overlap.")
}

# Context-only neighbours can be hidden while genes used by the active overlay
# remain visible. Keep nodes in the widget so toggling preserves their positions.
apply_context_visibility <- function(data, graph, overlay_symbols = character(), show_context = FALSE) {
  context_ids <- graph$nodes$gene_id[graph$nodes$type == "Network neighbor" &
    !graph$nodes$gene_name %in% overlay_symbols]
  data$nodes$hidden <- !show_context & data$nodes$id %in% context_ids
  data$nodes$physics <- !data$nodes$hidden
  hidden_ids <- data$nodes$id[data$nodes$hidden]
  data$edges$hidden <- data$edges$from %in% hidden_ids | data$edges$to %in% hidden_ids
  data$edges$physics <- !data$edges$hidden
  data
}

normalize_gene_id <- function(x) {
  str_remove(as.character(x), "\\.[0-9]+$")
}

get_drug_target_ids <- function(drug_id) {
  row <- drug_network |>
    filter(drugId == !!drug_id) |>
    slice(1)

  if (nrow(row) == 0) {
    return(character())
  }

  extract_gene_ids(row$hasEdgesTo) |>
    normalize_gene_id() |>
    unique()
}

get_gene_display_names <- function(gene_ids) {
  gene_ids <- unique(normalize_gene_id(gene_ids))

  if (length(gene_ids) == 0) {
    return(tibble(gene_id = character(), gene_name = character()))
  }

  from_gene_data <- gene_data |>
    transmute(
      gene_id = gene_id_clean,
      gene_name = as.character(gene_name)
    ) |>
    filter(!is.na(gene_id), gene_id != "") |>
    distinct(gene_id, .keep_all = TRUE)

  from_tx2gene <- tx2gene_names |>
    transmute(
      gene_id,
      gene_name = as.character(gene_name)
    ) |>
    distinct(gene_id, .keep_all = TRUE)

  tibble(gene_id = gene_ids) |>
    left_join(from_gene_data, by = "gene_id") |>
    left_join(
      from_tx2gene |> rename(gene_name_tx = gene_name),
      by = "gene_id"
    ) |>
    mutate(
      gene_name = coalesce(
        na_if(gene_name, ""),
        na_if(gene_name_tx, ""),
        gene_id
      )
    ) |>
    select(gene_id, gene_name)
}

get_drug_interaction_subgraph <- function(
    drug_id,
    max_neighbors = interaction_network_max_neighbors
) {
  target_ids <- get_drug_target_ids(drug_id)

  if (length(target_ids) == 0) {
    return(list(
      nodes = tibble(),
      edges = tibble(),
      targets = character()
    ))
  }

  # Keep only interaction edges touching at least one drug target. This is the
  # one-hop neighborhood used on the main screen.
  touching_edges <- interaction_edges |>
    filter(source %in% target_ids | target %in% target_ids)

  neighbor_ids <- unique(c(
    touching_edges$source,
    touching_edges$target
  ))

  neighbor_ids <- setdiff(
    neighbor_ids,
    target_ids
  )

  if (length(neighbor_ids) > max_neighbors) {
    degree_counts <- bind_rows(
      touching_edges |>
        filter(source %in% target_ids, target %in% neighbor_ids) |>
        count(target, name = "degree") |>
        transmute(gene_id = target, degree),
      touching_edges |>
        filter(target %in% target_ids, source %in% neighbor_ids) |>
        count(source, name = "degree") |>
        transmute(gene_id = source, degree)
    ) |>
      group_by(gene_id) |>
      summarise(degree = sum(degree), .groups = "drop") |>
      arrange(desc(degree), gene_id)

    neighbor_ids <- head(degree_counts$gene_id, max_neighbors)
  }

  kept_ids <- unique(c(target_ids, neighbor_ids))

  edges <- touching_edges |>
    filter(source %in% kept_ids, target %in% kept_ids) |>
    mutate(
      source = normalize_gene_id(source),
      target = normalize_gene_id(target)
    ) |>
    distinct(source, target)

  names <- get_gene_display_names(kept_ids)

  target_degree <- bind_rows(
    edges |>
      count(source, name = "degree") |>
      transmute(gene_id = source, degree),
    edges |>
      count(target, name = "degree") |>
      transmute(gene_id = target, degree)
  ) |>
    group_by(gene_id) |>
    summarise(degree = sum(degree), .groups = "drop")

  nodes <- names |>
    mutate(
      type = if_else(gene_id %in% target_ids, "Target", "Network neighbor")
    ) |>
    left_join(target_degree, by = "gene_id") |>
    mutate(degree = coalesce(degree, 0L))

  list(
    nodes = nodes,
    edges = edges,
    targets = target_ids
  )
}
