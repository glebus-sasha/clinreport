# Drug / STRING network visualization

prepare_vis_network <- function(subgraph, drug, pathway_genes, significant_ids) {
  gene_nodes <- subgraph$nodes |>
    mutate(
      label = if_else(is.na(gene_name) | gene_name == "", gene_id, gene_name),
      significant = gene_id %in% significant_ids,
      pathway_member = label %in% pathway_genes,
      pathway_hit = significant & pathway_member,
      group = if_else(type == "Target", "Target", "STRING neighbor"),
      shape = "dot",
      size = case_when(
        pathway_hit ~ 27,
        pathway_member ~ 17,
        type == "Target" ~ 22,
        TRUE ~ 11
      ),
      color = case_when(
        pathway_hit ~ "#d97706",
        pathway_member ~ "#7c3aed",
        type == "Target" ~ "#2563eb",
        TRUE ~ "#94a3b8"
      ),
      borderWidth = if_else(significant, 4, 1),
      title = paste0(
        "<b>", htmlEscape(label), "</b><br>",
        "Ensembl: ", htmlEscape(gene_id), "<br>",
        if_else(type == "Target", "Drug target", "STRING neighbor"),
        "<br>DE status: ", if_else(significant, "significant", "not significant"),
        "<br>Selected Hallmark pathway: ", if_else(pathway_member, "member", "not a member")
      )
    ) |>
    transmute(id = gene_id, label, group, shape, size, color, borderWidth, title)

  nodes <- bind_rows(
    tibble(
      id = "__DRUG__",
      label = drug$label,
      group = "Drug",
      shape = "diamond",
      size = 30,
      color = "#17202a",
      title = paste0("<b>", htmlEscape(drug$label), "</b><br>Drug")
    ),
    gene_nodes
  ) |>
    distinct(id, .keep_all = TRUE)

  nodes <- nodes |>
    mutate(
      x = if_else(id == "__DRUG__", -430, NA_real_),
      y = if_else(id == "__DRUG__", 0, NA_real_),
      fixed = id == "__DRUG__"
    )

  edges <- bind_rows(
    tibble(
      from = "__DRUG__",
      to = subgraph$targets,
      title = "Drug target",
      color = "#2563eb",
      width = 3,
      dashes = FALSE
    ),
    subgraph$edges |>
      transmute(
        from = source,
        to = target,
        title = "STRING interaction",
        color = "#cbd5e1",
        width = 1,
        dashes = TRUE
      )
  ) |>
    filter(from %in% nodes$id, to %in% nodes$id) |>
    distinct(from, to, .keep_all = TRUE)

  list(nodes = nodes, edges = edges)
}

register_drug_network_outputs <- function(
    output,
    selected_drug,
    subgraph,
    selected_pathways,
    significant_ids
) {
  output$drug_network_summary <- renderUI({
    drug <- selected_drug()

    if (is.null(drug) || nrow(drug) == 0) {
      return(
        div(
          class = "network-empty-message",
          "Select a drug in the table to display its STRING subgraph."
        )
      )
    }

    graph <- subgraph()
    div(
      class = "network-summary",
      span(class = "network-summary-drug", drug$label),
      span(paste0(
        " · ", length(graph$targets), " target genes · ",
        nrow(graph$nodes), " genes in subgraph · ",
        nrow(graph$edges), " STRING edges"
      ))
    )
  })

  output$drug_network_graph <- renderVisNetwork({
    drug <- selected_drug()
    req(drug)

    graph <- subgraph()
    req(nrow(graph$nodes) > 0)

    data <- prepare_vis_network(
      graph,
      drug,
      pathway_genes = character(),
      significant_ids = significant_ids
    )
    req(nrow(data$edges) > 0)

    visNetwork(data$nodes, data$edges, width = "100%", height = "500px") |>
      visNodes(
        borderWidth = 1,
        font = list(size = 14, color = "#17202a", face = "Arial")
      ) |>
      visEdges(
        smooth = FALSE,
        color = list(color = "#cbd5e1", highlight = "#475569")
      ) |>
      visOptions(
        highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
        collapse = FALSE
      ) |>
      visInteraction(
        hover = TRUE,
        navigationButtons = FALSE,
        dragNodes = TRUE,
        zoomView = TRUE
      ) |>
      visPhysics(
        solver = "forceAtlas2Based",
        forceAtlas2Based = list(
          gravitationalConstant = -50,
          centralGravity = 0.01,
          springLength = 120,
          springConstant = 0.05,
          damping = 0.4
        ),
        stabilization = list(enabled = TRUE, iterations = 300)
      )
  })

  observeEvent(selected_pathways(), {
    drug <- selected_drug()
    req(drug)

    graph <- subgraph()
    selected <- selected_pathways()
    pathway_genes <- hallmark_pathways |>
      filter(pathway %in% selected) |>
      pull(genes) |>
      unlist(use.names = FALSE) |>
      unique()

    nodes <- prepare_vis_network(
      graph,
      drug,
      pathway_genes = pathway_genes,
      significant_ids = significant_ids
    )$nodes |>
      select(id, size, color, borderWidth, title)

    visNetworkProxy("drug_network_graph") |>
      visUpdateNodes(nodes)
  }, ignoreInit = TRUE)
}
