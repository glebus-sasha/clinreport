# Drug / STRING network visualization

prepare_vis_network <- function(subgraph, drug, pathway_gene_colors, significant_ids, color_by_direction) {
  gene_nodes <- subgraph$nodes |>
    left_join(
      gene_data |>
        transmute(gene_id = gene_id_clean, log2FoldChange) |>
        distinct(gene_id, .keep_all = TRUE),
      by = "gene_id"
    ) |>
    mutate(
      label = if_else(is.na(gene_name) | gene_name == "", gene_id, gene_name),
      pathway_color = unname(pathway_gene_colors[label]),
      pathway_member = !is.na(pathway_color),
      significant = gene_id %in% significant_ids,
      group = if_else(type == "Target", "Target", "STRING neighbor"),
      shape = "dot",
      size = case_when(
        type == "Target" ~ 22,
        pathway_member ~ 17,
        TRUE ~ 11
      ),
      color = if (color_by_direction) {
        case_when(
          is.na(log2FoldChange) ~ hallmark_direction_colors[["unknown"]],
          log2FoldChange >= 0 ~ hallmark_direction_colors[["up"]],
          TRUE ~ hallmark_direction_colors[["down"]]
        )
      } else {
        case_when(
          pathway_member ~ pathway_color,
          type == "Target" ~ "#2563eb",
          TRUE ~ "#94a3b8"
        )
      },
      borderWidth = if_else(significant, 4, 1),
      title = paste0(
        "<b>", htmlEscape(label), "</b><br>",
        "Ensembl: ", htmlEscape(gene_id), "<br>",
        if_else(type == "Target", "Drug target", "STRING neighbor"),
        "<br>DE status: ", if_else(significant, "significant", "not significant"),
        "<br>Expression: ", case_when(
          is.na(log2FoldChange) ~ "unknown",
          log2FoldChange >= 0 ~ "up",
          TRUE ~ "down"
        ),
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
    significant_ids,
    selected_network_gene,
    color_by_direction
) {
  extra_pathway_node_ids <- reactiveVal(character())
  extra_pathway_edge_ids <- reactiveVal(character())
  rendered_pathways <- reactiveVal(character())

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
      pathway_gene_colors = setNames(character(), character()),
      significant_ids = significant_ids,
      # Colour changes are applied through the proxy below.  Keeping this
      # reactive value isolated prevents visNetwork from rebuilding itself.
      color_by_direction = isolate(color_by_direction())
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
      ) |>
      visEvents(
        select = JS(
          "function(properties) {
             Shiny.setInputValue('selected_network_node', properties.nodes[0] || null, {priority: 'event'});
           }"
        ),
        deselectNode = JS(
          "function() {
             Shiny.setInputValue('selected_network_node', null, {priority: 'event'});
           }"
        )
      )
  })

  observeEvent(list(selected_pathways(), color_by_direction()), {
    drug <- selected_drug()
    req(drug)

    graph <- subgraph()
    selected <- selected_pathways()
    structure_changed <- !identical(selected, rendered_pathways())
    selected_sets <- hallmark_pathways |>
      filter(pathway %in% selected) |>
      mutate(
        pathway_color = if (color_by_direction()) {
          pathway_direction_color(nes)
        } else {
          pathway_identity_color(pathway)
        }
      )
    pathway_gene_colors <- setNames(character(), character())
    if (nrow(selected_sets) > 0) {
      pathway_gene_colors <- unlist(
        lapply(seq_len(nrow(selected_sets)), function(i) {
          setNames(
            rep(selected_sets$pathway_color[i], length(selected_sets$genes[[i]])),
            selected_sets$genes[[i]]
          )
        }),
        use.names = TRUE
      )
      pathway_gene_colors <- pathway_gene_colors[!duplicated(names(pathway_gene_colors))]
    }

    base_data <- prepare_vis_network(
      graph,
      drug,
      pathway_gene_colors = pathway_gene_colors,
      significant_ids = significant_ids,
      color_by_direction = color_by_direction()
    )

    nodes <- base_data$nodes |>
      select(id, size, color, borderWidth, title)

    proxy <- visNetworkProxy("drug_network_graph") |>
      visUpdateNodes(nodes)

    # Changing the colour mode must not reconstruct the network: the proxy
    # updates only visual attributes for the nodes and links that are present.
    # Structural work is reserved for adding/removing selected pathways.
    if (structure_changed) {
      old_extra_ids <- extra_pathway_node_ids()
      old_extra_edge_ids <- extra_pathway_edge_ids()
      if (length(old_extra_edge_ids) > 0) {
        proxy <- proxy |> visRemoveEdges(old_extra_edge_ids)
      }
      if (length(old_extra_ids) > 0) {
        proxy <- proxy |> visRemoveNodes(old_extra_ids)
      }
    }

    if (length(selected) == 0) {
      extra_pathway_node_ids(character())
      extra_pathway_edge_ids(character())
      rendered_pathways(selected)
      return(invisible())
    }

    base_symbols <- base_data$nodes$label
    extra_symbols <- setdiff(unique(unlist(selected_sets$genes)), base_symbols)

    gene_attributes <- gene_data |>
      transmute(
        label = gene_name,
        significant = gene_id_clean %in% significant_ids,
        log2FoldChange
      ) |>
      distinct(label, .keep_all = TRUE)

    extra_genes <- tibble(label = extra_symbols) |>
      left_join(gene_attributes, by = "label") |>
      mutate(pathway_color = unname(pathway_gene_colors[label])) |>
      mutate(
        significant = coalesce(significant, FALSE),
        id = paste0("__PATHWAY_GENE__", label),
        group = "Selected pathway gene",
        shape = "dot",
        size = 10,
        color = if (color_by_direction()) {
          case_when(
            is.na(log2FoldChange) ~ hallmark_direction_colors[["unknown"]],
            log2FoldChange >= 0 ~ hallmark_direction_colors[["up"]],
            TRUE ~ hallmark_direction_colors[["down"]]
          )
        } else {
          pathway_color
        },
        borderWidth = if_else(significant, 3, 1),
        title = paste0(
          "<b>", htmlEscape(label), "</b><br>",
          "Selected Hallmark pathway gene<br>DE status: ",
          if_else(significant, "significant", "not significant"),
          "<br>Expression: ", case_when(
            is.na(log2FoldChange) ~ "unknown",
            log2FoldChange >= 0 ~ "up",
            TRUE ~ "down"
          )
        )
      ) |>
      select(id, label, group, shape, size, color, borderWidth, title)

    hubs <- selected_sets |>
      transmute(
        id = paste0("__PATHWAY_SET__", pathway),
        label = sub("^HALLMARK_", "", pathway),
        group = "Selected pathway",
        shape = "box",
        size = 18,
        color = pathway_color,
        borderWidth = 2,
        title = paste0("<b>", htmlEscape(label), "</b><br>Hallmark pathway")
      )

    if (nrow(extra_genes) + nrow(hubs) > 0) {
      proxy <- proxy |> visUpdateNodes(bind_rows(extra_genes, hubs))
    }

    pathway_edges <- lapply(seq_len(nrow(selected_sets)), function(i) {
      pathway <- selected_sets$pathway[i]
      symbols <- selected_sets$genes[[i]]
      target_ids <- if_else(
        symbols %in% base_symbols,
        base_data$nodes$id[match(symbols, base_symbols)],
        paste0("__PATHWAY_GENE__", symbols)
      )

      tibble(
        id = paste0("__PATHWAY_EDGE__", pathway, "__", symbols),
        from = paste0("__PATHWAY_SET__", pathway),
        to = target_ids,
        title = "Hallmark pathway membership",
        color = selected_sets$pathway_color[i],
        width = 1,
        dashes = TRUE
      )
    }) |>
      bind_rows()

    # Build the induced STRING network for every gene in the selected Hallmark
    # sets.  These are genuine protein/gene associations, unlike the dashed
    # membership edges above.
    gene_id_lookup <- bind_rows(
      gene_data |>
        transmute(gene_id = gene_id_clean, label = gene_name),
      tx2gene_names |>
        transmute(gene_id, label = gene_name)
    ) |>
      filter(!is.na(gene_id), gene_id != "", !is.na(label), label != "") |>
      distinct(label, .keep_all = TRUE)

    pathway_node_lookup <- tibble(label = unique(unlist(selected_sets$genes))) |>
      left_join(gene_id_lookup, by = "label") |>
      left_join(
        base_data$nodes |>
          transmute(label, base_node_id = id),
        by = "label"
      ) |>
      mutate(
        node_id = coalesce(base_node_id, paste0("__PATHWAY_GENE__", label))
      ) |>
      filter(!is.na(gene_id)) |>
      distinct(gene_id, .keep_all = TRUE)

    pathway_string_edges <- string_edges |>
      filter(
        source %in% pathway_node_lookup$gene_id,
        target %in% pathway_node_lookup$gene_id
      ) |>
      left_join(
        pathway_node_lookup |>
          transmute(source = gene_id, from = node_id),
        by = "source"
      ) |>
      left_join(
        pathway_node_lookup |>
          transmute(target = gene_id, to = node_id),
        by = "target"
      ) |>
      transmute(
        id = paste0("__PATHWAY_STRING__", source, "__", target),
        from,
        to,
        title = "STRING interaction between selected pathway genes",
        color = "rgba(15, 118, 110, 0.48)",
        width = 1.6,
        dashes = FALSE
      ) |>
      distinct(id, .keep_all = TRUE)

    dynamic_edges <- bind_rows(pathway_edges, pathway_string_edges)
    if (nrow(dynamic_edges) > 0) {
      proxy |> visUpdateEdges(dynamic_edges)
    }

    extra_pathway_node_ids(c(extra_genes$id, hubs$id))
    extra_pathway_edge_ids(dynamic_edges$id)
    rendered_pathways(selected)
  }, ignoreInit = TRUE)
}
