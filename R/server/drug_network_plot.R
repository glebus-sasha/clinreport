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
      label = if_else(is.na(gene_name) | gene_name == "", gene_id, gene_name)
    ) |>
    left_join(
      gene_data |>
        transmute(label = gene_name, log2FoldChange_by_name = log2FoldChange) |>
        distinct(label, .keep_all = TRUE),
      by = "label"
    ) |>
    mutate(log2FoldChange = coalesce(log2FoldChange, log2FoldChange_by_name)) |>
    left_join(
      wes_gene_summary |>
        rename(label = gene_symbol),
      by = "label"
    ) |>
    mutate(
      pathway_color = unname(pathway_gene_colors[label]),
      pathway_member = !is.na(pathway_color),
      significant = gene_id %in% significant_ids,
      mutation_count = coalesce(mutation_count, 0L),
      mutation_classes = coalesce(mutation_classes, "none"),
      mutated = mutation_count > 0,
      group = if_else(type == "Network neighbor", paste(interaction_network_name, "neighbor"), type),
      shape = "dot",
      size = case_when(
        type == "Target" ~ 22,
        pathway_member ~ 17,
        TRUE ~ 11
      ),
      color = if (color_by_direction) {
        case_when(
          is.na(log2FoldChange) ~ pathway_direction_colors[["unknown"]],
          log2FoldChange >= 0 ~ pathway_direction_colors[["up"]],
          TRUE ~ pathway_direction_colors[["down"]]
        )
      } else {
        case_when(
          pathway_member ~ pathway_color,
          type == "Target" ~ "#2563eb",
          TRUE ~ "#94a3b8"
        )
      },
      borderWidth = case_when(
        significant ~ 4,
        TRUE ~ 1
      ),
      title = paste0(
        "<b>", htmlEscape(label), "</b><br>",
        "Ensembl: ", htmlEscape(gene_id), "<br>",
        if_else(type == "Target", "Drug target", if_else(type == "Network neighbor", paste(interaction_network_name, "neighbor"), type)),
        "<br>DE status: ", if_else(significant, "significant", "not significant"),
        "<br>Expression: ", case_when(
          is.na(log2FoldChange) ~ "unknown",
          log2FoldChange >= 0 ~ "up",
          TRUE ~ "down"
        ),
        "<br>Selected ", pathway_collection_name, " pathway: ", if_else(pathway_member, "member", "not a member"),
        "<br>WES variants: ", mutation_count, " (", htmlEscape(mutation_classes), ")"
      )
    ) |>
    transmute(id = gene_id, label, group, shape, size, color, borderWidth, mutated, title)

  nodes <- bind_rows(
    tibble(
      id = "__DRUG__",
      label = drug$label,
      group = "Drug",
      shape = "diamond",
      size = 30,
      color = "#17202a",
      mutated = FALSE,
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
        title = paste(interaction_network_name, "interaction"),
        color = "#cbd5e1",
        width = 1,
        dashes = TRUE
      )
  ) |>
    filter(from %in% nodes$id, to %in% nodes$id) |>
    distinct(from, to, .keep_all = TRUE)

  edges$id <- paste("__BASE__", edges$from, edges$to, sep = "::")
  edges$arrows.to.enabled <- FALSE
  edges$arrows.to.type <- "arrow"
  edges$smooth.enabled <- FALSE
  edges$smooth.type <- "curvedCW"
  edges$smooth.roundness <- 0.12
  list(nodes = nodes, edges = edges)
}

register_drug_network_outputs <- function(
    output,
    input,
    selected_drug,
    subgraph,
    selected_pathways,
    significant_ids,
    selected_network_gene,
    color_by_direction,
    selected_tfs,
    overlay_mode,
    target_connected_only,
    show_context_genes
) {
  extra_pathway_node_ids <- reactiveVal(character())
  extra_pathway_edge_ids <- reactiveVal(character())
  rendered_pathways <- reactiveVal(character())
  graph_generation <- reactiveVal(0L)
  visible_wes_count <- reactiveVal(0L)

  output$drug_network_summary <- renderUI({
    drug <- selected_drug()

    if (is.null(drug) || nrow(drug) == 0) {
      return(
        div(
          class = "network-empty-message",
          paste("Select a drug in the table to display its", interaction_network_name, "subgraph.")
        )
      )
    }

    graph <- subgraph()
    selected <- if (overlay_mode() == "tf") character() else selected_pathways()
    target_nodes <- graph$nodes |>
      filter(type == "Target")
    neighbor_count <- sum(graph$nodes$type != "Target")
    target_symbols <- target_nodes |>
      pull(gene_name) |>
      unique() |>
      discard(is.na)
    target_tf_count <- if (has_tf_analysis) sum(target_symbols %in% tf_results$TF) else 0L
    wes_target_count <- if (has_wes) sum(target_symbols %in% wes_gene_summary$gene_symbol) else 0L
    pathway_count <- sum(vapply(pathway_sets$genes,
      function(gs) any(target_symbols %in% gs), logical(1)))
    div(
      class = "network-summary",
      span(class = "network-summary-drug", drug$label),
      span(paste0(" · ", nrow(target_nodes), " drug-linked targets · ",
        target_tf_count, " target TFs · ", wes_target_count, " WES targets · ",
        pathway_count, " connected ", pathway_collection_name, " pathways")),
      if (overlay_mode() == "tf" && has_tf_analysis)
        span(class = "network-summary-tf",
          paste0(" · ", n_distinct(get_target_tf_links(target_symbols)$tf),
            " linked TF regulators · ", length(selected_tfs()), " selected TFs")),
      if (show_context_genes()) span(paste0(" · ", neighbor_count, " context genes")),
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
    data <- apply_context_visibility(data, graph, show_context = isolate(show_context_genes()))
    visible_wes_count(count_visible_wes_genes(data$nodes))

    # Proxy messages are safe only after this particular widget has mounted.
    # A generation also rejects late ready events from the previous drug.
    generation <- isolate(graph_generation()) + 1L
    graph_generation(generation)
    extra_pathway_node_ids(character())
    extra_pathway_edge_ids(character())
    rendered_pathways(character())

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
        ),
        afterDrawing = JS(
          "function(ctx) {
             var nodes = this.body.nodes;
             Object.keys(nodes).forEach(function(id) {
               var node = nodes[id];
               if (!node || !node.options || node.options.hidden || !node.options.mutated) return;
               ctx.save();
               ctx.beginPath();
               ctx.arc(node.x, node.y, 3.4, 0, 2 * Math.PI);
               ctx.fillStyle = '#be185d';
               ctx.fill();
               ctx.lineWidth = 1.2;
               ctx.strokeStyle = '#ffffff';
               ctx.stroke();
               ctx.restore();
             });
           }"
        )
      ) |>
      htmlwidgets::onRender(JS(sprintf(
        "function(el, x) {
           Shiny.setInputValue('drug_network_ready', %d, {priority: 'event'});
         }", generation
      )))
  })

  observeEvent(list(input$drug_network_ready, graph_generation(), selected_pathways(), selected_tfs(), overlay_mode(), color_by_direction(), target_connected_only(), show_context_genes(), subgraph()), {
    req(graph_generation() > 0L,
        identical(as.integer(input$drug_network_ready), graph_generation()))
    drug <- selected_drug()
    req(drug)

    graph <- subgraph()
    selected <- selected_pathways()
    key <- c(drug$drugId, overlay_mode(), as.character(target_connected_only()),
             if (overlay_mode() == "tf") selected_tfs() else selected)
    structure_changed <- !identical(key, rendered_pathways())
    if (overlay_mode() == "tf") {
      data <- build_tf_network(graph, drug, selected_tfs(), significant_ids, color_by_direction(),
                               target_connected_only())
      overlay_ids <- unique(c(data$tf_edges$from, data$tf_edges$to,
                             data$nodes$id[data$nodes$shape == "triangle"]))
      data <- apply_context_visibility(data, graph,
        data$nodes$label[data$nodes$id %in% overlay_ids], show_context_genes())
      visible_wes_count(count_visible_wes_genes(data$nodes))
      proxy <- visNetworkProxy("drug_network_graph")
      if (structure_changed && length(extra_pathway_edge_ids())) proxy <- visRemoveEdges(proxy, extra_pathway_edge_ids())
      if (structure_changed && length(extra_pathway_node_ids())) proxy <- visRemoveNodes(proxy, extra_pathway_node_ids())
      proxy <- visUpdateNodes(proxy, data$nodes)
      proxy <- visUpdateEdges(proxy, data$edges)
      extra_pathway_node_ids(data$extra_ids)
      extra_pathway_edge_ids(data$tf_edges$id)
      rendered_pathways(key)
      return(invisible())
    }
    selected_sets <- pathway_sets |>
      filter(pathway %in% selected) |>
      mutate(
        pathway_color = if (color_by_direction()) {
          pathway_direction_color(nes)
        } else {
          pathway_identity_color(pathway)
        }
      )
    if (target_connected_only()) {
      selected_sets$genes <- lapply(selected_sets$genes, intersect, y = graph$nodes$gene_name)
    }
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
    base_data <- apply_context_visibility(base_data, graph,
      unique(unlist(selected_sets$genes)), show_context_genes())

    nodes <- base_data$nodes |>
      select(id, shape, size, color, borderWidth, mutated, title, hidden, physics)

    proxy <- visNetworkProxy("drug_network_graph") |>
      visUpdateNodes(nodes) |>
      visUpdateEdges(base_data$edges)

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
      visible_wes_count(count_visible_wes_genes(base_data$nodes))
      extra_pathway_node_ids(character())
      extra_pathway_edge_ids(character())
      rendered_pathways(key)
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
      left_join(wes_gene_summary, by = c("label" = "gene_symbol")) |>
      mutate(pathway_color = unname(pathway_gene_colors[label])) |>
      mutate(
        significant = coalesce(significant, FALSE),
        mutation_count = coalesce(mutation_count, 0L),
        mutation_classes = coalesce(mutation_classes, "none"),
        mutated = mutation_count > 0,
        id = paste0("__PATHWAY_GENE__", label),
        group = "Selected pathway gene",
        shape = "dot",
        size = 10,
        color = if (color_by_direction()) {
          case_when(
            is.na(log2FoldChange) ~ pathway_direction_colors[["unknown"]],
            log2FoldChange >= 0 ~ pathway_direction_colors[["up"]],
            TRUE ~ pathway_direction_colors[["down"]]
          )
        } else {
          pathway_color
        },
        borderWidth = case_when(
          significant ~ 3,
          TRUE ~ 1
        ),
        title = paste0(
          "<b>", htmlEscape(label), "</b><br>",
          paste0("Selected ", pathway_collection_name, " pathway gene<br>DE status: "),
          if_else(significant, "significant", "not significant"),
          "<br>Expression: ", case_when(
            is.na(log2FoldChange) ~ "unknown",
            log2FoldChange >= 0 ~ "up",
            TRUE ~ "down"
          ),
          "<br>WES variants: ", mutation_count, " (", htmlEscape(mutation_classes), ")"
        )
      ) |>
      select(id, label, group, shape, size, color, borderWidth, mutated, title)

    hubs <- selected_sets |>
      transmute(
        id = paste0("__PATHWAY_SET__", pathway),
        label = display_pathway_name(pathway),
        group = "Selected pathway",
        shape = "box",
        size = 18,
        color = pathway_color,
        borderWidth = 2,
        mutated = FALSE,
        title = paste0("<b>", htmlEscape(label), "</b><br>", pathway_collection_name, " pathway")
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
        title = paste(pathway_collection_name, "pathway membership"),
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

    pathway_interaction_edges <- interaction_edges |>
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
        id = paste0("__PATHWAY_INTERACTION__", source, "__", target),
        from,
        to,
        title = paste(interaction_network_name, "interaction between selected pathway genes"),
        color = "rgba(15, 118, 110, 0.48)",
        width = 1.6,
        dashes = FALSE
      ) |>
      distinct(id, .keep_all = TRUE)

    dynamic_edges <- bind_rows(pathway_edges, pathway_interaction_edges)
    if (nrow(dynamic_edges) > 0) {
      proxy |> visUpdateEdges(dynamic_edges)
    }

    extra_pathway_node_ids(c(extra_genes$id, hubs$id))
    visible_wes_count(count_visible_wes_genes(bind_rows(base_data$nodes, extra_genes)))
    extra_pathway_edge_ids(dynamic_edges$id)
    rendered_pathways(key)
  }, ignoreInit = FALSE)
}
