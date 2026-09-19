# Hallmark GSEA pathway controls

register_gsea_outputs <- function(
    output,
    input,
    selected_pathways,
    significant_ids,
    subgraph,
    selected_network_gene,
    selected_drug,
    drug_target_symbols
) {
  sort_pathways <- function(pathway_table) {
    if (identical(input$gsea_sort, "fdr")) {
      pathway_table |>
        arrange(fdr_q, desc(abs(nes)), pathway)
    } else {
      pathway_table |>
        arrange(desc(abs(nes)), fdr_q, pathway)
    }
  }

  output$gsea_pathway_tiles <- renderUI({
    active <- selected_pathways()
    focused_gene <- selected_network_gene()
    pathway_table <- sort_pathways(hallmark_pathways)

    lapply(seq_len(nrow(pathway_table)), function(i) {
      pathway <- pathway_table[i, ]
      selected_class <- if (pathway$pathway %in% active) " active" else ""
      direction_class <- if (is.na(pathway$nes) || pathway$nes >= 0) " up" else " down"
      click <- sprintf(
        "Shiny.setInputValue('selected_gsea_pathway', '%s', {priority: 'event'})",
        pathway$pathway
      )

      div(
        class = paste0("gsea-pathway-tile", selected_class, direction_class),
        onclick = click,
        div(class = "gsea-pathway-name", sub("^HALLMARK_", "", pathway$pathway)),
        div(
          class = "gsea-pathway-metrics",
          span(ifelse(is.na(pathway$nes), "NES —", sprintf("NES %+.2f", pathway$nes))),
          span(ifelse(is.na(pathway$fdr_q), "FDR —", sprintf("FDR %.2g", pathway$fdr_q)))
        )
      )
    })
  })

  output$gsea_matching_pathways <- renderUI({
    focused_gene <- selected_network_gene()
    drug <- selected_drug()
    target_symbols <- drug_target_symbols()

    make_pathway_tiles <- function(pathway_table) {
      lapply(seq_len(nrow(pathway_table)), function(i) {
        pathway <- pathway_table[i, ]
        selected_class <- if (pathway$pathway %in% selected_pathways()) " active" else ""
        direction_class <- if (is.na(pathway$nes) || pathway$nes >= 0) " up" else " down"
        click <- sprintf(
          "Shiny.setInputValue('selected_gsea_pathway', '%s', {priority: 'event'})",
          pathway$pathway
        )

        div(
          class = paste0("gsea-pathway-tile", selected_class, direction_class),
          onclick = click,
          div(class = "gsea-pathway-name", sub("^HALLMARK_", "", pathway$pathway)),
          div(
            class = "gsea-pathway-metrics",
            span(sprintf("NES %+.2f", pathway$nes)),
            span(sprintf("FDR %.2g", pathway$fdr_q))
          )
        )
      })
    }

    drug_pathways <- hallmark_pathways |>
      filter(vapply(genes, function(gene_set) {
        any(target_symbols %in% gene_set)
      }, logical(1))) |>
      sort_pathways()

    gene_pathways <- if (!is.null(focused_gene)) {
      hallmark_pathways |>
        filter(vapply(genes, function(gene_set) {
          focused_gene$symbol %in% gene_set
        }, logical(1))) |>
        sort_pathways()
    } else {
      hallmark_pathways[0, ]
    }

    req(nrow(drug_pathways) > 0 || nrow(gene_pathways) > 0)

    div(
      class = "gsea-matching-pathways",
      if (nrow(drug_pathways) > 0) tagList(
        div(
          class = "gsea-matching-title",
          paste0("Pathways containing targets of ", drug$label)
        ),
        make_pathway_tiles(drug_pathways)
      ),
      if (nrow(gene_pathways) > 0) tagList(
        div(
          class = "gsea-matching-title",
          paste0("Pathways containing ", focused_gene$symbol)
        ),
        make_pathway_tiles(gene_pathways)
      )
    )
  })

  output$gsea_pathway_details <- renderUI({
    selected <- selected_pathways()
    req(length(selected) > 0)
    focused_gene <- selected_network_gene()

    network_symbols <- subgraph()$nodes$gene_name

    lapply(selected, function(pathway_name) {
      pathway <- hallmark_pathways |>
        filter(pathway == pathway_name) |>
        slice(1)

      genes <- pathway$genes[[1]]
      gene_status <- tibble(gene_name = genes) |>
        left_join(
          gene_data |>
            transmute(
              gene_name,
              significant = gene_id_clean %in% significant_ids
            ) |>
            distinct(gene_name, .keep_all = TRUE),
          by = "gene_name"
        ) |>
        mutate(
          significant = coalesce(significant, FALSE),
          in_network = gene_name %in% network_symbols,
          chip_class = case_when(
            significant & in_network ~ "deg-network",
            significant ~ "deg-only",
            in_network ~ "network-only",
            TRUE ~ "pathway-only"
          )
        )

      div(
        class = "gsea-details",
        div(
          class = "gsea-details-heading",
          sub("^HALLMARK_", "", pathway$pathway)
        ),
        div(
          class = "gsea-details-stats",
          span(sprintf("NES %+.2f", pathway$nes)),
          span(sprintf("FDR %.3g", pathway$fdr_q)),
          span(paste0("set size ", pathway$size)),
          span(paste0(sum(gene_status$significant), " DE-significant"))
        ),
        tags$a(
          href = pathway$msigdb_url,
          target = "_blank",
          rel = "noopener noreferrer",
          class = "gsea-msigdb-link",
          "Open in MSigDB ↗"
        ),
        div(
          class = "gsea-gene-list",
          lapply(seq_len(nrow(gene_status)), function(i) {
            span(
              class = paste("gsea-gene-chip", gene_status$chip_class[i]),
              onclick = sprintf(
                "Shiny.setInputValue('selected_network_gene_symbol', '%s', {priority: 'event'})",
                gene_status$gene_name[i]
              ),
              title = paste(
                gene_status$gene_name[i],
                if (gene_status$significant[i]) "DE-significant" else "not DE-significant",
                if (gene_status$in_network[i]) "in drug network" else "not in drug network",
                sep = " · "
              ),
              gene_status$gene_name[i]
            )
          })
        )
      )
    })
  })
}
