# Hallmark GSEA pathway controls

register_gsea_outputs <- function(output, selected_pathways, significant_ids, subgraph) {
  output$gsea_pathway_tiles <- renderUI({
    active <- selected_pathways()

    lapply(seq_len(nrow(hallmark_pathways)), function(i) {
      pathway <- hallmark_pathways[i, ]
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

  output$gsea_pathway_details <- renderUI({
    selected <- selected_pathways()
    req(length(selected) > 0)

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
