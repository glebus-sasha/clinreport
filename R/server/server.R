# Shiny server orchestration

# 22. SERVER
# ============================================================

server <- function(
    input,
    output,
    session
) {
  
  
  # ==========================================================
  # DRUG TABLE
  # ==========================================================
  
  output$drug_table <- renderDT({
    
    table_data <- drug_network |>
      transmute(
        
        `Drug` = label,
        
        `Status` = status,
        
        `Score` = score
      )
    
    
    datatable(
      
      table_data,
      
      selection = "single",
      
      rownames = FALSE,
      
      filter = "none",
      
      escape = FALSE,
      
      options = list(
        
        pageLength = 18,
        
        lengthChange = FALSE,
        
        searching = TRUE,
        
        scrollX = TRUE,
        
        autoWidth = FALSE,
        
        dom = "ftip",
        
        order = list(
          list(
            2,
            "desc"
          )
        ),
        
        columnDefs = list(
          
          list(
            width = "55%",
            targets = 0
          ),
          
          list(
            width = "25%",
            targets = 1
          ),
          
          list(
            width = "20%",
            targets = 2
          )
        )
      )
    ) |>
      formatStyle(
        "Drug",
        fontWeight = "600"
      ) |>
      formatStyle(
        "Status",
        color = "#66707a"
      ) |>
      formatRound(
        "Score",
        digits = 3
      )
  })
  
  
  # ==========================================================
  # SELECTED DRUG
  # ==========================================================
  
  selected_drug <- reactive({
    selected_row <- input$drug_table_rows_selected

    if (is.null(selected_row) || length(selected_row) == 0) {
      return(NULL)
    }

    drug_network[
      selected_row,
      ,
      drop = FALSE
    ]
  })
  
  
  # ==========================================================
  # CURRENT DRUG GENES
  # ==========================================================
  
  current_drug_genes <- reactive({
    
    req(
      selected_drug()
    )
    
    drug_id <- selected_drug()$drugId
    
    get_drug_genes(
      drug_id
    )
  })
  
  
  # ==========================================================
  # DRUG / INTERACTION-NETWORK SUBGRAPH
  # ==========================================================

  drug_interaction_subgraph <- reactive({
    req(selected_drug())

    get_drug_interaction_subgraph(
      selected_drug()$drugId,
      max_neighbors = interaction_network_max_neighbors
    )
  })


  # ==========================================================
  # VIEW MODE
  # ==========================================================
  
  view_mode <- reactiveVal(
    "drug"
  )
  
  selected_gene <- reactiveVal(
    NULL
  )

  selected_pathways <- reactiveVal(character())
  selected_tfs <- reactiveVal(character())
  target_connected_only <- reactive(!identical(input$target_connected_only, FALSE))
  show_context_genes <- reactive(isTRUE(input$show_context_genes))
  overlay_mode <- reactive(if (identical(input$overlay_mode, 'tf')) 'tf' else 'pathways')
  selected_network_gene <- reactiveVal(NULL)

  set_network_gene_focus <- function(gene_symbol) {
    gene_lookup <- bind_rows(
      gene_data |>
        transmute(gene_id_clean, gene_name),
      tx2gene_names |>
        transmute(gene_id_clean = gene_id, gene_name)
    ) |>
      filter(!is.na(gene_name), gene_name == gene_symbol) |>
      distinct(gene_name, .keep_all = TRUE)

    if (nrow(gene_lookup) == 0) {
      selected_network_gene(NULL)
      return(invisible())
    }

    selected_network_gene(list(
      id = gene_lookup$gene_id_clean[1],
      symbol = gene_symbol
    ))
  }

  significant_gene_ids <- gene_data |>
    filter(
      !is.na(padj),
      padj < padj_cutoff,
      !is.na(log2FoldChange),
      abs(log2FoldChange) >= log2fc_cutoff
    ) |>
    pull(gene_id_clean) |>
    unique()

  observeEvent(input$selected_gsea_pathway, {
    pathway <- input$selected_gsea_pathway
    selected <- selected_pathways()

    selected_pathways(
      if (pathway %in% selected) {
        setdiff(selected, pathway)
      } else {
        c(selected, pathway)
      }
    )
  })

  current_drug_target_symbols <- reactive({
    drug <- selected_drug()
    req(drug)

    get_gene_display_names(get_drug_target_ids(drug$drugId)) |>
      pull(gene_name) |>
      discard(is.na) |>
      unique()
  })

  observeEvent(selected_drug()$drugId, {
    targets <- current_drug_target_symbols()
    selected_pathways(character())
    selected_tfs(unique(get_target_tf_links(targets)$tf))
  }, priority = 10)

  observeEvent(input$select_all_pathways, {
    targets <- current_drug_target_symbols()
    selected_pathways(pathway_sets$pathway[vapply(pathway_sets$genes,
      function(genes) any(targets %in% genes), logical(1))])
  })

  register_tf_outputs(input, output, session, selected_tfs,
    target_connected_only = target_connected_only,
    drug_target_symbols = current_drug_target_symbols,
    selected_drug = selected_drug, selected_network_gene = selected_network_gene,
    significant_ids = significant_gene_ids)

  drug_table_proxy <- dataTableProxy("drug_table")

  observeEvent(input$selected_network_node, {
    node_id <- input$selected_network_node

    if (is.null(node_id) || node_id == "__DRUG__" || startsWith(node_id, "__PATHWAY_SET__")) {
      selected_network_gene(NULL)
      return()
    }

    gene_symbol <- if (startsWith(node_id, "__TF_GENE__")) {
      sub('^__TF_GENE__', '', node_id)
    } else if (startsWith(node_id, "__PATHWAY_GENE__")) {
      sub("^__PATHWAY_GENE__", "", node_id)
    } else {
      gene_data |>
        filter(gene_id_clean == node_id) |>
        pull(gene_name) |>
        first(default = NA_character_)
    }

    if (is.na(gene_symbol) || !nzchar(gene_symbol)) {
      selected_network_gene(NULL)
      return()
    }

    set_network_gene_focus(gene_symbol)
    selected_gene(gene_data |>
      filter(gene_name == gene_symbol) |>
      slice(1))
  })

  output$focused_drugs <- renderUI({
    gene <- selected_network_gene()
    req(!is.null(gene))

    related_drugs <- drug_network |>
      filter(vapply(
        drugId,
        function(drug_id) gene$id %in% get_drug_target_ids(drug_id),
        logical(1)
      )) |>
      arrange(desc(score), label)

    div(
      class = "focused-drugs-panel",
      div(class = "focused-drugs-title", paste0("Drugs targeting ", gene$symbol)),
      if (nrow(related_drugs) > 0) {
        div(
          class = "focused-drugs-list",
          lapply(seq_len(nrow(related_drugs)), function(i) {
            div(
              class = "focused-drug-chip",
              onclick = sprintf(
                "Shiny.setInputValue('selected_focused_drug_id', '%s', {priority: 'event'})",
                related_drugs$drugId[i]
              ),
              span(class = "focused-drug-name", related_drugs$label[i]),
              span(class = "focused-drug-status", related_drugs$status[i]),
              span(class = "focused-drug-score", sprintf("%.3f", related_drugs$score[i]))
            )
          })
        )
      } else {
        div(class = "focused-drugs-empty", "No direct drug targets in this dataset.")
      }
    )
  })

  observeEvent(input$clear_gsea_pathways, {
    selected_pathways(character())
    selected_network_gene(NULL)
  })

  observeEvent(input$selected_focused_drug_id, {
    row_index <- match(input$selected_focused_drug_id, drug_network$drugId)
    req(!is.na(row_index))
    selectRows(drug_table_proxy, row_index)
  })

  observeEvent(input$selected_network_gene_symbol, {
    set_network_gene_focus(input$selected_network_gene_symbol)
  })

  output$wes_variants_table <- renderDT({
    variant_table <- wes_variants |>
      mutate(
        classification_rank = case_when(
          classification %in% c("Oncogenic", "Likely Oncogenic", "Pathogenic", "Likely Pathogenic") ~ 1L,
          classification == "VUS" ~ 2L,
          TRUE ~ 3L
        )
      ) |>
      arrange(classification_rank, desc(impact), gene_symbol, position)

    escape_html <- function(x) as.character(htmltools::htmlEscape(x))
    shiny_button <- function(label, input_id, value, class_name) {
      paste0(
        "<button type='button' class='", class_name, "' onclick='",
        "Shiny.setInputValue(\"", input_id, "\", ",
        jsonlite::toJSON(value, auto_unbox = TRUE),
        ", {priority: \"event\"});'>", escape_html(label), "</button>"
      )
    }

    table_data <- lapply(seq_len(nrow(variant_table)), function(i) {
      variant <- variant_table[i, ]
      gene_url <- if (!is.na(variant$gene_id) && nzchar(variant$gene_id)) {
        ensembl_url(variant$gene_id)
      } else {
        ncbi_gene_search_url(variant$gene_symbol)
      }
      related_pathways <- pathway_sets |>
        filter(vapply(genes, function(gene_set) variant$gene_symbol %in% gene_set, logical(1)))
      pathway_html <- if (nrow(related_pathways) == 0) {
        "—"
      } else {
        paste(vapply(seq_len(nrow(related_pathways)), function(j) {
          pathway <- related_pathways[j, ]
          shiny_button(
            display_pathway_name(pathway$pathway),
            "selected_gsea_pathway",
            pathway$pathway,
            "wes-table-pathway"
          )
        }, character(1)), collapse = " ")
      }
      has_expression <- any(gene_data$gene_name == variant$gene_symbol, na.rm = TRUE)
      links <- c(
        if (has_expression) shiny_button("Gene profile", "selected_wes_expression_gene_symbol", variant$gene_symbol, "wes-table-profile"),
        paste0("<a href='", escape_html(gene_url), "' target='_blank' rel='noopener noreferrer'>Ensembl ↗</a>"),
        if (!is.null(clinvar_allele_url(variant$clinvar_allele_id))) {
          paste0("<a href='", escape_html(clinvar_allele_url(variant$clinvar_allele_id)), "' target='_blank' rel='noopener noreferrer'>ClinVar ↗</a>")
        },
        if (!is.null(dbsnp_url(variant$dbsnp_rsid))) {
          paste0("<a href='", escape_html(dbsnp_url(variant$dbsnp_rsid)), "' target='_blank' rel='noopener noreferrer'>dbSNP ↗</a>")
        }
      )

      data.frame(
        Gene = shiny_button(variant$gene_symbol, "selected_wes_gene_symbol", variant$gene_symbol, "wes-table-gene"),
        Variant = escape_html(variant$variant_label),
        Annotation = escape_html(if_else(is.na(variant$protein_change), variant$consequence, variant$protein_change)),
        Impact = paste0("<span class='wes-table-impact'>", escape_html(variant$impact), "</span>"),
        Classification = paste0("<span class='wes-table-classification'>", escape_html(variant$classification), "</span>"),
        VAF = if_else(is.na(variant$vaf), "—", sprintf("%.1f%%", variant$vaf * 100)),
        Depth = if_else(is.na(variant$depth), "—", as.character(variant$depth)),
        Pathways = pathway_html,
        Links = paste(links, collapse = " · "),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }) |>
      bind_rows()

    datatable(
      table_data,
      rownames = FALSE,
      escape = FALSE,
      filter = "top",
      class = "compact stripe hover",
      options = list(
        pageLength = 15,
        lengthMenu = list(c(15, 25, 50), c("15", "25", "50")),
        scrollX = TRUE,
        autoWidth = FALSE,
        order = list(list(4, "asc"), list(0, "asc")),
        columnDefs = list(
          list(targets = c(7, 8), orderable = FALSE)
        )
      )
    )
  })

  observeEvent(input$selected_wes_gene_symbol, {
    set_network_gene_focus(input$selected_wes_gene_symbol)
    focused_gene <- selected_network_gene()
    if (!is.null(focused_gene)) {
      session$sendCustomMessage(
        "focus-network-gene",
        list(id = focused_gene$id, symbol = focused_gene$symbol)
      )
    }
  })

  observeEvent(input$selected_wes_expression_gene_symbol, {
    expression_gene <- gene_data |>
      filter(gene_name == input$selected_wes_expression_gene_symbol) |>
      slice(1)
    req(nrow(expression_gene) > 0)

    selected_gene(expression_gene)
  })
  
  
  # ==========================================================
  # RIGHT PANEL
  # ==========================================================
  
  output$right_panel <- renderUI({
    drug <- selected_drug()

    if (is.null(drug) || nrow(drug) == 0) {
      return(render_input_overview())
    }
    
    
    drug_id <- drug$drugId
    
    drug_name <- drug$label
    
    
    data <- load_drug(
      drug_id
    )
    
    
    genes <- current_drug_genes()
    
    
    has_json <- any(
      map_int(
        data,
        length
      ) > 0
    )
    
    
    if (!has_json) {
      
      return(
        
        tagList(
          
          div(
            class = "drug-header",
            
            div(
              class = "drug-title",
              drug_name
            ),
            
            div(
              class = "drug-id",
              drug_id
            )
          ),
          
          div(
            class = "empty-message",
            
            "Detailed information is not available ",
            "for this drug."
          ),

          div(
            class = "network-gsea-layout",
            div(
              class = "drug-network-panel top-network-panel",
              div(
                class = "network-visualization-header",
                div(class = "network-visualization-title", network_display_name),
                div(
                  class = "network-visualization-subtitle",
                  paste("Drug targets, differential expression, and", pathway_collection_name, "pathway coverage")
                )
              ),
              uiOutput("drug_network_summary"),
              visNetworkOutput("drug_network_graph", height = "500px")
            ),
            div(
              class = "gsea-panel",
              render_network_overlay_panel(),
              conditionalPanel("input.overlay_mode !== 'tf'", uiOutput("gsea_pathway_details")),
              tabsetPanel(render_tf_tab())
            )
          )
        )
      )
    }
    
    
    render_drug_page(
      
      drug_id,
      
      drug_name,
      
      data,
      
      genes,

      network_display_name
    )
  })
  
  
  # ==========================================================
  # NEW DRUG SELECTED
  # ==========================================================
  
  observeEvent(
    
    input$drug_table_rows_selected,
    
    {
      
      view_mode(
        "drug"
      )
      
      selected_gene(
        NULL
      )

      selected_network_gene(NULL)
      
    },
    
    ignoreInit = TRUE
  )

  output$inline_gene_panel <- renderUI({
    gene <- selected_gene()
    if (is.null(gene)) return(NULL)
    div(class = "inline-gene-panel", render_gene_page(gene, inline = TRUE))
  })
  observeEvent(input$close_gene_profile, selected_gene(NULL))
  
  
  # ==========================================================
  # GENE CLICK
  # ==========================================================
  
  observeEvent(
    
    input$selected_gene_index,
    
    {
      
      req(
        selected_drug()
      )
      
      
      genes <- current_drug_genes()
      
      
      i <- as.integer(
        input$selected_gene_index
      )
      
      
      if (
        is.na(i) ||
        i < 1 ||
        i > nrow(genes)
      ) {
        
        return()
      }
      
      
      selected_gene(
        genes[
          i,
          ,
          drop = FALSE
        ]
      )
      
      
    }
  )
  
  
  # ==========================================================
  # BACK TO DRUG
  # ==========================================================
  
  observeEvent(
    
    input$back_to_drug,
    
    {
      
      selected_gene(
        NULL
      )
      
    }
  )
  
  
  color_by_direction <- reactive(isTRUE(input$gsea_color_by_direction))

  register_drug_network_outputs(
    output = output,
    input = input,
    selected_drug = selected_drug,
    subgraph = drug_interaction_subgraph,
    selected_pathways = selected_pathways,
    significant_ids = significant_gene_ids,
    selected_network_gene = selected_network_gene,
    color_by_direction = color_by_direction,
    selected_tfs = selected_tfs, overlay_mode = overlay_mode,
    target_connected_only = target_connected_only,
    show_context_genes = show_context_genes
  )

  register_gsea_outputs(
    output = output,
    input = input,
    selected_pathways = selected_pathways,
    significant_ids = significant_gene_ids,
    subgraph = drug_interaction_subgraph,
    selected_network_gene = selected_network_gene,
    selected_drug = selected_drug,
    drug_target_symbols = current_drug_target_symbols
  )


  # ==========================================================
  # VOLCANO
  # ==========================================================
  
  # Volcano plot output is registered separately.
  register_volcano_output(
    output = output,
    selected_gene = selected_gene,
    gene_data = gene_data,
    padj_cutoff = padj_cutoff,
    log2fc_cutoff = log2fc_cutoff
  )
}
