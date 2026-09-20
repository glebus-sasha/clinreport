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

  drug_table_proxy <- dataTableProxy("drug_table")

  observeEvent(input$selected_network_node, {
    node_id <- input$selected_network_node

    if (is.null(node_id) || node_id == "__DRUG__" || startsWith(node_id, "__PATHWAY_SET__")) {
      selected_network_gene(NULL)
      return()
    }

    gene_symbol <- if (startsWith(node_id, "__PATHWAY_GENE__")) {
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
  
  
  # ==========================================================
  # RIGHT PANEL
  # ==========================================================
  
  output$right_panel <- renderUI({
    drug <- selected_drug()

    if (is.null(drug) || nrow(drug) == 0) {
      return(render_input_overview())
    }
    
    
    if (
      view_mode() == "gene"
    ) {
      
      gene <- selected_gene()
      
      req(
        gene
      )
      
      return(
        render_gene_page(
          gene
        )
      )
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
              div(class = "gsea-panel-title", paste(pathway_collection_name, "pathways")),
              div(class = "gsea-panel-subtitle", paste("All tested", pathway_collection_name, "gene sets · click to highlight DE genes")),
              div(class = "gsea-pathway-list", uiOutput("gsea_pathway_tiles")),
              uiOutput("gsea_pathway_details")
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

      selected_pathways(character())
      selected_network_gene(NULL)
      
    },
    
    ignoreInit = TRUE
  )
  
  
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
      
      
      view_mode(
        "gene"
      )
      
    }
  )
  
  
  # ==========================================================
  # BACK TO DRUG
  # ==========================================================
  
  observeEvent(
    
    input$back_to_drug,
    
    {
      
      view_mode(
        "drug"
      )
      
      selected_gene(
        NULL
      )
      
    }
  )
  
  
  color_by_direction <- reactive(isTRUE(input$gsea_color_by_direction))

  register_drug_network_outputs(
    output = output,
    selected_drug = selected_drug,
    subgraph = drug_interaction_subgraph,
    selected_pathways = selected_pathways,
    significant_ids = significant_gene_ids,
    selected_network_gene = selected_network_gene,
    color_by_direction = color_by_direction
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
