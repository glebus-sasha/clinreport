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
    
    req(
      input$drug_table_rows_selected
    )
    
    drug_network[
      input$drug_table_rows_selected,
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
  # DRUG / STRING SUBGRAPH
  # ==========================================================

  drug_string_subgraph <- reactive({
    req(selected_drug())

    get_drug_string_subgraph(
      selected_drug()$drugId,
      max_neighbors = network_max_neighbors
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
  
  
  # ==========================================================
  # RIGHT PANEL
  # ==========================================================
  
  output$right_panel <- renderUI({
    
    req(
      selected_drug()
    )
    
    
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
    
    
    drug <- selected_drug()
    
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
                  "Selected drug, direct targets, and one-hop STRING neighbors"
                )
              ),
              uiOutput("drug_network_summary"),
              visNetworkOutput("drug_network_graph", height = "500px")
            ),
            div(
              class = "gsea-panel",
              div(class = "gsea-panel-title", "Hallmark pathways"),
              div(class = "gsea-panel-subtitle", "All tested gene sets · click to highlight DE genes"),
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
  
  
  register_drug_network_outputs(
    output = output,
    selected_drug = selected_drug,
    subgraph = drug_string_subgraph,
    selected_pathways = selected_pathways,
    significant_ids = significant_gene_ids
  )

  register_gsea_outputs(
    output = output,
    selected_pathways = selected_pathways,
    significant_ids = significant_gene_ids,
    subgraph = drug_string_subgraph
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
