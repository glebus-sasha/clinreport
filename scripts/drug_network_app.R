library(shiny)
library(tidyverse)
library(visNetwork)
library(glue)

source("load_data.R")
source("network.R")


# =========================================================
# Paths
# =========================================================

RAW_DIR <- "../raw"

DRUG_FILE <- file.path(
  RAW_DIR,
  "all_samples_string_human_links_v12_0_min900_Ensembl_diamond_trustrank.csv"
)

PROTEIN_FILE <- file.path(
  RAW_DIR,
  "/network/nedrex.reviewed_proteins_exp.Ensembl.edges.tsv"
)

GENE_FILE <- file.path(
  RAW_DIR,
  "tx2gene.tsv"
)


# =========================================================
# Load data
# =========================================================

message("Loading gene symbols...")

gene_symbols <- load_gene_symbols(
  GENE_FILE
)

message(
  "Gene symbols: ",
  format(nrow(gene_symbols), big.mark = ",")
)


message("Loading protein network...")

protein_edges <- load_protein_edges(
  PROTEIN_FILE
)

message(
  "Protein edges: ",
  format(nrow(protein_edges), big.mark = ",")
)


message("Loading drug-target data...")

drug_targets <- load_drug_targets(
  DRUG_FILE
)

message(
  "Drug-target edges: ",
  format(nrow(drug_targets), big.mark = ",")
)


message("Building network...")

network <- build_network(
  drug_targets = drug_targets,
  protein_edges = protein_edges,
  gene_symbols = gene_symbols
)


message(
  "Network nodes: ",
  format(nrow(network$nodes), big.mark = ",")
)

message(
  "Network edges: ",
  format(nrow(network$edges), big.mark = ",")
)


# =========================================================
# UI
# =========================================================

ui <- fluidPage(
  
  tags$head(
    
    tags$style(
      HTML("
        body {
          overflow: hidden;
        }

        .sidebar {
          height: calc(100vh - 80px);
          overflow-y: auto;
        }

        #network {
          height: calc(100vh - 80px) !important;
        }

        .info-box {
          padding: 10px;
          background: #f7f7f7;
          border-radius: 6px;
          margin-top: 10px;
        }
      ")
    )
  ),
  
  titlePanel(
    "Drug–Target Network Explorer"
  ),
  
  sidebarLayout(
    
    sidebarPanel(
      class = "sidebar",
      
      selectizeInput(
        inputId = "drug",
        label = "Drug",
        choices = NULL,
        multiple = FALSE,
        options = list(
          placeholder = "Search for a drug..."
        )
      ),
      
      sliderInput(
        inputId = "hops",
        label = "Network depth",
        min = 0,
        max = 2,
        value = 1,
        step = 1
      ),
      
      checkboxInput(
        inputId = "show_proteins",
        label = "Show protein–protein network",
        value = FALSE
      ),
      
      sliderInput(
        inputId = "max_nodes",
        label = "Maximum nodes",
        min = 50,
        max = 2000,
        value = 500,
        step = 50
      ),
      
      hr(),
      
      h4("Selected drug"),
      
      div(
        class = "info-box",
        
        verbatimTextOutput(
          "drug_info"
        )
      ),
      
      hr(),
      
      h4("Network statistics"),
      
      div(
        class = "info-box",
        
        verbatimTextOutput(
          "network_info"
        )
      )
    ),
    
    mainPanel(
      
      visNetworkOutput(
        outputId = "network",
        width = "100%",
        height = "100%"
      )
    )
  )
)


# =========================================================
# Server
# =========================================================

server <- function(
    input,
    output,
    session
) {
  
  # -------------------------------------------------------
  # Drug selector
  # -------------------------------------------------------
  
  observe({
    
    drugs <- network$nodes |>
      filter(type == "drug") |>
      arrange(label)
    
    choices <- setNames(
      drugs$id,
      drugs$label
    )
    
    updateSelectizeInput(
      session = session,
      inputId = "drug",
      choices = choices,
      selected = drugs$id[[1]],
      server = TRUE
    )
  })
  
  
  # -------------------------------------------------------
  # Selected drug
  # -------------------------------------------------------
  
  selected_drug <- reactive({
    
    req(input$drug)
    
    network$nodes |>
      filter(
        id == input$drug,
        type == "drug"
      )
  })
  
  
  # -------------------------------------------------------
  # Selected subnetwork
  # -------------------------------------------------------
  
  selected_network <- reactive({
    
    req(input$drug)
    
    get_drug_network(
      network = network,
      drug_id = input$drug,
      hops = input$hops,
      show_proteins = input$show_proteins,
      max_nodes = input$max_nodes
    )
  })
  
  
  # -------------------------------------------------------
  # Drug information
  # -------------------------------------------------------
  
  output$drug_info <- renderText({
    
    drug <- selected_drug()
    
    req(nrow(drug) == 1)
    
    glue(
      "Name:       {drug$label}\n",
      "DrugBank:   {drug$id}\n",
      "Status:     {drug$status}\n",
      "Score:      {round(drug$score, 4)}"
    )
  })
  
  
  # -------------------------------------------------------
  # Network information
  # -------------------------------------------------------
  
  output$network_info <- renderText({
    
    graph <- selected_network()
    
    req(nrow(graph$nodes) > 0)
    
    n_drugs <- graph$nodes |>
      filter(type == "drug") |>
      nrow()
    
    n_targets <- graph$nodes |>
      filter(type == "target") |>
      nrow()
    
    n_proteins <- graph$nodes |>
      filter(type == "protein") |>
      nrow()
    
    glue(
      "Nodes:       {nrow(graph$nodes)}\n",
      "Edges:       {nrow(graph$edges)}\n",
      "Drugs:       {n_drugs}\n",
      "Targets:     {n_targets}\n",
      "Proteins:    {n_proteins}"
    )
  })
  
  
  # -------------------------------------------------------
  # Network visualization
  # -------------------------------------------------------
  
  output$network <- renderVisNetwork({
    
    graph <- selected_network()
    
    req(nrow(graph$nodes) > 0)
    
    # -----------------------------------------------------
    # Nodes
    # -----------------------------------------------------
    
    nodes <- graph$nodes |>
      mutate(
        
        color = case_when(
          
          type == "drug" ~ "#E74C3C",
          
          type == "target" ~ "#3498DB",
          
          type == "protein" ~ "#95A5A6",
          
          TRUE ~ "#BDC3C7"
        ),
        
        shape = case_when(
          
          type == "drug" ~ "diamond",
          
          type == "target" ~ "dot",
          
          type == "protein" ~ "dot",
          
          TRUE ~ "dot"
        ),
        
        size = case_when(
          
          type == "drug" ~ 35,
          
          type == "target" ~ 22,
          
          type == "protein" ~ 10,
          
          TRUE ~ 10
        ),
        
        title = case_when(
          
          type == "drug" ~ paste0(
            "<b>",
            label,
            "</b><br>",
            "DrugBank: ",
            id,
            "<br>",
            "Status: ",
            status,
            "<br>",
            "Score: ",
            round(score, 4)
          ),
          
          TRUE ~ paste0(
            "<b>",
            label,
            "</b><br>",
            "Ensembl: ",
            id
          )
        )
      ) |>
      select(
        id,
        label,
        group = type,
        color,
        shape,
        size,
        title
      )
    
    
    # -----------------------------------------------------
    # Edges
    # -----------------------------------------------------
    
    edges <- graph$edges |>
      transmute(
        from = source,
        to = target,
        title = type
      )
    
    
    # -----------------------------------------------------
    # Visualization
    # -----------------------------------------------------
    
    visNetwork(
      nodes = nodes,
      edges = edges
    ) |>
      
      visNodes(
        borderWidth = 1,
        font = list(
          size = 16
        )
      ) |>
      
      visEdges(
        smooth = FALSE,
        color = list(
          color = "#BBBBBB",
          highlight = "#333333"
        )
      ) |>
      
      visOptions(
        highlightNearest = list(
          enabled = TRUE,
          degree = 1,
          hover = TRUE
        ),
        nodesIdSelection = list(
          enabled = TRUE,
          useLabels = TRUE
        )
      ) |>
      
      visInteraction(
        hover = TRUE,
        navigationButtons = TRUE,
        keyboard = TRUE
      ) |>
      
      visPhysics(
        enabled = TRUE,
        solver = "forceAtlas2Based",
        
        forceAtlas2Based = list(
          gravitationalConstant = -100,
          centralGravity = 0.02,
          springLength = 120,
          springConstant = 0.08,
          damping = 0.4,
          avoidOverlap = 1
        ),
        
        stabilization = list(
          enabled = TRUE,
          iterations = 300
        )
      )
  })
}


# =========================================================
# Run application
# =========================================================

shinyApp(
  ui = ui,
  server = server
)