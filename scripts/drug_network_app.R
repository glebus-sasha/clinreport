library(shiny)
library(tidyverse)
library(igraph)
library(visNetwork)
library(plotly)
library(networkD3)
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


# ============================================================
# Load data
# ============================================================

message("Loading gene symbols...")
gene_symbols <- load_gene_symbols(GENE_FILE)

message("Loading protein network...")
protein_edges <- load_protein_edges(PROTEIN_FILE)

message("Loading drug-target data...")
drug_targets <- load_drug_targets(DRUG_FILE)

message("Building network...")

network <- build_network(
  drug_targets = drug_targets,
  protein_edges = protein_edges,
  gene_symbols = gene_symbols
)

message(
  glue(
    "Loaded: ",
    nrow(drug_targets), " drug-target relations, ",
    nrow(protein_edges), " protein-protein edges, ",
    nrow(network$nodes), " nodes, ",
    nrow(network$edges), " edges."
  )
)

# ============================================================
# Helper functions
# ============================================================

get_drug_choices <- function() {
  drug_targets |>
    distinct(
      drug_id,
      drug_label
    ) |>
    arrange(drug_label) |>
    mutate(
      display = paste0(
        drug_label,
        " [",
        drug_id,
        "]"
      )
    )
}

drug_choices <- get_drug_choices()

# ------------------------------------------------------------
# Get selected subgraph
# ------------------------------------------------------------

get_selected_subgraph <- function(
    drug_id,
    hops,
    show_proteins,
    max_nodes
) {
  
  get_drug_network(
    network = network,
    drug_id = drug_id,
    hops = hops,
    show_proteins = show_proteins,
    max_nodes = max_nodes
  )
}

# ------------------------------------------------------------
# Node metadata
# ------------------------------------------------------------

add_node_metadata <- function(nodes) {
  
  nodes |>
    left_join(
      gene_symbols,
      by = c("id" = "gene_id")
    ) |>
    mutate(
      gene_symbol = coalesce(
        gene_symbol,
        label
      ),
      display_label = case_when(
        type == "drug" ~ label,
        TRUE ~ coalesce(
          gene_symbol,
          id
        )
      )
    )
}

# ------------------------------------------------------------
# Distance from drug
# ------------------------------------------------------------

calculate_distances <- function(
    nodes,
    edges,
    drug_id
) {
  
  if (nrow(nodes) == 0) {
    return(
      nodes |>
        mutate(distance = NA_integer_)
    )
  }
  
  graph <- graph_from_data_frame(
    d = edges,
    vertices = nodes,
    directed = FALSE
  )
  
  distances <- igraph::distances(
    graph,
    v = drug_id,
    to = V(graph)$name
  )
  
  tibble(
    id = V(graph)$name,
    distance = as.numeric(distances[1, ])
  ) |>
    left_join(
      nodes,
      by = "id"
    )
}

# ============================================================
# Network preparation
# ============================================================

prepare_visual_network <- function(
    subgraph
) {
  
  if (nrow(subgraph$nodes) == 0) {
    return(
      list(
        nodes = tibble(),
        edges = tibble()
      )
    )
  }
  
  nodes <- add_node_metadata(
    subgraph$nodes
  )
  
  nodes <- calculate_distances(
    nodes = nodes,
    edges = subgraph$edges,
    drug_id = nodes |>
      filter(type == "drug") |>
      pull(id) |>
      first()
  )
  
  nodes <- nodes |>
    mutate(
      shape = case_when(
        type == "drug" ~ "diamond",
        type == "target" ~ "dot",
        type == "protein" ~ "dot",
        TRUE ~ "dot"
      ),
      color = case_when(
        type == "drug" ~ "#D62728",
        type == "target" ~ "#1F77B4",
        type == "protein" ~ "#999999",
        TRUE ~ "#777777"
      ),
      size = case_when(
        type == "drug" ~ 30,
        type == "target" ~ 18,
        type == "protein" ~ 10,
        TRUE ~ 12
      ),
      title = case_when(
        type == "drug" ~ glue(
          "<b>{label}</b><br>",
          "DrugBank: {id}<br>",
          "Status: {status}<br>",
          "Score: {round(score, 4)}"
        ),
        TRUE ~ glue(
          "<b>{display_label}</b><br>",
          "Ensembl: {id}<br>",
          "Type: {type}"
        )
      )
    )
  
  edges <- subgraph$edges |>
    mutate(
      arrows = "to",
      color = case_when(
        type == "drug-target" ~ "#D62728",
        type == "protein-protein" ~ "#BBBBBB",
        TRUE ~ "#999999"
      ),
      width = case_when(
        type == "drug-target" ~ 3,
        TRUE ~ 1
      ),
      dashes = case_when(
        type == "protein-protein" ~ TRUE,
        TRUE ~ FALSE
      )
    )
  
  list(
    nodes = nodes,
    edges = edges
  )
}

# ============================================================
# Radial layout
# ============================================================

make_radial_layout <- function(
    subgraph,
    drug_id
) {
  
  prepared <- prepare_visual_network(
    subgraph
  )
  
  nodes <- prepared$nodes
  edges <- prepared$edges
  
  if (nrow(nodes) == 0) {
    return(
      list(
        nodes = nodes,
        edges = edges
      )
    )
  }
  
  target_nodes <- nodes |>
    filter(type == "target") |>
    arrange(label)
  
  protein_nodes <- nodes |>
    filter(type == "protein") |>
    arrange(label)
  
  drug_nodes <- nodes |>
    filter(type == "drug")
  
  coordinates <- tibble(
    id = character(),
    x = numeric(),
    y = numeric()
  )
  
  if (nrow(drug_nodes) > 0) {
    
    coordinates <- bind_rows(
      coordinates,
      tibble(
        id = drug_nodes$id,
        x = 0,
        y = 0
      )
    )
  }
  
  if (nrow(target_nodes) > 0) {
    
    angles <- seq(
      0,
      2 * pi,
      length.out = nrow(target_nodes) + 1
    )[-(nrow(target_nodes) + 1)]
    
    coordinates <- bind_rows(
      coordinates,
      tibble(
        id = target_nodes$id,
        x = cos(angles) * 400,
        y = sin(angles) * 400
      )
    )
  }
  
  if (nrow(protein_nodes) > 0) {
    
    angles <- seq(
      0,
      2 * pi,
      length.out = nrow(protein_nodes) + 1
    )[-(nrow(protein_nodes) + 1)]
    
    coordinates <- bind_rows(
      coordinates,
      tibble(
        id = protein_nodes$id,
        x = cos(angles) * 800,
        y = sin(angles) * 800
      )
    )
  }
  
  nodes <- nodes |>
    left_join(
      coordinates,
      by = "id"
    )
  
  list(
    nodes = nodes,
    edges = edges
  )
}

# ============================================================
# 3D coordinates
# ============================================================

make_3d_layout <- function(
    subgraph,
    drug_id
) {
  
  prepared <- prepare_visual_network(
    subgraph
  )
  
  nodes <- prepared$nodes
  edges <- prepared$edges
  
  if (nrow(nodes) == 0) {
    return(
      list(
        nodes = nodes,
        edges = edges
      )
    )
  }
  
  graph <- graph_from_data_frame(
    d = edges,
    vertices = nodes,
    directed = FALSE
  )
  
  set.seed(123)
  
  coords <- layout_with_fr(
    graph,
    dim = 3,
    niter = 300
  )
  
  coordinates <- tibble(
    id = V(graph)$name,
    x = coords[, 1],
    y = coords[, 2],
    z = coords[, 3]
  )
  
  nodes <- nodes |>
    left_join(
      coordinates,
      by = "id"
    )
  
  list(
    nodes = nodes,
    edges = edges
  )
}

# ============================================================
# Matrix data
# ============================================================

make_matrix_data <- function(
    drug_ids
) {
  
  selected <- drug_targets |>
    filter(
      drug_id %in% drug_ids
    ) |>
    left_join(
      gene_symbols,
      by = c("target" = "gene_id")
    ) |>
    mutate(
      gene_symbol = coalesce(
        gene_symbol,
        target
      )
    )
  
  if (nrow(selected) == 0) {
    return(
      list(
        matrix = matrix(
          numeric(),
          nrow = 0,
          ncol = 0
        ),
        drugs = character(),
        targets = character()
      )
    )
  }
  
  matrix_data <- selected |>
    group_by(
      drug_id,
      drug_label,
      gene_symbol
    ) |>
    summarise(
      score = max(
        score,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  drugs <- matrix_data |>
    distinct(
      drug_id,
      drug_label
    ) |>
    arrange(drug_label)
  
  targets <- matrix_data |>
    distinct(gene_symbol) |>
    arrange(gene_symbol) |>
    pull(gene_symbol)
  
  mat <- matrix(
    0,
    nrow = nrow(drugs),
    ncol = length(targets),
    dimnames = list(
      drugs$drug_label,
      targets
    )
  )
  
  for (i in seq_len(nrow(matrix_data))) {
    
    row_id <- matrix_data$drug_label[i]
    col_id <- matrix_data$gene_symbol[i]
    
    mat[row_id, col_id] <-
      matrix_data$score[i]
  }
  
  list(
    matrix = mat,
    drugs = drugs,
    targets = targets
  )
}

# ============================================================
# Sankey
# ============================================================

make_sankey_data <- function(
    subgraph
) {
  
  if (nrow(subgraph$nodes) == 0) {
    return(
      list(
        nodes = tibble(),
        links = tibble()
      )
    )
  }
  
  nodes <- add_node_metadata(
    subgraph$nodes
  )
  
  edges <- subgraph$edges
  
  edges <- edges |>
    left_join(
      nodes |>
        select(
          source = id,
          source_label = display_label,
          source_type = type
        ),
      by = "source"
    ) |>
    left_join(
      nodes |>
        select(
          target = id,
          target_label = display_label,
          target_type = type
        ),
      by = "target"
    )
  
  links <- edges |>
    transmute(
      source_label,
      target_label,
      value = 1,
      source_type,
      target_type
    )
  
  sankey_nodes <- tibble(
    name = unique(
      c(
        links$source_label,
        links$target_label
      )
    )
  )
  
  links <- links |>
    mutate(
      source = match(
        source_label,
        sankey_nodes$name
      ) - 1,
      target = match(
        target_label,
        sankey_nodes$name
      ) - 1
    )
  
  list(
    nodes = sankey_nodes,
    links = links
  )
}

# ============================================================
# Community detection
# ============================================================

make_communities <- function(
    subgraph
) {
  
  protein_nodes <- subgraph$nodes |>
    filter(type == "protein")
  
  protein_edges <- subgraph$edges |>
    filter(type == "protein-protein")
  
  if (
    nrow(protein_nodes) < 2 ||
    nrow(protein_edges) == 0
  ) {
    
    return(
      tibble()
    )
  }
  
  graph <- graph_from_data_frame(
    d = protein_edges,
    vertices = protein_nodes,
    directed = FALSE
  )
  
  communities <- cluster_louvain(
    graph
  )
  
  tibble(
    id = V(graph)$name,
    community = as.character(
      membership(communities)
    )
  ) |>
    left_join(
      protein_nodes,
      by = "id"
    ) |>
    arrange(
      community,
      label
    )
}

# ============================================================
# UI
# ============================================================

ui <- fluidPage(
  
  tags$head(
    
    tags$style(
      HTML(
        "
        body {
          font-family: -apple-system, BlinkMacSystemFont,
                       'Segoe UI', sans-serif;
        }

        .app-title {
          font-size: 28px;
          font-weight: 700;
          margin-bottom: 4px;
        }

        .app-subtitle {
          color: #666;
          margin-bottom: 20px;
        }

        .stat-card {
          background: #f7f7f7;
          border-radius: 10px;
          padding: 14px;
          margin-bottom: 10px;
        }

        .stat-number {
          font-size: 24px;
          font-weight: 700;
        }

        .stat-label {
          color: #666;
          font-size: 13px;
        }

        .legend-box {
          padding: 10px;
          background: #fafafa;
          border-radius: 8px;
          margin-top: 10px;
        }

        .drug-color {
          color: #D62728;
          font-weight: 700;
        }

        .target-color {
          color: #1F77B4;
          font-weight: 700;
        }

        .protein-color {
          color: #777777;
          font-weight: 700;
        }
        "
      )
    )
  ),
  
  div(
    class = "app-title",
    "Drug–Target–Protein Network Explorer"
  ),
  
  div(
    class = "app-subtitle",
    "Interactive exploration of drug targets and protein neighborhoods"
  ),
  
  sidebarLayout(
    
    sidebarPanel(
      
      width = 3,
      
      selectizeInput(
        inputId = "drug",
        label = "Drug",
        choices = NULL,
        selected = NULL
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
        label = "Show protein network",
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
      
      h4("Compare drugs"),
      
      selectizeInput(
        inputId = "compare_drugs",
        label = "Drugs for matrix comparison",
        choices = NULL,
        multiple = TRUE,
        options = list(
          maxItems = 10,
          placeholder = "Select up to 10 drugs"
        )
      ),
      
      hr(),
      
      h4("Legend"),
      
      div(
        class = "legend-box",
        
        p(
          span(
            class = "drug-color",
            "◆ Drug"
          )
        ),
        
        p(
          span(
            class = "target-color",
            "● Target"
          )
        ),
        
        p(
          span(
            class = "protein-color",
            "● Protein"
          )
        )
      ),
      
      hr(),
      
      uiOutput("drug_info"),
      
      uiOutput("network_info")
    ),
    
    mainPanel(
      
      width = 9,
      
      tabsetPanel(
        
        id = "main_tabs",
        
        # ----------------------------------------------------
        # Network
        # ----------------------------------------------------
        
        tabPanel(
          "Network",
          
          br(),
          
          visNetworkOutput(
            "network_plot",
            height = "750px"
          )
        ),
        
        # ----------------------------------------------------
        # Radial
        # ----------------------------------------------------
        
        tabPanel(
          "Radial",
          
          br(),
          
          p(
            "Drug at the center, direct targets on the first ring, ",
            "and protein neighbors on the outer ring."
          ),
          
          visNetworkOutput(
            "radial_plot",
            height = "750px"
          )
        ),
        
        # ----------------------------------------------------
        # 3D
        # ----------------------------------------------------
        
        tabPanel(
          "3D",
          
          br(),
          
          p(
            "Rotate, zoom and inspect the selected network in 3D."
          ),
          
          plotlyOutput(
            "plot_3d",
            height = "750px"
          )
        ),
        
        # ----------------------------------------------------
        # Matrix
        # ----------------------------------------------------
        
        tabPanel(
          "Matrix",
          
          br(),
          
          p(
            "Drug × target matrix. Cell intensity represents ",
            "the drug-target score."
          ),
          
          plotlyOutput(
            "matrix_plot",
            height = "700px"
          )
        ),
        
        # ----------------------------------------------------
        # Sankey
        # ----------------------------------------------------
        
        tabPanel(
          "Sankey",
          
          br(),
          
          p(
            "Flow representation of drug → target → protein ",
            "relationships."
          ),
          
          sankeyNetworkOutput(
            "sankey_plot",
            height = "750px"
          )
        ),
        
        # ----------------------------------------------------
        # Communities
        # ----------------------------------------------------
        
        tabPanel(
          "Communities",
          
          br(),
          
          p(
            "Protein communities detected with Louvain clustering."
          ),
          
          visNetworkOutput(
            "community_plot",
            height = "750px"
          )
        ),
        
        # ----------------------------------------------------
        # Compare
        # ----------------------------------------------------
        
        tabPanel(
          "Compare",
          
          br(),
          
          h4("Shared targets"),
          
          plotlyOutput(
            "compare_plot",
            height = "600px"
          )
        ),
        
        # ----------------------------------------------------
        # Statistics
        # ----------------------------------------------------
        
        tabPanel(
          "Statistics",
          
          br(),
          
          fluidRow(
            
            column(
              width = 4,
              uiOutput("stat_nodes")
            ),
            
            column(
              width = 4,
              uiOutput("stat_edges")
            ),
            
            column(
              width = 4,
              uiOutput("stat_targets")
            )
          ),
          
          br(),
          
          fluidRow(
            
            column(
              width = 6,
              
              h4("Node counts"),
              
              plotlyOutput(
                "node_type_plot",
                height = "400px"
              )
            ),
            
            column(
              width = 6,
              
              h4("Top connected nodes"),
              
              plotlyOutput(
                "degree_plot",
                height = "400px"
              )
            )
          )
        )
      )
    )
  )
)

# ============================================================
# Server
# ============================================================

server <- function(
    input,
    output,
    session
) {
  
  # ----------------------------------------------------------
  # Drug selectors
  # ----------------------------------------------------------
  
  updateSelectizeInput(
    session,
    "drug",
    choices = setNames(
      drug_choices$drug_id,
      drug_choices$display
    ),
    selected = drug_choices$drug_id[1]
  )
  
  updateSelectizeInput(
    session,
    "compare_drugs",
    choices = setNames(
      drug_choices$drug_id,
      drug_choices$display
    )
  )
  
  # ----------------------------------------------------------
  # Selected drug
  # ----------------------------------------------------------
  
  selected_drug <- reactive({
    req(input$drug)
    
    input$drug
  })
  
  # ----------------------------------------------------------
  # Selected subgraph
  # ----------------------------------------------------------
  
  selected_subgraph <- reactive({
    
    req(
      selected_drug(),
      input$hops,
      input$max_nodes
    )
    
    get_selected_subgraph(
      drug_id = selected_drug(),
      hops = input$hops,
      show_proteins = input$show_proteins,
      max_nodes = input$max_nodes
    )
  })
  
  # ----------------------------------------------------------
  # Prepared visual network
  # ----------------------------------------------------------
  
  prepared_network <- reactive({
    
    prepare_visual_network(
      selected_subgraph()
    )
  })
  
  # ==========================================================
  # Drug information
  # ==========================================================
  
  output$drug_info <- renderUI({
    
    req(selected_drug())
    
    info <- drug_targets |>
      filter(
        drug_id == selected_drug()
      ) |>
      distinct(
        drug_id,
        drug_label,
        status,
        score
      ) |>
      slice(1)
    
    if (nrow(info) == 0) {
      return(NULL)
    }
    
    div(
      
      h4(info$drug_label),
      
      p(
        strong("DrugBank: "),
        info$drug_id
      ),
      
      p(
        strong("Status: "),
        info$status
      ),
      
      p(
        strong("Score: "),
        round(info$score, 4)
      )
    )
  })
  
  # ==========================================================
  # Network information
  # ==========================================================
  
  output$network_info <- renderUI({
    
    subgraph <- selected_subgraph()
    
    if (nrow(subgraph$nodes) == 0) {
      return(
        div(
          class = "stat-card",
          "No network data."
        )
      )
    }
    
    n_drugs <- subgraph$nodes |>
      filter(type == "drug") |>
      nrow()
    
    n_targets <- subgraph$nodes |>
      filter(type == "target") |>
      nrow()
    
    n_proteins <- subgraph$nodes |>
      filter(type == "protein") |>
      nrow()
    
    n_edges <- nrow(
      subgraph$edges
    )
    
    div(
      
      h4("Current network"),
      
      p(
        strong("Drugs: "),
        n_drugs
      ),
      
      p(
        strong("Targets: "),
        n_targets
      ),
      
      p(
        strong("Proteins: "),
        n_proteins
      ),
      
      p(
        strong("Edges: "),
        n_edges
      )
    )
  })
  
  # ==========================================================
  # NETWORK
  # ==========================================================
  
  output$network_plot <- renderVisNetwork({
    
    prepared <- prepared_network()
    
    req(
      nrow(prepared$nodes) > 0
    )
    
    visNetwork(
      prepared$nodes,
      prepared$edges
    ) |>
      
      visNodes(
        borderWidth = 1
      ) |>
      
      visEdges(
        smooth = FALSE
      ) |>
      
      visOptions(
        highlightNearest = list(
          enabled = TRUE,
          degree = 1,
          hover = TRUE
        ),
        selectedBy = "type",
        collapse = FALSE
      ) |>
      
      visInteraction(
        hover = TRUE,
        navigationButtons = TRUE,
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
        stabilization = list(
          enabled = TRUE,
          iterations = 300
        )
      ) |>
      
      visLayout(
        randomSeed = 42
      )
  })
  
  # ==========================================================
  # RADIAL
  # ==========================================================
  
  output$radial_plot <- renderVisNetwork({
    
    radial <- make_radial_layout(
      subgraph = selected_subgraph(),
      drug_id = selected_drug()
    )
    
    req(
      nrow(radial$nodes) > 0
    )
    
    visNetwork(
      radial$nodes,
      radial$edges
    ) |>
      
      visNodes(
        fixed = list(
          x = TRUE,
          y = TRUE
        )
      ) |>
      
      visEdges(
        smooth = FALSE
      ) |>
      
      visOptions(
        highlightNearest = list(
          enabled = TRUE,
          degree = 1,
          hover = TRUE
        )
      ) |>
      
      visInteraction(
        hover = TRUE,
        navigationButtons = TRUE
      ) |>
      
      visPhysics(
        enabled = FALSE
      ) |>
      
      visLayout(
        randomSeed = 42
      )
  })
  
  # ==========================================================
  # 3D
  # ==========================================================
  
  output$plot_3d <- renderPlotly({
    
    layout <- make_3d_layout(
      subgraph = selected_subgraph(),
      drug_id = selected_drug()
    )
    
    nodes <- layout$nodes
    edges <- layout$edges
    
    req(
      nrow(nodes) > 0
    )
    
    edge_plot <- tibble(
      x = numeric(),
      y = numeric(),
      z = numeric()
    )
    
    for (i in seq_len(nrow(edges))) {
      
      source_node <- nodes |>
        filter(id == edges$source[i])
      
      target_node <- nodes |>
        filter(id == edges$target[i])
      
      if (
        nrow(source_node) == 0 ||
        nrow(target_node) == 0
      ) {
        next
      }
      
      edge_plot <- bind_rows(
        edge_plot,
        tibble(
          x = c(
            source_node$x,
            target_node$x,
            NA
          ),
          y = c(
            source_node$y,
            target_node$y,
            NA
          ),
          z = c(
            source_node$z,
            target_node$z,
            NA
          )
        )
      )
    }
    
    plot <- plot_ly()
    
    if (nrow(edge_plot) > 0) {
      
      plot <- plot |>
        add_trace(
          data = edge_plot,
          x = ~x,
          y = ~y,
          z = ~z,
          type = "scatter3d",
          mode = "lines",
          line = list(
            width = 2
          ),
          hoverinfo = "none",
          showlegend = FALSE
        )
    }
    
    plot |>
      add_trace(
        data = nodes,
        x = ~x,
        y = ~y,
        z = ~z,
        type = "scatter3d",
        mode = "markers",
        marker = list(
          size = ~size,
          color = ~color,
          opacity = 0.85
        ),
        text = ~title,
        hoverinfo = "text",
        showlegend = FALSE
      ) |>
      
      layout(
        scene = list(
          xaxis = list(
            title = ""
          ),
          yaxis = list(
            title = ""
          ),
          zaxis = list(
            title = ""
          )
        ),
        margin = list(
          l = 0,
          r = 0,
          b = 0,
          t = 0
        )
      )
  })
  
  # ==========================================================
  # MATRIX
  # ==========================================================
  
  output$matrix_plot <- renderPlotly({
    
    drug_ids <- input$compare_drugs
    
    if (
      is.null(drug_ids) ||
      length(drug_ids) == 0
    ) {
      drug_ids <- selected_drug()
    }
    
    matrix_data <- make_matrix_data(
      drug_ids
    )
    
    mat <- matrix_data$matrix
    
    req(
      length(mat) > 0
    )
    
    plot_ly(
      x = colnames(mat),
      y = rownames(mat),
      z = mat,
      type = "heatmap",
      hovertemplate = paste(
        "Drug: %{y}<br>",
        "Target: %{x}<br>",
        "Score: %{z:.4f}",
        "<extra></extra>"
      )
    ) |>
      
      layout(
        xaxis = list(
          title = "Target",
          tickangle = -45
        ),
        yaxis = list(
          title = "Drug"
        )
      )
  })
  
  # ==========================================================
  # SANKEY
  # ==========================================================
  
  output$sankey_plot <- renderSankeyNetwork({
    
    sankey <- make_sankey_data(
      selected_subgraph()
    )
    
    req(
      nrow(sankey$nodes) > 0,
      nrow(sankey$links) > 0
    )
    
    sankeyNetwork(
      Links = sankey$links,
      Nodes = sankey$nodes,
      Source = "source",
      Target = "target",
      Value = "value",
      NodeID = "name",
      fontSize = 13,
      nodeWidth = 25,
      sinksRight = TRUE,
      iterations = 0
    )
  })
  
  # ==========================================================
  # COMMUNITIES
  # ==========================================================
  
  output$community_plot <- renderVisNetwork({
    
    subgraph <- selected_subgraph()
    
    communities <- make_communities(
      subgraph
    )
    
    req(
      nrow(communities) > 0
    )
    
    protein_nodes <- subgraph$nodes |>
      filter(type == "protein") |>
      left_join(
        communities |>
          select(
            id,
            community
          ),
        by = "id"
      ) |>
      mutate(
        label = coalesce(
          gene_symbol,
          label
        ),
        group = paste0(
          "community_",
          community
        ),
        shape = "dot",
        size = 12,
        title = glue(
          "<b>{label}</b><br>",
          "Ensembl: {id}<br>",
          "Community: {community}"
        )
      )
    
    protein_edges <- subgraph$edges |>
      filter(
        type == "protein-protein",
        source %in% protein_nodes$id,
        target %in% protein_nodes$id
      ) |>
      mutate(
        color = "#BBBBBB",
        width = 1
      )
    
    visNetwork(
      protein_nodes,
      protein_edges
    ) |>
      
      visNodes(
        borderWidth = 1
      ) |>
      
      visEdges(
        smooth = FALSE
      ) |>
      
      visOptions(
        highlightNearest = list(
          enabled = TRUE,
          degree = 1,
          hover = TRUE
        ),
        groups = TRUE
      ) |>
      
      visInteraction(
        hover = TRUE,
        navigationButtons = TRUE
      ) |>
      
      visPhysics(
        solver = "forceAtlas2Based",
        forceAtlas2Based = list(
          gravitationalConstant = -40,
          springLength = 80,
          springConstant = 0.05,
          damping = 0.4
        ),
        stabilization = list(
          enabled = TRUE,
          iterations = 300
        )
      )
  })
  
  # ==========================================================
  # COMPARE
  # ==========================================================
  
  output$compare_plot <- renderPlotly({
    
    drug_ids <- input$compare_drugs
    
    if (
      is.null(drug_ids) ||
      length(drug_ids) < 2
    ) {
      
      return(
        plotly_empty(
          type = "scatter",
          mode = "text"
        ) |>
          layout(
            annotations = list(
              list(
                text = "Select at least two drugs.",
                x = 0.5,
                y = 0.5,
                showarrow = FALSE
              )
            )
          )
      )
    }
    
    selected <- drug_targets |>
      filter(
        drug_id %in% drug_ids
      ) |>
      left_join(
        gene_symbols,
        by = c("target" = "gene_id")
      ) |>
      mutate(
        target_label = coalesce(
          gene_symbol,
          target
        )
      )
    
    target_counts <- selected |>
      distinct(
        drug_id,
        target_label
      ) |>
      count(
        target_label,
        name = "n_drugs"
      ) |>
      arrange(
        desc(n_drugs),
        target_label
      )
    
    selected <- selected |>
      left_join(
        target_counts,
        by = "target_label"
      )
    
    plot_ly(
      selected,
      x = ~drug_id,
      y = ~target_label,
      z = ~score,
      type = "heatmap",
      hovertemplate = paste(
        "Drug: %{x}<br>",
        "Target: %{y}<br>",
        "Score: %{z:.4f}",
        "<extra></extra>"
      )
    ) |>
      
      layout(
        xaxis = list(
          title = "Drug"
        ),
        yaxis = list(
          title = "Target"
        )
      )
  })
  
  # ==========================================================
  # STATISTICS
  # ==========================================================
  
  output$stat_nodes <- renderUI({
    
    n <- nrow(
      selected_subgraph()$nodes
    )
    
    div(
      class = "stat-card",
      div(
        class = "stat-number",
        n
      ),
      div(
        class = "stat-label",
        "Nodes"
      )
    )
  })
  
  output$stat_edges <- renderUI({
    
    n <- nrow(
      selected_subgraph()$edges
    )
    
    div(
      class = "stat-card",
      div(
        class = "stat-number",
        n
      ),
      div(
        class = "stat-label",
        "Edges"
      )
    )
  })
  
  output$stat_targets <- renderUI({
    
    n <- selected_subgraph()$nodes |>
      filter(type == "target") |>
      nrow()
    
    div(
      class = "stat-card",
      div(
        class = "stat-number",
        n
      ),
      div(
        class = "stat-label",
        "Targets"
      )
    )
  })
  
  # ==========================================================
  # NODE TYPE PLOT
  # ==========================================================
  
  output$node_type_plot <- renderPlotly({
    
    data <- selected_subgraph()$nodes |>
      count(
        type
      )
    
    plot_ly(
      data,
      x = ~type,
      y = ~n,
      type = "bar",
      text = ~n,
      textposition = "auto"
    ) |>
      
      layout(
        xaxis = list(
          title = ""
        ),
        yaxis = list(
          title = "Number of nodes"
        )
      )
  })
  
  # ==========================================================
  # DEGREE PLOT
  # ==========================================================
  
  output$degree_plot <- renderPlotly({
    
    subgraph <- selected_subgraph()
    
    req(
      nrow(subgraph$nodes) > 0
    )
    
    graph <- graph_from_data_frame(
      d = subgraph$edges,
      vertices = subgraph$nodes,
      directed = FALSE
    )
    
    degree_data <- tibble(
      id = V(graph)$name,
      degree = degree(graph)
    ) |>
      left_join(
        subgraph$nodes |>
          select(
            id,
            label,
            type
          ),
        by = "id"
      ) |>
      arrange(
        desc(degree)
      ) |>
      slice_head(
        n = 20
      ) |>
      mutate(
        label = coalesce(
          label,
          id
        )
      )
    
    plot_ly(
      degree_data,
      x = ~degree,
      y = ~reorder(label, degree),
      type = "bar",
      orientation = "h",
      text = ~degree,
      textposition = "auto"
    ) |>
      
      layout(
        xaxis = list(
          title = "Degree"
        ),
        yaxis = list(
          title = ""
        )
      )
  })
}

# ============================================================
# Run application
# ============================================================

shinyApp(
  ui = ui,
  server = server
)