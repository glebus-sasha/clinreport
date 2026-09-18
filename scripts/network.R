library(tidyverse)
library(igraph)


# ---------------------------------------------------------
# Drug nodes
# ---------------------------------------------------------

build_drug_nodes <- function(drug_targets) {
  
  drug_targets |>
    distinct(
      drug_id,
      drug_label,
      status,
      score
    ) |>
    transmute(
      id = drug_id,
      label = drug_label,
      type = "drug",
      status,
      score
    )
}


# ---------------------------------------------------------
# Direct target nodes
# ---------------------------------------------------------

build_target_nodes <- function(
    drug_targets,
    gene_symbols
) {
  
  drug_targets |>
    distinct(target) |>
    rename(id = target) |>
    left_join(
      gene_symbols,
      by = c("id" = "gene_id")
    ) |>
    transmute(
      id,
      label = coalesce(gene_symbol, id),
      type = "target",
      status = NA_character_,
      score = NA_real_
    )
}


# ---------------------------------------------------------
# Protein nodes
# ---------------------------------------------------------

build_protein_nodes <- function(
    protein_edges,
    gene_symbols
) {
  
  bind_rows(
    protein_edges |>
      select(id = source),
    
    protein_edges |>
      select(id = target)
  ) |>
    distinct() |>
    left_join(
      gene_symbols,
      by = c("id" = "gene_id")
    ) |>
    transmute(
      id,
      label = coalesce(gene_symbol, id),
      type = "protein",
      status = NA_character_,
      score = NA_real_
    )
}


# ---------------------------------------------------------
# Drug-target edges
# ---------------------------------------------------------

build_drug_target_edges <- function(drug_targets) {
  
  drug_targets |>
    transmute(
      source = drug_id,
      target,
      type = "drug-target"
    ) |>
    distinct()
}


# ---------------------------------------------------------
# Protein-protein edges
# ---------------------------------------------------------

build_protein_edges <- function(protein_edges) {
  
  protein_edges |>
    transmute(
      source,
      target,
      type = "protein-protein"
    ) |>
    distinct()
}


# ---------------------------------------------------------
# Complete network
# ---------------------------------------------------------

build_network <- function(
    drug_targets,
    protein_edges,
    gene_symbols
) {
  
  drug_nodes <- build_drug_nodes(
    drug_targets
  )
  
  target_nodes <- build_target_nodes(
    drug_targets,
    gene_symbols
  )
  
  protein_nodes <- build_protein_nodes(
    protein_edges,
    gene_symbols
  )
  
  nodes <- bind_rows(
    drug_nodes,
    target_nodes,
    protein_nodes
  ) |>
    distinct(id, .keep_all = TRUE)
  
  edges <- bind_rows(
    build_drug_target_edges(
      drug_targets
    ),
    
    build_protein_edges(
      protein_edges
    )
  ) |>
    distinct(
      source,
      target,
      .keep_all = TRUE
    )
  
  list(
    nodes = nodes,
    edges = edges
  )
}


# ---------------------------------------------------------
# Neighborhood
# ---------------------------------------------------------

get_neighborhood <- function(
    network,
    node_id,
    hops = 1
) {
  
  graph <- igraph::graph_from_data_frame(
    d = network$edges,
    vertices = network$nodes,
    directed = FALSE
  )
  
  if (!node_id %in% igraph::V(graph)$name) {
    
    return(
      list(
        nodes = tibble(),
        edges = tibble()
      )
    )
  }
  
  selected <- igraph::ego(
    graph,
    order = hops,
    nodes = node_id,
    mode = "all"
  )[[1]]
  
  selected_ids <- igraph::V(graph)$name[selected]
  
  nodes <- network$nodes |>
    filter(id %in% selected_ids)
  
  edges <- network$edges |>
    filter(
      source %in% selected_ids,
      target %in% selected_ids
    )
  
  list(
    nodes = nodes,
    edges = edges
  )
}


# ---------------------------------------------------------
# Drug neighborhood
# ---------------------------------------------------------

get_drug_network <- function(
    network,
    drug_id,
    hops = 1,
    show_proteins = FALSE,
    max_nodes = 500
) {
  
  subgraph <- get_neighborhood(
    network = network,
    node_id = drug_id,
    hops = hops
  )
  
  if (nrow(subgraph$nodes) == 0) {
    return(subgraph)
  }
  
  if (!show_proteins) {
    
    allowed_ids <- subgraph$nodes |>
      filter(
        type %in% c(
          "drug",
          "target"
        )
      ) |>
      pull(id)
    
    subgraph$nodes <- subgraph$nodes |>
      filter(id %in% allowed_ids)
    
    subgraph$edges <- subgraph$edges |>
      filter(
        source %in% allowed_ids,
        target %in% allowed_ids
      )
  }
  
  # -------------------------------------------------------
  # Limit graph size
  # -------------------------------------------------------
  
  if (nrow(subgraph$nodes) > max_nodes) {
    
    drug_ids <- subgraph$nodes |>
      filter(type == "drug") |>
      pull(id)
    
    target_ids <- subgraph$nodes |>
      filter(type == "target") |>
      pull(id)
    
    remaining <- max(
      0,
      max_nodes - length(drug_ids)
    )
    
    target_ids <- target_ids |>
      head(remaining)
    
    keep_ids <- c(
      drug_ids,
      target_ids
    )
    
    subgraph$nodes <- subgraph$nodes |>
      filter(id %in% keep_ids)
    
    subgraph$edges <- subgraph$edges |>
      filter(
        source %in% keep_ids,
        target %in% keep_ids
      )
  }
  
  subgraph
}