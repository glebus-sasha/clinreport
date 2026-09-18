library(tidyverse)

# ---------------------------------------------------------
# Helpers
# ---------------------------------------------------------

clean_ensembl_id <- function(x) {
  x |>
    as.character() |>
    stringr::str_remove("\\..*$") |>
    stringr::str_trim()
}


# ---------------------------------------------------------
# Gene symbols
# ---------------------------------------------------------

load_gene_symbols <- function(path) {
  
  readr::read_tsv(
    file = path,
    show_col_types = FALSE,
    progress = FALSE
  ) |>
    transmute(
      gene_id = clean_ensembl_id(gene_id),
      gene_symbol = gene_name
    ) |>
    filter(
      !is.na(gene_id),
      gene_id != "",
      !is.na(gene_symbol),
      gene_symbol != ""
    ) |>
    distinct(gene_id, .keep_all = TRUE)
}


# ---------------------------------------------------------
# Protein-protein edges
# ---------------------------------------------------------

load_protein_edges <- function(path) {
  
  readr::read_tsv(
    file = path,
    col_names = c("source", "target"),
    show_col_types = FALSE,
    progress = FALSE
  ) |>
    transmute(
      source = clean_ensembl_id(source),
      target = clean_ensembl_id(target)
    ) |>
    filter(
      !is.na(source),
      !is.na(target),
      source != "",
      target != "",
      source != target
    ) |>
    distinct()
}


# ---------------------------------------------------------
# Parse hasEdgesTo
# ---------------------------------------------------------

parse_target_list <- function(x) {
  
  if (is.na(x) || x == "") {
    return(character())
  }
  
  x |>
    stringr::str_remove("^\\[") |>
    stringr::str_remove("\\]$") |>
    stringr::str_remove_all("'") |>
    stringr::str_remove_all("\"") |>
    stringr::str_split(",\\s*") |>
    purrr::pluck(1) |>
    stringr::str_trim()
}


# ---------------------------------------------------------
# Drug-target network
# ---------------------------------------------------------

load_drug_targets <- function(path) {
  
  raw <- readr::read_csv(
    file = path,
    show_col_types = FALSE,
    progress = FALSE,
    quote = "\""
  )
  
  raw |>
    transmute(
      drug_id = as.character(drugId),
      drug_label = as.character(label),
      status = as.character(status),
      drugstone_type = as.character(drugstoneType),
      score = as.numeric(score),
      targets = purrr::map(
        hasEdgesTo,
        parse_target_list
      )
    ) |>
    tidyr::unnest_longer(
      targets,
      values_to = "target"
    ) |>
    transmute(
      drug_id,
      drug_label,
      status,
      drugstone_type,
      score,
      target = clean_ensembl_id(target)
    ) |>
    filter(
      !is.na(drug_id),
      drug_id != "",
      !is.na(target),
      target != ""
    ) |>
    distinct()
}