# Drug-to-gene business logic

# 11. GENES FOR DRUG
# ============================================================

extract_gene_ids <- function(x) {
  
  if (is.null(x)) {
    return(character())
  }
  
  x <- as.character(x)
  
  if (
    length(x) == 0 ||
    all(is.na(x))
  ) {
    return(character())
  }
  
  ids <- str_extract_all(
    x,
    "ENSG[0-9]+"
  )[[1]]
  
  unique(ids)
}


get_drug_genes <- function(
    drug_id
) {
  
  row <- drug_network |>
    filter(
      drugId == !!drug_id
    ) |>
    slice(1)
  
  if (nrow(row) == 0) {
    return(tibble())
  }
  
  gene_ids <- extract_gene_ids(
    row$hasEdgesTo
  )
  
  if (length(gene_ids) == 0) {
    return(tibble())
  }
  
  gene_data |>
    filter(
      gene_id_clean %in% gene_ids
    ) |>
    select(
      gene_id,
      gene_id_clean,
      gene_name,
      baseMean,
      log2FoldChange,
      lfcSE,
      pvalue,
      padj
    ) |>
    distinct(
      gene_id_clean,
      .keep_all = TRUE
    ) |>
    mutate(
      
      significant_padj =
        !is.na(padj) &
        padj < padj_cutoff,
      
      significant_fc =
        !is.na(log2FoldChange) &
        abs(log2FoldChange) >= log2fc_cutoff,
      
      significant_both =
        significant_padj &
        significant_fc,
      
      direction =
        case_when(
          is.na(log2FoldChange) ~ "Neutral",
          log2FoldChange > 0 ~ "Up",
          log2FoldChange < 0 ~ "Down",
          TRUE ~ "Neutral"
        )
    ) |>
    arrange(
      is.na(padj),
      padj
    )
}


