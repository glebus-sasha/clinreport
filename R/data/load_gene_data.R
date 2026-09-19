# Gene expression loading and validation

# 4. LOAD GENE DATA
# ============================================================

gene_data <- read.delim(
  gene_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_gene_columns <- c(
  "gene_id",
  "gene_name",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "pvalue",
  "padj"
)

missing_gene_columns <- setdiff(
  required_gene_columns,
  names(gene_data)
)

if (length(missing_gene_columns) > 0) {
  stop(
    "Missing columns in gene file: ",
    paste(missing_gene_columns, collapse = ", ")
  )
}

gene_data <- gene_data |>
  mutate(
    baseMean = suppressWarnings(as.numeric(baseMean)),
    log2FoldChange = suppressWarnings(as.numeric(log2FoldChange)),
    lfcSE = suppressWarnings(as.numeric(lfcSE)),
    pvalue = suppressWarnings(as.numeric(pvalue)),
    padj = suppressWarnings(as.numeric(padj))
  ) |>
  mutate(
    gene_id_clean = str_remove(
      as.character(gene_id),
      "\\.[0-9]+$"
    )
  )


