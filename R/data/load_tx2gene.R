# Transcript-to-gene mapping

tx2gene <- read.delim(
  tx2gene_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_tx2gene_columns <- c("transcript_id", "gene_id", "gene_name")
missing_tx2gene_columns <- setdiff(required_tx2gene_columns, names(tx2gene))

if (length(missing_tx2gene_columns) > 0) {
  stop(
    "Missing columns in tx2gene file: ",
    paste(missing_tx2gene_columns, collapse = ", ")
  )
}

tx2gene <- tx2gene |>
  transmute(
    transcript_id = as.character(transcript_id),
    gene_id = str_remove(as.character(gene_id), "\\.[0-9]+$"),
    gene_name = as.character(gene_name)
  ) |>
  filter(
    !is.na(gene_id),
    gene_id != ""
  ) |>
  distinct(gene_id, gene_name)

# One display name per gene. gene_data remains the preferred source when available.
tx2gene_names <- tx2gene |>
  filter(!is.na(gene_name), gene_name != "") |>
  group_by(gene_id) |>
  summarise(
    gene_name = first(gene_name),
    .groups = "drop"
  )
