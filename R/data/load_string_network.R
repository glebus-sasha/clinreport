# STRING network loading

# The file contains an edge list with Ensembl gene identifiers.
interaction_edges <- read.delim(
  interaction_network_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = c("character", "character")
)

required_string_columns <- c("source", "target")
missing_string_columns <- setdiff(required_string_columns, names(interaction_edges))

if (length(missing_string_columns) > 0) {
  stop(
    "Missing columns in interaction network edge file: ",
    paste(missing_string_columns, collapse = ", ")
  )
}

interaction_edges <- interaction_edges |>
  transmute(
    source = str_remove(as.character(source), "\\.[0-9]+$"),
    target = str_remove(as.character(target), "\\.[0-9]+$")
  ) |>
  filter(
    !is.na(source),
    !is.na(target),
    source != "",
    target != "",
    source != target
  ) |>
  distinct()
