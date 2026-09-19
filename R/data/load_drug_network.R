# Drug network loading and validation

# ============================================================
# 3. LOAD DRUG NETWORK
# ============================================================

drug_network <- read.csv(
  drug_network_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (
  ncol(drug_network) > 0 &&
  names(drug_network)[1] == ""
) {
  drug_network <- drug_network[, -1, drop = FALSE]
}

required_drug_columns <- c(
  "drugId",
  "label",
  "status",
  "drugstoneType",
  "score",
  "hasEdgesTo",
  "isResult",
  "isConnector"
)

missing_drug_columns <- setdiff(
  required_drug_columns,
  names(drug_network)
)

if (length(missing_drug_columns) > 0) {
  stop(
    "Missing columns in drug network CSV: ",
    paste(missing_drug_columns, collapse = ", ")
  )
}


