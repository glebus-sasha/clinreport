# PubChem description records may include an upstream attribution (e.g. ChEBI).
get_drug_description <- function(data) {
  records <- get_nested(data$pubchem$pubchem_description, "InformationList", "Information")
  if (!is.list(records)) return(NULL)
  for (record in records) {
    if (!is.list(record)) next
    description <- safe_text(record$Description, "")
    if (!nzchar(description)) next
    source_url <- safe_text(record$DescriptionURL, "")
    if (!grepl("^https://", source_url, ignore.case = TRUE)) source_url <- ""
    return(list(text = description,
      source = safe_text(record$DescriptionSourceName, "PubChem"),
      url = source_url))
  }
  NULL
}
