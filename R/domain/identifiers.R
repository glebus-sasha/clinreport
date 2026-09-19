# External identifiers

# ============================================================
# 8. IDENTIFIERS
# ============================================================

get_chembl_id <- function(data) {
  
  candidates <- list(
    get_nested(
      data$chembl$chembl_id,
      "chembl_id"
    ),
    
    get_nested(
      data$chembl$chembl_molecule,
      "chembl_id"
    )
  )
  
  candidates <- map_chr(
    candidates,
    ~ safe_text(.x, "")
  )
  
  candidates[
    nzchar(candidates)
  ][1] %||% NULL
}


get_pubchem_cid <- function(data) {
  
  x <- get_nested(
    data$pubchem$pubchem_cids,
    "IdentifierList",
    "CID"
  )
  
  value <- safe_text(
    x,
    ""
  )
  
  if (!nzchar(value)) {
    return(NULL)
  }
  
  str_split(
    value,
    ";"
  )[[1]][1]
}


`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x[1])) y else x
}


