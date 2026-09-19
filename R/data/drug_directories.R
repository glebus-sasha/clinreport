# Discovery of per-source drug directories

# 5. FIND JSON DIRECTORIES
# ============================================================

sources <- c(
  "chembl",
  "opentargets",
  "pubchem"
)

find_drug_dirs <- function(source) {
  
  source_dir <- file.path(
    clinreport_dir,
    source
  )
  
  if (!dir.exists(source_dir)) {
    return(tibble())
  }
  
  dirs <- list.dirs(
    source_dir,
    recursive = FALSE,
    full.names = TRUE
  )
  
  if (length(dirs) == 0) {
    return(tibble())
  }
  
  tibble(
    source = source,
    path = dirs,
    folder = basename(dirs),
    drug_id = str_extract(
      basename(dirs),
      "DB[0-9]+$"
    )
  ) |>
    filter(!is.na(drug_id))
}

drug_dirs <- map_dfr(
  sources,
  find_drug_dirs
)


