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

# Accept both the development layout (clinreport/<source>/<drug>/) and the
# release layout where source directories may be nested one level deeper.
# Nextflow stages the annotations directory as a single path, preserving this
# hierarchy instead of flattening it with a wildcard stageAs pattern.
drug_dirs <- map_dfr(sources, find_drug_dirs)
if (nrow(drug_dirs) == 0 && nzchar(clinreport_dir) && dir.exists(clinreport_dir)) {
  candidates <- list.dirs(clinreport_dir, recursive = TRUE, full.names = TRUE)
  # Release exports may name directories as *_DB00142_chembl (source suffix)
  # instead of nesting them under a chembl/ directory.
  candidates <- candidates[grepl("DB[0-9]+(_(chembl|opentargets|opentargits|pubchem))?$", basename(candidates))]
  if (length(candidates)) {
    drug_dirs <- tibble(
      path = candidates,
      folder = basename(candidates),
      drug_id = str_extract(basename(candidates), "DB[0-9]+"),
      source = vapply(candidates, function(p) {
        bits <- strsplit(normalizePath(p, winslash = "/"), "/", fixed = TRUE)[[1]]
        hit <- bits[bits %in% sources]
        if (length(hit)) hit[[length(hit)]] else {
          suffix <- str_extract(basename(p), "(chembl|opentargets|opentargits|pubchem)$")
          if (identical(suffix, "opentargits")) suffix <- "opentargets"
          ifelse(is.na(suffix), "unknown", suffix)
        }
      }, character(1))
    ) |>
      select(source, path, folder, drug_id)
  }
}
