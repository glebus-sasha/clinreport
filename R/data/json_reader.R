# JSON loading layer

# ============================================================
# 6. JSON READER
# ============================================================

read_json_safe <- function(path) {
  
  if (!file.exists(path)) {
    return(NULL)
  }
  
  tryCatch(
    
    fromJSON(
      path,
      simplifyVector = FALSE
    ),
    
    error = function(e) {
      
      message(
        "JSON error: ",
        path,
        " | ",
        e$message
      )
      
      NULL
    }
  )
}


read_drug_source <- function(
    source,
    drug_id
) {
  
  x <- drug_dirs |>
    filter(
      source == !!source,
      drug_id == !!drug_id
    )
  
  if (nrow(x) == 0) {
    return(list())
  }
  
  json_files <- list.files(
    x$path[1],
    pattern = "\\.json$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(json_files) == 0) {
    return(list())
  }
  
  result <- list()
  
  for (file in json_files) {
    
    object_name <- tools::file_path_sans_ext(
      basename(file)
    )
    
    result[[object_name]] <- read_json_safe(
      file
    )
  }
  
  result
}


load_drug <- function(drug_id) {
  
  list(
    chembl = read_drug_source(
      "chembl",
      drug_id
    ),
    
    opentargets = read_drug_source(
      "opentargets",
      drug_id
    ),
    
    pubchem = read_drug_source(
      "pubchem",
      drug_id
    )
  )
}


