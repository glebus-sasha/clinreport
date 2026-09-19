# Drug name resolution

# ============================================================
# 9. DRUG NAME
# ============================================================

get_drug_name_from_json <- function(
    drug_id
) {
  
  chembl <- read_drug_source(
    "chembl",
    drug_id
  )
  
  if ("chembl_id" %in% names(chembl)) {
    
    x <- chembl$chembl_id
    
    if (!is.null(x$drug_name)) {
      
      return(
        as.character(
          x$drug_name[[1]]
        )
      )
    }
  }
  
  if ("chembl_molecule" %in% names(chembl)) {
    
    x <- chembl$chembl_molecule
    
    if (!is.null(x$pref_name)) {
      
      return(
        as.character(
          x$pref_name[[1]]
        )
      )
    }
  }
  
  drug_id
}


