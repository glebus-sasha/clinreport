# External source buttons

# ============================================================
# 12. EXTERNAL SOURCE BUTTONS
# ============================================================

source_buttons <- function(
    drug_id,
    chembl_id,
    pubchem_cid
) {
  
  buttons <- list()
  
  if (!is.null(drug_id)) {
    
    buttons <- append(
      buttons,
      list(
        tags$a(
          href = drugbank_url(drug_id),
          target = "_blank",
          rel = "noopener noreferrer",
          class = "source-button",
          "DrugBank"
        )
      )
    )
  }
  
  if (!is.null(chembl_id)) {
    
    buttons <- append(
      buttons,
      list(
        tags$a(
          href = chembl_url(chembl_id),
          target = "_blank",
          rel = "noopener noreferrer",
          class = "source-button",
          "ChEMBL"
        )
      )
    )
  }
  
  if (!is.null(pubchem_cid)) {
    
    buttons <- append(
      buttons,
      list(
        tags$a(
          href = pubchem_url(pubchem_cid),
          target = "_blank",
          rel = "noopener noreferrer",
          class = "source-button",
          "PubChem"
        )
      )
    )
  }
  
  if (length(buttons) == 0) {
    return(NULL)
  }
  
  div(
    class = "source-buttons",
    buttons
  )
}


