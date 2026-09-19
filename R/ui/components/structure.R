# 2D structure component

# 14. 2D STRUCTURE
# ============================================================

render_structure <- function(
    pubchem_cid
) {
  
  if (
    is.null(pubchem_cid) ||
    !nzchar(pubchem_cid)
  ) {
    
    return(
      div(
        class = "structure-empty",
        "2D structure is not available from PubChem for this record."
      )
    )
  }
  
  image_url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
    pubchem_cid,
    "/PNG"
  )
  
  div(
    class = "structure-container",
    
    tags$img(
      src = image_url,
      class = "chemical-structure",
      alt = "2D chemical structure"
    ),
    
    div(
      class = "structure-caption",
      
      "PubChem CID ",
      
      tags$a(
        href = pubchem_url(pubchem_cid),
        target = "_blank",
        rel = "noopener noreferrer",
        class = "external-link",
        pubchem_cid
      )
    )
  )
}


