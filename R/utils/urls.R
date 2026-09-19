# External URL builders

# 10. EXTERNAL LINKS
# ============================================================

pubchem_url <- function(cid) {
  
  if (
    is.null(cid) ||
    !nzchar(cid)
  ) {
    return(NULL)
  }
  
  paste0(
    "https://pubchem.ncbi.nlm.nih.gov/compound/",
    cid
  )
}


chembl_url <- function(chembl_id) {
  
  if (
    is.null(chembl_id) ||
    !nzchar(chembl_id)
  ) {
    return(NULL)
  }
  
  paste0(
    "https://www.ebi.ac.uk/chembl/explore/compound/",
    chembl_id
  )
}


drugbank_url <- function(drug_id) {
  
  if (
    is.null(drug_id) ||
    !nzchar(drug_id)
  ) {
    return(NULL)
  }
  
  paste0(
    "https://go.drugbank.com/drugs/",
    drug_id
  )
}


ensembl_url <- function(gene_id) {
  
  if (
    is.null(gene_id) ||
    !nzchar(gene_id)
  ) {
    return(NULL)
  }
  
  paste0(
    "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=",
    gene_id
  )
}


