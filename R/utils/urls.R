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


ncbi_gene_search_url <- function(gene_symbol) {
  if (is.null(gene_symbol) || !nzchar(gene_symbol)) {
    return(NULL)
  }

  paste0(
    "https://www.ncbi.nlm.nih.gov/gene/?term=",
    utils::URLencode(paste0(gene_symbol, "[sym]"), reserved = TRUE)
  )
}


clinvar_allele_url <- function(allele_id) {
  if (is.null(allele_id) || is.na(allele_id) || !nzchar(allele_id)) {
    return(NULL)
  }

  paste0(
    "https://www.ncbi.nlm.nih.gov/clinvar/?term=",
    utils::URLencode(paste0(allele_id, "[alleleid]"), reserved = TRUE)
  )
}


dbsnp_url <- function(rsid) {
  if (is.null(rsid) || is.na(rsid) || !nzchar(rsid)) {
    return(NULL)
  }

  paste0("https://www.ncbi.nlm.nih.gov/snp/", rsid)
}

