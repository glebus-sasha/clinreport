# PCGR-annotated whole-exome variants

extract_vcf_info <- function(info, key) {
  prefix <- paste0(key, "=")
  vapply(
    strsplit(info, ";", fixed = TRUE),
    function(fields) {
      match <- fields[startsWith(fields, prefix)]
      if (length(match) == 0) {
        return(NA_character_)
      }
      sub(paste0("^", key, "="), "", match[1])
    },
    character(1)
  )
}

first_alt_value <- function(x) {
  vapply(
    strsplit(as.character(x), ",", fixed = TRUE),
    function(values) values[1],
    character(1)
  )
}

format_mutation_classes <- function(classes) {
  counts <- sort(table(classes), decreasing = TRUE)
  paste0(names(counts), " (", as.integer(counts), ")", collapse = ", ")
}

wes_connection <- gzfile(wes_vcf_file, "rt")
wes_vcf <- read.delim(
  wes_connection,
  header = FALSE,
  comment.char = "#",
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE,
  fill = TRUE
)
close(wes_connection)

if (ncol(wes_vcf) < 8) {
  stop("WES VCF must contain the eight mandatory VCF columns.")
}

names(wes_vcf)[seq_len(8)] <- c(
  "chrom", "position", "variant_id", "ref", "alt", "quality", "filter", "info"
)

wes_variants <- wes_vcf |>
  transmute(
    chrom = as.character(chrom),
    position = as.integer(position),
    ref = as.character(ref),
    alt = as.character(alt),
    filter = as.character(filter),
    gene_symbol = na_if(extract_vcf_info(info, "SYMBOL"), ""),
    gene_id = str_remove(extract_vcf_info(info, "ENSEMBL_GENE_ID"), "\\.[0-9]+$"),
    consequence = na_if(extract_vcf_info(info, "Consequence"), ""),
    impact = na_if(extract_vcf_info(info, "IMPACT"), ""),
    protein_change = coalesce(
      na_if(extract_vcf_info(info, "HGVSp_short"), ""),
      na_if(extract_vcf_info(info, "ALTERATION"), "")
    ),
    oncogenicity = na_if(extract_vcf_info(info, "ONCOGENICITY"), ""),
    clinvar_classification = na_if(extract_vcf_info(info, "CLINVAR_CLASSIFICATION"), ""),
    clinvar_allele_id = na_if(extract_vcf_info(info, "CLINVAR_ALLELE_ID"), ""),
    dbsnp_rsid = na_if(extract_vcf_info(info, "DBSNP_RSID"), ""),
    vaf = suppressWarnings(as.numeric(first_alt_value(extract_vcf_info(info, "TVAF")))),
    depth = suppressWarnings(as.integer(extract_vcf_info(info, "TDP")))
  ) |>
  mutate(
    classification = coalesce(oncogenicity, clinvar_classification, "Unclassified"),
    classification = str_replace_all(classification, "_", " "),
    variant_label = paste0(chrom, ":", position, " ", ref, ">", alt)
  ) |>
  filter(!is.na(gene_symbol), gene_symbol != "")

wes_gene_summary <- wes_variants |>
  group_by(gene_symbol) |>
  summarise(
    mutation_count = n(),
    mutation_classes = format_mutation_classes(classification),
    .groups = "drop"
  )
