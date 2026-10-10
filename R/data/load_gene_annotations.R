# Gene annotations may be supplied separately or staged inside the Nextflow
# annotation directory. Retain only fields used by the gene card.
annotation_path <- gene_annotations_file
if (!nzchar(annotation_path) && nzchar(clinreport_dir) && dir.exists(clinreport_dir)) {
  candidates <- list.files(clinreport_dir, pattern = "(^|_)gene_annotations\\.json$",
    recursive = TRUE, full.names = TRUE)
  if (length(candidates) > 1L) {
    stop("Multiple gene annotation JSON files found; supply --gene-annotations-file explicitly")
  }
  if (length(candidates) == 1L) annotation_path <- candidates[[1L]]
}

gene_annotations <- list()
gene_annotation_source <- "Gene annotations"
if (nzchar(annotation_path)) {
  annotations_json <- jsonlite::fromJSON(annotation_path, simplifyVector = FALSE)
  if (!is.list(annotations_json$genes) || is.null(names(annotations_json$genes))) {
    stop("Gene annotations JSON must contain a named 'genes' object: ", annotation_path)
  }
  gene_annotations <- lapply(annotations_json$genes, function(record) {
    list(name = record$name, summary = record$summary)
  })
  if (is.character(annotations_json$metadata$source) &&
      length(annotations_json$metadata$source) == 1L &&
      nzchar(annotations_json$metadata$source)) {
    gene_annotation_source <- annotations_json$metadata$source
  }
  rm(annotations_json)
}
