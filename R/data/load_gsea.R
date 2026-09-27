# GSEA reports and gene sets, grouped by the supplied GMT collection list.

read_gsea_report <- function(path, direction, collection) {
  report <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)

  # GSEA appends an unnamed trailing field to these reports. Extract only the
  # stable named columns before using dplyr.
  tibble(
    pathway = as.character(report[["NAME"]]),
    collection = collection,
    direction = direction,
    size = as.integer(report[["SIZE"]]),
    nes = as.numeric(report[["NES"]]),
    nominal_p = as.numeric(report[["NOM p-val"]]),
    fdr_q = as.numeric(report[["FDR q-val"]]),
    leading_edge = as.character(report[["LEADING EDGE"]])
  )
}

read_pathway_gmt <- function(path, collection) {
  lines <- readLines(path, warn = FALSE)
  fields <- strsplit(lines, "\t", fixed = TRUE)

  tibble(
    pathway = vapply(fields, `[`, character(1), 1),
    collection = collection,
    msigdb_url = vapply(fields, `[`, character(1), 2),
    genes = lapply(fields, function(x) unique(x[-c(1, 2)]))
  )
}

gsea_files <- if (exists("gsea_files") && length(gsea_files)) gsea_files else
  c(gsea_carcinoma_report_file, gsea_normal_report_file)
gmt_files <- if (exists("gmt_files") && length(gmt_files)) gmt_files else pathway_gmt_file
collection_list <- if (exists("gmt_list_file") && nzchar(gmt_list_file)) {
  read.csv(gmt_list_file, stringsAsFactors = FALSE)
} else {
  data.frame(name = pathway_collection_name, url = basename(gmt_files))
}
if (!all(c("name", "url") %in% names(collection_list)) || !nrow(collection_list)) {
  stop("GMT collection list must contain name and url columns")
}
collection_list$file <- basename(sub("[?#].*$", "", collection_list$url))
if (anyNA(collection_list$name) || any(!nzchar(collection_list$name)) ||
    anyDuplicated(collection_list$name) || anyDuplicated(collection_list$file) ||
    !setequal(collection_list$file, basename(gmt_files))) {
  stop("GMT files do not match the collection list")
}
pathway_collections <- collection_list$name
display_collection_name <- function(x) gsub("_", " ", x, fixed = TRUE)
gmt_paths <- setNames(gmt_files, basename(gmt_files))
report_names <- basename(gsea_files)
used_reports <- character()
pathway_sets <- bind_rows(lapply(seq_len(nrow(collection_list)), function(i) {
  gmt_path <- gmt_paths[[collection_list$file[i]]]
  stem <- tools::file_path_sans_ext(collection_list$file[i])
  matches <- grepl(paste0(stem, ".gsea_report_for_"), report_names, fixed = TRUE)
  reports <- gsea_files[matches]
  if (!length(reports)) stop("No GSEA reports for collection: ", collection_list$name[i])
  used_reports <<- c(used_reports, reports)
  results <- bind_rows(lapply(reports, function(path) {
    report_name <- basename(path)
    group <- sub("^.*\\.gsea_report_for_", "", report_name)
    group <- sub("\\.tsv$", "", group, ignore.case = TRUE)
    if (!nzchar(group)) stop("GSEA report has no group name: ", report_name)
    direction <- tools::toTitleCase(gsub("[_-]", " ", group))
    read_gsea_report(path, direction, collection_list$name[i])
  })) |>
    arrange(fdr_q, desc(abs(nes))) |>
    distinct(pathway, .keep_all = TRUE)
  gene_sets <- read_pathway_gmt(gmt_path, collection_list$name[i])
  missing_sets <- setdiff(results$pathway, gene_sets$pathway)
  if (length(missing_sets)) {
    stop("GSEA pathways absent from ", basename(gmt_path), ": ",
      paste(head(missing_sets, 5), collapse = ", "))
  }
  gene_sets |>
    inner_join(results, by = c("pathway", "collection"))
})) |>
  arrange(collection, desc(abs(nes)), pathway)
if (!setequal(used_reports, gsea_files)) stop("GSEA reports without a matching GMT collection: ",
  paste(basename(setdiff(gsea_files, used_reports)), collapse = ", "))

pathways_for_collection <- function(collection) {
  if (length(collection) != 1L || is.na(collection) || !collection %in% pathway_collections) {
    stop("Unknown GSEA collection: ", collection)
  }
  pathway_sets[pathway_sets$collection == collection, , drop = FALSE]
}

# A pathway keeps this identity colour everywhere in the default view: its
# tile, hub, membership links and member genes.  Direction colours are kept
# separate because they are an optional, different visual encoding.
# Golden-angle hues keep neighbouring tiles visually distinct after sorting.
pathway_set_palette <- grDevices::hcl(
  h = ((seq_len(nrow(pathway_sets)) - 1) * 137.508) %% 360,
  c = 78,
  l = 48
)
pathway_set_colors <- setNames(
  pathway_set_palette,
  pathway_sets$pathway
)
pathway_direction_colors <- c(up = "#dc2626", down = "#2563eb", unknown = "#94a3b8")

pathway_identity_color <- function(pathways) {
  unname(pathway_set_colors[as.character(pathways)])
}

pathway_direction_color <- function(nes) {
  ifelse(
    is.na(nes), pathway_direction_colors[["unknown"]],
    ifelse(nes >= 0, pathway_direction_colors[["up"]], pathway_direction_colors[["down"]])
  )
}
