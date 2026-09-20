# Hallmark GSEA reports and gene sets

read_gsea_report <- function(path, direction) {
  report <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)

  # GSEA appends an unnamed trailing field to these reports. Extract only the
  # stable named columns before using dplyr.
  tibble(
    pathway = as.character(report[["NAME"]]),
    direction = direction,
    size = as.integer(report[["SIZE"]]),
    nes = as.numeric(report[["NES"]]),
    nominal_p = as.numeric(report[["NOM p-val"]]),
    fdr_q = as.numeric(report[["FDR q-val"]]),
    leading_edge = as.character(report[["LEADING EDGE"]])
  )
}

read_hallmark_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  fields <- strsplit(lines, "\t", fixed = TRUE)

  tibble(
    pathway = vapply(fields, `[`, character(1), 1),
    msigdb_url = vapply(fields, `[`, character(1), 2),
    genes = lapply(fields, function(x) unique(x[-c(1, 2)]))
  )
}

gsea_results <- bind_rows(
  read_gsea_report(gsea_carcinoma_report_file, "Carcinoma"),
  read_gsea_report(gsea_normal_report_file, "Normal")
) |>
  distinct(pathway, .keep_all = TRUE)

hallmark_pathways <- read_hallmark_gmt(hallmark_gmt_file) |>
  left_join(gsea_results, by = "pathway") |>
  arrange(desc(abs(nes)), pathway)

# A pathway keeps this identity colour everywhere in the default view: its
# tile, hub, membership links and member genes.  Direction colours are kept
# separate because they are an optional, different visual encoding.
# Golden-angle hues keep neighbouring tiles visually distinct even after the
# table is re-sorted.  Each of the 50 Hallmark sets therefore has its own
# stable, non-grey identity colour.
hallmark_pathway_palette <- grDevices::hcl(
  h = ((seq_len(nrow(hallmark_pathways)) - 1) * 137.508) %% 360,
  c = 78,
  l = 48
)
hallmark_pathway_colors <- setNames(
  hallmark_pathway_palette,
  hallmark_pathways$pathway
)
hallmark_direction_colors <- c(up = "#dc2626", down = "#2563eb", unknown = "#94a3b8")

pathway_identity_color <- function(pathways) {
  unname(hallmark_pathway_colors[as.character(pathways)])
}

pathway_direction_color <- function(nes) {
  ifelse(
    is.na(nes), hallmark_direction_colors[["unknown"]],
    ifelse(nes >= 0, hallmark_direction_colors[["up"]], hallmark_direction_colors[["down"]])
  )
}
