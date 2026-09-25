# WES variant table shell

render_wes_variants <- function() {
  if (!has_wes) {
    return(div(class = "empty-message", "WES variants were not supplied for this report."))
  }
  if (nrow(wes_variants) == 0) {
    return(div(class = "empty-message", "No gene-level variants were found in the WES VCF."))
  }

  vus_count <- sum(wes_variants$classification == "VUS")

  div(
    class = "wes-variants-panel",
    div(
      class = "targets-header",
      div(
        div(class = "targets-title", "WES variants"),
        div(
          class = "targets-criteria",
          tags$span(class = "criteria-chip", "PASS"),
          tags$span(class = "criteria-chip", "PCGR annotated"),
          tags$span(class = "criteria-chip", paste0(vus_count, " VUS"))
        )
      ),
      div(
        class = "targets-count",
        strong(nrow(wes_variants)),
        " variants · ",
        n_distinct(wes_variants$gene_symbol),
        " genes"
      )
    ),
    div(
      class = "wes-variants-note",
      "Variants are shown for genes present in the selected drug network."
    ),
    div(class = "wes-variants-table", DT::DTOutput("wes_variants_table"))
  )
}
