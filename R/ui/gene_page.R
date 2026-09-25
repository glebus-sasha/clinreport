# Gene page composition

# 20. GENE PAGE
# ============================================================

render_gene_page <- function(
    gene,
    inline = FALSE
) {
  
  gene_name <- ifelse(
    is.na(gene$gene_name) ||
      !nzchar(gene$gene_name),
    "Unknown gene",
    gene$gene_name
  )

  gene_variants <- wes_variants |>
    filter(gene_symbol == gene_name)
  mutation_count <- nrow(gene_variants)
  drug_count <- sum(vapply(drug_network$drugId, function(drug_id)
    gene$gene_id_clean %in% get_drug_target_ids(drug_id), logical(1)))
  
  
  div(
    class = if (inline) "gene-page gene-page-inline" else "gene-page",
    if (inline) actionButton("close_gene_profile", "Close gene profile", class = "back-button"),
    
    div(
      class = "gene-header",
      
      div(
        class = "gene-title-large",
        gene_name
      ),
      
      div(
        class = "gene-id",
        
        tags$a(
          href = ensembl_url(
            gene$gene_id_clean
          ),
          target = "_blank",
          rel = "noopener noreferrer",
          class = "external-link",
          gene$gene_id_clean
        )
      )
    ),
    
    
    div(
      class = "gene-summary",
      
      div(
        class = "summary-card",
        
        div(
          class = "summary-label",
          "log2 Fold Change"
        ),
        
        div(
          class = "summary-value",
          
          ifelse(
            is.na(gene$log2FoldChange),
            "—",
            sprintf(
              "%+.3f",
              gene$log2FoldChange
            )
          )
        ),
        
        div(
          class = "summary-note",
          
          ifelse(
            is.na(gene$log2FoldChange),
            "Not available",
            ifelse(
              gene$log2FoldChange > 0,
              "Up-regulated",
              ifelse(
                gene$log2FoldChange < 0,
                "Down-regulated",
                "No change"
              )
            )
          )
        )
      ),
      
      
      div(
        class = "summary-card",
        
        div(
          class = "summary-label",
          "Adjusted p-value"
        ),
        
        div(
          class = "summary-value",
          
          ifelse(
            is.na(gene$padj),
            "—",
            format.pval(
              gene$padj,
              digits = 3,
              eps = 1e-300
            )
          )
        ),
        
        div(
          class = "summary-note",
          
          ifelse(
            !is.na(gene$padj) &&
              gene$padj < padj_cutoff,
            "Below threshold",
            "Not below threshold"
          )
        )
      ),
      
      
      div(
        class = "summary-card",
        
        div(
          class = "summary-label",
          "baseMean"
        ),
        
        div(
          class = "summary-value",
          
          ifelse(
            is.na(gene$baseMean),
            "—",
            format(
              round(
                gene$baseMean,
                2
              ),
              big.mark = ","
            )
          )
        ),
        
        div(
          class = "summary-note",
          "Mean normalized expression"
        )
      ),

      div(
        class = "summary-card",
        div(class = "summary-label", "Drug-linked modules"),
        div(class = "summary-value", drug_count),
        div(class = "summary-note", "Drugs linked to this gene")
      ),

      div(
        class = "summary-card",

        div(
          class = "summary-label",
          "WES variants"
        ),

        div(
          class = "summary-value",
          if (has_wes) mutation_count else "—"
        ),

        div(
          class = "summary-note",
          if (!has_wes) {
            "Not supplied"
          } else if (mutation_count == 0) {
            "No PASS variants"
          } else {
            wes_gene_summary |>
              filter(gene_symbol == gene_name) |>
              pull(mutation_classes) |>
              first(default = "PCGR-annotated")
          }
        )
      )
    ),

    if (has_wes && mutation_count > 0) div(
      class = "gene-wes-variants",
      div(class = "gene-wes-variants-title", "WES variant calls"),
      div(
        class = "gene-wes-variants-list",
        lapply(seq_len(mutation_count), function(i) {
          variant <- gene_variants[i, ]
          protein_text <- if (is.na(variant$protein_change)) variant$consequence else variant$protein_change
          clinvar_url <- clinvar_allele_url(variant$clinvar_allele_id)
          rsid_url <- dbsnp_url(variant$dbsnp_rsid)
          div(
            class = "gene-wes-variant",
            span(class = "gene-wes-coordinate", variant$variant_label),
            span(class = "gene-wes-protein", protein_text),
            span(class = "gene-wes-classification", variant$classification),
            if (!is.na(variant$vaf)) span(sprintf("VAF %.1f%%", variant$vaf * 100)),
            if (!is.null(clinvar_url)) tags$a(href = clinvar_url, target = "_blank", rel = "noopener noreferrer", "ClinVar ↗"),
            if (!is.null(rsid_url)) tags$a(href = rsid_url, target = "_blank", rel = "noopener noreferrer", "dbSNP ↗")
          )
        })
      )
    ),
    
    
    div(
      class = "gene-filter-status",
      
      if (
        !is.na(gene$padj) &&
        gene$padj < padj_cutoff
      ) {
        
        span(
          class = "gene-filter-badge pass",
          paste0(
            "padj < ",
            padj_cutoff
          )
        )
        
      } else {
        
        span(
          class = "gene-filter-badge fail",
          paste0(
            "padj ≥ ",
            padj_cutoff
          )
        )
      },
      
      
      if (
        !is.na(gene$log2FoldChange) &&
        abs(gene$log2FoldChange) >= log2fc_cutoff
      ) {
        
        span(
          class = "gene-filter-badge pass",
          paste0(
            "|log2FC| ≥ ",
            log2fc_cutoff
          )
        )
        
      } else {
        
        span(
          class = "gene-filter-badge fail",
          paste0(
            "|log2FC| < ",
            log2fc_cutoff
          )
        )
      }
    ),
    
    
    div(
      class = "section-header",
      
      div(
        class = "section-title",
        "Differential expression"
      ),
      
      div(
        class = "section-subtitle",
        "Genome-wide differential-expression context"
      )
    ),
    
    plotlyOutput(
      "gene_volcano",
      height = "165px"
    )
  )
}
