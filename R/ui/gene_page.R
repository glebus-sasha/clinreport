# Gene page composition

# 20. GENE PAGE
# ============================================================

render_gene_page <- function(
    gene
) {
  
  gene_name <- ifelse(
    is.na(gene$gene_name) ||
      !nzchar(gene$gene_name),
    "Unknown gene",
    gene$gene_name
  )
  
  
  div(
    
    actionButton(
      "back_to_drug",
      "← Back to drug",
      class = "back-button"
    ),
    
    
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
        "Gene expression overview"
      ),
      
      div(
        class = "section-subtitle",
        "Genome-wide differential expression context"
      )
    ),
    
    plotlyOutput(
      "gene_volcano",
      height = "600px"
    )
  )
}


