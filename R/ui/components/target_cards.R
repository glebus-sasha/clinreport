# Target card components

# 18. TARGET CARDS
# ============================================================

render_target_cards <- function(
    genes
) {
  
  if (nrow(genes) == 0) {
    
    return(
      div(
        class = "empty-message",
        "No target genes found."
      )
    )
  }
  
  
  map(
    seq_len(nrow(genes)),
    function(i) {
      
      g <- genes[i, ]
      
      
      padj_pass <- isTRUE(
        !is.na(g$padj) &&
          g$padj < padj_cutoff
      )
      
      fc_pass <- isTRUE(
        !is.na(g$log2FoldChange) &&
          abs(g$log2FoldChange) >= log2fc_cutoff
      )
      
      both_pass <- padj_pass && fc_pass
      
      
      direction_class <-
        
        if (
          is.na(g$log2FoldChange)
        ) {
          
          "neutral"
          
        } else if (
          g$log2FoldChange > 0
        ) {
          
          "up"
          
        } else if (
          g$log2FoldChange < 0
        ) {
          
          "down"
          
        } else {
          
          "neutral"
        }
      
      
      cutoff_class <-
        
        if (both_pass) {
          
          "passes-both"
          
        } else if (padj_pass) {
          
          "passes-padj"
          
        } else if (fc_pass) {
          
          "passes-fc"
          
        } else {
          
          "passes-neither"
        }
      
      
      card_class <- paste(
        "gene-card",
        direction_class,
        cutoff_class
      )
      
      
      onclick_js <- sprintf(
        paste0(
          "Shiny.setInputValue(",
          "'selected_gene_index', ",
          "%d, ",
          "{priority: 'event'}",
          ")"
        ),
        i
      )
      
      
      gene_name <- ifelse(
        is.na(g$gene_name) ||
          !nzchar(g$gene_name),
        "Unknown",
        g$gene_name
      )
      
      
      fc_text <- ifelse(
        is.na(g$log2FoldChange),
        "—",
        sprintf(
          "%+.2f",
          g$log2FoldChange
        )
      )
      
      
      padj_text <- ifelse(
        is.na(g$padj),
        "—",
        format.pval(
          g$padj,
          digits = 2,
          eps = 1e-300
        )
      )
      
      
      div(
        
        class = card_class,
        
        onclick = onclick_js,
        
        
        div(
          class = "gene-main",
          
          div(
            class = "gene-symbol",
            gene_name
          ),
          
          div(
            class = "gene-id-small",
            
            tags$a(
              href = ensembl_url(
                g$gene_id_clean
              ),
              target = "_blank",
              rel = "noopener noreferrer",
              class = "gene-external-link",
              g$gene_id_clean
            )
          )
        ),
        
        
        div(
          class = "gene-expression",
          
          div(
            class = "expression-direction",
            
            span(
              class = paste0(
                "direction-arrow ",
                direction_class
              ),
              
              ifelse(
                direction_class == "up",
                "↑",
                ifelse(
                  direction_class == "down",
                  "↓",
                  "—"
                )
              )
            ),
            
            span(
              class = "expression-value",
              fc_text
            )
          )
        ),
        
        
        div(
          class = "gene-padj",
          
          div(
            class = "gene-stat-label",
            "padj"
          ),
          
          div(
            class = "gene-stat-value",
            padj_text
          )
        ),
        
        
        div(
          class = "gene-status",
          
          if (both_pass) {
            
            span(
              class = "cutoff-badge both",
              "meets criteria"
            )
            
          } else if (padj_pass) {
            
            span(
              class = "cutoff-badge padj-only",
              "padj"
            )
            
          } else if (fc_pass) {
            
            span(
              class = "cutoff-badge fc-only",
              "fold change"
            )
            
          } else {
            
            span(
              class = "cutoff-badge neither",
              "below criteria"
            )
          }
        ),
        
        
        div(
          class = "gene-arrow",
          "›"
        )
      )
    }
  )
}


