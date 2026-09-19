# Volcano plot registration

register_volcano_output <- function(
    output,
    selected_gene,
    gene_data,
    padj_cutoff,
    log2fc_cutoff
) {
    output$gene_volcano <- renderPlotly({
      
      gene <- selected_gene()
      
      req(
        gene
      )
      
      
      selected_id <-
        gene$gene_id_clean
      
      
      plot_data <- gene_data |>
        mutate(
          
          neg_log10_padj = ifelse(
            is.na(padj),
            NA_real_,
            -log10(
              pmax(
                padj,
                .Machine$double.xmin
              )
            )
          ),
          
          selected =
            gene_id_clean == selected_id,
          
          pass_padj =
            !is.na(padj) &
            padj < padj_cutoff,
          
          pass_fc =
            !is.na(log2FoldChange) &
            abs(log2FoldChange) >= log2fc_cutoff,
          
          pass_both =
            pass_padj &
            pass_fc
        )
      
      
      p <- plot_ly()
      
      
      # --------------------------------------------------------
      # BACKGROUND
      # --------------------------------------------------------
      
      background_data <- plot_data |>
        filter(
          !selected,
          !pass_both
        )
      
      
      if (nrow(background_data) > 0) {
        
        p <- p |>
          
          add_trace(
            
            data = background_data,
            
            x = ~log2FoldChange,
            
            y = ~neg_log10_padj,
            
            type = "scatter",
            
            mode = "markers",
            
            text = ~paste0(
              "<b>",
              ifelse(
                is.na(gene_name),
                "Unknown",
                gene_name
              ),
              "</b><br>",
              "Gene ID: ",
              gene_id_clean,
              "<br>log2FC: ",
              ifelse(
                is.na(log2FoldChange),
                "NA",
                round(
                  log2FoldChange,
                  3
                )
              ),
              "<br>padj: ",
              ifelse(
                is.na(padj),
                "NA",
                signif(
                  padj,
                  3
                )
              )
            ),
            
            hoverinfo = "text",
            
            marker = list(
              size = 6,
              opacity = 0.22
            ),
            
            name = "Below criteria",
            
            inherit = FALSE
          )
      }
      
      
      # --------------------------------------------------------
      # SIGNIFICANT
      # --------------------------------------------------------
      
      significant_data <- plot_data |>
        filter(
          !selected,
          pass_both
        )
      
      
      if (nrow(significant_data) > 0) {
        
        p <- p |>
          
          add_trace(
            
            data = significant_data,
            
            x = ~log2FoldChange,
            
            y = ~neg_log10_padj,
            
            type = "scatter",
            
            mode = "markers",
            
            text = ~paste0(
              "<b>",
              ifelse(
                is.na(gene_name),
                "Unknown",
                gene_name
              ),
              "</b><br>",
              "Gene ID: ",
              gene_id_clean,
              "<br>log2FC: ",
              round(
                log2FoldChange,
                3
              ),
              "<br>padj: ",
              signif(
                padj,
                3
              )
            ),
            
            hoverinfo = "text",
            
            marker = list(
              size = 7,
              opacity = 0.70
            ),
            
            name = "Meets criteria",
            
            inherit = FALSE
          )
      }
      
      
      # --------------------------------------------------------
      # CUTOFFS
      # --------------------------------------------------------
      
      padj_y <- -log10(
        max(
          padj_cutoff,
          .Machine$double.xmin
        )
      )
      
      
      p <- p |>
        
        layout(
          
          title = list(
            text = "Differential expression",
            font = list(
              size = 16
            )
          ),
          
          xaxis = list(
            title = "log2 Fold Change",
            zeroline = TRUE
          ),
          
          yaxis = list(
            title = "-log10 adjusted p-value"
          ),
          
          hovermode = "closest",
          
          margin = list(
            l = 65,
            r = 25,
            b = 55,
            t = 55
          ),
          
          shapes = list(
            
            list(
              type = "line",
              x0 = -Inf,
              x1 = Inf,
              y0 = padj_y,
              y1 = padj_y,
              line = list(
                dash = "dash",
                width = 1
              )
            ),
            
            list(
              type = "line",
              x0 = -log2fc_cutoff,
              x1 = -log2fc_cutoff,
              y0 = 0,
              y1 = Inf,
              line = list(
                dash = "dash",
                width = 1
              )
            ),
            
            list(
              type = "line",
              x0 = log2fc_cutoff,
              x1 = log2fc_cutoff,
              y0 = 0,
              y1 = Inf,
              line = list(
                dash = "dash",
                width = 1
              )
            )
          )
        )
      
      
      # --------------------------------------------------------
      # SELECTED GENE
      # --------------------------------------------------------
      
      selected_row <- plot_data |>
        filter(
          selected
        )
      
      
      if (nrow(selected_row) > 0) {
        
        selected_padj <- selected_row$padj[1]
        
        selected_log2fc <-
          selected_row$log2FoldChange[1]
        
        selected_name <-
          selected_row$gene_name[1]
        
        
        selected_padj_text <-
          
          if (
            is.na(selected_padj)
          ) {
            
            "NA"
            
          } else {
            
            signif(
              selected_padj,
              3
            )
          }
        
        
        selected_fc_text <-
          
          if (
            is.na(selected_log2fc)
          ) {
            
            "NA"
            
          } else {
            
            round(
              selected_log2fc,
              3
            )
          }
        
        
        p <- p |>
          
          add_trace(
            
            data = selected_row,
            
            x = ~log2FoldChange,
            
            y = ~neg_log10_padj,
            
            type = "scatter",
            
            mode = "markers+text",
            
            text = ~gene_name,
            
            textposition = "top center",
            
            marker = list(
              size = 16,
              symbol = "diamond",
              line = list(
                width = 2
              )
            ),
            
            hoverinfo = "text",
            
            hovertext = paste0(
              "<b>",
              selected_name,
              "</b><br>",
              "Gene ID: ",
              selected_id,
              "<br>log2FC: ",
              selected_fc_text,
              "<br>padj: ",
              selected_padj_text,
              "<br><b>Selected gene</b>"
            ),
            
            name = "Selected gene",
            
            inherit = FALSE
          )
      }
      
      
      p
    })
}
