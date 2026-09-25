# Drug page composition

# 19. DRUG PAGE
# ============================================================

render_drug_page <- function(
    drug_id,
    drug_name,
    data,
    genes,
    network_name
) {
  
  chembl_id <- get_chembl_id(data)
  pubchem_cid <- get_pubchem_cid(data)

  evidence_badge <- function(label, value, note, status = "neutral", gene_list = character()) {
    div(class = paste("drug-module-evidence", status),
      div(class = "drug-module-evidence-label", label),
      div(class = "drug-module-evidence-value", value),
      div(class = "drug-module-evidence-note", note),
      if (length(gene_list) > 0) tags$details(
        class = "drug-module-gene-list",
        tags$summary("Show genes"),
        div(paste(sort(unique(gene_list)), collapse = " · "))))
  }

  deg_count <- sum(genes$significant_both %in% TRUE, na.rm = TRUE)
  wes_symbols <- if (has_wes) wes_gene_summary$gene_symbol else character()
  tf_symbols <- if (has_tf_analysis) unique(tf_results$TF) else character()
  wes_count <- sum(genes$gene_name %in% wes_symbols, na.rm = TRUE)
  tf_count <- sum(genes$gene_name %in% tf_symbols, na.rm = TRUE)
  deg_genes <- genes$gene_name[genes$significant_both %in% TRUE]
  tf_genes <- genes$gene_name[genes$gene_name %in% tf_symbols]
  wes_genes <- genes$gene_name[genes$gene_name %in% wes_symbols]
  expanded_genes <- genes$gene_name[!(genes$significant_both %in% TRUE |
    genes$gene_name %in% wes_symbols | genes$gene_name %in% tf_symbols)]
  expanded_count <- sum(!(genes$significant_both %in% TRUE |
    genes$gene_name %in% wes_symbols |
    genes$gene_name %in% tf_symbols), na.rm = TRUE)
  
  
  tagList(
    div(class = "drug-module-shell",
    div(class = "drug-module-shell-title", drug_name),
    tags$details(
    class = "drug-module",
    open = "open",
    tags$summary(class = "drug-module-toggle", "Drug module"),

    div(class = "drug-module-summary",
      div(class = "drug-module-summary-heading",
        div(class = "drug-module-summary-subtitle",
          paste(nrow(genes), "drug-linked target genes in the supplied network"))),
      div(class = "drug-module-evidence-grid",
        evidence_badge("DEG", deg_count, paste("of", nrow(genes), "targets"),
          if (deg_count > 0) "positive" else "neutral", deg_genes),
        if (isTRUE(as.logical(use_tf_activity)) && has_tf_analysis)
          evidence_badge("TF activity", tf_count, "targets that are TFs", "positive", tf_genes),
        if (isTRUE(as.logical(use_wes)) && has_wes)
          evidence_badge(wes_display_name, wes_count, "targets with variants", "positive", wes_genes),
        evidence_badge(paste(interaction_network_name, "expansion"), expanded_count,
          "drug-linked targets from network expansion", "expanded", expanded_genes)),
    div(
      class = "drug-network-panel combined-network-panel",
      div(
        class = "network-content-grid",
        div(
          class = "network-graph-column",
        div(
          class = "network-visualization-header",
          div(
            class = "network-header-copy",
            div(
              class = "network-visualization-title",
              network_name
            ),
            div(
              class = "network-visualization-subtitle",
              "Drug-linked targets and analysis context"
            )
          )
        ),
        uiOutput("drug_network_summary"),
        div(
          class = "gsea-sort network-color-toggle",
          tags$label(
            class = "radio-inline",
            tags$input(
              id = "gsea_color_by_direction",
              type = "checkbox",
              class = "shiny-input-checkbox"
            ),
            "Expression direction"
          ),
          tags$label(
            class = "radio-inline",
            title = "Show drug targets, their direct network neighbours and TF regulators. Disable to show all genes in selected overlays.",
            tags$input(
              id = "target_connected_only",
              type = "checkbox",
              class = "shiny-input-checkbox",
              checked = "checked"
            ),
            "Target-linked genes"
          ),
          tags$label(
            class = "radio-inline",
            title = "Show network neighbours even when they are not part of the selected pathways or TF regulons.",
            tags$input(id = "show_context_genes", type = "checkbox", class = "shiny-input-checkbox"),
            "Context genes"
          )
        ),
        div(class = "drug-module-legend",
          span(class = "legend-title", "Network legend"),
          span(class = "legend-item legend-target", "● Drug → target"),
          span(class = "legend-item legend-string", "┄ STRING interaction"),
          if (has_wes) span(class = "legend-item legend-wes", "· WES variant"),
          conditionalPanel("input.gsea_color_by_direction === true",
            span(class = "legend-item legend-expression", "↑ up · ↓ down")),
          conditionalPanel("input.overlay_mode !== 'tf'",
            span(class = "legend-item legend-pathway", "● Pathway membership")),
          if (has_tf_analysis) conditionalPanel("input.overlay_mode === 'tf'",
            tagList(
              span(class = "legend-item legend-tf", "▲ TF regulation"),
              span(class = "legend-item legend-regulation", "→ activation · ⊣ repression")))),
        visNetworkOutput("drug_network_graph", height = "500px")
        ),
        div(
          class = "gsea-panel embedded-gsea-panel",
        tags$button(
          id = "drug_network_toggle",
          type = "button",
          class = "network-panel-collapse",
          title = "Collapse network analysis panel",
          "›"
        ),
        div(class = "network-panel-body",
          uiOutput("inline_gene_panel"),
          render_network_overlay_panel())
        )
      ),
      conditionalPanel("input.overlay_mode !== 'tf'", div(class = "gsea-details-strip", uiOutput("gsea_pathway_details")))
    ),
    )
    ),
    
    div(
      class = "drug-information-sections",
      tags$details(
        class = "drug-information-section",
        open = "open",
        tags$summary(class = "drug-information-section-title", "Drug information"),
        tabsetPanel(
          id = "drug_information_tabs",
          selected = "Drug overview",
          if (has_drug_details) {
            tabPanel("Drug overview", render_general(drug_id, data))
          } else {
            tabPanel("Drug overview", div(
              class = "empty-message",
              "Drug annotation files were not supplied. Pharmacology, safety and structure details are unavailable."
            ))
          },
          if (has_drug_details) tabPanel("Mechanism & indications", render_pharmacology(data)),
          if (has_drug_details) tabPanel("Safety", render_safety(data))
        )
      ),
      tags$details(
        class = "drug-information-section",
        open = "open",
        tags$summary(class = "drug-information-section-title", "Molecular evidence"),
        tabsetPanel(
          id = "molecular_evidence_tabs",
          selected = "Direct targets",
          tabPanel(
            "Direct targets",
        
        div(
          class = "targets-header",
          
          div(
            
            div(
              class = "targets-title",
              "Target genes"
            ),
            
            div(
              class = "targets-criteria",
              
              tags$span(
                class = "criteria-chip",
                paste0(
                  "padj < ",
                  padj_cutoff
                )
              ),
              
              tags$span(
                class = "criteria-chip",
                paste0(
                  "|log2FC| ≥ ",
                  log2fc_cutoff
                )
              )
            )
          ),
          
          div(
            class = "targets-count",
            
            strong(
              nrow(genes)
            ),
            
            " ",
            
            ifelse(
              nrow(genes) == 1,
              "gene",
              "genes"
            )
          )
        ),
        
        div(
          class = "targets-list",
          
          render_target_cards(
            genes
          )
        )
          ),
          if (has_wes) tabPanel("WES evidence", render_wes_variants()),
          render_tf_tab()
        )
      )
    )
  )
  )
}
