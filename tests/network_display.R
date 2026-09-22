# Pure display rules: no Shiny installation or patient data required.
source('R/domain/drug_network.R')
source('R/domain/tf_network.R')
graph <- list(nodes = data.frame(
  gene_id = c('target', 'member', 'context'),
  gene_name = c('T', 'M', 'C'),
  type = c('Target', 'Network neighbor', 'Network neighbor')))
data <- list(
  nodes = data.frame(id = c('__DRUG__', graph$nodes$gene_id, 'tf')),
  edges = data.frame(from = c('__DRUG__', 'target', 'target', 'tf'),
                     to = c('target', 'member', 'context', 'target')))
hidden <- apply_context_visibility(data, graph)
stopifnot(identical(hidden$nodes$hidden, c(FALSE, FALSE, TRUE, TRUE, FALSE)),
          identical(hidden$edges$hidden, c(FALSE, TRUE, TRUE, FALSE)),
          identical(hidden$nodes$physics, !hidden$nodes$hidden),
          identical(hidden$edges$physics, !hidden$edges$hidden))
overlay <- apply_context_visibility(data, graph, 'M')
stopifnot(identical(overlay$nodes$hidden, c(FALSE, FALSE, FALSE, TRUE, FALSE)),
          identical(overlay$edges$hidden, c(FALSE, FALSE, TRUE, FALSE)))
visible <- apply_context_visibility(data, graph, show_context = TRUE)
stopifnot(!any(visible$nodes$hidden), !any(visible$edges$hidden),
          identical(visible$nodes$id, data$nodes$id))
tf_results <- data.frame(TF = c('STAT3', 'MYC', 'TP53'))
colors <- tf_identity_color(tf_results$TF)
stopifnot(length(unique(colors)) == 3L, !anyNA(colors),
          identical(tf_identity_color('MYC'), colors[2]))
tf_results <- tf_results[3:1, , drop = FALSE]
stopifnot(identical(tf_identity_color(c('STAT3', 'MYC', 'TP53')), colors))
cat('PASS: context visibility, overlay members, edge visibility and stable TF colors\n')

# A gene can carry multiple evidence flags; the remainder is not total minus
# the sum of the three counts. Imported regulators are not automatically targets.
summary <- target_evidence_summary(c('a', 'b', 'c'), c('A', 'B', 'C'),
  'a', 'A', c('B', 'OUTSIDE_REGULATOR'))
stopifnot(grepl('1 with significant expression changes', summary, fixed = TRUE),
  grepl('1 with WES variants', summary, fixed = TRUE),
  grepl('1 TFs in the imported analysis', summary, fixed = TRUE),
  grepl('1 without these flags', summary, fixed = TRUE))
missing <- target_evidence_summary('a', 'A', character())
stopifnot(grepl('WES not supplied', missing, fixed = TRUE),
  grepl('TF analysis not supplied', missing, fixed = TRUE))
cat('PASS: target evidence overlap, remainder, TF membership and missing inputs\n')

# Hidden WES-positive neighbours must not count, but visible overlay genes do.
data$nodes$label <- c('Drug', 'T', 'M', 'C', 'TF')
data$nodes$mutated <- c(FALSE, FALSE, TRUE, TRUE, TRUE)
stopifnot(count_visible_wes_genes(apply_context_visibility(data, graph)$nodes) == 1L,
  count_visible_wes_genes(apply_context_visibility(data, graph, 'M')$nodes) == 2L,
  count_visible_wes_genes(apply_context_visibility(data, graph, show_context = TRUE)$nodes) == 3L)
data$nodes$mutated <- c(FALSE, FALSE, FALSE, TRUE, FALSE)
stopifnot(count_visible_wes_genes(apply_context_visibility(data, graph)$nodes) == 0L)
cat('PASS: WES count follows visible context and overlay genes\n')
