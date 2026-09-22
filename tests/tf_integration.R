invisible(source('app.R'))
stopifnot(is.null(dorothea_resource$error), nrow(tf_regulons)>0)
g <- get_drug_interaction_subgraph(drug_network$drugId[1])
base <- prepare_vis_network(g,drug_network[1,],setNames(character(),character()),character(),FALSE)
for (selected in list(character(),tf_results$TF[1],tf_results$TF)) {
 d <- build_tf_network(g,drug_network[1,],selected,character(),TRUE)
 stopifnot(!anyDuplicated(d$nodes$id),!anyDuplicated(d$edges$id),
   all(d$edges$from %in% d$nodes$id),all(d$edges$to %in% d$nodes$id),
   all(base$edges$id %in% d$edges$id),all(selected %in% d$nodes$label))
 cat(length(selected),'TFs:',nrow(d$tf_edges),'regulatory edges\n')
}
shiny::testServer(server, {
 session$setInputs(drug_table_rows_selected=1, overlay_mode='tf')
 targets <- current_drug_target_symbols()
 expected_pathways <- pathway_sets$pathway[vapply(pathway_sets$genes,
   function(genes) any(targets %in% genes), logical(1))]
 stopifnot(target_connected_only(), setequal(selected_pathways(), expected_pathways),
   setequal(selected_tfs(), get_target_tf_links(targets)$tf))
 session$setInputs(select_all_tfs=1, select_all_pathways=1)
 stopifnot(setequal(selected_pathways(), expected_pathways),
   setequal(selected_tfs(), get_target_tf_links(targets)$tf), !show_context_genes())
 session$setInputs(show_context_genes=TRUE)
 stopifnot(show_context_genes())
 invisible(output$drug_network_graph)
 session$setInputs(clear_tfs=1); stopifnot(length(selected_tfs())==0)
 session$setInputs(selected_tf_index=1); stopifnot(identical(selected_tfs(),tf_results$TF[1]))
 session$setInputs(overlay_mode='pathways',selected_gsea_pathway=pathway_sets$pathway[1])
 invisible(output$drug_network_graph)
 session$setInputs(overlay_mode='tf',gsea_color_by_direction=TRUE,drug_table_rows_selected=2)
 stopifnot(setequal(selected_tfs(), get_target_tf_links(current_drug_target_symbols())$tf))
 session$setInputs(target_connected_only=FALSE)
 stopifnot(!target_connected_only())
 session$setInputs(target_connected_only=TRUE)
 invisible(output$drug_network_graph); invisible(output$tf_table)
})
filtered <- build_tf_network(g, drug_network[1, ], tf_results$TF, character(), FALSE, TRUE)
target_symbols <- g$nodes$gene_name[g$nodes$gene_id %in% g$targets]
allowed <- unique(c(g$nodes$gene_name, get_target_tf_links(target_symbols)$tf))
stopifnot(all(filtered$nodes$label[filtered$nodes$id != '__DRUG__'] %in% allowed),
  all(base$edges$id %in% filtered$edges$id),
  all(filtered$tf_edges$to %in% g$nodes$gene_id),
  all(filtered$edges$from %in% filtered$nodes$id),
  all(filtered$edges$to %in% filtered$nodes$id))
cat('PASS: graph identities, preserved base edges, TF selection and overlay switching\n')

# A browser mounts the base widget before it can receive overlay proxy updates.
shiny::testServer(server, {
 network_messages <- list()
 original_send <- session$sendCustomMessage
 session$sendCustomMessage <- function(type, message) {
   if (grepl('drug_network_graph', jsonlite::toJSON(message, auto_unbox = TRUE), fixed = TRUE)) {
     network_messages[[length(network_messages) + 1L]] <<- message
   }
   original_send(type, message)
 }
 ready_generation <- function() {
   widget <- as.character(output$drug_network_graph)
   token <- regmatches(widget, regexpr("drug_network_ready', [0-9]+", widget))
   stopifnot(length(token) == 1L, nzchar(token))
   as.integer(sub("drug_network_ready', ", '', token, fixed = TRUE))
 }
 session$setInputs(drug_table_rows_selected = 1, overlay_mode = 'pathways')
 generation <- ready_generation()
 defaults <- selected_pathways()
 stopifnot(length(network_messages) == 0L)
 session$setInputs(drug_network_ready = generation)
 stopifnot(length(network_messages) > 0L, identical(defaults, selected_pathways()))

 session$setInputs(drug_table_rows_selected = 2)
 next_generation <- ready_generation()
 stopifnot(next_generation > generation)
 network_messages <- list()
 session$setInputs(drug_network_ready = generation)
 stopifnot(length(network_messages) == 0L)
 session$setInputs(drug_network_ready = next_generation)
 stopifnot(length(network_messages) > 0L)
})
cat('PASS: initial overlays wait for widget readiness; stale ready events are ignored\n')
