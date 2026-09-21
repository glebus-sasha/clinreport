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
 session$setInputs(drug_table_rows_selected=1, overlay_mode='tf',select_all_tfs=1)
 stopifnot(length(selected_tfs())==nrow(tf_results)); invisible(output$drug_network_graph)
 session$setInputs(clear_tfs=1); stopifnot(length(selected_tfs())==0)
 session$setInputs(selected_tf_index=1); stopifnot(identical(selected_tfs(),tf_results$TF[1]))
 session$setInputs(overlay_mode='pathways',selected_gsea_pathway=pathway_sets$pathway[1])
 invisible(output$drug_network_graph)
 session$setInputs(overlay_mode='tf',gsea_color_by_direction=TRUE,drug_table_rows_selected=2)
 invisible(output$drug_network_graph); invisible(output$tf_table)
})
cat('PASS: graph identities, preserved base edges, TF selection and overlay switching\n')
