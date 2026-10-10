# Reproducible PDF report generation shared by the Shiny UI and Nextflow.

report_dictionary <- list(
  en = list(
    title = "ClinReport - final drug list", report_title = "Report title", patient = "Patient / sample:", status = "List status",
    automatic = "Automatically generated", expert = "Reviewed and changed by an expert",
    automatic_note = "No expert changes were made to the composition or order.",
    expert_note = "The composition or order differs from the automatic list for the selected Top N.",
    drug_count = "Number of drugs:", graph_count = "Network graphs:", generated = "Generated:",
    expert_label = "Expert:", comment_label = "Expert comment",
    disclaimer = "This report supports expert assessment and does not replace a clinical decision.",
    final_list = "Final list", score_note = "priority score is not a clinical probability of efficacy",
    molecular = "Molecular rationale", molecular_note = "Summary of evidence available in the current analysis run",
    clinical = "Clinical and pharmacological context",
    clinical_note = "Summary from available annotations; missing records do not mean missing clinical evidence",
    network = "Network context:", network_note = "One-drug graph; neighbours limited to",
    network_missing = "Network data are not available for this drug", drug_target = "Drug target",
    network_neighbour = "Network neighbour", history = "Expert correction history",
    added = "Added:", removed = "Removed:", moved = "Drug order changed.",
    targets = "Targets", de = "DE", wes = "WES", tf = "TF", drug = "Drug", id = "ID",
    status_header = "Status", original_rank = "Original rank", change = "Change",
    rank = "Rank", approval = "First approval", phase = "Phase", indications = "Indications",
    warnings = "Warnings", fallback_title = "ClinReport - packaging mode",
    fallback_text = "The full drug list and network graphs will be available when the prepared application is launched.",
    fallback_defaults = "Defaults:", fallback_reason = "Technical reason:",
    available = "available", not_provided = "not provided"
  ),
  ru = list(
    title = "ClinReport - итоговый список лекарств", report_title = "Название отчёта", patient = "Пациент / образец:", status = "Статус списка",
    automatic = "Сформирован автоматически", expert = "Проверен и изменён экспертом",
    automatic_note = "Экспертская корректировка состава и порядка не внесена.",
    expert_note = "Состав или порядок лекарств отличается от автоматического списка для выбранного Top N.",
    drug_count = "Количество лекарств:", graph_count = "Графиков сетевого контекста:", generated = "Сформирован:",
    expert_label = "Эксперт:", comment_label = "Комментарий эксперта",
    disclaimer = "Отчёт предназначен для поддержки экспертной оценки и не заменяет клиническое решение.",
    final_list = "Финальный список", score_note = "priority score не является клинической вероятностью эффективности",
    molecular = "Молекулярное обоснование", molecular_note = "Сводка по данным, доступным в текущем запуске анализа",
    clinical = "Клинический и фармакологический контекст",
    clinical_note = "Сводная информация из доступных аннотаций; отсутствие записи не означает отсутствие клинических данных",
    network = "Сетевой контекст:", network_note = "Показан граф для одного препарата; соседей ограничено",
    network_missing = "Сетевые данные для препарата недоступны", drug_target = "Мишень лекарства",
    network_neighbour = "Сетевой сосед", history = "История экспертской корректировки",
    added = "Добавлены:", removed = "Исключены:", moved = "Изменён порядок лекарств.",
    targets = "Мишени", de = "DE", wes = "WES", tf = "TF", drug = "Лекарство", id = "ID",
    status_header = "Статус", original_rank = "Исх. ранг", change = "Изменение",
    rank = "Ранг", approval = "Первое одобрение", phase = "Фаза", indications = "Показания",
    warnings = "Предупреждения", fallback_title = "ClinReport - режим упаковки",
    fallback_text = "Полный список лекарств и сетевые графики будут доступны при запуске подготовленного приложения.",
    fallback_defaults = "Параметры по умолчанию:", fallback_reason = "Техническая причина:",
    available = "доступен", not_provided = "не предоставлен"
  )
)

report_tr <- function(key, language = "en") {
  language <- if (language %in% names(report_dictionary)) language else "en"
  value <- report_dictionary[[language]][[key]]
  if (is.null(value)) key else value
}

report_default_order <- function(drugs = drug_network) {
  drugs |>
    transmute(
      drug_id = as.character(drugId),
      label = as.character(label),
      score = suppressWarnings(as.numeric(score))
    ) |>
    arrange(desc(score), label, drug_id) |>
    pull(drug_id)
}

report_scalar <- function(x, default = "-") {
  if (is.null(x) || length(x) == 0) {
    return(default)
  }
  value <- as.character(x[[1]])[1]
  if (is.na(value) || !nzchar(value)) return(default)
  value
}

report_wrap <- function(x, width = 92) {
  x <- report_scalar(x, "")
  if (!nzchar(x)) return(character())
  strwrap(x, width = width, simplify = TRUE)
}

report_canvas <- function() {
  par(mar = c(0.25, 0.25, 0.25, 0.25))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), asp = 1)
}

report_text <- function(x, y, text, width = 92, cex = 0.82, col = "#17202a", font = 1,
                        line_height = 0.032, adj = 0) {
  lines <- report_wrap(text, width)
  if (!length(lines)) return(y)
  text(x, y - seq_along(lines) * line_height, lines, adj = adj, cex = cex, col = col, font = font)
  y - length(lines) * line_height
}

report_header <- function(title, subtitle = NULL) {
  text(0.06, 0.95, title, adj = 0, cex = 1.35, font = 2, col = "#0f766e")
  if (!is.null(subtitle) && nzchar(subtitle)) {
    text(0.06, 0.915, subtitle, adj = 0, cex = 0.78, col = "#64748b")
  }
  segments(0.06, 0.885, 0.94, 0.885, col = "#cbd5e1")
}

report_draw_table <- function(headers, rows, widths, y_top = 0.84, row_height = 0.045,
                              font_size = 0.63, header_size = 0.65) {
  if (!length(rows)) return(invisible(y_top))
  x <- c(0.06, 0.06 + cumsum(widths) * 0.88)
  rect(x[1], y_top - row_height, x[length(x)], y_top, col = "#e6f2f0", border = NA)
  for (i in seq_along(headers)) {
    text((x[i] + x[i + 1]) / 2, y_top - row_height * 0.68, headers[[i]],
      cex = header_size, font = 2, col = "#134e4a")
  }
  for (r in seq_len(nrow(rows))) {
    y1 <- y_top - row_height * (r + 1)
    y2 <- y_top - row_height * r
    if (r %% 2 == 0) rect(x[1], y1, x[length(x)], y2, col = "#f8fafc", border = NA)
    for (i in seq_along(headers)) {
      value <- report_scalar(rows[[i]][r], "")
      text((x[i] + x[i + 1]) / 2, (y1 + y2) / 2, value,
        cex = font_size, col = "#17202a")
    }
  }
  segments(x[1], y_top, x[length(x)], y_top, col = "#94a3b8")
  segments(x[1], y_top - row_height * (nrow(rows) + 1), x[length(x)],
    y_top - row_height * (nrow(rows) + 1), col = "#94a3b8")
  invisible(y_top - row_height * (nrow(rows) + 1))
}

report_target_metrics <- function(drug_id) {
  genes <- tryCatch(get_drug_genes(drug_id), error = function(e) data.frame())
  target_symbols <- if (nrow(genes) && "gene_name" %in% names(genes)) {
    unique(as.character(genes$gene_name))
  } else {
    character()
  }
  deg <- if (nrow(genes) && "significant_both" %in% names(genes)) {
    sum(genes$significant_both %in% TRUE, na.rm = TRUE)
  } else {
    0L
  }
  wes_symbols <- if (exists("has_wes") && isTRUE(has_wes) && exists("wes_gene_summary") &&
      "gene_symbol" %in% names(wes_gene_summary)) wes_gene_summary$gene_symbol else character()
  tf_symbols <- if (exists("has_tf_analysis") && isTRUE(has_tf_analysis) && exists("tf_results") &&
      "TF" %in% names(tf_results)) unique(tf_results$TF) else character()
  data.frame(
    targets = nrow(genes),
    deg = deg,
    wes = sum(target_symbols %in% wes_symbols, na.rm = TRUE),
    tf = sum(target_symbols %in% tf_symbols, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

report_clinical_context <- function(drug_id) {
  data <- tryCatch(load_drug(drug_id), error = function(e) list())
  has_data <- length(data) > 0
  molecule <- data$chembl$chembl_molecule
  indications <- data$chembl$chembl_indications$drug_indications
  warnings <- data$chembl$chembl_warnings$drug_warnings
  indication_text <- if (length(indications)) {
    paste(unique(vapply(indications, function(x) report_scalar(x$efo_term, ""), character(1))), collapse = "; ")
  } else "-"
  mechanism <- data$chembl$chembl_mechanisms$mechanisms
  mechanism_text <- if (length(mechanism)) {
    paste(unique(vapply(mechanism, function(x) report_scalar(x$mechanism_of_action, ""), character(1))), collapse = "; ")
  } else "-"
  data.frame(
    approval = report_scalar(molecule$first_approval),
    phase = report_scalar(molecule$max_phase),
    indications = indication_text,
    mechanism = mechanism_text,
    warnings = if (!has_data) "-" else if (length(warnings)) as.character(length(warnings)) else "0",
    stringsAsFactors = FALSE
  )
}

report_build_context <- function(final_ids, baseline_ids, top_n, graph_top_n,
                                 expert_name = "", expert_comment = "", report_title = "") {
  final_ids <- unique(as.character(final_ids))
  baseline_ids <- unique(as.character(baseline_ids))
  current_positions <- setNames(seq_along(final_ids), final_ids)
  baseline_positions <- setNames(seq_along(baseline_ids), baseline_ids)
  added <- setdiff(final_ids, baseline_ids)
  removed <- setdiff(baseline_ids, final_ids)
  same_set <- setequal(final_ids, baseline_ids)
  moved <- if (same_set && !identical(final_ids, baseline_ids)) final_ids else character()
  expert_adjusted <- !identical(final_ids, baseline_ids)
  list(
    final_ids = final_ids,
    baseline_ids = baseline_ids,
    top_n = as.integer(top_n),
    graph_top_n = min(as.integer(graph_top_n), length(final_ids)),
    expert_adjusted = expert_adjusted,
    expert_name = report_scalar(expert_name, ""),
    expert_comment = report_scalar(expert_comment, ""),
    report_title = report_scalar(report_title, ""),
    added = added,
    removed = removed,
    moved = moved,
    current_positions = current_positions,
    baseline_positions = baseline_positions,
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
}

report_draw_network <- function(drug_id, language = "en") {
  drug_row <- drug_network[match(drug_id, drug_network$drugId), , drop = FALSE]
  drug_label <- if (nrow(drug_row)) report_scalar(drug_row$label, drug_id) else drug_id
  graph <- tryCatch(
    get_drug_interaction_subgraph(drug_id, max_neighbors = min(interaction_network_max_neighbors, 24)),
    error = function(e) list(nodes = data.frame(), edges = data.frame(), targets = character())
  )
  report_canvas()
  report_header(paste(report_tr("network", language), drug_label),
    paste(report_tr("network_note", language), min(interaction_network_max_neighbors, 24)))
  if (!nrow(graph$nodes)) {
    text(0.5, 0.52, report_tr("network_missing", language), cex = 1, col = "#64748b")
    return(invisible())
  }
  nodes <- graph$nodes
  nodes$gene_id <- as.character(nodes$gene_id)
  targets <- nodes$gene_id[nodes$type == "Target"]
  neighbors <- nodes$gene_id[nodes$type != "Target"]
  neighbors <- head(neighbors, 24)
  ids <- c(targets, neighbors)
  angles <- if (length(targets)) seq(0, 2 * pi, length.out = length(targets) + 1)[-1] else numeric()
  coords <- data.frame(id = ids, x = NA_real_, y = NA_real_, stringsAsFactors = FALSE)
  if (length(targets)) {
    coords$x[seq_along(targets)] <- 0.5 + 0.22 * cos(angles)
    coords$y[seq_along(targets)] <- 0.53 + 0.22 * sin(angles)
  }
  if (length(neighbors)) {
    angles2 <- seq(0, 2 * pi, length.out = length(neighbors) + 1)[-1]
    idx <- length(targets) + seq_along(neighbors)
    coords$x[idx] <- 0.5 + 0.38 * cos(angles2)
    coords$y[idx] <- 0.53 + 0.38 * sin(angles2)
  }
  if (nrow(graph$edges)) {
    for (i in seq_len(nrow(graph$edges))) {
      a <- coords[match(graph$edges$source[i], coords$id), ]
      b <- coords[match(graph$edges$target[i], coords$id), ]
      if (nrow(a) && nrow(b) && all(is.finite(c(a$x, a$y, b$x, b$y)))) {
        segments(a$x, a$y, b$x, b$y, col = "#cbd5e1")
      }
    }
  }
  points(0.5, 0.53, pch = 23, bg = "#0f766e", col = "#134e4a", cex = 2.5)
  text(0.5, 0.47, drug_label, cex = 0.76, font = 2, col = "#134e4a")
  if (nrow(coords)) {
    is_target <- coords$id %in% targets
    points(coords$x[is_target], coords$y[is_target], pch = 21, bg = "#dbeafe", col = "#2563eb", cex = 1.55)
    points(coords$x[!is_target], coords$y[!is_target], pch = 21, bg = "#f1f5f9", col = "#94a3b8", cex = 1.2)
    labels <- nodes$gene_name[match(coords$id, nodes$gene_id)]
    labels <- ifelse(is.na(labels) | !nzchar(as.character(labels)), coords$id, as.character(labels))
    labels <- substr(labels, 1, 16)
    text(coords$x, coords$y - 0.035, labels, cex = 0.52, col = "#334155")
  }
  legend("bottom", legend = c(report_tr("drug_target", language), report_tr("network_neighbour", language)), pch = c(21, 21),
    pt.bg = c("#dbeafe", "#f1f5f9"), col = c("#2563eb", "#94a3b8"), bty = "n", cex = 0.7)
  invisible()
}

report_open_device <- function(path) {
  tryCatch(
    grDevices::cairo_pdf(path, width = 8.27, height = 11.69, onefile = TRUE, family = "DejaVu Sans"),
    error = function(e) grDevices::pdf(path, width = 8.27, height = 11.69, onefile = TRUE, useDingbats = FALSE)
  )
}

generate_clinreport_fallback_pdf <- function(path, patient = "-", reason = "", top_n = 15,
                                             graph_top_n = 3, state_json_path = NULL,
                                             language = "en", report_title = "") {
  language <- if (language %in% names(report_dictionary)) language else "en"
  report_title <- report_scalar(report_title,
    paste(patient, format(Sys.time(), "%Y-%m-%d %H-%M-%S")))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  report_open_device(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  report_canvas()
  report_header(report_title, paste(report_tr("fallback_title", language), "-", report_tr("patient", language), patient))
  text(0.06, 0.80, report_tr("fallback_title", language), adj = 0, cex = 1.0, font = 2, col = "#0f766e")
  fallback_tail <- if (language == "ru") "и графиков." else "and graphs."
  report_text(0.06, 0.73,
    paste(report_tr("fallback_text", language), report_tr("fallback_defaults", language),
      paste("Top", top_n), paste("Top", graph_top_n), fallback_tail), cex = 0.82)
  if (nzchar(reason)) report_text(0.06, 0.58, paste(report_tr("fallback_reason", language), reason), cex = 0.68, col = "#64748b")
  if (!is.null(state_json_path)) {
    if (requireNamespace("jsonlite", quietly = TRUE)) {
      jsonlite::write_json(list(
        generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
        patient_id = patient,
        report_title = report_title,
        language = language,
        top_n = as.integer(top_n),
        graph_top_n = as.integer(graph_top_n),
        expert_adjusted = FALSE,
        baseline_ids = character(),
        final_ids = character()
      ), state_json_path, auto_unbox = TRUE, pretty = TRUE)
    } else {
      writeLines(c(
        "{",
        paste0('  "patient_id": "', gsub('"', '\\\"', patient, fixed = TRUE), '",'),
        paste0('  "top_n": ', as.integer(top_n), ','),
        paste0('  "graph_top_n": ', as.integer(graph_top_n), ','),
        '  "expert_adjusted": false,',
        '  "baseline_ids": [],',
        '  "final_ids": []',
        "}"
      ), state_json_path, useBytes = TRUE)
    }
  }
  invisible(path)
}

generate_clinreport_pdf <- function(path, final_ids, baseline_ids = final_ids, top_n = length(final_ids),
                                    graph_top_n = 3, expert_name = "", expert_comment = "",
                                    state_json_path = NULL, language = "en", report_title = "") {
  if (!length(final_ids)) stop("The final report list must contain at least one drug")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  language <- if (language %in% names(report_dictionary)) language else "en"
  patient <- if (exists("patient_id")) report_scalar(patient_id, "-") else "-"
  report_title <- report_scalar(report_title,
    paste(patient, format(Sys.time(), "%Y-%m-%d %H-%M-%S")))
  context <- report_build_context(final_ids, baseline_ids, top_n, graph_top_n,
    expert_name, expert_comment, report_title)
  report_open_device(path)
  on.exit(grDevices::dev.off(), add = TRUE)

  report_canvas()
  report_header(context$report_title, paste(report_tr("title", language), "-", report_tr("patient", language), patient))
  text(0.06, 0.82, report_tr("status", language), adj = 0, cex = 0.78, font = 2, col = "#475569")
  status <- if (context$expert_adjusted) report_tr("expert", language) else report_tr("automatic", language)
  status_note <- if (context$expert_adjusted) {
    report_tr("expert_note", language)
  } else {
    report_tr("automatic_note", language)
  }
  rect(0.06, 0.69, 0.94, 0.79, col = if (context$expert_adjusted) "#fff7ed" else "#ecfdf5",
    border = if (context$expert_adjusted) "#f59e0b" else "#10b981")
  text(0.08, 0.75, status, adj = 0, cex = 1.0, font = 2,
    col = if (context$expert_adjusted) "#9a3412" else "#047857")
  text(0.08, 0.715, report_wrap(status_note, 105), adj = 0, cex = 0.72, col = "#475569")
  y <- 0.62
  y <- report_text(0.06, y, paste(report_tr("drug_count", language), length(context$final_ids)), cex = 0.82, font = 2)
  y <- report_text(0.06, y - 0.02, paste(report_tr("graph_count", language), context$graph_top_n), cex = 0.82)
  y <- report_text(0.06, y - 0.02, paste(report_tr("generated", language), context$generated_at), cex = 0.82)
  if (nzchar(context$expert_name)) y <- report_text(0.06, y - 0.02, paste(report_tr("expert_label", language), context$expert_name), cex = 0.82)
  if (nzchar(context$expert_comment)) {
    text(0.06, y - 0.02, report_tr("comment_label", language), adj = 0, cex = 0.82, font = 2, col = "#475569")
    report_text(0.06, y - 0.055, context$expert_comment, cex = 0.78)
  }
  text(0.06, 0.18, report_tr("disclaimer", language),
    adj = 0, cex = 0.7, col = "#64748b")

  drug_rows <- drug_network[match(context$final_ids, drug_network$drugId), , drop = FALSE]
  baseline_rank <- match(context$final_ids, context$baseline_ids)
  change <- vapply(seq_along(context$final_ids), function(i) {
    id <- context$final_ids[i]
    if (id %in% context$added) report_tr("added", language) else if (id %in% context$moved) report_tr("moved", language) else ""
  }, character(1))
  table_rows <- data.frame(
    rank = seq_along(context$final_ids),
    drug = as.character(drug_rows$label),
    id = as.character(drug_rows$drugId),
    status = ifelse(nzchar(as.character(drug_rows$status)), as.character(drug_rows$status), "-"),
    score = formatC(as.numeric(drug_rows$score), format = "f", digits = 3),
    original_rank = ifelse(is.na(baseline_rank), "-", baseline_rank),
    change = change,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  report_canvas()
  report_header(report_tr("final_list", language), paste("Top", context$top_n, "-", report_tr("score_note", language)))
  report_draw_table(c(report_tr("rank", language), report_tr("drug", language), report_tr("id", language),
    report_tr("status_header", language), "Score", report_tr("original_rank", language), report_tr("change", language)),
    table_rows, c(0.055, 0.27, 0.14, 0.12, 0.10, 0.10, 0.215),
    row_height = min(0.045, 0.75 / max(1, nrow(table_rows))), font_size = 0.58)

  evidence_rows <- lapply(context$final_ids, function(id) {
    m <- report_target_metrics(id)
    row <- drug_network[match(id, drug_network$drugId), , drop = FALSE]
    data.frame(
      drug = report_scalar(row$label, id),
      targets = m$targets,
      de = m$deg,
      wes = m$wes,
      tf = m$tf,
      stringsAsFactors = FALSE
    )
  }) |> bind_rows()
  report_canvas()
  report_header(report_tr("molecular", language), report_tr("molecular_note", language))
  report_draw_table(c(report_tr("drug", language), report_tr("targets", language), report_tr("de", language), report_tr("wes", language), report_tr("tf", language)), evidence_rows,
    c(0.42, 0.14, 0.14, 0.14, 0.16),
    row_height = min(0.045, 0.75 / max(1, nrow(evidence_rows))), font_size = 0.65)
  y <- 0.14
  report_text(0.06, y, paste("DE:", "padj <", padj_cutoff, "; |log2FC| >=", log2fc_cutoff), cex = 0.72)
  report_text(0.06, y - 0.05,
    paste("WES:", if (exists("has_wes") && isTRUE(has_wes)) report_tr("available", language) else report_tr("not_provided", language),
      "- TF:", if (exists("has_tf_analysis") && isTRUE(has_tf_analysis)) report_tr("available", language) else report_tr("not_provided", language)), cex = 0.72)

  clinical_rows <- lapply(context$final_ids, function(id) {
    row <- drug_network[match(id, drug_network$drugId), , drop = FALSE]
    clinical <- report_clinical_context(id)
    data.frame(
      drug = report_scalar(row$label, id),
      approval = clinical$approval,
      phase = clinical$phase,
      indications = substr(clinical$indications, 1, 42),
      warnings = clinical$warnings,
      stringsAsFactors = FALSE
    )
  }) |> bind_rows()
  report_canvas()
  report_header(report_tr("clinical", language), report_tr("clinical_note", language))
  report_draw_table(c(report_tr("drug", language), report_tr("approval", language), report_tr("phase", language),
    report_tr("indications", language), report_tr("warnings", language)), clinical_rows,
    c(0.29, 0.16, 0.10, 0.31, 0.14),
    row_height = min(0.045, 0.75 / max(1, nrow(clinical_rows))), font_size = 0.48)

  context_ids <- head(context$final_ids, context$graph_top_n)
  for (id in context_ids) report_draw_network(id, language)

  if (length(context$added) || length(context$removed) || length(context$moved)) {
    report_canvas()
    report_header(report_tr("history", language))
    y <- 0.82
    if (length(context$added)) y <- report_text(0.08, y, paste(report_tr("added", language), paste(context$added, collapse = ", ")), cex = 0.78)
    if (length(context$removed)) y <- report_text(0.08, y - 0.02, paste(report_tr("removed", language), paste(context$removed, collapse = ", ")), cex = 0.78)
    if (length(context$moved)) y <- report_text(0.08, y - 0.02, report_tr("moved", language), cex = 0.78)
    if (nzchar(context$expert_comment)) report_text(0.08, y - 0.02, paste(report_tr("comment_label", language), context$expert_comment), cex = 0.78)
  }

  if (!is.null(state_json_path)) {
    jsonlite::write_json(
      list(
        generated_at = context$generated_at,
        patient_id = patient,
        report_title = context$report_title,
        language = language,
        top_n = context$top_n,
        graph_top_n = context$graph_top_n,
        expert_adjusted = context$expert_adjusted,
        expert_name = context$expert_name,
        expert_comment = context$expert_comment,
        baseline_ids = context$baseline_ids,
        final_ids = context$final_ids
      ),
      state_json_path, auto_unbox = TRUE, pretty = TRUE, null = "null"
    )
  }
  invisible(path)
}
