# Reproducible PDF report generation shared by the Shiny UI and Nextflow.

report_dictionary <- list(
  en = list(
    title = "O2T ClinReport - Drug List", report_title = "Report title", patient = "Patient:", status = "List status",
    automatic = "Automatically generated", expert = "Reviewed and changed by an expert",
    automatic_note = "No expert changes were made to the drug list composition or order.",
    expert_note = "The drug list composition or order differs from the automatic list for the selected Top N.",
    drug_count = "Number of drugs:", graph_count = "Number of drug modules:", generated = "Generated:",
    expert_label = "Expert:", comment_label = "Expert comment",
    disclaimer = "This report is informational and intended to support, not replace, review by a qualified healthcare professional. It must not be used as the sole basis for diagnosis or treatment decisions.",
    final_list = "Drug list", score_note = "priority score is not a clinical probability of efficacy",
    molecular = "Molecular rationale", molecular_note = "Summary of evidence available in the current analysis run",
    clinical = "Clinical and pharmacological context", clinical_continued = "Clinical and pharmacological context - continued",
    clinical_note = "Summary from available annotations; absence of an indication, warning, or other annotation does not mean there is no clinical evidence or risk",
    network = "Drug module:", network_note = "Drug targets only; colours show expression direction",
    network_missing = "Network data are not available for this drug", drug_target = "Drug - target",
    network_neighbour = "Network neighbour", direct_targets = "Direct targets", context_genes = "Context genes",
    de_up = "DE up", de_down = "DE down", expression_unknown = "Expression unknown", de_significant = "DE significant", de_pair = "DE up/down", wes_marker = "WES variant", tf_marker = "TF link", tf_target = "TF target",
    string_interaction = "STRING interaction",
    markers = "Markers:", targets_shown = "targets shown:", context_shown = "context genes shown:",
    history = "Expert correction history",
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
    title = "O2T ClinReport - Список лекарств", report_title = "Название отчёта", patient = "Пациент:", status = "Статус списка",
    automatic = "Сформирован автоматически", expert = "Проверен и изменён экспертом",
    automatic_note = "Экспертская корректировка состава и порядка списка лекарств не внесена.",
    expert_note = "Состав или порядок списка лекарств отличается от автоматического списка для выбранного Top N.",
    drug_count = "Количество лекарств:", graph_count = "Количество модулей препаратов:", generated = "Сформирован:",
    expert_label = "Эксперт:", comment_label = "Комментарий эксперта",
    disclaimer = "Отчёт носит информационный характер и предназначен для поддержки, а не замены оценки квалифицированного медицинского специалиста. Он не должен использоваться как единственное основание для диагностики или выбора лечения.",
    final_list = "Список лекарств", score_note = "priority score не является клинической вероятностью эффективности",
    molecular = "Молекулярное обоснование", molecular_note = "Сводка по данным, доступным в текущем запуске анализа",
    clinical = "Клинический и фармакологический контекст", clinical_continued = "Клинический и фармакологический контекст - продолжение",
    clinical_note = "Аннотации из доступных источников; отсутствие показания или предупреждения не означает отсутствия данных или рисков",
    network = "Модуль препарата:", network_note = "Только мишени лекарства; цвет показывает направление экспрессии",
    network_missing = "Сетевые данные для препарата недоступны", drug_target = "Лекарство - мишень",
    network_neighbour = "Сетевой сосед", direct_targets = "Прямые мишени", context_genes = "Контекстные гены",
    de_up = "DE вверх", de_down = "DE вниз", expression_unknown = "Направление экспрессии неизвестно", de_significant = "Значимое DE", de_pair = "DE вверх/вниз", wes_marker = "Вариант WES", tf_marker = "Связь с TF", tf_target = "TF-таргет",
    string_interaction = "Взаимодействие STRING",
    markers = "Маркеры:", targets_shown = "показано мишеней:", context_shown = "показано контекстных генов:",
    history = "История экспертской корректировки",
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

report_tr <- function(key, language = "ru") {
  language <- if (language %in% names(report_dictionary)) language else "en"
  value <- report_dictionary[[language]][[key]]
  if (is.null(value)) return(key)
  if (is.character(value)) Encoding(value) <- "UTF-8"
  value
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
  subtitle_lines <- character()
  if (!is.null(subtitle) && nzchar(subtitle)) {
    subtitle_lines <- report_wrap(subtitle, width = 112)
    text(0.06, 0.915 - (seq_along(subtitle_lines) - 1L) * 0.022, subtitle_lines,
      adj = 0, cex = 0.72, col = "#64748b")
  }
  separator_y <- 0.885 - max(0, length(subtitle_lines) - 1L) * 0.022
  segments(0.06, separator_y, 0.94, separator_y, col = "#cbd5e1")
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

report_draw_wrapped_table <- function(headers, rows, widths, wrap_cols = integer(), y_top = 0.84,
                                      base_row_height = 0.045, font_size = 0.48,
                                      header_size = 0.65, wrap_widths = NULL) {
  if (!length(rows) || !nrow(rows)) return(invisible(y_top))
  x <- c(0.06, 0.06 + cumsum(widths) * 0.88)
  if (is.null(wrap_widths)) wrap_widths <- pmax(10, round(widths * 105))
  wrap_widths <- rep(wrap_widths, length.out = length(headers))
  header_cells <- lapply(seq_along(headers), function(i) {
    if (i %in% wrap_cols) strwrap(headers[[i]], width = wrap_widths[[i]], simplify = TRUE) else headers[[i]]
  })
  header_height <- max(base_row_height, 0.014 + max(vapply(header_cells, length, integer(1))) * 0.026)
  rect(x[1], y_top - header_height, x[length(x)], y_top, col = "#e6f2f0", border = NA)
  for (i in seq_along(headers)) {
    if (i %in% wrap_cols) {
      lines <- header_cells[[i]]
      text(x[i] + 0.006, y_top - 0.022 - (seq_along(lines) - 1L) * 0.026,
        lines, adj = 0, cex = header_size, font = 2, col = "#134e4a")
    } else {
      text((x[i] + x[i + 1]) / 2, y_top - header_height * 0.68, header_cells[[i]],
        cex = header_size, font = 2, col = "#134e4a")
    }
  }
  y <- y_top - header_height
  for (r in seq_len(nrow(rows))) {
    cells <- lapply(seq_along(headers), function(i) {
      value <- report_scalar(rows[[i]][r], "")
      if (!(i %in% wrap_cols)) return(value)
      if (!nzchar(value)) return("-")
      strwrap(value, width = wrap_widths[[i]], simplify = TRUE)
    })
    line_count <- max(1L, vapply(cells, length, integer(1)))
    rh <- max(base_row_height, 0.014 + line_count * 0.026)
    y1 <- y - rh
    y2 <- y
    if (r %% 2 == 0) rect(x[1], y1, x[length(x)], y2, col = "#f8fafc", border = NA)
    for (i in seq_along(headers)) {
      if (i %in% wrap_cols) {
        lines <- cells[[i]]
        text(x[i] + 0.006, y2 - 0.022 - (seq_along(lines) - 1L) * 0.026,
          lines, adj = 0, cex = font_size, col = "#17202a")
      } else {
        text((x[i] + x[i + 1]) / 2, (y1 + y2) / 2, cells[[i]],
          cex = font_size, col = "#17202a")
      }
    }
    y <- y1
  }
  segments(x[-c(1, length(x))], y_top, x[-c(1, length(x))], y, col = "#e2e8f0")
  segments(x[1], y_top, x[length(x)], y_top, col = "#94a3b8")
  segments(x[1], y, x[length(x)], y, col = "#94a3b8")
  invisible(y)
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

report_network_annotations <- function(nodes) {
  nodes <- as.data.frame(nodes, stringsAsFactors = FALSE)
  gene_ids <- as.character(nodes$gene_id)
  labels <- as.character(nodes$gene_name)
  labels[is.na(labels) | !nzchar(labels)] <- gene_ids[is.na(labels) | !nzchar(labels)]
  log_fc <- rep(NA_real_, length(gene_ids))
  significant <- rep(FALSE, length(gene_ids))
  if (exists("gene_data") && is.data.frame(gene_data) && nrow(gene_data)) {
    id_col <- if ("gene_id_clean" %in% names(gene_data)) "gene_id_clean" else NULL
    name_col <- if ("gene_name" %in% names(gene_data)) "gene_name" else NULL
    fc_col <- if ("log2FoldChange" %in% names(gene_data)) "log2FoldChange" else NULL
    padj_col <- if ("padj" %in% names(gene_data)) "padj" else NULL
    sig_col <- if ("significant_both" %in% names(gene_data)) "significant_both" else NULL
    id_match <- if (!is.null(id_col)) match(gene_ids, as.character(gene_data[[id_col]])) else rep(NA_integer_, length(gene_ids))
    name_match <- if (!is.null(name_col)) match(labels, as.character(gene_data[[name_col]])) else rep(NA_integer_, length(gene_ids))
    id_match[is.na(id_match)] <- name_match[is.na(id_match)]
    if (!is.null(fc_col)) log_fc <- suppressWarnings(as.numeric(gene_data[[fc_col]][id_match]))
    if (!is.null(padj_col) && exists("padj_cutoff") && exists("log2fc_cutoff")) {
      padj <- suppressWarnings(as.numeric(gene_data[[padj_col]][id_match]))
      significant <- !is.na(padj) & padj < padj_cutoff &
        !is.na(log_fc) & abs(log_fc) >= log2fc_cutoff
    } else if (!is.null(sig_col)) {
      significant <- as.logical(gene_data[[sig_col]][id_match])
    }
    significant[is.na(significant)] <- FALSE
  }
  wes_symbols <- if (exists("has_wes") && isTRUE(has_wes) && exists("wes_gene_summary") &&
      "gene_symbol" %in% names(wes_gene_summary)) as.character(wes_gene_summary$gene_symbol) else character()
  tf_symbols <- if (exists("has_tf_analysis") && isTRUE(has_tf_analysis) && exists("tf_results") &&
      "TF" %in% names(tf_results)) as.character(tf_results$TF) else character()
  data.frame(
    gene_id = gene_ids,
    label = labels,
    type = as.character(nodes$type),
    log_fc = log_fc,
    significant = significant,
    wes = labels %in% wes_symbols,
    tf = labels %in% tf_symbols,
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

report_draw_network <- function(drug_id, language = "ru") {
  drug_row <- drug_network[match(drug_id, drug_network$drugId), , drop = FALSE]
  drug_label <- if (nrow(drug_row)) report_scalar(drug_row$label, drug_id) else drug_id
  graph <- tryCatch(
    get_drug_interaction_subgraph(drug_id, max_neighbors = min(interaction_network_max_neighbors, 16)),
    error = function(e) list(nodes = data.frame(), edges = data.frame(), targets = character())
  )
  report_canvas()
  report_header(paste(report_tr("network", language), drug_label))
  if (!nrow(graph$nodes)) {
    text(0.5, 0.52, report_tr("network_missing", language), cex = 1, col = "#64748b")
    return(invisible())
  }
  annotations <- report_network_annotations(graph$nodes)
  annotations$degree <- vapply(annotations$gene_id, function(id) {
    if (!nrow(graph$edges)) return(0L)
    sum(as.character(graph$edges$source) == id | as.character(graph$edges$target) == id, na.rm = TRUE)
  }, integer(1))
  evidence_order <- with(annotations, order(-as.integer(significant), -as.integer(wes),
    -as.integer(tf), -degree, label, na.last = TRUE))
  annotations <- annotations[evidence_order, , drop = FALSE]
  # Match the report's default direction-coloured app view: hide context
  # genes, but keep every direct drug target.
  target_nodes <- annotations[annotations$type == "Target", , drop = FALSE]
  target_ids <- target_nodes$gene_id
  drug_x <- 0.12
  drug_y <- 0.52
  node_radius <- max(0.014, min(0.035, 0.12 / sqrt(max(1, nrow(target_nodes)))))
  label_cex <- max(0.34, min(0.46, 0.95 / sqrt(max(1, nrow(target_nodes)))))
  # Use a deterministic force-directed layout so the PDF reads as a network,
  # rather than as a circular chart.  The grid is only the initial state;
  # STRING-linked targets are pulled together and close nodes are separated.
  network_positions <- function(ids, edges, radius) {
    n <- length(ids)
    if (!n) return(data.frame(gene_id = character(), x = numeric(), y = numeric()))
    if (n == 1) return(data.frame(gene_id = ids, x = 0.56, y = 0.52))
    columns <- max(3, ceiling(sqrt(n * 1.35)))
    rows <- ceiling(n / columns)
    grid_x <- seq(0.28, 0.88, length.out = columns)
    grid_y <- if (rows == 1) 0.52 else seq(0.25, 0.79, length.out = rows)
    index <- seq_len(n)
    positions <- cbind(
      grid_x[((index - 1) %% columns) + 1],
      grid_y[((index - 1) %/% columns) + 1]
    )
    colnames(positions) <- c("x", "y")
    adjacency <- matrix(FALSE, nrow = n, ncol = n)
    if (nrow(edges)) {
      from <- match(as.character(edges$source), ids)
      to <- match(as.character(edges$target), ids)
      keep <- which(!is.na(from) & !is.na(to) & from != to)
      if (length(keep)) {
        adjacency[cbind(from[keep], to[keep])] <- TRUE
        adjacency[cbind(to[keep], from[keep])] <- TRUE
      }
    }
    preferred <- max(0.12, min(0.18, 0.34 / sqrt(n) + 0.08))
    minimum <- max(0.032, radius * 2 + 0.012)
    for (iteration in seq_len(90)) {
      force <- matrix(0, nrow = n, ncol = 2)
      for (i in seq_len(n - 1)) {
        for (j in (i + 1):n) {
          delta <- positions[i, ] - positions[j, ]
          distance <- sqrt(sum(delta * delta))
          distance <- max(distance, 0.001)
          if (distance < preferred * 1.4) {
            push <- 0.0025 * (preferred * 1.4 - distance) / distance
            force[i, ] <- force[i, ] + delta * push
            force[j, ] <- force[j, ] - delta * push
          }
          if (adjacency[i, j]) {
            pull <- 0.0018 * (distance - preferred) / distance
            force[i, ] <- force[i, ] - delta * pull
            force[j, ] <- force[j, ] + delta * pull
          }
        }
      }
      force[, 1] <- force[, 1] + (0.58 - positions[, 1]) * 0.0008
      force[, 2] <- force[, 2] + (0.52 - positions[, 2]) * 0.0008
      positions <- positions + pmax(-0.025, pmin(0.025, force))
      positions[, 1] <- pmax(0.24, pmin(0.90, positions[, 1]))
      positions[, 2] <- pmax(0.22, pmin(0.82, positions[, 2]))
    }
    for (iteration in seq_len(40)) {
      for (i in seq_len(n - 1)) {
        for (j in (i + 1):n) {
          delta <- positions[i, ] - positions[j, ]
          distance <- sqrt(sum(delta * delta))
          if (distance > 0 && distance < minimum) {
            shift <- delta / distance * (minimum - distance) / 2
            positions[i, ] <- positions[i, ] + shift
            positions[j, ] <- positions[j, ] - shift
          }
        }
      }
      positions[, 1] <- pmax(0.24, pmin(0.90, positions[, 1]))
      positions[, 2] <- pmax(0.22, pmin(0.82, positions[, 2]))
    }
    data.frame(gene_id = ids, x = positions[, "x"], y = positions[, "y"], stringsAsFactors = FALSE)
  }
  target_edges <- if (nrow(graph$edges)) graph$edges else data.frame(source = character(), target = character())
  # Do not let DE significance determine the visual position of a node.  The
  # evidence ordering is useful for selecting and annotating nodes, but the
  # layout starts from degree/name order so significant genes are not placed
  # in one artificial band.
  layout_ids <- target_nodes$gene_id[order(-target_nodes$degree, target_nodes$label)]
  coords <- network_positions(layout_ids, target_edges, node_radius)
  if (nrow(coords) && nrow(graph$edges)) {
    for (i in seq_len(nrow(graph$edges))) {
      source_id <- as.character(graph$edges$source[i])
      target_id <- as.character(graph$edges$target[i])
      if (!(source_id %in% coords$gene_id && target_id %in% coords$gene_id)) next
      a <- coords[match(source_id, coords$gene_id), ]
      b <- coords[match(target_id, coords$gene_id), ]
      segments(a$x, a$y, b$x, b$y, col = "#cbd5e1", lwd = 0.8, lty = 2)
    }
  }
  if (length(target_ids)) {
    for (i in seq_along(target_ids)) {
      target <- coords[match(target_ids[i], coords$gene_id), ]
      segments(drug_x + 0.045, drug_y, target$x - node_radius, target$y,
        col = "#2563eb", lwd = 1.8)
    }
  }
  # Match the app's direction-coloured graph: charcoal diamond for the drug
  # and filled circles whose colours encode log2 fold-change.
  diamond_x <- drug_x + c(0, 0.045, 0, -0.045)
  diamond_y <- drug_y + c(0.055, 0, -0.055, 0)
  polygon(diamond_x, diamond_y, col = "#17202a", border = "#17202a")
  text(drug_x, drug_y, substr(drug_label, 1, 18), cex = 0.54, font = 2, col = "white")
  draw_node <- function(row, x, y) {
    fill <- if (is.na(row$log_fc) || !is.finite(row$log_fc)) "#94a3b8" else if (row$log_fc >= 0) "#dc2626" else "#2563eb"
    border_width <- if (isTRUE(row$significant)) 4 else 1
    # Keep the border in the same colour family as the node, but darken it so
    # the significant-node cartouche remains visible in a static PDF.
    border <- switch(fill,
      "#dc2626" = "#991b1b",
      "#2563eb" = "#1d4ed8",
      "#94a3b8" = "#64748b",
      fill
    )
    if (isTRUE(row$tf)) {
      triangle_angles <- pi / 2 + c(0, 2 * pi / 3, 4 * pi / 3)
      polygon(x + node_radius * cos(triangle_angles), y + node_radius * sin(triangle_angles),
        col = fill, border = border, lwd = border_width)
    } else {
      symbols(x, y, circles = node_radius, inches = FALSE, add = TRUE,
        bg = fill, fg = border, lwd = border_width)
    }
    if (isTRUE(row$wes)) {
      points(x, y, pch = 21, cex = 0.58, bg = "#be185d", col = "white", lwd = 1)
    }
    text(x, y - node_radius - 0.014, substr(row$label, 1, 14), cex = label_cex, col = "#17202a")
  }
  if (nrow(target_nodes)) for (i in seq_len(nrow(target_nodes))) {
    target <- coords[match(target_ids[i], coords$gene_id), ]
    draw_node(target_nodes[i, ], target$x, target$y)
  }
  legend("bottom", legend = c(report_tr("drug_target", language), report_tr("string_interaction", language),
    report_tr("de_up", language), report_tr("de_down", language), report_tr("expression_unknown", language),
    report_tr("de_significant", language), report_tr("wes_marker", language), report_tr("tf_target", language)),
    pch = c(NA, NA, 21, 21, 21, 21, 21, 24),
    pt.bg = c(NA, NA, "#dc2626", "#2563eb", "#94a3b8", "#94a3b8", "#be185d", "#94a3b8"),
    col = c("#2563eb", "#cbd5e1", "#dc2626", "#2563eb", "#94a3b8", "#17202a", "white", "#17202a"),
    lty = c(1, 2, 0, 0, 0, 0, 0, 0), lwd = c(2, 1, 1, 1, 1, 3, 1, 1),
    bty = "n", cex = 0.50, ncol = 4)
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
                                             language = "ru", report_title = "") {
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
                                    state_json_path = NULL, language = "ru", report_title = "") {
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
  report_text(0.06, 0.21, report_tr("disclaimer", language), width = 125,
    cex = 0.62, col = "#64748b", line_height = 0.024)

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
      indications = clinical$indications,
      warnings = clinical$warnings,
      stringsAsFactors = FALSE
    )
  }) |> bind_rows()
  if (nrow(clinical_rows)) {
    clinical_page_size <- 5L
    for (page_start in seq(1L, nrow(clinical_rows), by = clinical_page_size)) {
      page_end <- min(nrow(clinical_rows), page_start + clinical_page_size - 1L)
      page_rows <- clinical_rows[page_start:page_end, , drop = FALSE]
      report_canvas()
      report_header(
        if (page_start == 1L) report_tr("clinical", language) else report_tr("clinical_continued", language),
        if (page_start == 1L) report_tr("clinical_note", language) else NULL
      )
      report_draw_wrapped_table(
        c(report_tr("drug", language), report_tr("approval", language), report_tr("phase", language),
          report_tr("indications", language), report_tr("warnings", language)),
        page_rows,
        c(0.23, 0.14, 0.09, 0.35, 0.19),
        wrap_cols = c(4L, 5L),
        base_row_height = 0.052,
        font_size = 0.46,
        header_size = 0.58,
        wrap_widths = c(26, 12, 8, 39, 16)
      )
    }
  }

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
