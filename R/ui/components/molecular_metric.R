# Molecular profile metric component

# ============================================================
# 13. MOLECULAR PROFILE
# ============================================================

molecular_metric <- function(
    label,
    value,
    max_value,
    suffix = ""
) {
  
  numeric_value <- suppressWarnings(
    as.numeric(value)
  )
  
  pct <- if (
    !is.na(numeric_value) &&
    max_value > 0
  ) {
    min(
      100,
      max(
        0,
        numeric_value / max_value * 100
      )
    )
  } else {
    0
  }
  
  div(
    class = "molecular-metric",
    
    div(
      class = "metric-top",
      
      span(
        class = "metric-label",
        label
      ),
      
      span(
        class = "metric-value",
        ifelse(
          is.na(numeric_value),
          "—",
          paste0(
            format(
              numeric_value,
              trim = TRUE,
              scientific = FALSE
            ),
            suffix
          )
        )
      )
    ),
    
    div(
      class = "metric-track",
      
      div(
        class = "metric-fill",
        style = paste0(
          "width:",
          pct,
          "%"
        )
      )
    )
  )
}


# ============================================================
