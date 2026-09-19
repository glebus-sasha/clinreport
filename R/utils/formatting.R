# Formatting helpers

safe_text <- function(x, fallback = "—") {
  
  value <- value_to_text(x)
  
  if (
    is.null(value) ||
    !nzchar(trimws(value))
  ) {
    return(fallback)
  }
  
  value
}


field_row <- function(
    label,
    value,
    mono = FALSE,
    link = NULL
) {
  
  value <- safe_text(value)
  
  value_ui <- if (!is.null(link)) {
    
    tags$a(
      href = link,
      target = "_blank",
      rel = "noopener noreferrer",
      class = "external-link",
      value
    )
    
  } else if (mono) {
    
    tags$span(
      class = "mono-value",
      value
    )
    
  } else {
    
    htmlEscape(value)
  }
  
  div(
    class = "field-row",
    
    div(
      class = "field-label",
      label
    ),
    
    div(
      class = "field-value",
      value_ui
    )
  )
}


