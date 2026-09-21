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

