# Generic nested-value helpers

# 7. HELPERS
# ============================================================

get_nested <- function(
    x,
    ...
) {
  
  path <- list(...)
  
  value <- x
  
  for (p in path) {
    
    if (
      !is.list(value) ||
      is.null(value[[p]])
    ) {
      return(NULL)
    }
    
    value <- value[[p]]
  }
  
  value
}


value_to_text <- function(x) {
  
  if (is.null(x)) {
    return(NULL)
  }
  
  if (!is.list(x)) {
    
    x <- as.character(x)
    
    x <- x[
      !is.na(x)
    ]
    
    if (length(x) == 0) {
      return(NULL)
    }
    
    return(
      paste(
        x,
        collapse = "; "
      )
    )
  }
  
  if (
    length(x) > 0 &&
    all(
      map_lgl(
        x,
        ~ !is.list(.x)
      )
    )
  ) {
    
    values <- map_chr(
      x,
      ~ as.character(.x)
    )
    
    values <- values[
      !is.na(values)
    ]
    
    if (length(values) == 0) {
      return(NULL)
    }
    
    return(
      paste(
        values,
        collapse = "; "
      )
    )
  }
  
  NULL
}


