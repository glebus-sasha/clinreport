# Safety tab renderer

# 17. SAFETY
# ============================================================

render_safety <- function(
    data
) {
  
  warnings <- data$chembl$chembl_warnings
  molecule <- data$chembl$chembl_molecule
  
  
  warning_cards <- list()
  
  if (
    !is.null(warnings) &&
    !is.null(warnings$drug_warnings)
  ) {
    
    warning_cards <- map(
      warnings$drug_warnings,
      function(x) {
        
        div(
          class = "warning-card",
          
          div(
            class = "warning-icon",
            "!"
          ),
          
          div(
            class = "warning-content",
            
            div(
              class = "warning-class",
              safe_text(
                x$warning_class,
                "Safety warning"
              )
            ),
            
            div(
              class = "warning-meta",
              
              div(
                class = "warning-type",
                safe_text(
                  x$warning_type,
                  "Unknown warning type"
                )
              ),
              
              div(
                class = "warning-country",
                safe_text(
                  x$warning_country,
                  "Unknown country"
                )
              )
            )
          )
        )
      }
    )
  }
  
  
  synonyms <- character()
  
  if (
    !is.null(molecule) &&
    !is.null(molecule$molecule_synonyms)
  ) {
    
    synonyms <- map_chr(
      molecule$molecule_synonyms,
      function(x) {
        safe_text(
          x$molecule_synonym,
          ""
        )
      }
    )
    
    synonyms <- synonyms[
      nzchar(synonyms)
    ]
    
    synonyms <- unique(
      synonyms
    )
  }
  
  
  synonym_chips <- map(
    synonyms,
    function(x) {
      
      span(
        class = "synonym-chip",
        x
      )
    }
  )
  
  
  tagList(
    
    div(
      class = "section-header",
      
      div(
        class = "section-title",
        "Safety warnings"
      ),
      
      div(
        class = "section-subtitle",
        "Regulatory and safety-related records"
      )
    ),
    
    if (length(warning_cards) > 0) {
      
      div(
        class = "warning-list",
        warning_cards
      )
      
    } else {
      
      div(
        class = "empty-message",
        "No warning records available."
      )
    },
    
    
    div(
      class = "section-header safety-section",
      
      div(
        class = "section-title",
        "Synonyms"
      ),
      
      div(
        class = "section-subtitle",
        paste(
          length(synonyms),
          "unique names"
        )
      )
    ),
    
    if (length(synonym_chips) > 0) {
      
      div(
        class = "synonym-cloud",
        synonym_chips
      )
      
    } else {
      
      div(
        class = "empty-message",
        "No synonyms available."
      )
    }
  )
}


# ============================================================
