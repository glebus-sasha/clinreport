# Pharmacology tab renderer

# 16. PHARMACOLOGY
# ============================================================

render_pharmacology <- function(
    data
) {
  
  indications <- data$chembl$chembl_indications
  mechanisms <- data$chembl$chembl_mechanisms
  
  
  indication_cards <- list()
  
  if (
    !is.null(indications) &&
    !is.null(indications$drug_indications)
  ) {
    
    indication_cards <- map(
      indications$drug_indications,
      function(x) {
        
        phase <- safe_text(
          x$max_phase_for_ind,
          ""
        )
        
        div(
          class = "indication-card",
          
          div(
            class = "indication-top",
            
            div(
              class = "indication-name",
              safe_text(
                x$efo_term,
                "Unknown indication"
              )
            ),
            
            if (nzchar(phase)) {
              
              div(
                class = "phase-badge",
                paste0(
                  "Phase ",
                  phase
                )
              )
            }
          ),
          
          div(
            class = "indication-meta",
            
            div(
              class = "meta-item",
              
              span(
                class = "meta-label",
                "EFO"
              ),
              
              span(
                class = "meta-value mono-value",
                safe_text(
                  x$efo_id
                )
              )
            ),
            
            div(
              class = "meta-item",
              
              span(
                class = "meta-label",
                "MeSH"
              ),
              
              span(
                class = "meta-value",
                safe_text(
                  x$mesh_heading
                )
              )
            )
          )
        )
      }
    )
  }
  
  
  mechanism_cards <- list()
  
  if (
    !is.null(mechanisms) &&
    !is.null(mechanisms$mechanisms)
  ) {
    
    mechanism_cards <- map(
      mechanisms$mechanisms,
      function(x) {
        
        target_id <- safe_text(
          x$target_chembl_id,
          ""
        )
        
        target_link <- if (
          nzchar(target_id)
        ) {
          chembl_url(target_id)
        } else {
          NULL
        }
        
        div(
          class = "mechanism-card",
          
          div(
            class = "mechanism-main",
            
            div(
              class = "mechanism-title",
              safe_text(
                x$mechanism_of_action,
                "Unknown mechanism"
              )
            ),
            
            div(
              class = "mechanism-target",
              
              span(
                class = "target-label",
                "Target"
              ),
              
              if (!is.null(target_link)) {
                
                tags$a(
                  href = target_link,
                  target = "_blank",
                  rel = "noopener noreferrer",
                  class = "target-link",
                  target_id
                )
                
              } else {
                target_id
              }
            )
          ),
          
          div(
            class = "mechanism-meta",
            
            div(
              class = "mechanism-pill",
              safe_text(
                x$action_type,
                "Unknown"
              )
            ),
            
            div(
              class = "mechanism-detail",
              
              span(
                class = "mechanism-detail-label",
                "Direct interaction"
              ),
              
              span(
                class = "mechanism-detail-value",
                safe_text(
                  x$direct_interaction
                )
              )
            )
          )
        )
      }
    )
  }
  
  
  tagList(
    
    div(
      class = "section-header",
      
      div(
        class = "section-title",
        "Indications"
      ),
      
      div(
        class = "section-subtitle",
        "Disease associations reported by ChEMBL"
      )
    ),
    
    if (length(indication_cards) > 0) {
      
      div(
        class = "indications-grid",
        indication_cards
      )
      
    } else {
      
      div(
        class = "empty-message",
        "No indication records available."
      )
    },
    
    
    div(
      class = "section-header pharmacology-section",
      
      div(
        class = "section-title",
        "Mechanisms of action"
      ),
      
      div(
        class = "section-subtitle",
        "Target-level pharmacological relationships"
      )
    ),
    
    if (length(mechanism_cards) > 0) {
      
      div(
        class = "mechanism-list",
        mechanism_cards
      )
      
    } else {
      
      div(
        class = "empty-message",
        "No mechanism records available."
      )
    }
  )
}


