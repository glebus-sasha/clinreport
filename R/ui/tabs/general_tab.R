# General tab renderer

# ============================================================
# 15. GENERAL TAB
# ============================================================

render_general <- function(
    drug_id,
    data
) {
  
  chembl_id <- get_chembl_id(data)
  pubchem_cid <- get_pubchem_cid(data)
  
  molecule <- data$chembl$chembl_molecule
  
  molecular_weight <- suppressWarnings(
    as.numeric(
      safe_text(
        get_nested(
          molecule,
          "molecule_properties",
          "full_mwt"
        ),
        NA
      )
    )
  )
  
  alogp <- suppressWarnings(
    as.numeric(
      safe_text(
        get_nested(
          molecule,
          "molecule_properties",
          "alogp"
        ),
        NA
      )
    )
  )
  
  hba <- suppressWarnings(
    as.numeric(
      safe_text(
        get_nested(
          molecule,
          "molecule_properties",
          "hba"
        ),
        NA
      )
    )
  )
  
  hbd <- suppressWarnings(
    as.numeric(
      safe_text(
        get_nested(
          molecule,
          "molecule_properties",
          "hbd"
        ),
        NA
      )
    )
  )
  
  psa <- suppressWarnings(
    as.numeric(
      safe_text(
        get_nested(
          molecule,
          "molecule_properties",
          "psa"
        ),
        NA
      )
    )
  )
  
  rtb <- suppressWarnings(
    as.numeric(
      safe_text(
        get_nested(
          molecule,
          "molecule_properties",
          "rtb"
        ),
        NA
      )
    )
  )
  
  
  tagList(
    
    # --------------------------------------------------------
    # IDENTIFIERS
    # --------------------------------------------------------
    
    div(
      class = "section-header",
      
      div(
        class = "section-title",
        "Identifiers"
      ),
      
      div(
        class = "section-subtitle",
        "Cross-references to external resources"
      )
    ),
    
    div(
      class = "identifier-grid",
      
      div(
        class = "identifier-card",
        
        div(
          class = "identifier-label",
          "DrugBank ID"
        ),
        
        div(
          class = "identifier-value mono-value",
          drug_id
        )
      ),
      
      div(
        class = "identifier-card",
        
        div(
          class = "identifier-label",
          "ChEMBL ID"
        ),
        
        div(
          class = "identifier-value",
          
          if (!is.null(chembl_id)) {
            
            tags$a(
              href = chembl_url(chembl_id),
              target = "_blank",
              rel = "noopener noreferrer",
              class = "identifier-link",
              chembl_id
            )
            
          } else {
            "—"
          }
        )
      ),
      
      div(
        class = "identifier-card",
        
        div(
          class = "identifier-label",
          "PubChem CID"
        ),
        
        div(
          class = "identifier-value",
          
          if (!is.null(pubchem_cid)) {
            
            tags$a(
              href = pubchem_url(pubchem_cid),
              target = "_blank",
              rel = "noopener noreferrer",
              class = "identifier-link",
              pubchem_cid
            )
            
          } else {
            "—"
          }
        )
      ),
      
      div(
        class = "identifier-card wide",
        
        div(
          class = "identifier-label",
          "Preferred name"
        ),
        
        div(
          class = "identifier-value preferred-name",
          safe_text(
            get_nested(
              molecule,
              "pref_name"
            )
          )
        )
      )
    ),
    
    source_buttons(
      drug_id,
      chembl_id,
      pubchem_cid
    ),
    
    
    # --------------------------------------------------------
    # MOLECULAR PROPERTIES
    # --------------------------------------------------------
    
    div(
      class = "section-header",
      
      div(
        class = "section-title",
        "Molecular properties"
      ),
      
      div(
        class = "section-subtitle",
        "Compact physicochemical profile"
      )
    ),
    
    div(
      class = "molecular-grid",
      
      molecular_metric(
        "Molecular weight",
        molecular_weight,
        800,
        " Da"
      ),
      
      molecular_metric(
        "AlogP",
        alogp,
        8
      ),
      
      molecular_metric(
        "H-bond acceptors",
        hba,
        12
      ),
      
      molecular_metric(
        "H-bond donors",
        hbd,
        8
      ),
      
      molecular_metric(
        "Polar surface area",
        psa,
        160,
        " Å²"
      ),
      
      molecular_metric(
        "Rotatable bonds",
        rtb,
        12
      )
    ),
    
    div(
      class = "formula-card",
      
      div(
        class = "formula-label",
        "Molecular formula"
      ),
      
      div(
        class = "formula-value mono-value",
        safe_text(
          get_nested(
            molecule,
            "molecule_properties",
            "full_molformula"
          )
        )
      )
    ),
    
    
    # --------------------------------------------------------
    # DEVELOPMENT
    # --------------------------------------------------------
    
    div(
      class = "section-header",
      
      div(
        class = "section-title",
        "Development"
      )
    ),
    
    div(
      class = "development-grid",
      
      div(
        class = "development-card",
        
        div(
          class = "development-label",
          "First approval"
        ),
        
        div(
          class = "development-value",
          safe_text(
            get_nested(
              molecule,
              "first_approval"
            )
          )
        )
      ),
      
      div(
        class = "development-card",
        
        div(
          class = "development-label",
          "Maximum phase"
        ),
        
        div(
          class = "phase-badge large",
          paste0(
            "Phase ",
            safe_text(
              get_nested(
                molecule,
                "max_phase"
              )
            )
          )
        )
      )
    ),
    
    
    # --------------------------------------------------------
    # STRUCTURES
    # --------------------------------------------------------
    
    div(
      class = "section-header",
      
      div(
        class = "section-title",
        "Structures"
      ),
      
      div(
        class = "section-subtitle",
        "2D representation and machine-readable identifiers"
      )
    ),
    
    div(
      class = "structure-layout",
      
      div(
        class = "structure-visual",
        
        render_structure(
          pubchem_cid
        )
      ),
      
      div(
        class = "structure-data",
        
        div(
          class = "structure-field",
          
          div(
            class = "structure-field-label",
            "Canonical SMILES"
          ),
          
          div(
            class = "structure-code",
            safe_text(
              get_nested(
                molecule,
                "molecule_structures",
                "canonical_smiles"
              )
            )
          )
        ),
        
        div(
          class = "structure-field",
          
          div(
            class = "structure-field-label",
            "Standard InChI"
          ),
          
          div(
            class = "structure-code",
            safe_text(
              get_nested(
                molecule,
                "molecule_structures",
                "standard_inchi"
              )
            )
          )
        ),
        
        div(
          class = "structure-field",
          
          div(
            class = "structure-field-label",
            "InChI Key"
          ),
          
          div(
            class = "structure-code",
            safe_text(
              get_nested(
                molecule,
                "molecule_structures",
                "standard_inchi_key"
              )
            )
          )
        )
      )
    )
  )
}


