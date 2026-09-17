library(shiny)
library(DT)
library(jsonlite)
library(dplyr)
library(purrr)
library(stringr)
library(htmltools)


# ============================================================
# 1. PATHS
# ============================================================

root_dir <- "C:/projects/clinreport/raw/clinreport"

sources <- c(
  "chembl",
  "opentargets",
  "pubchem"
)


# ============================================================
# 2. FIND DRUG DIRECTORIES
# ============================================================

find_drug_dirs <- function(source) {
  
  source_dir <- file.path(
    root_dir,
    source
  )
  
  dirs <- list.dirs(
    source_dir,
    recursive = FALSE,
    full.names = TRUE
  )
  
  if (length(dirs) == 0) {
    return(tibble())
  }
  
  tibble(
    source = source,
    path = dirs,
    folder = basename(dirs),
    drug_id = str_extract(
      basename(dirs),
      "DB[0-9]+$"
    )
  ) |>
    filter(!is.na(drug_id))
}


drug_dirs <- map_dfr(
  sources,
  find_drug_dirs
)


# ============================================================
# 3. READ JSON
# ============================================================

read_json_safe <- function(path) {
  
  if (!file.exists(path)) {
    return(NULL)
  }
  
  tryCatch(
    
    fromJSON(
      path,
      simplifyVector = FALSE
    ),
    
    error = function(e) {
      
      message(
        "JSON ERROR: ",
        path,
        "\n",
        e$message
      )
      
      NULL
    }
  )
}


# ============================================================
# 4. READ ALL JSON FOR ONE SOURCE / DRUG
# ============================================================

read_drug_source <- function(
    source,
    drug_id
) {
  
  x <- drug_dirs |>
    filter(
      source == !!source,
      drug_id == !!drug_id
    )
  
  if (nrow(x) == 0) {
    return(list())
  }
  
  json_files <- list.files(
    x$path[1],
    pattern = "\\.json$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(json_files) == 0) {
    return(list())
  }
  
  result <- list()
  
  for (file in json_files) {
    
    object_name <- tools::file_path_sans_ext(
      basename(file)
    )
    
    result[[object_name]] <- read_json_safe(
      file
    )
  }
  
  result
}


# ============================================================
# 5. LOAD ONE DRUG
# ============================================================

load_drug <- function(drug_id) {
  
  list(
    
    chembl = read_drug_source(
      "chembl",
      drug_id
    ),
    
    opentargets = read_drug_source(
      "opentargets",
      drug_id
    ),
    
    pubchem = read_drug_source(
      "pubchem",
      drug_id
    )
  )
}


# ============================================================
# 6. BASIC DRUG INFORMATION
# ============================================================

get_drug_name <- function(drug_id) {
  
  chembl <- read_drug_source(
    "chembl",
    drug_id
  )
  
  if ("chembl_id" %in% names(chembl)) {
    
    x <- chembl$chembl_id
    
    if (!is.null(x$drug_name)) {
      
      return(
        as.character(
          x$drug_name[[1]]
        )
      )
    }
  }
  
  
  if ("chembl_molecule" %in% names(chembl)) {
    
    x <- chembl$chembl_molecule
    
    if (!is.null(x$pref_name)) {
      
      return(
        as.character(
          x$pref_name[[1]]
        )
      )
    }
  }
  
  
  opentargets <- read_drug_source(
    "opentargets",
    drug_id
  )
  
  if ("opentargets_id" %in% names(opentargets)) {
    
    x <- opentargets$opentargets_id
    
    if (!is.null(x$drug_name)) {
      
      return(
        as.character(
          x$drug_name[[1]]
        )
      )
    }
  }
  
  
  drug_id
}


get_chembl_id <- function(drug_id) {
  
  chembl <- read_drug_source(
    "chembl",
    drug_id
  )
  
  if ("chembl_id" %in% names(chembl)) {
    
    x <- chembl$chembl_id
    
    if (!is.null(x$chembl_id)) {
      
      return(
        as.character(
          x$chembl_id[[1]]
        )
      )
    }
  }
  
  NA_character_
}


get_pubchem_cid <- function(drug_id) {
  
  pubchem <- read_drug_source(
    "pubchem",
    drug_id
  )
  
  if ("pubchem_cids" %in% names(pubchem)) {
    
    x <- pubchem$pubchem_cids
    
    if (
      !is.null(x$IdentifierList) &&
      !is.null(x$IdentifierList$CID)
    ) {
      
      return(
        as.character(
          x$IdentifierList$CID[[1]]
        )
      )
    }
  }
  
  NA_character_
}


# ============================================================
# 7. BUILD DRUG TABLE
# ============================================================

drug_ids <- sort(
  unique(drug_dirs$drug_id)
)


drug_table <- tibble(
  drug_id = drug_ids
) |>
  mutate(
    
    drug_name = map_chr(
      drug_id,
      get_drug_name
    ),
    
    chembl_id = map_chr(
      drug_id,
      get_chembl_id
    ),
    
    pubchem_cid = map_chr(
      drug_id,
      get_pubchem_cid
    )
  )


# ============================================================
# 8. HELPERS
# ============================================================

pretty_name <- function(x) {
  
  x |>
    str_replace_all(
      "_",
      " "
    ) |>
    str_replace_all(
      "(?i)chembl",
      "ChEMBL"
    ) |>
    str_replace_all(
      "(?i)pubchem",
      "PubChem"
    ) |>
    str_replace_all(
      "(?i)opentargets",
      "Open Targets"
    ) |>
    str_replace_all(
      "(?i)inchi",
      "InChI"
    ) |>
    str_replace_all(
      "(?i)smiles",
      "SMILES"
    ) |>
    str_to_sentence()
}


is_empty_value <- function(x) {
  
  if (is.null(x)) {
    return(TRUE)
  }
  
  if (length(x) == 0) {
    return(TRUE)
  }
  
  if (!is.list(x) && all(is.na(x))) {
    return(TRUE)
  }
  
  FALSE
}


# ============================================================
# 9. CONVERT JSON VALUE TO DISPLAYABLE TEXT
# ============================================================

value_to_text <- function(x) {
  
  if (is.null(x)) {
    return(NULL)
  }
  
  
  # atomic value
  
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
  
  
  # list of atomic values
  
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


# ============================================================
# 10. RENDER A SIMPLE FIELD
# ============================================================

field_row <- function(
    label,
    value
) {
  
  if (is.null(value)) {
    return(NULL)
  }
  
  value <- value_to_text(value)
  
  if (
    is.null(value) ||
    !nzchar(trimws(value))
  ) {
    return(NULL)
  }
  
  div(
    
    class = "field-row",
    
    div(
      class = "field-label",
      label
    ),
    
    div(
      class = "field-value",
      htmlEscape(value)
    )
  )
}


# ============================================================
# 11. GET VALUE FROM NESTED OBJECT
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


# ============================================================
# 12. CHOOSE FIRST AVAILABLE VALUE
# ============================================================

first_value <- function(...) {
  
  values <- list(...)
  
  for (x in values) {
    
    if (!is_empty_value(x)) {
      
      txt <- value_to_text(x)
      
      if (
        !is.null(txt) &&
        nzchar(trimws(txt))
      ) {
        return(txt)
      }
    }
  }
  
  NULL
}


# ============================================================
# 13. GENERAL INFORMATION
# ============================================================

render_general <- function(
    drug_id,
    data
) {
  
  chembl_id <- data$chembl$chembl_id
  
  molecule <- data$chembl$chembl_molecule
  
  pubchem_properties <-
    data$pubchem$pubchem_properties
  
  properties <- NULL
  
  if (
    !is.null(pubchem_properties) &&
    !is.null(pubchem_properties$PropertyTable) &&
    !is.null(pubchem_properties$PropertyTable$Properties)
  ) {
    
    properties <-
      pubchem_properties$PropertyTable$Properties[[1]]
  }
  
  
  tagList(
    
    div(
      class = "section-title",
      "Identifiers"
    ),
    
    field_row(
      "DrugBank ID",
      drug_id
    ),
    
    field_row(
      "ChEMBL ID",
      get_nested(
        chembl_id,
        "chembl_id"
      )
    ),
    
    field_row(
      "PubChem CID",
      get_nested(
        data$pubchem$pubchem_cids,
        "IdentifierList",
        "CID"
      )
    ),
    
    field_row(
      "Preferred name",
      get_nested(
        molecule,
        "pref_name"
      )
    ),
    
    
    div(
      class = "section-title",
      "Molecular properties"
    ),
    
    field_row(
      "Molecular formula",
      get_nested(
        molecule,
        "molecule_properties",
        "full_molformula"
      )
    ),
    
    field_row(
      "Molecular weight",
      get_nested(
        molecule,
        "molecule_properties",
        "full_mwt"
      )
    ),
    
    field_row(
      "Freebase molecular weight",
      get_nested(
        molecule,
        "molecule_properties",
        "mw_freebase"
      )
    ),
    
    field_row(
      "AlogP",
      get_nested(
        molecule,
        "molecule_properties",
        "alogp"
      )
    ),
    
    field_row(
      "H-bond acceptors",
      get_nested(
        molecule,
        "molecule_properties",
        "hba"
      )
    ),
    
    field_row(
      "H-bond donors",
      get_nested(
        molecule,
        "molecule_properties",
        "hbd"
      )
    ),
    
    field_row(
      "Polar surface area",
      get_nested(
        molecule,
        "molecule_properties",
        "psa"
      )
    ),
    
    field_row(
      "Rotatable bonds",
      get_nested(
        molecule,
        "molecule_properties",
        "rtb"
      )
    ),
    
    field_row(
      "Aromatic rings",
      get_nested(
        molecule,
        "molecule_properties",
        "aromatic_rings"
      )
    ),
    
    
    div(
      class = "section-title",
      "Development"
    ),
    
    field_row(
      "First approval",
      get_nested(
        molecule,
        "first_approval"
      )
    ),
    
    field_row(
      "Maximum phase",
      get_nested(
        molecule,
        "max_phase"
      )
    ),
    
    field_row(
      "First in class",
      get_nested(
        molecule,
        "first_in_class"
      )
    ),
    
    
    div(
      class = "section-title",
      "Structures"
    ),
    
    field_row(
      "Canonical SMILES",
      get_nested(
        molecule,
        "molecule_structures",
        "canonical_smiles"
      )
    ),
    
    field_row(
      "Standard InChI",
      get_nested(
        molecule,
        "molecule_structures",
        "standard_inchi"
      )
    ),
    
    field_row(
      "InChI Key",
      get_nested(
        molecule,
        "molecule_structures",
        "standard_inchi_key"
      )
    )
  )
}


# ============================================================
# 14. PHARMACOLOGY
# ============================================================

render_pharmacology <- function(
    data
) {
  
  indications <- data$chembl$chembl_indications
  
  mechanisms <- data$chembl$chembl_mechanisms
  
  opentargets <- data$opentargets
  
  
  tagList(
    
    div(
      class = "section-title",
      "Indications"
    ),
    
    if (
      !is.null(indications) &&
      !is.null(indications$drug_indications)
    ) {
      
      map(
        indications$drug_indications,
        function(x) {
          
          tagList(
            
            field_row(
              "EFO term",
              x$efo_term
            ),
            
            field_row(
              "EFO ID",
              x$efo_id
            ),
            
            field_row(
              "Mesh heading",
              x$mesh_heading
            ),
            
            field_row(
              "Max phase",
              x$max_phase_for_ind
            )
          )
        }
      )
    },
    
    
    div(
      class = "section-title",
      "Mechanisms of action"
    ),
    
    if (
      !is.null(mechanisms) &&
      !is.null(mechanisms$mechanisms)
    ) {
      
      map(
        mechanisms$mechanisms,
        function(x) {
          
          tagList(
            
            field_row(
              "Mechanism",
              x$mechanism_of_action
            ),
            
            field_row(
              "Action type",
              x$action_type
            ),
            
            field_row(
              "Molecular mechanism",
              x$molecular_mechanism
            ),
            
            field_row(
              "Target",
              x$target_chembl_id
            ),
            
            field_row(
              "Direct interaction",
              x$direct_interaction
            )
          )
        }
      )
    },
    
    
    div(
      class = "section-title",
      "Open Targets"
    ),
    
    if (
      length(opentargets) > 0
    ) {
      
      map(
        opentargets,
        function(x) {
          
          if (
            !is.null(x$search_data) &&
            !is.null(x$search_data$search_hits)
          ) {
            
            map(
              x$search_data$search_hits,
              function(hit) {
                
                tagList(
                  
                  field_row(
                    "Target",
                    hit$id
                  ),
                  
                  field_row(
                    "Name",
                    hit$name
                  ),
                  
                  field_row(
                    "Entity",
                    hit$entity
                  ),
                  
                  field_row(
                    "Description",
                    hit$description
                  )
                )
              }
            )
          }
        }
      )
    }
  )
}


# ============================================================
# 15. SAFETY / SOURCES
# ============================================================

render_safety <- function(
    data
) {
  
  warnings <- data$chembl$chembl_warnings
  
  molecule <- data$chembl$chembl_molecule
  
  
  tagList(
    
    div(
      class = "section-title",
      "Warnings"
    ),
    
    if (
      !is.null(warnings) &&
      !is.null(warnings$drug_warnings)
    ) {
      
      map(
        warnings$drug_warnings,
        function(x) {
          
          tagList(
            
            field_row(
              "Warning class",
              x$warning_class
            ),
            
            field_row(
              "Warning type",
              x$warning_type
            ),
            
            field_row(
              "Country",
              x$warning_country
            ),
            
            field_row(
              "EFO ID",
              x$efo_id_for_warning_class
            )
          )
        }
      )
    },
    
    
    div(
      class = "section-title",
      "Synonyms"
    ),
    
    if (
      !is.null(molecule) &&
      !is.null(molecule$molecule_synonyms)
    ) {
      
      map(
        molecule$molecule_synonyms,
        function(x) {
          
          field_row(
            "Synonym",
            x$molecule_synonym
          )
        }
      )
    },
    
    
    div(
      class = "section-title",
      "PubChem references"
    ),
    
    if (
      !is.null(data$pubchem$pubchem_pugview) &&
      !is.null(
        data$pubchem$pubchem_pugview$Record
      )
    ) {
      
      record <-
        data$pubchem$pubchem_pugview$Record
      
      field_row(
        "Record title",
        record$RecordTitle
      )
      
      field_row(
        "Record number",
        record$RecordNumber
      )
    }
  )
}


# ============================================================
# 16. UI
# ============================================================

ui <- fluidPage(
  
  tags$head(
    
    tags$style(
      HTML(
        "
        body {
          background: #f4f5f7;
          font-family: -apple-system, BlinkMacSystemFont,
                       'Segoe UI', sans-serif;
        }

        .app-title {
          font-size: 28px;
          font-weight: 600;
          margin: 20px 0;
        }

        .main-layout {
          display: flex;
          gap: 20px;
          align-items: stretch;
        }

        /* LEFT */

        .drug-list {
          width: 43%;
          background: white;
          border-radius: 10px;
          padding: 15px;
          box-shadow: 0 2px 8px rgba(0,0,0,0.08);
        }

        /* RIGHT */

        .drug-details {
          flex: 1;
          background: white;
          border-radius: 10px;
          padding: 25px;
          min-height: 650px;
          box-shadow: 0 2px 8px rgba(0,0,0,0.08);
        }

        .drug-header {
          margin-bottom: 15px;
        }

        .drug-title {
          font-size: 28px;
          font-weight: 600;
          line-height: 1.2;
        }

        .drug-id {
          color: #777;
          font-size: 14px;
          margin-top: 4px;
        }

        .section-title {
          font-size: 18px;
          font-weight: 600;
          margin-top: 25px;
          margin-bottom: 10px;
          padding-bottom: 7px;
          border-bottom: 1px solid #ddd;
        }

        .field-row {
          display: flex;
          padding: 8px 0;
          border-bottom: 1px solid #f0f0f0;
        }

        .field-label {
          width: 32%;
          min-width: 170px;
          color: #666;
          font-weight: 500;
        }

        .field-value {
          width: 68%;
          word-break: break-word;
          white-space: pre-wrap;
        }

        .empty-message {
          color: #999;
          padding: 40px;
          text-align: center;
        }

        .nav-tabs {
          margin-top: 20px;
        }

        .tab-content {
          padding-top: 5px;
        }

        table.dataTable tbody tr.selected {
          background-color: #e8f0fe !important;
        }

        @media (max-width: 1000px) {

          .main-layout {
            flex-direction: column;
          }

          .drug-list {
            width: 100%;
          }

          .drug-details {
            width: 100%;
          }
        }
        "
      )
    )
  ),
  
  
  div(
    class = "app-title",
    "ClinReport — Drug Reference"
  ),
  
  
  div(
    class = "main-layout",
    
    div(
      class = "drug-list",
      
      DTOutput(
        "drug_table"
      )
    ),
    
    
    div(
      class = "drug-details",
      
      uiOutput(
        "drug_details"
      )
    )
  )
)


# ============================================================
# 17. SERVER
# ============================================================

server <- function(
    input,
    output,
    session
) {
  
  
  # ----------------------------------------------------------
  # TABLE
  # ----------------------------------------------------------
  
  output$drug_table <- renderDT({
    
    datatable(
      
      drug_table |>
        select(
          drug_id,
          drug_name,
          chembl_id,
          pubchem_cid
        ),
      
      selection = "single",
      
      rownames = FALSE,
      
      filter = "top",
      
      colnames = c(
        "Drug ID",
        "Drug name",
        "ChEMBL ID",
        "PubChem CID"
      ),
      
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        autoWidth = TRUE
      )
    )
  })
  
  
  # ----------------------------------------------------------
  # SELECTED DRUG
  # ----------------------------------------------------------
  
  selected_drug <- reactive({
    
    req(
      input$drug_table_rows_selected
    )
    
    drug_table[
      input$drug_table_rows_selected,
      ,
      drop = FALSE
    ]
  })
  
  
  # ----------------------------------------------------------
  # DETAILS
  # ----------------------------------------------------------
  
  output$drug_details <- renderUI({
    
    req(
      selected_drug()
    )
    
    drug <- selected_drug()
    
    drug_id <- drug$drug_id
    
    data <- load_drug(
      drug_id
    )
    
    
    tagList(
      
      div(
        class = "drug-header",
        
        div(
          class = "drug-title",
          drug$drug_name
        ),
        
        div(
          class = "drug-id",
          drug_id
        )
      ),
      
      
      tabsetPanel(
        
        tabPanel(
          "Основная информация",
          
          render_general(
            drug_id,
            data
          )
        ),
        
        
        tabPanel(
          "Показания и механизм",
          
          render_pharmacology(
            data
          )
        ),
        
        
        tabPanel(
          "Безопасность и источники",
          
          render_safety(
            data
          )
        )
      )
    )
  })
}


# ============================================================
# 18. RUN
# ============================================================

shinyApp(
  ui = ui,
  server = server
)