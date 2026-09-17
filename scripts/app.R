library(shiny)
library(DT)
library(jsonlite)
library(dplyr)
library(purrr)
library(stringr)
library(htmltools)
library(plotly)


# ============================================================
# 1. INPUT PARAMETERS
# ============================================================

drug_network_file <- file.path(
  "C:/projects/clinreport/raw",
  "all_samples_string_human_links_v12_0_min900_Ensembl_diamond_trustrank.csv"
)

clinreport_dir <- file.path(
  "C:/projects/clinreport/raw",
  "clinreport"
)

gene_file <- file.path(
  "C:/projects/clinreport/raw",
  "carcinoma_vs_normal_gene_names_added.tsv"
)

padj_cutoff <- 0.05
log2fc_cutoff <- 1.0


# ============================================================
# 2. VALIDATION
# ============================================================

if (!file.exists(drug_network_file)) {
  stop("Drug network file does not exist: ", drug_network_file)
}

if (!dir.exists(clinreport_dir)) {
  stop("ClinReport directory does not exist: ", clinreport_dir)
}

if (!file.exists(gene_file)) {
  stop("Gene expression file does not exist: ", gene_file)
}

if (
  !is.numeric(padj_cutoff) ||
  length(padj_cutoff) != 1 ||
  is.na(padj_cutoff) ||
  padj_cutoff <= 0 ||
  padj_cutoff >= 1
) {
  stop("padj_cutoff must be a single number between 0 and 1.")
}

if (
  !is.numeric(log2fc_cutoff) ||
  length(log2fc_cutoff) != 1 ||
  is.na(log2fc_cutoff) ||
  log2fc_cutoff < 0
) {
  stop("log2fc_cutoff must be a single non-negative number.")
}


# ============================================================
# 3. LOAD DRUG NETWORK
# ============================================================

drug_network <- read.csv(
  drug_network_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (
  ncol(drug_network) > 0 &&
  names(drug_network)[1] == ""
) {
  drug_network <- drug_network[, -1, drop = FALSE]
}

required_drug_columns <- c(
  "drugId",
  "label",
  "status",
  "drugstoneType",
  "score",
  "hasEdgesTo",
  "isResult",
  "isConnector"
)

missing_drug_columns <- setdiff(
  required_drug_columns,
  names(drug_network)
)

if (length(missing_drug_columns) > 0) {
  stop(
    "Missing columns in drug network CSV: ",
    paste(missing_drug_columns, collapse = ", ")
  )
}


# ============================================================
# 4. LOAD GENE DATA
# ============================================================

gene_data <- read.delim(
  gene_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_gene_columns <- c(
  "gene_id",
  "gene_name",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "pvalue",
  "padj"
)

missing_gene_columns <- setdiff(
  required_gene_columns,
  names(gene_data)
)

if (length(missing_gene_columns) > 0) {
  stop(
    "Missing columns in gene file: ",
    paste(missing_gene_columns, collapse = ", ")
  )
}

gene_data <- gene_data |>
  mutate(
    baseMean = suppressWarnings(as.numeric(baseMean)),
    log2FoldChange = suppressWarnings(as.numeric(log2FoldChange)),
    lfcSE = suppressWarnings(as.numeric(lfcSE)),
    pvalue = suppressWarnings(as.numeric(pvalue)),
    padj = suppressWarnings(as.numeric(padj))
  ) |>
  mutate(
    gene_id_clean = str_remove(
      as.character(gene_id),
      "\\.[0-9]+$"
    )
  )


# ============================================================
# 5. FIND JSON DIRECTORIES
# ============================================================

sources <- c(
  "chembl",
  "opentargets",
  "pubchem"
)

find_drug_dirs <- function(source) {
  
  source_dir <- file.path(
    clinreport_dir,
    source
  )
  
  if (!dir.exists(source_dir)) {
    return(tibble())
  }
  
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
# 6. JSON READER
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
        "JSON error: ",
        path,
        " | ",
        e$message
      )
      
      NULL
    }
  )
}


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


# ============================================================
# 8. IDENTIFIERS
# ============================================================

get_chembl_id <- function(data) {
  
  candidates <- list(
    get_nested(
      data$chembl$chembl_id,
      "chembl_id"
    ),
    
    get_nested(
      data$chembl$chembl_molecule,
      "chembl_id"
    )
  )
  
  candidates <- map_chr(
    candidates,
    ~ safe_text(.x, "")
  )
  
  candidates[
    nzchar(candidates)
  ][1] %||% NULL
}


get_pubchem_cid <- function(data) {
  
  x <- get_nested(
    data$pubchem$pubchem_cids,
    "IdentifierList",
    "CID"
  )
  
  value <- safe_text(
    x,
    ""
  )
  
  if (!nzchar(value)) {
    return(NULL)
  }
  
  str_split(
    value,
    ";"
  )[[1]][1]
}


`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x[1])) y else x
}


# ============================================================
# 9. DRUG NAME
# ============================================================

get_drug_name_from_json <- function(
    drug_id
) {
  
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
  
  drug_id
}


# ============================================================
# 10. EXTERNAL LINKS
# ============================================================

pubchem_url <- function(cid) {
  
  if (
    is.null(cid) ||
    !nzchar(cid)
  ) {
    return(NULL)
  }
  
  paste0(
    "https://pubchem.ncbi.nlm.nih.gov/compound/",
    cid
  )
}


chembl_url <- function(chembl_id) {
  
  if (
    is.null(chembl_id) ||
    !nzchar(chembl_id)
  ) {
    return(NULL)
  }
  
  paste0(
    "https://www.ebi.ac.uk/chembl/explore/compound/",
    chembl_id
  )
}


drugbank_url <- function(drug_id) {
  
  if (
    is.null(drug_id) ||
    !nzchar(drug_id)
  ) {
    return(NULL)
  }
  
  paste0(
    "https://go.drugbank.com/drugs/",
    drug_id
  )
}


ensembl_url <- function(gene_id) {
  
  if (
    is.null(gene_id) ||
    !nzchar(gene_id)
  ) {
    return(NULL)
  }
  
  paste0(
    "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=",
    gene_id
  )
}


# ============================================================
# 11. GENES FOR DRUG
# ============================================================

extract_gene_ids <- function(x) {
  
  if (is.null(x)) {
    return(character())
  }
  
  x <- as.character(x)
  
  if (
    length(x) == 0 ||
    all(is.na(x))
  ) {
    return(character())
  }
  
  ids <- str_extract_all(
    x,
    "ENSG[0-9]+"
  )[[1]]
  
  unique(ids)
}


get_drug_genes <- function(
    drug_id
) {
  
  row <- drug_network |>
    filter(
      drugId == !!drug_id
    ) |>
    slice(1)
  
  if (nrow(row) == 0) {
    return(tibble())
  }
  
  gene_ids <- extract_gene_ids(
    row$hasEdgesTo
  )
  
  if (length(gene_ids) == 0) {
    return(tibble())
  }
  
  gene_data |>
    filter(
      gene_id_clean %in% gene_ids
    ) |>
    select(
      gene_id,
      gene_id_clean,
      gene_name,
      baseMean,
      log2FoldChange,
      lfcSE,
      pvalue,
      padj
    ) |>
    distinct(
      gene_id_clean,
      .keep_all = TRUE
    ) |>
    mutate(
      
      significant_padj =
        !is.na(padj) &
        padj < padj_cutoff,
      
      significant_fc =
        !is.na(log2FoldChange) &
        abs(log2FoldChange) >= log2fc_cutoff,
      
      significant_both =
        significant_padj &
        significant_fc,
      
      direction =
        case_when(
          is.na(log2FoldChange) ~ "Neutral",
          log2FoldChange > 0 ~ "Up",
          log2FoldChange < 0 ~ "Down",
          TRUE ~ "Neutral"
        )
    ) |>
    arrange(
      is.na(padj),
      padj
    )
}


# ============================================================
# 12. EXTERNAL SOURCE BUTTONS
# ============================================================

source_buttons <- function(
    drug_id,
    chembl_id,
    pubchem_cid
) {
  
  buttons <- list()
  
  if (!is.null(drug_id)) {
    
    buttons <- append(
      buttons,
      list(
        tags$a(
          href = drugbank_url(drug_id),
          target = "_blank",
          rel = "noopener noreferrer",
          class = "source-button",
          "DrugBank"
        )
      )
    )
  }
  
  if (!is.null(chembl_id)) {
    
    buttons <- append(
      buttons,
      list(
        tags$a(
          href = chembl_url(chembl_id),
          target = "_blank",
          rel = "noopener noreferrer",
          class = "source-button",
          "ChEMBL"
        )
      )
    )
  }
  
  if (!is.null(pubchem_cid)) {
    
    buttons <- append(
      buttons,
      list(
        tags$a(
          href = pubchem_url(pubchem_cid),
          target = "_blank",
          rel = "noopener noreferrer",
          class = "source-button",
          "PubChem"
        )
      )
    )
  }
  
  if (length(buttons) == 0) {
    return(NULL)
  }
  
  div(
    class = "source-buttons",
    buttons
  )
}


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
# 14. 2D STRUCTURE
# ============================================================

render_structure <- function(
    pubchem_cid
) {
  
  if (
    is.null(pubchem_cid) ||
    !nzchar(pubchem_cid)
  ) {
    
    return(
      div(
        class = "structure-empty",
        "2D structure is not available from PubChem for this record."
      )
    )
  }
  
  image_url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
    pubchem_cid,
    "/PNG"
  )
  
  div(
    class = "structure-container",
    
    tags$img(
      src = image_url,
      class = "chemical-structure",
      alt = "2D chemical structure"
    ),
    
    div(
      class = "structure-caption",
      
      "PubChem CID ",
      
      tags$a(
        href = pubchem_url(pubchem_cid),
        target = "_blank",
        rel = "noopener noreferrer",
        class = "external-link",
        pubchem_cid
      )
    )
  )
}


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


# ============================================================
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


# ============================================================
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
# 18. TARGET CARDS
# ============================================================

render_target_cards <- function(
    genes
) {
  
  if (nrow(genes) == 0) {
    
    return(
      div(
        class = "empty-message",
        "No target genes found."
      )
    )
  }
  
  
  map(
    seq_len(nrow(genes)),
    function(i) {
      
      g <- genes[i, ]
      
      
      padj_pass <- isTRUE(
        !is.na(g$padj) &&
          g$padj < padj_cutoff
      )
      
      fc_pass <- isTRUE(
        !is.na(g$log2FoldChange) &&
          abs(g$log2FoldChange) >= log2fc_cutoff
      )
      
      both_pass <- padj_pass && fc_pass
      
      
      direction_class <-
        
        if (
          is.na(g$log2FoldChange)
        ) {
          
          "neutral"
          
        } else if (
          g$log2FoldChange > 0
        ) {
          
          "up"
          
        } else if (
          g$log2FoldChange < 0
        ) {
          
          "down"
          
        } else {
          
          "neutral"
        }
      
      
      cutoff_class <-
        
        if (both_pass) {
          
          "passes-both"
          
        } else if (padj_pass) {
          
          "passes-padj"
          
        } else if (fc_pass) {
          
          "passes-fc"
          
        } else {
          
          "passes-neither"
        }
      
      
      card_class <- paste(
        "gene-card",
        direction_class,
        cutoff_class
      )
      
      
      onclick_js <- sprintf(
        paste0(
          "Shiny.setInputValue(",
          "'selected_gene_index', ",
          "%d, ",
          "{priority: 'event'}",
          ")"
        ),
        i
      )
      
      
      gene_name <- ifelse(
        is.na(g$gene_name) ||
          !nzchar(g$gene_name),
        "Unknown",
        g$gene_name
      )
      
      
      fc_text <- ifelse(
        is.na(g$log2FoldChange),
        "—",
        sprintf(
          "%+.2f",
          g$log2FoldChange
        )
      )
      
      
      padj_text <- ifelse(
        is.na(g$padj),
        "—",
        format.pval(
          g$padj,
          digits = 2,
          eps = 1e-300
        )
      )
      
      
      div(
        
        class = card_class,
        
        onclick = onclick_js,
        
        
        div(
          class = "gene-main",
          
          div(
            class = "gene-symbol",
            gene_name
          ),
          
          div(
            class = "gene-id-small",
            
            tags$a(
              href = ensembl_url(
                g$gene_id_clean
              ),
              target = "_blank",
              rel = "noopener noreferrer",
              class = "gene-external-link",
              g$gene_id_clean
            )
          )
        ),
        
        
        div(
          class = "gene-expression",
          
          div(
            class = "expression-direction",
            
            span(
              class = paste0(
                "direction-arrow ",
                direction_class
              ),
              
              ifelse(
                direction_class == "up",
                "↑",
                ifelse(
                  direction_class == "down",
                  "↓",
                  "—"
                )
              )
            ),
            
            span(
              class = "expression-value",
              fc_text
            )
          )
        ),
        
        
        div(
          class = "gene-padj",
          
          div(
            class = "gene-stat-label",
            "padj"
          ),
          
          div(
            class = "gene-stat-value",
            padj_text
          )
        ),
        
        
        div(
          class = "gene-status",
          
          if (both_pass) {
            
            span(
              class = "cutoff-badge both",
              "meets criteria"
            )
            
          } else if (padj_pass) {
            
            span(
              class = "cutoff-badge padj-only",
              "padj"
            )
            
          } else if (fc_pass) {
            
            span(
              class = "cutoff-badge fc-only",
              "fold change"
            )
            
          } else {
            
            span(
              class = "cutoff-badge neither",
              "below criteria"
            )
          }
        ),
        
        
        div(
          class = "gene-arrow",
          "›"
        )
      )
    }
  )
}


# ============================================================
# 19. DRUG PAGE
# ============================================================

render_drug_page <- function(
    drug_id,
    drug_name,
    data,
    genes
) {
  
  chembl_id <- get_chembl_id(data)
  pubchem_cid <- get_pubchem_cid(data)
  
  
  tagList(
    
    div(
      class = "drug-header",
      
      div(
        class = "drug-title-row",
        
        div(
          
          div(
            class = "drug-title",
            drug_name
          ),
          
          div(
            class = "drug-id",
            drug_id
          )
        ),
        
        source_buttons(
          drug_id,
          chembl_id,
          pubchem_cid
        )
      )
    ),
    
    
    tabsetPanel(
      
      id = "drug_tabs",
      
      
      tabPanel(
        "General",
        
        render_general(
          drug_id,
          data
        )
      ),
      
      
      tabPanel(
        "Indications & Mechanism",
        
        render_pharmacology(
          data
        )
      ),
      
      
      tabPanel(
        "Safety",
        
        render_safety(
          data
        )
      ),
      
      
      tabPanel(
        "Targets",
        
        div(
          class = "targets-header",
          
          div(
            
            div(
              class = "targets-title",
              "Target genes"
            ),
            
            div(
              class = "targets-criteria",
              
              tags$span(
                class = "criteria-chip",
                paste0(
                  "padj < ",
                  padj_cutoff
                )
              ),
              
              tags$span(
                class = "criteria-chip",
                paste0(
                  "|log2FC| ≥ ",
                  log2fc_cutoff
                )
              )
            )
          ),
          
          div(
            class = "targets-count",
            
            strong(
              nrow(genes)
            ),
            
            " ",
            
            ifelse(
              nrow(genes) == 1,
              "gene",
              "genes"
            )
          )
        ),
        
        div(
          class = "targets-list",
          
          render_target_cards(
            genes
          )
        )
      )
    )
  )
}


# ============================================================
# 20. GENE PAGE
# ============================================================

render_gene_page <- function(
    gene
) {
  
  gene_name <- ifelse(
    is.na(gene$gene_name) ||
      !nzchar(gene$gene_name),
    "Unknown gene",
    gene$gene_name
  )
  
  
  div(
    
    actionButton(
      "back_to_drug",
      "← Back to drug",
      class = "back-button"
    ),
    
    
    div(
      class = "gene-header",
      
      div(
        class = "gene-title-large",
        gene_name
      ),
      
      div(
        class = "gene-id",
        
        tags$a(
          href = ensembl_url(
            gene$gene_id_clean
          ),
          target = "_blank",
          rel = "noopener noreferrer",
          class = "external-link",
          gene$gene_id_clean
        )
      )
    ),
    
    
    div(
      class = "gene-summary",
      
      div(
        class = "summary-card",
        
        div(
          class = "summary-label",
          "log2 Fold Change"
        ),
        
        div(
          class = "summary-value",
          
          ifelse(
            is.na(gene$log2FoldChange),
            "—",
            sprintf(
              "%+.3f",
              gene$log2FoldChange
            )
          )
        ),
        
        div(
          class = "summary-note",
          
          ifelse(
            is.na(gene$log2FoldChange),
            "Not available",
            ifelse(
              gene$log2FoldChange > 0,
              "Up-regulated",
              ifelse(
                gene$log2FoldChange < 0,
                "Down-regulated",
                "No change"
              )
            )
          )
        )
      ),
      
      
      div(
        class = "summary-card",
        
        div(
          class = "summary-label",
          "Adjusted p-value"
        ),
        
        div(
          class = "summary-value",
          
          ifelse(
            is.na(gene$padj),
            "—",
            format.pval(
              gene$padj,
              digits = 3,
              eps = 1e-300
            )
          )
        ),
        
        div(
          class = "summary-note",
          
          ifelse(
            !is.na(gene$padj) &&
              gene$padj < padj_cutoff,
            "Below threshold",
            "Not below threshold"
          )
        )
      ),
      
      
      div(
        class = "summary-card",
        
        div(
          class = "summary-label",
          "baseMean"
        ),
        
        div(
          class = "summary-value",
          
          ifelse(
            is.na(gene$baseMean),
            "—",
            format(
              round(
                gene$baseMean,
                2
              ),
              big.mark = ","
            )
          )
        ),
        
        div(
          class = "summary-note",
          "Mean normalized expression"
        )
      )
    ),
    
    
    div(
      class = "gene-filter-status",
      
      if (
        !is.na(gene$padj) &&
        gene$padj < padj_cutoff
      ) {
        
        span(
          class = "gene-filter-badge pass",
          paste0(
            "padj < ",
            padj_cutoff
          )
        )
        
      } else {
        
        span(
          class = "gene-filter-badge fail",
          paste0(
            "padj ≥ ",
            padj_cutoff
          )
        )
      },
      
      
      if (
        !is.na(gene$log2FoldChange) &&
        abs(gene$log2FoldChange) >= log2fc_cutoff
      ) {
        
        span(
          class = "gene-filter-badge pass",
          paste0(
            "|log2FC| ≥ ",
            log2fc_cutoff
          )
        )
        
      } else {
        
        span(
          class = "gene-filter-badge fail",
          paste0(
            "|log2FC| < ",
            log2fc_cutoff
          )
        )
      }
    ),
    
    
    div(
      class = "section-header",
      
      div(
        class = "section-title",
        "Gene expression overview"
      ),
      
      div(
        class = "section-subtitle",
        "Genome-wide differential expression context"
      )
    ),
    
    plotlyOutput(
      "gene_volcano",
      height = "600px"
    )
  )
}


# ============================================================
# 21. UI
# ============================================================

ui <- fluidPage(
  
  tags$head(
    
    tags$style(
      
      HTML(
        "
        /* ==================================================
           GLOBAL
           ================================================== */

        * {
          box-sizing: border-box;
        }

        body {
          background: #f3f5f7;
          font-family:
            -apple-system,
            BlinkMacSystemFont,
            'Segoe UI',
            sans-serif;
          color: #20242a;
          margin: 0;
        }

        .container-fluid {
          padding-left: 26px;
          padding-right: 26px;
        }

        .app-title {
          font-size: 28px;
          font-weight: 650;
          letter-spacing: -0.5px;
          margin: 24px 0 5px 0;
          color: #171b20;
        }

        .app-subtitle {
          color: #8a9199;
          font-size: 13px;
          margin-bottom: 20px;
        }


        /* ==================================================
           MAIN LAYOUT
           ================================================== */

        .main-layout {
          display: grid;
          grid-template-columns: minmax(380px, 43%) minmax(0, 57%);
          gap: 16px;
          align-items: stretch;
        }


        /* ==================================================
           LEFT / NETWORK ZONE
           ================================================== */

        .drug-list {
          background: #f8fafb;
          border: 1px solid #e2e6ea;
          border-radius: 14px;
          padding: 14px;
          box-shadow:
            0 2px 8px rgba(20, 30, 40, 0.035);
          min-height: 760px;
        }

        .network-header {
          padding: 7px 6px 12px 6px;
        }

        .network-title {
          font-size: 16px;
          font-weight: 650;
          color: #30363d;
        }

        .network-subtitle {
          color: #9299a1;
          font-size: 11px;
          margin-top: 3px;
        }


        /* ==================================================
           RIGHT / DETAIL ZONE
           ================================================== */

        .drug-details {
          background: white;
          border: 1px solid #e4e7ea;
          border-radius: 14px;
          padding: 28px 30px;
          min-height: 760px;
          box-shadow:
            0 2px 10px rgba(20, 30, 40, 0.045);
          overflow: hidden;
        }


        /* ==================================================
           DRUG HEADER
           ================================================== */

        .drug-header {
          margin-bottom: 4px;
        }

        .drug-title-row {
          display: flex;
          justify-content: space-between;
          align-items: flex-start;
          gap: 20px;
        }

        .drug-title {
          font-size: 29px;
          font-weight: 650;
          line-height: 1.15;
          letter-spacing: -0.5px;
          color: #171b20;
        }

        .drug-id {
          color: #9aa1a8;
          font-size: 12px;
          margin-top: 6px;
          font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
        }


        /* ==================================================
           SOURCE BUTTONS
           ================================================== */

        .source-buttons {
          display: flex;
          gap: 6px;
          flex-wrap: wrap;
          justify-content: flex-end;
          margin-top: 2px;
        }

        .source-button {
          display: inline-flex;
          align-items: center;
          padding: 6px 9px;
          border: 1px solid #dfe3e7;
          border-radius: 7px;
          background: #fafbfc;
          color: #616971 !important;
          font-size: 10px;
          font-weight: 600;
          text-decoration: none !important;
          transition:
            background 0.12s ease,
            border-color 0.12s ease;
        }

        .source-button:hover {
          background: white;
          border-color: #c8ced5;
        }


        /* ==================================================
           TABS
           ================================================== */

        .nav-tabs {
          margin-top: 22px;
          border-bottom: 1px solid #e6e8eb;
        }

        .nav-tabs > li > a {
          color: #7d858d;
          border: none !important;
          padding: 10px 14px;
          font-size: 13px;
          transition:
            color 0.12s ease;
        }

        .nav-tabs > li > a:hover {
          background: transparent !important;
          color: #394149;
        }

        .nav-tabs > li.active > a,
        .nav-tabs > li.active > a:hover,
        .nav-tabs > li.active > a:focus {
          color: #252b31;
          background: transparent !important;
          border: none !important;
          border-bottom: 2px solid #3d4650 !important;
          font-weight: 650;
        }

        .tab-content {
          padding-top: 8px;
        }


        /* ==================================================
           SECTION HEADERS
           ================================================== */

        .section-header {
          margin-top: 23px;
          margin-bottom: 12px;
        }

        .section-title {
          font-size: 16px;
          font-weight: 650;
          color: #30363d;
        }

        .section-subtitle {
          color: #9aa1a8;
          font-size: 11px;
          margin-top: 3px;
        }


        /* ==================================================
           IDENTIFIERS
           ================================================== */

        .identifier-grid {
          display: grid;
          grid-template-columns:
            minmax(130px, 1fr)
            minmax(130px, 1fr)
            minmax(130px, 1fr);
          gap: 8px;
        }

        .identifier-card {
          background: #fafbfc;
          border: 1px solid #e7e9ec;
          border-radius: 9px;
          padding: 11px 12px;
          min-width: 0;
        }

        .identifier-card.wide {
          grid-column: 1 / -1;
        }

        .identifier-label {
          color: #969da5;
          font-size: 10px;
          text-transform: uppercase;
          letter-spacing: 0.35px;
          margin-bottom: 5px;
        }

        .identifier-value {
          color: #343a41;
          font-size: 13px;
          font-weight: 550;
          overflow-wrap: anywhere;
        }

        .preferred-name {
          font-size: 14px;
        }

        .identifier-link,
        .external-link,
        .target-link,
        .gene-external-link {
          color: #526c84 !important;
          text-decoration: none !important;
        }

        .identifier-link:hover,
        .external-link:hover,
        .target-link:hover,
        .gene-external-link:hover {
          text-decoration: underline !important;
        }


        /* ==================================================
           MOLECULAR PROFILE
           ================================================== */

        .molecular-grid {
          display: grid;
          grid-template-columns: repeat(2, minmax(0, 1fr));
          gap: 9px;
        }

        .molecular-metric {
          background: #fafbfc;
          border: 1px solid #e7e9ec;
          border-radius: 9px;
          padding: 11px 12px;
        }

        .metric-top {
          display: flex;
          align-items: center;
          justify-content: space-between;
          gap: 10px;
          margin-bottom: 8px;
        }

        .metric-label {
          color: #7e868e;
          font-size: 11px;
        }

        .metric-value {
          color: #30363d;
          font-size: 13px;
          font-weight: 650;
          white-space: nowrap;
        }

        .metric-track {
          height: 4px;
          background: #e9ecef;
          border-radius: 5px;
          overflow: hidden;
        }

        .metric-fill {
          height: 100%;
          background: #7e93a8;
          border-radius: 5px;
          min-width: 2px;
        }

        .formula-card {
          margin-top: 9px;
          padding: 11px 12px;
          border: 1px solid #e7e9ec;
          border-radius: 9px;
          background: #fcfcfd;
        }

        .formula-label {
          color: #969da5;
          font-size: 10px;
          text-transform: uppercase;
          letter-spacing: 0.35px;
          margin-bottom: 5px;
        }

        .formula-value {
          font-size: 14px;
          color: #343a41;
        }


        /* ==================================================
           DEVELOPMENT
           ================================================== */

        .development-grid {
          display: grid;
          grid-template-columns: 1fr 1fr;
          gap: 9px;
        }

        .development-card {
          border: 1px solid #e7e9ec;
          background: #fafbfc;
          border-radius: 9px;
          padding: 12px;
        }

        .development-label {
          color: #9299a1;
          font-size: 10px;
          margin-bottom: 6px;
        }

        .development-value {
          font-size: 18px;
          font-weight: 650;
          color: #343a41;
        }

        .phase-badge {
          display: inline-flex;
          align-items: center;
          padding: 4px 7px;
          border-radius: 5px;
          background: #edf1f4;
          color: #596773;
          font-size: 10px;
          font-weight: 650;
          white-space: nowrap;
        }

        .phase-badge.large {
          font-size: 13px;
          padding: 6px 9px;
        }


        /* ==================================================
           STRUCTURES
           ================================================== */

        .structure-layout {
          display: grid;
          grid-template-columns: minmax(250px, 38%) minmax(0, 62%);
          gap: 13px;
          align-items: stretch;
        }

        .structure-visual {
          min-height: 280px;
          border: 1px solid #e5e8eb;
          border-radius: 10px;
          background: #fbfcfd;
          display: flex;
          align-items: center;
          justify-content: center;
        }

        .structure-container {
          width: 100%;
          height: 100%;
          min-height: 280px;
          display: flex;
          flex-direction: column;
          align-items: center;
          justify-content: center;
          padding: 14px;
        }

        .chemical-structure {
          max-width: 100%;
          max-height: 240px;
          object-fit: contain;
        }

        .structure-caption {
          margin-top: 7px;
          color: #9aa1a8;
          font-size: 10px;
        }

        .structure-data {
          display: flex;
          flex-direction: column;
          gap: 9px;
        }

        .structure-field {
          border: 1px solid #e7e9ec;
          background: #fafbfc;
          border-radius: 9px;
          padding: 11px 12px;
        }

        .structure-field-label {
          color: #969da5;
          font-size: 10px;
          text-transform: uppercase;
          letter-spacing: 0.35px;
          margin-bottom: 6px;
        }

        .structure-code {
          font-family:
            ui-monospace,
            SFMono-Regular,
            Menlo,
            Monaco,
            Consolas,
            monospace;
          font-size: 11px;
          line-height: 1.5;
          color: #454c54;
          overflow-wrap: anywhere;
          word-break: break-word;
        }

        .mono-value {
          font-family:
            ui-monospace,
            SFMono-Regular,
            Menlo,
            Monaco,
            Consolas,
            monospace;
        }

        .structure-empty {
          color: #9aa1a8;
          text-align: center;
          padding: 25px;
          font-size: 12px;
        }


        /* ==================================================
           INDICATIONS
           ================================================== */

        .indications-grid {
          display: grid;
          grid-template-columns: repeat(2, minmax(0, 1fr));
          gap: 9px;
        }

        .indication-card {
          border: 1px solid #e5e8eb;
          border-radius: 10px;
          background: #fafbfc;
          padding: 12px;
        }

        .indication-top {
          display: flex;
          align-items: flex-start;
          justify-content: space-between;
          gap: 10px;
        }

        .indication-name {
          font-size: 14px;
          font-weight: 650;
          color: #353c43;
          line-height: 1.35;
        }

        .indication-meta {
          display: flex;
          flex-direction: column;
          gap: 5px;
          margin-top: 10px;
          padding-top: 9px;
          border-top: 1px solid #e8eaed;
        }

        .meta-item {
          display: flex;
          gap: 8px;
          font-size: 10px;
        }

        .meta-label {
          width: 35px;
          color: #9aa1a8;
          text-transform: uppercase;
        }

        .meta-value {
          color: #606870;
          overflow-wrap: anywhere;
        }


        /* ==================================================
           MECHANISMS
           ================================================== */

        .pharmacology-section {
          margin-top: 28px;
        }

        .mechanism-list {
          display: flex;
          flex-direction: column;
          gap: 8px;
        }

        .mechanism-card {
          display: flex;
          justify-content: space-between;
          align-items: center;
          gap: 15px;
          padding: 13px 14px;
          border: 1px solid #e5e8eb;
          border-radius: 10px;
          background: #fafbfc;
        }

        .mechanism-main {
          min-width: 0;
          flex: 1;
        }

        .mechanism-title {
          font-size: 14px;
          font-weight: 600;
          color: #343b42;
        }

        .mechanism-target {
          margin-top: 6px;
          font-size: 11px;
          color: #8d959d;
        }

        .target-label {
          margin-right: 6px;
          color: #a0a6ad;
        }

        .mechanism-meta {
          display: flex;
          align-items: center;
          gap: 10px;
          flex-shrink: 0;
        }

        .mechanism-pill {
          padding: 5px 8px;
          border-radius: 5px;
          background: #edf1f4;
          color: #596773;
          font-size: 9px;
          font-weight: 650;
        }

        .mechanism-detail {
          text-align: right;
        }

        .mechanism-detail-label {
          display: block;
          color: #a0a6ad;
          font-size: 9px;
        }

        .mechanism-detail-value {
          display: block;
          margin-top: 2px;
          color: #555e67;
          font-size: 11px;
          font-weight: 600;
        }


        /* ==================================================
           SAFETY
           ================================================== */

        .safety-section {
          margin-top: 28px;
        }

        .warning-list {
          display: flex;
          flex-direction: column;
          gap: 9px;
        }

        .warning-card {
          display: flex;
          align-items: center;
          gap: 12px;
          padding: 13px 14px;
          border: 1px solid #eadfdf;
          border-left: 4px solid #b76d6d;
          border-radius: 9px;
          background: #fdfafa;
        }

        .warning-icon {
          width: 26px;
          height: 26px;
          flex: 0 0 26px;
          display: flex;
          align-items: center;
          justify-content: center;
          border-radius: 50%;
          background: #f0dddd;
          color: #9d5555;
          font-weight: 700;
          font-size: 13px;
        }

        .warning-content {
          min-width: 0;
        }

        .warning-class {
          color: #4b3636;
          font-size: 14px;
          font-weight: 650;
        }

        .warning-meta {
          display: flex;
          gap: 8px;
          margin-top: 5px;
        }

        .warning-type,
        .warning-country {
          display: inline-flex;
          padding: 4px 6px;
          border-radius: 4px;
          background: #f2eded;
          color: #766666;
          font-size: 9px;
          font-weight: 600;
        }

        .synonym-cloud {
          display: flex;
          flex-wrap: wrap;
          gap: 6px;
        }

        .synonym-chip {
          display: inline-flex;
          padding: 6px 8px;
          border-radius: 6px;
          border: 1px solid #e1e4e7;
          background: #fafbfc;
          color: #626a72;
          font-size: 10px;
        }


        /* ==================================================
           TARGETS
           ================================================== */

        .targets-header {
          display: flex;
          align-items: flex-end;
          justify-content: space-between;
          margin-top: 10px;
          margin-bottom: 16px;
        }

        .targets-title {
          font-size: 20px;
          font-weight: 650;
          color: #30363d;
        }

        .targets-criteria {
          display: flex;
          gap: 5px;
          margin-top: 7px;
        }

        .criteria-chip {
          padding: 4px 6px;
          border-radius: 5px;
          background: #f0f2f4;
          color: #7c858e;
          font-size: 9px;
          font-weight: 600;
        }

        .targets-count {
          color: #8f979f;
          font-size: 12px;
          white-space: nowrap;
        }

        .targets-list {
          display: flex;
          flex-direction: column;
          gap: 7px;
        }

        .gene-card {
          display: flex;
          align-items: center;
          padding: 11px 13px;
          background: #fafbfc;
          border: 1px solid #e5e8eb;
          border-radius: 9px;
          cursor: pointer;
          transition:
            transform 0.12s ease,
            box-shadow 0.12s ease,
            background 0.12s ease,
            border-color 0.12s ease;
        }

        .gene-card:hover {
          transform: translateY(-1px);
          box-shadow:
            0 4px 12px rgba(20,30,40,0.06);
          background: white;
          border-color: #d5dbe0;
        }

        .gene-card.up {
          border-left: 4px solid #7698b8;
        }

        .gene-card.down {
          border-left: 4px solid #b77a7a;
        }

        .gene-card.neutral {
          border-left: 4px solid #aeb5bc;
        }

        .gene-card.passes-both {
          background: #f8fafc;
        }

        .gene-card.passes-neither {
          opacity: 0.62;
        }

        .gene-main {
          flex: 1;
          min-width: 0;
        }

        .gene-symbol {
          font-size: 15px;
          font-weight: 650;
          color: #30363d;
        }

        .gene-id-small {
          font-size: 9px;
          margin-top: 3px;
          overflow: hidden;
          text-overflow: ellipsis;
          white-space: nowrap;
        }

        .gene-expression {
          width: 90px;
          text-align: right;
        }

        .expression-direction {
          display: flex;
          justify-content: flex-end;
          align-items: center;
          gap: 5px;
        }

        .direction-arrow {
          font-size: 15px;
          font-weight: 700;
        }

        .direction-arrow.up {
          color: #5d83a5;
        }

        .direction-arrow.down {
          color: #a96767;
        }

        .direction-arrow.neutral {
          color: #9da4ab;
        }

        .expression-value {
          font-size: 12px;
          font-weight: 650;
          color: #4a525a;
        }

        .gene-padj {
          width: 90px;
          text-align: right;
          margin-left: 10px;
        }

        .gene-stat-label {
          color: #a0a6ad;
          font-size: 9px;
          margin-bottom: 2px;
        }

        .gene-stat-value {
          color: #4a525a;
          font-size: 11px;
          font-weight: 650;
        }

        .gene-status {
          width: 105px;
          text-align: center;
          margin-left: 10px;
        }

        .cutoff-badge {
          display: inline-block;
          padding: 4px 7px;
          border-radius: 5px;
          font-size: 9px;
          font-weight: 650;
          white-space: nowrap;
        }

        .cutoff-badge.both {
          background: #eaf0f5;
          color: #506b83;
        }

        .cutoff-badge.padj-only {
          background: #eef3f6;
          color: #5b778b;
        }

        .cutoff-badge.fc-only {
          background: #f1eff4;
          color: #6e647d;
        }

        .cutoff-badge.neither {
          background: #eef0f2;
          color: #858c93;
        }

        .gene-arrow {
          color: #b1b7bd;
          font-size: 21px;
          width: 17px;
          text-align: center;
          transition:
            transform 0.12s ease;
        }

        .gene-card:hover .gene-arrow {
          transform: translateX(2px);
        }


        /* ==================================================
           GENE PAGE
           ================================================== */

        .back-button {
          margin-bottom: 20px;
          border-radius: 7px;
          color: #59616a;
          border-color: #dfe3e7;
          background: white;
          font-size: 12px;
        }

        .back-button:hover {
          background: #f7f8f9;
        }

        .gene-title-large {
          font-size: 30px;
          font-weight: 650;
          letter-spacing: -0.4px;
          color: #20252b;
        }

        .gene-id {
          color: #8d959d;
          font-size: 11px;
          margin-top: 5px;
        }

        .gene-summary {
          display: grid;
          grid-template-columns: repeat(3, 1fr);
          gap: 9px;
          margin: 15px 0 13px 0;
        }

        .summary-card {
          background: #fafbfc;
          border: 1px solid #e5e7eb;
          border-radius: 9px;
          padding: 14px;
        }

        .summary-label {
          color: #858d95;
          font-size: 10px;
          margin-bottom: 5px;
        }

        .summary-value {
          font-size: 20px;
          font-weight: 650;
          color: #30363d;
        }

        .summary-note {
          color: #9ba2a9;
          font-size: 9px;
          margin-top: 4px;
        }

        .gene-filter-status {
          display: flex;
          gap: 6px;
          margin: 5px 0 18px 0;
        }

        .gene-filter-badge {
          display: inline-block;
          padding: 4px 8px;
          border-radius: 5px;
          font-size: 9px;
          font-weight: 650;
        }

        .gene-filter-badge.pass {
          background: #edf3ef;
          color: #58725f;
        }

        .gene-filter-badge.fail {
          background: #f0f1f3;
          color: #858c93;
        }


        /* ==================================================
           DATATABLE / INPUT TABLE
           ================================================== */

        .drug-list .dataTables_wrapper {
          padding: 0;
        }

        .drug-list .dataTables_filter {
          float: none;
          margin: 4px 0 10px 0;
        }

        .drug-list .dataTables_filter label {
          width: 100%;
          color: #9299a1;
          font-size: 10px;
        }

        .drug-list .dataTables_filter input {
          width: 100%;
          margin-left: 0 !important;
          margin-top: 5px;
          border: 1px solid #dfe3e7 !important;
          border-radius: 7px !important;
          padding: 7px 9px !important;
          font-size: 12px !important;
          box-shadow: none !important;
        }

        table.dataTable {
          width: 100% !important;
          font-size: 12px;
          border-collapse: separate !important;
          border-spacing: 0 4px !important;
        }

        table.dataTable thead th {
          color: #8b939b;
          font-weight: 650;
          font-size: 10px;
          text-transform: uppercase;
          letter-spacing: 0.3px;
          border-bottom: none !important;
          padding: 8px 8px !important;
          white-space: nowrap;
        }

        table.dataTable tbody tr {
          background: white;
        }

        table.dataTable tbody td {
          vertical-align: middle;
          border-top: 1px solid #e8eaed;
          border-bottom: 1px solid #e8eaed;
          padding: 9px 8px !important;
          color: #4c545c;
        }

        table.dataTable tbody td:first-child {
          border-left: 1px solid #e8eaed;
          border-radius: 7px 0 0 7px;
        }

        table.dataTable tbody td:last-child {
          border-right: 1px solid #e8eaed;
          border-radius: 0 7px 7px 0;
        }

        table.dataTable tbody tr:hover td {
          background: #f8fafc !important;
        }

        table.dataTable tbody tr.selected td {
          background: #eef3f7 !important;
          border-color: #dce5ec;
        }

        .drug-table-name {
          font-weight: 600;
          color: #343b42;
        }

        .status-badge {
          display: inline-block;
          padding: 4px 6px;
          border-radius: 5px;
          background: #eef1f3;
          color: #66707a;
          font-size: 9px;
          font-weight: 650;
        }

        .score-value {
          font-family:
            ui-monospace,
            SFMono-Regular,
            Menlo,
            monospace;
          font-size: 10px;
          color: #626b74;
        }


        /* ==================================================
           EMPTY
           ================================================== */

        .empty-message {
          color: #9aa1a8;
          padding: 38px 20px;
          text-align: center;
          font-size: 12px;
          border: 1px dashed #e1e4e7;
          border-radius: 9px;
          background: #fbfcfd;
        }


        /* ==================================================
           RESPONSIVE
           ================================================== */

        @media (max-width: 1100px) {

          .main-layout {
            grid-template-columns: 1fr;
          }

          .drug-list,
          .drug-details {
            min-height: auto;
          }

          .structure-layout {
            grid-template-columns: 1fr;
          }
        }

        @media (max-width: 750px) {

          .container-fluid {
            padding-left: 12px;
            padding-right: 12px;
          }

          .drug-details {
            padding: 20px;
          }

          .identifier-grid,
          .molecular-grid,
          .indications-grid,
          .development-grid,
          .gene-summary {
            grid-template-columns: 1fr;
          }

          .identifier-card.wide {
            grid-column: auto;
          }

          .drug-title-row {
            flex-direction: column;
          }

          .source-buttons {
            justify-content: flex-start;
          }

          .mechanism-card {
            flex-direction: column;
            align-items: flex-start;
          }

          .mechanism-meta {
            width: 100%;
            justify-content: space-between;
          }

          .gene-status {
            display: none;
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
    class = "app-subtitle",
    "Drug network · pharmacology · safety · molecular properties · target genes"
  ),
  
  
  div(
    class = "main-layout",
    
    
    # ========================================================
    # LEFT ZONE
    # ========================================================
    
    div(
      class = "drug-list",
      
      div(
        class = "network-header",
        
        div(
          class = "network-title",
          "Drug network"
        ),
        
        div(
          class = "network-subtitle",
          paste(
            nrow(drug_network),
            "records · select a drug to inspect details"
          )
        )
      ),
      
      DTOutput(
        "drug_table"
      )
    ),
    
    
    # ========================================================
    # RIGHT ZONE
    # ========================================================
    
    div(
      class = "drug-details",
      
      uiOutput(
        "right_panel"
      )
    )
  )
)


# ============================================================
# 22. SERVER
# ============================================================

server <- function(
    input,
    output,
    session
) {
  
  
  # ==========================================================
  # DRUG TABLE
  # ==========================================================
  
  output$drug_table <- renderDT({
    
    table_data <- drug_network |>
      transmute(
        
        `Drug` = label,
        
        `Status` = status,
        
        `Score` = score
      )
    
    
    datatable(
      
      table_data,
      
      selection = "single",
      
      rownames = FALSE,
      
      filter = "none",
      
      escape = FALSE,
      
      options = list(
        
        pageLength = 18,
        
        lengthChange = FALSE,
        
        searching = TRUE,
        
        scrollX = TRUE,
        
        autoWidth = FALSE,
        
        dom = "ftip",
        
        order = list(
          list(
            2,
            "desc"
          )
        ),
        
        columnDefs = list(
          
          list(
            width = "55%",
            targets = 0
          ),
          
          list(
            width = "25%",
            targets = 1
          ),
          
          list(
            width = "20%",
            targets = 2
          )
        )
      )
    ) |>
      formatStyle(
        "Drug",
        fontWeight = "600"
      ) |>
      formatStyle(
        "Status",
        color = "#66707a"
      ) |>
      formatRound(
        "Score",
        digits = 3
      )
  })
  
  
  # ==========================================================
  # SELECTED DRUG
  # ==========================================================
  
  selected_drug <- reactive({
    
    req(
      input$drug_table_rows_selected
    )
    
    drug_network[
      input$drug_table_rows_selected,
      ,
      drop = FALSE
    ]
  })
  
  
  # ==========================================================
  # CURRENT DRUG GENES
  # ==========================================================
  
  current_drug_genes <- reactive({
    
    req(
      selected_drug()
    )
    
    drug_id <- selected_drug()$drugId
    
    get_drug_genes(
      drug_id
    )
  })
  
  
  # ==========================================================
  # VIEW MODE
  # ==========================================================
  
  view_mode <- reactiveVal(
    "drug"
  )
  
  selected_gene <- reactiveVal(
    NULL
  )
  
  
  # ==========================================================
  # RIGHT PANEL
  # ==========================================================
  
  output$right_panel <- renderUI({
    
    req(
      selected_drug()
    )
    
    
    if (
      view_mode() == "gene"
    ) {
      
      gene <- selected_gene()
      
      req(
        gene
      )
      
      return(
        render_gene_page(
          gene
        )
      )
    }
    
    
    drug <- selected_drug()
    
    drug_id <- drug$drugId
    
    drug_name <- drug$label
    
    
    data <- load_drug(
      drug_id
    )
    
    
    genes <- current_drug_genes()
    
    
    has_json <- any(
      map_int(
        data,
        length
      ) > 0
    )
    
    
    if (!has_json) {
      
      return(
        
        tagList(
          
          div(
            class = "drug-header",
            
            div(
              class = "drug-title",
              drug_name
            ),
            
            div(
              class = "drug-id",
              drug_id
            )
          ),
          
          div(
            class = "empty-message",
            
            "Detailed information is not available ",
            "for this drug."
          )
        )
      )
    }
    
    
    render_drug_page(
      
      drug_id,
      
      drug_name,
      
      data,
      
      genes
    )
  })
  
  
  # ==========================================================
  # NEW DRUG SELECTED
  # ==========================================================
  
  observeEvent(
    
    input$drug_table_rows_selected,
    
    {
      
      view_mode(
        "drug"
      )
      
      selected_gene(
        NULL
      )
      
    },
    
    ignoreInit = TRUE
  )
  
  
  # ==========================================================
  # GENE CLICK
  # ==========================================================
  
  observeEvent(
    
    input$selected_gene_index,
    
    {
      
      req(
        selected_drug()
      )
      
      
      genes <- current_drug_genes()
      
      
      i <- as.integer(
        input$selected_gene_index
      )
      
      
      if (
        is.na(i) ||
        i < 1 ||
        i > nrow(genes)
      ) {
        
        return()
      }
      
      
      selected_gene(
        genes[
          i,
          ,
          drop = FALSE
        ]
      )
      
      
      view_mode(
        "gene"
      )
      
    }
  )
  
  
  # ==========================================================
  # BACK TO DRUG
  # ==========================================================
  
  observeEvent(
    
    input$back_to_drug,
    
    {
      
      view_mode(
        "drug"
      )
      
      selected_gene(
        NULL
      )
      
    }
  )
  
  
  # ==========================================================
  # VOLCANO
  # ==========================================================
  
  output$gene_volcano <- renderPlotly({
    
    gene <- selected_gene()
    
    req(
      gene
    )
    
    
    selected_id <-
      gene$gene_id_clean
    
    
    plot_data <- gene_data |>
      mutate(
        
        neg_log10_padj = ifelse(
          is.na(padj),
          NA_real_,
          -log10(
            pmax(
              padj,
              .Machine$double.xmin
            )
          )
        ),
        
        selected =
          gene_id_clean == selected_id,
        
        pass_padj =
          !is.na(padj) &
          padj < padj_cutoff,
        
        pass_fc =
          !is.na(log2FoldChange) &
          abs(log2FoldChange) >= log2fc_cutoff,
        
        pass_both =
          pass_padj &
          pass_fc
      )
    
    
    p <- plot_ly()
    
    
    # --------------------------------------------------------
    # BACKGROUND
    # --------------------------------------------------------
    
    background_data <- plot_data |>
      filter(
        !selected,
        !pass_both
      )
    
    
    if (nrow(background_data) > 0) {
      
      p <- p |>
        
        add_trace(
          
          data = background_data,
          
          x = ~log2FoldChange,
          
          y = ~neg_log10_padj,
          
          type = "scatter",
          
          mode = "markers",
          
          text = ~paste0(
            "<b>",
            ifelse(
              is.na(gene_name),
              "Unknown",
              gene_name
            ),
            "</b><br>",
            "Gene ID: ",
            gene_id_clean,
            "<br>log2FC: ",
            ifelse(
              is.na(log2FoldChange),
              "NA",
              round(
                log2FoldChange,
                3
              )
            ),
            "<br>padj: ",
            ifelse(
              is.na(padj),
              "NA",
              signif(
                padj,
                3
              )
            )
          ),
          
          hoverinfo = "text",
          
          marker = list(
            size = 6,
            opacity = 0.22
          ),
          
          name = "Below criteria",
          
          inherit = FALSE
        )
    }
    
    
    # --------------------------------------------------------
    # SIGNIFICANT
    # --------------------------------------------------------
    
    significant_data <- plot_data |>
      filter(
        !selected,
        pass_both
      )
    
    
    if (nrow(significant_data) > 0) {
      
      p <- p |>
        
        add_trace(
          
          data = significant_data,
          
          x = ~log2FoldChange,
          
          y = ~neg_log10_padj,
          
          type = "scatter",
          
          mode = "markers",
          
          text = ~paste0(
            "<b>",
            ifelse(
              is.na(gene_name),
              "Unknown",
              gene_name
            ),
            "</b><br>",
            "Gene ID: ",
            gene_id_clean,
            "<br>log2FC: ",
            round(
              log2FoldChange,
              3
            ),
            "<br>padj: ",
            signif(
              padj,
              3
            )
          ),
          
          hoverinfo = "text",
          
          marker = list(
            size = 7,
            opacity = 0.70
          ),
          
          name = "Meets criteria",
          
          inherit = FALSE
        )
    }
    
    
    # --------------------------------------------------------
    # CUTOFFS
    # --------------------------------------------------------
    
    padj_y <- -log10(
      max(
        padj_cutoff,
        .Machine$double.xmin
      )
    )
    
    
    p <- p |>
      
      layout(
        
        title = list(
          text = "Differential expression",
          font = list(
            size = 16
          )
        ),
        
        xaxis = list(
          title = "log2 Fold Change",
          zeroline = TRUE
        ),
        
        yaxis = list(
          title = "-log10 adjusted p-value"
        ),
        
        hovermode = "closest",
        
        margin = list(
          l = 65,
          r = 25,
          b = 55,
          t = 55
        ),
        
        shapes = list(
          
          list(
            type = "line",
            x0 = -Inf,
            x1 = Inf,
            y0 = padj_y,
            y1 = padj_y,
            line = list(
              dash = "dash",
              width = 1
            )
          ),
          
          list(
            type = "line",
            x0 = -log2fc_cutoff,
            x1 = -log2fc_cutoff,
            y0 = 0,
            y1 = Inf,
            line = list(
              dash = "dash",
              width = 1
            )
          ),
          
          list(
            type = "line",
            x0 = log2fc_cutoff,
            x1 = log2fc_cutoff,
            y0 = 0,
            y1 = Inf,
            line = list(
              dash = "dash",
              width = 1
            )
          )
        )
      )
    
    
    # --------------------------------------------------------
    # SELECTED GENE
    # --------------------------------------------------------
    
    selected_row <- plot_data |>
      filter(
        selected
      )
    
    
    if (nrow(selected_row) > 0) {
      
      selected_padj <- selected_row$padj[1]
      
      selected_log2fc <-
        selected_row$log2FoldChange[1]
      
      selected_name <-
        selected_row$gene_name[1]
      
      
      selected_padj_text <-
        
        if (
          is.na(selected_padj)
        ) {
          
          "NA"
          
        } else {
          
          signif(
            selected_padj,
            3
          )
        }
      
      
      selected_fc_text <-
        
        if (
          is.na(selected_log2fc)
        ) {
          
          "NA"
          
        } else {
          
          round(
            selected_log2fc,
            3
          )
        }
      
      
      p <- p |>
        
        add_trace(
          
          data = selected_row,
          
          x = ~log2FoldChange,
          
          y = ~neg_log10_padj,
          
          type = "scatter",
          
          mode = "markers+text",
          
          text = ~gene_name,
          
          textposition = "top center",
          
          marker = list(
            size = 16,
            symbol = "diamond",
            line = list(
              width = 2
            )
          ),
          
          hoverinfo = "text",
          
          hovertext = paste0(
            "<b>",
            selected_name,
            "</b><br>",
            "Gene ID: ",
            selected_id,
            "<br>log2FC: ",
            selected_fc_text,
            "<br>padj: ",
            selected_padj_text,
            "<br><b>Selected gene</b>"
          ),
          
          name = "Selected gene",
          
          inherit = FALSE
        )
    }
    
    
    p
  })
}


# ============================================================
# 23. RUN
# ============================================================

shinyApp(
  ui = ui,
  server = server
)