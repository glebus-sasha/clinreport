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

# Adjusted p-value cutoff
padj_cutoff <- 0.05

# Absolute log2 fold-change cutoff
log2fc_cutoff <- 1.0


# ============================================================
# 2. VALIDATION
# ============================================================

if (!file.exists(drug_network_file)) {
  stop(
    "Drug network file does not exist: ",
    drug_network_file
  )
}

if (!dir.exists(clinreport_dir)) {
  stop(
    "ClinReport directory does not exist: ",
    clinreport_dir
  )
}

if (!file.exists(gene_file)) {
  stop(
    "Gene expression file does not exist: ",
    gene_file
  )
}

if (
  !is.numeric(padj_cutoff) ||
  length(padj_cutoff) != 1 ||
  is.na(padj_cutoff) ||
  padj_cutoff <= 0 ||
  padj_cutoff >= 1
) {
  stop(
    "padj_cutoff must be a single number between 0 and 1."
  )
}

if (
  !is.numeric(log2fc_cutoff) ||
  length(log2fc_cutoff) != 1 ||
  is.na(log2fc_cutoff) ||
  log2fc_cutoff < 0
) {
  stop(
    "log2fc_cutoff must be a single non-negative number."
  )
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
    paste(
      missing_drug_columns,
      collapse = ", "
    )
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
    paste(
      missing_gene_columns,
      collapse = ", "
    )
  )
}


gene_data <- gene_data |>
  mutate(
    baseMean = suppressWarnings(
      as.numeric(baseMean)
    ),
    
    log2FoldChange = suppressWarnings(
      as.numeric(log2FoldChange)
    ),
    
    lfcSE = suppressWarnings(
      as.numeric(lfcSE)
    ),
    
    pvalue = suppressWarnings(
      as.numeric(pvalue)
    ),
    
    padj = suppressWarnings(
      as.numeric(padj)
    )
  )


gene_data <- gene_data |>
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
    filter(
      !is.na(drug_id)
    )
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


field_row <- function(
    label,
    value
) {
  
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
# 8. DRUG NAME
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
# 9. GENES FOR DRUG
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
        significant_fc
    ) |>
    arrange(
      is.na(padj),
      padj
    )
}


# ============================================================
# 10. DRUG — GENERAL
# ============================================================

render_general <- function(
    drug_id,
    data
) {
  
  chembl_id <- data$chembl$chembl_id
  
  molecule <- data$chembl$chembl_molecule
  
  
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
# 11. DRUG — INDICATIONS & MECHANISM
# ============================================================

render_pharmacology <- function(
    data
) {
  
  indications <-
    data$chembl$chembl_indications
  
  mechanisms <-
    data$chembl$chembl_mechanisms
  
  
  tagList(
    
    div(
      class = "section-title",
      "Indications"
    ),
    
    if (
      !is.null(indications) &&
      !is.null(
        indications$drug_indications
      )
    ) {
      
      map(
        indications$drug_indications,
        function(x) {
          
          div(
            class = "info-card",
            
            field_row(
              "EFO term",
              x$efo_term
            ),
            
            field_row(
              "EFO ID",
              x$efo_id
            ),
            
            field_row(
              "MeSH heading",
              x$mesh_heading
            ),
            
            field_row(
              "Maximum phase",
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
      !is.null(
        mechanisms$mechanisms
      )
    ) {
      
      map(
        mechanisms$mechanisms,
        function(x) {
          
          div(
            class = "info-card",
            
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
    }
  )
}


# ============================================================
# 12. DRUG — SAFETY
# ============================================================

render_safety <- function(
    data
) {
  
  warnings <-
    data$chembl$chembl_warnings
  
  molecule <-
    data$chembl$chembl_molecule
  
  
  tagList(
    
    div(
      class = "section-title",
      "Warnings"
    ),
    
    if (
      !is.null(warnings) &&
      !is.null(
        warnings$drug_warnings
      )
    ) {
      
      map(
        warnings$drug_warnings,
        function(x) {
          
          div(
            class = "info-card",
            
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
      !is.null(
        molecule$molecule_synonyms
      )
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
    }
  )
}


# ============================================================
# 13. TARGET CARDS
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
      
      
      div(
        
        class = card_class,
        
        onclick = onclick_js,
        
        
        div(
          class = "gene-main",
          
          div(
            class = "gene-symbol",
            
            ifelse(
              is.na(g$gene_name) ||
                !nzchar(g$gene_name),
              "Unknown",
              g$gene_name
            )
          ),
          
          div(
            class = "gene-id-small",
            g$gene_id
          )
        ),
        
        
        div(
          class = "gene-stats",
          
          div(
            class = "gene-stat",
            
            div(
              class = "gene-stat-label",
              "log2FC"
            ),
            
            div(
              class = "gene-stat-value",
              
              ifelse(
                is.na(
                  g$log2FoldChange
                ),
                "—",
                sprintf(
                  "%+.2f",
                  g$log2FoldChange
                )
              )
            )
          ),
          
          
          div(
            class = "gene-stat",
            
            div(
              class = "gene-stat-label",
              "padj"
            ),
            
            div(
              class = "gene-stat-value",
              
              ifelse(
                is.na(g$padj),
                "—",
                format.pval(
                  g$padj,
                  digits = 2,
                  eps = 1e-300
                )
              )
            )
          )
        ),
        
        
        div(
          class = "gene-cutoff-status",
          
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
# 14. DRUG PAGE
# ============================================================

render_drug_page <- function(
    drug_id,
    drug_name,
    data,
    genes
) {
  
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
              paste0(
                "Criteria · padj < ",
                padj_cutoff,
                " · |log2FC| ≥ ",
                log2fc_cutoff
              )
            )
          ),
          
          div(
            class = "targets-count",
            
            paste(
              nrow(genes),
              ifelse(
                nrow(genes) == 1,
                "gene",
                "genes"
              )
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
# 15. GENE PAGE
# ============================================================

render_gene_page <- function(
    gene
) {
  
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
        gene$gene_name
      ),
      
      div(
        class = "gene-id",
        gene$gene_id
      )
    ),
    
    
    div(
      class = "section-title",
      "Differential expression"
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
      class = "section-title",
      "Gene expression overview"
    ),
    
    plotlyOutput(
      "gene_volcano",
      height = "600px"
    )
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
        /* ==================================================
           GLOBAL
           ================================================== */

        body {
          background: #f4f5f7;
          font-family:
            -apple-system,
            BlinkMacSystemFont,
            'Segoe UI',
            sans-serif;
          color: #1f2937;
        }

        .app-title {
          font-size: 28px;
          font-weight: 600;
          letter-spacing: -0.3px;
          margin: 22px 0;
          color: #20252b;
        }

        .main-layout {
          display: flex;
          gap: 20px;
          align-items: stretch;
        }


        /* ==================================================
           LEFT PANEL
           ================================================== */

        .drug-list {
          width: 43%;
          background: white;
          border-radius: 12px;
          padding: 15px;
          box-shadow:
            0 2px 8px rgba(0,0,0,0.06);
        }


        /* ==================================================
           RIGHT PANEL
           ================================================== */

        .drug-details {
          flex: 1;
          background: white;
          border-radius: 12px;
          padding: 27px;
          min-height: 700px;
          box-shadow:
            0 2px 8px rgba(0,0,0,0.06);
        }


        /* ==================================================
           DRUG HEADER
           ================================================== */

        .drug-header {
          margin-bottom: 12px;
        }

        .drug-title {
          font-size: 30px;
          font-weight: 600;
          line-height: 1.2;
          letter-spacing: -0.4px;
          color: #20252b;
        }

        .drug-id {
          color: #9aa0a6;
          font-size: 13px;
          margin-top: 5px;
        }


        /* ==================================================
           TABS
           ================================================== */

        .nav-tabs {
          margin-top: 20px;
          border-bottom: 1px solid #e5e7eb;
        }

        .nav-tabs > li > a {
          color: #6b7280;
          border: none !important;
          padding: 10px 14px;
          font-size: 14px;
          transition:
            color 0.15s ease;
        }

        .nav-tabs > li > a:hover {
          background: transparent !important;
          color: #374151;
        }

        .nav-tabs > li.active > a,
        .nav-tabs > li.active > a:hover,
        .nav-tabs > li.active > a:focus {
          color: #1f2937;
          background: transparent !important;
          border: none !important;
          border-bottom: 2px solid #374151 !important;
          font-weight: 600;
        }

        .tab-content {
          padding-top: 10px;
        }


        /* ==================================================
           SECTION HEADERS
           ================================================== */

        .section-title {
          font-size: 17px;
          font-weight: 600;
          margin-top: 26px;
          margin-bottom: 12px;
          padding-bottom: 8px;
          border-bottom: 1px solid #eceef1;
          color: #30343b;
        }


        /* ==================================================
           FIELDS
           ================================================== */

        .field-row {
          display: flex;
          padding: 9px 0;
          border-bottom: 1px solid #f2f3f5;
        }

        .field-label {
          width: 32%;
          min-width: 180px;
          color: #737980;
          font-weight: 500;
          font-size: 13px;
        }

        .field-value {
          width: 68%;
          word-break: break-word;
          white-space: pre-wrap;
          color: #30343b;
          font-size: 13px;
        }


        /* ==================================================
           INFO CARDS
           ================================================== */

        .info-card {
          background: #fafbfc;
          border: 1px solid #e7e9ec;
          border-radius: 8px;
          padding: 10px 15px;
          margin-bottom: 12px;
        }


        /* ==================================================
           TARGET HEADER
           ================================================== */

        .targets-header {
          display: flex;
          align-items: flex-end;
          justify-content: space-between;
          margin-top: 10px;
          margin-bottom: 17px;
        }

        .targets-title {
          font-size: 21px;
          font-weight: 600;
          color: #2f343a;
        }

        .targets-criteria {
          margin-top: 4px;
          color: #9aa0a6;
          font-size: 11px;
          letter-spacing: 0.1px;
        }

        .targets-count {
          color: #9aa0a6;
          font-size: 13px;
          white-space: nowrap;
        }


        /* ==================================================
           TARGET LIST
           ================================================== */

        .targets-list {
          display: flex;
          flex-direction: column;
          gap: 8px;
        }


        /* ==================================================
           GENE CARD
           ================================================== */

        .gene-card {
          display: flex;
          align-items: center;
          padding: 13px 15px;
          background: #fafbfc;
          border: 1px solid #e6e8eb;
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
            0 4px 12px rgba(0,0,0,0.07);
          background: white;
          border-color: #d9dde3;
        }


        /* ==================================================
           DIRECTION
           ================================================== */

        .gene-card.up {
          border-left: 4px solid #4f86c6;
        }

        .gene-card.down {
          border-left: 4px solid #d96b6b;
        }

        .gene-card.neutral {
          border-left: 4px solid #aeb4bb;
        }


        /* ==================================================
           CUTOFF STATE
           ================================================== */

        .gene-card.passes-both {
          background: #f7faff;
          border-top-color: #d9e6f5;
          border-right-color: #d9e6f5;
          border-bottom-color: #d9e6f5;
        }

        .gene-card.passes-padj {
          background: #fafcff;
        }

        .gene-card.passes-fc {
          background: #fbfaff;
        }

        .gene-card.passes-neither {
          opacity: 0.62;
        }

        .gene-card.passes-both .gene-symbol {
          font-weight: 700;
        }


        /* ==================================================
           GENE MAIN
           ================================================== */

        .gene-main {
          flex: 1;
          min-width: 0;
        }

        .gene-symbol {
          font-size: 16px;
          font-weight: 600;
          color: #30343b;
        }

        .gene-id-small {
          color: #9aa0a6;
          font-size: 10px;
          margin-top: 3px;
        }


        /* ==================================================
           GENE STATS
           ================================================== */

        .gene-stats {
          display: flex;
          gap: 22px;
          margin-right: 18px;
        }

        .gene-stat {
          min-width: 70px;
          text-align: right;
        }

        .gene-stat-label {
          font-size: 10px;
          color: #a1a6ad;
          margin-bottom: 2px;
        }

        .gene-stat-value {
          font-size: 13px;
          font-weight: 600;
          color: #3f444b;
        }


        /* ==================================================
           CUTOFF BADGES
           ================================================== */

        .gene-cutoff-status {
          min-width: 95px;
          text-align: center;
          margin-right: 10px;
        }

        .cutoff-badge {
          display: inline-block;
          padding: 4px 7px;
          border-radius: 5px;
          font-size: 9px;
          font-weight: 600;
          white-space: nowrap;
          letter-spacing: 0.1px;
        }

        .cutoff-badge.both {
          background: #e8f1fb;
          color: #416f9f;
        }

        .cutoff-badge.padj-only {
          background: #edf6fb;
          color: #39738e;
        }

        .cutoff-badge.fc-only {
          background: #f1edfa;
          color: #69539a;
        }

        .cutoff-badge.neither {
          background: #f0f1f3;
          color: #858b92;
        }


        .gene-arrow {
          color: #b3b8be;
          font-size: 22px;
          width: 18px;
          text-align: center;
          transition:
            color 0.12s ease,
            transform 0.12s ease;
        }

        .gene-card:hover .gene-arrow {
          color: #777e86;
          transform: translateX(2px);
        }


        /* ==================================================
           GENE PAGE
           ================================================== */

        .back-button {
          margin-bottom: 20px;
          border-radius: 7px;
          color: #555;
          border-color: #ddd;
          background: white;
        }

        .back-button:hover {
          background: #f7f7f7;
        }

        .gene-title-large {
          font-size: 30px;
          font-weight: 600;
          letter-spacing: -0.4px;
          color: #20252b;
        }

        .gene-summary {
          display: flex;
          gap: 12px;
          margin: 15px 0 14px 0;
        }

        .summary-card {
          flex: 1;
          background: #fafbfc;
          border: 1px solid #e5e7eb;
          border-radius: 9px;
          padding: 15px;
        }

        .summary-label {
          color: #858b92;
          font-size: 11px;
          margin-bottom: 5px;
        }

        .summary-value {
          font-size: 20px;
          font-weight: 600;
          color: #30343b;
        }

        .gene-filter-status {
          display: flex;
          gap: 7px;
          margin: 5px 0 20px 0;
        }

        .gene-filter-badge {
          display: inline-block;
          padding: 4px 8px;
          border-radius: 5px;
          font-size: 10px;
          font-weight: 600;
        }

        .gene-filter-badge.pass {
          background: #edf6ef;
          color: #4b7855;
        }

        .gene-filter-badge.fail {
          background: #f1f2f3;
          color: #858b92;
        }


        /* ==================================================
           EMPTY
           ================================================== */

        .empty-message {
          color: #9aa0a6;
          padding: 40px;
          text-align: center;
          font-size: 13px;
        }


        /* ==================================================
           DATATABLE
           ================================================== */

        table.dataTable {
          font-size: 13px;
        }

        table.dataTable thead th {
          color: #737980;
          font-weight: 600;
          font-size: 12px;
          border-bottom: 1px solid #ddd !important;
        }

        table.dataTable tbody td {
          vertical-align: middle;
        }

        table.dataTable tbody tr.selected {
          background-color: #eef3f9 !important;
        }

        table.dataTable tbody tr:hover {
          background-color: #f8fafc !important;
        }


        /* ==================================================
           RESPONSIVE
           ================================================== */

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

          .gene-stats {
            gap: 10px;
          }

          .gene-cutoff-status {
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
        "right_panel"
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
  
  
  # ==========================================================
  # DRUG TABLE
  #
  # Drug ID intentionally hidden from UI.
  # It remains in drug_network and is used internally.
  # ==========================================================
  
  output$drug_table <- renderDT({
    
    table_data <- drug_network |>
      select(
        label,
        status,
        score
      ) |>
      mutate(
        ` ` = row_number(),
        .before = 1
      )
    
    
    datatable(
      
      table_data,
      
      selection = "single",
      
      rownames = FALSE,
      
      filter = "top",
      
      colnames = c(
        "",
        "Drug name",
        "Status",
        "Score"
      ),
      
      options = list(
        pageLength = 20,
        scrollX = TRUE,
        autoWidth = TRUE,
        
        columnDefs = list(
          list(
            width = "45px",
            targets = 0,
            className = "dt-center"
          ),
          list(
            width = "auto",
            targets = 1
          ),
          list(
            width = "100px",
            targets = 2
          ),
          list(
            width = "90px",
            targets = 3
          )
        )
      )
    )
  })
  
  
  # ==========================================================
  # SELECTED DRUG
  #
  # ВАЖНО:
  #
  # DataTable теперь показывает другую таблицу без drugId,
  # но selected row index всё равно соответствует строке
  # исходного drug_network.
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
    # Background genes
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
            opacity = 0.25
          ),
          
          name = "Below criteria",
          
          inherit = FALSE
        )
    }
    
    
    # --------------------------------------------------------
    # Genes passing both criteria
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
            opacity = 0.72
          ),
          
          name = "Meets criteria",
          
          inherit = FALSE
        )
    }
    
    
    # --------------------------------------------------------
    # Cutoff lines
    # --------------------------------------------------------
    
    padj_y <- -log10(
      max(
        padj_cutoff,
        .Machine$double.xmin
      )
    )
    
    
    p <- p |>
      
      layout(
        
        title = "Differential expression",
        
        xaxis = list(
          title = "log2 Fold Change",
          zeroline = TRUE
        ),
        
        yaxis = list(
          title = "-log10 adjusted p-value"
        ),
        
        hovermode = "closest",
        
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
    # SELECTED GENE — ALWAYS ON TOP
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
# 18. RUN
# ============================================================

shinyApp(
  ui = ui,
  server = server
)
