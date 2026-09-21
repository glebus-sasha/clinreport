# Analysis inputs stay in raw; the public resource is fetched into memory only.
tf_file <- file.path("raw", "carcinoma_vs_normal_significant_tfs.tsv")
empty_tf_results <- function() tibble(TF = character(), logFC = double(), AveExpr = double(),
  t = double(), P.Value = double(), adj.P.Val = double())
tf_input <- tryCatch({
  x <- read.delim(tf_file, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(all(names(empty_tf_results()) %in% names(x)), !anyDuplicated(x$TF),
            all(!is.na(x$TF) & nzchar(x$TF)))
  for (key in setdiff(names(empty_tf_results()), "TF")) {
    if (!is.numeric(x[[key]])) stop("Non-numeric TF column: ", key)
  }
  list(data = as_tibble(x), error = NULL)
}, error = function(e) list(data = empty_tf_results(), error = conditionMessage(e)))
tf_results <- tf_input$data
dorothea_url <- "https://raw.githubusercontent.com/saezlab/dorothea/master/data/dorothea_hs.rda"
empty_tf_regulons <- function() tibble(tf = character(), target = character(), mor = double(), confidence = character())
download_dorothea <- function(resource_url = dorothea_url) {
  tryCatch({
    old <- options(timeout = 60)
    on.exit(options(old), add = TRUE)
    con <- url(resource_url, open = "rb")
    on.exit(close(con), add = TRUE)
    bytes <- readBin(con, what = "raw", n = 50L * 1024L * 1024L)
    decoded <- rawConnection(memDecompress(bytes, type = "unknown"))
    on.exit(close(decoded), add = TRUE)
    env <- new.env(parent = emptyenv())
    load(decoded, envir = env)
    x <- as_tibble(env$dorothea_hs)
    stopifnot(all(c("tf", "target", "mor", "confidence") %in% names(x)))
    x <- x |> filter(confidence %in% c("A", "B", "C"), mor %in% c(-1, 1),
                      !is.na(tf), !is.na(target)) |>
      select(tf, target, mor, confidence) |> distinct()
    if (!nrow(x)) stop("DoRothEA returned no A/B/C interactions")
    list(data = x, error = NULL)
  }, error = function(e) list(data = empty_tf_regulons(), error = conditionMessage(e)))
}
dorothea_resource <- download_dorothea()
tf_regulons <- dorothea_resource$data
