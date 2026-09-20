# Analysis inputs come exclusively from raw; the public reference is cached there.
tf_file <- file.path("raw", "carcinoma_vs_normal_significant_tfs.tsv")
tf_results <- if (file.exists(tf_file)) read.delim(tf_file, check.names = FALSE) else
  data.frame(TF = character(), logFC = numeric(), AveExpr = numeric(), t = numeric(),
             P.Value = numeric(), adj.P.Val = numeric())
stopifnot(all(c("TF", "logFC", "P.Value", "adj.P.Val") %in% names(tf_results)))
tf_results$TF <- trimws(as.character(tf_results$TF))
for (column in intersect(c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val"), names(tf_results))) {
  tf_results[[column]] <- suppressWarnings(as.numeric(tf_results[[column]]))
}
tf_results <- tf_results[!is.na(tf_results$TF) & nzchar(tf_results$TF) & !duplicated(tf_results$TF), ]
dorothea_revision <- "1461fb75e23e110c2281860526d4333920925282"
dorothea_url <- paste0("https://raw.githubusercontent.com/saezlab/dorothea/", dorothea_revision, "/data/dorothea_hs.rda")
dorothea_cache <- file.path("raw", "network", paste0("dorothea_hs_", dorothea_revision, ".rda"))
load_tf_reference <- function() {
  read_reference <- function(path) {
    env <- new.env(parent = emptyenv())
    load(path, envir = env)
    net <- as.data.frame(env$dorothea_hs)
    stopifnot(nrow(net) > 0, all(c("tf", "target", "mor", "confidence") %in% names(net)))
    net <- net[, c("tf", "target", "mor", "confidence")]
    net$mor <- as.numeric(net$mor)
    stopifnot(all(net$mor %in% c(-1, 1)))
    unique(net)
  }
  if (file.exists(dorothea_cache)) return(read_reference(dorothea_cache))
  dir.create(dirname(dorothea_cache), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = dirname(dorothea_cache), fileext = ".rda")
  on.exit(unlink(tmp), add = TRUE)
  download.file(dorothea_url, tmp, mode = "wb", quiet = TRUE)
  net <- read_reference(tmp)
  if (!file.rename(tmp, dorothea_cache)) stop("Could not save DoRothEA cache")
  net
}
