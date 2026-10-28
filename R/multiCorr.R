#' multiCorr
#' @description Calculate GCC, PCC, SCC, KCC correlations in one fast C++ parallel pass.
#' @param expression Numeric matrix (genes x samples).
#' @param edges Data frame with 'source' and 'target'.
#' @param methods Character vector of methods to return.
#' @param pernum Number of permutations for GCC significance.
#' @param useBestGCC Logical, whether to use the best GCC method.
#' @param asymmetricGCC Logical, whether to compute asymmetric GCC.
#' @param seed Integer seed for the GCC permutation test (results are reproducible).
#' @return Data frame with correlations and p-values.
#' @export
multiCorr <- function(expression, edges,
                      methods = c("GCC", "PCC", "SCC", "KCC"),
                      pernum = 0,
                      useBestGCC = FALSE, asymmetricGCC = FALSE,
                      seed = 12345) {

  if (is.data.frame(expression)) expression <- as.matrix(expression)
  edges <- data.table::as.data.table(edges)
  stopifnot(is.matrix(expression), all(c("source", "target") %in% colnames(edges)))

  gene_index <- setNames(seq_len(nrow(expression)), rownames(expression))
  src_idx <- gene_index[edges$source]
  tgt_idx <- gene_index[edges$target]

  # Strict filtering
  valid_edges <- !is.na(src_idx) & !is.na(tgt_idx) &
    src_idx >= 1 & src_idx <= nrow(expression) &
    tgt_idx >= 1 & tgt_idx <= nrow(expression)

  if (!any(valid_edges)) stop("No valid edges after filtering!, maybe your network genes don't match your genes in expression")

  # average ranks (ties handled), one row per gene
  ranks <- t(apply(expression, 1, rank, ties.method = "average"))
  storage.mode(ranks) <- "double"

  #wt <- 2 * (seq_len(ncol(expression))) - ncol(expression) - 1

  res_cpp <- multi_corr_necorr(
    expression, ranks,
    as.integer(src_idx[valid_edges] - 1),
    as.integer(tgt_idx[valid_edges] - 1),
    as.integer(pernum),
    useBestGCC,
    asymmetricGCC,
    rownames(edges)[valid_edges],  # Pass original rownames of valid edges
    as.numeric(seed),
    c("GCC", "PCC", "SCC", "KCC") %in% toupper(methods)   # only the requested coefficients are computed
  )
  # Convert to data.table and bind to edges
  cpp_dt <- safe_as_data_table(res_cpp)
  # Initialize with NA
  empty_results <- as.data.table(matrix(NA_real_, nrow = nrow(edges), ncol = ncol(cpp_dt)))
  setnames(empty_results, names(cpp_dt))
  # Fill only the valid edges
  empty_results[valid_edges, ] <- cpp_dt
  out <- data.table::as.data.table(cbind(edges, empty_results))

  # Filter for requested methods only
  keep_cols <- c("source", "target")
  for (m in methods) {
    keep_cols <- c(keep_cols, m, paste0(m, "_pvalue"))
    if (m == "GCC") keep_cols <- c(keep_cols, "GCC1", "GCC2")
  }
  return(out[, keep_cols, with = FALSE])
}
