#' Clean an expression matrix before correlation
#'
#' @description Handles non-finite values without distorting correlations.
#' \code{-Inf} (typically \code{log2(0)}) is set to the lowest finite value of
#' the matrix and \code{+Inf} to the highest, i.e. values are capped at the
#' observed range instead of being replaced by an arbitrary huge number.
#' Missing values are replaced by the gene's median. Genes with more than
#' \code{max_missing} non-finite values, or with zero variance, are removed
#' because no correlation can be computed for them.
#' @param mat Numeric matrix (genes x samples) with row names.
#' @param max_missing Maximum fraction of non-finite values tolerated per gene.
#' @return A list with the cleaned matrix \code{mat} and a named integer vector
#' \code{report} (counts of capped, imputed and removed entries).
#' @export
#' @examples
#' m <- matrix(c(1, 2, -Inf, 4, 5, 6, NA, 8), 2, byrow = TRUE,
#'             dimnames = list(c("g1", "g2"), paste0("s", 1:4)))
#' sanitize_expression(m, max_missing = 0.5)
sanitize_expression <- function(mat, max_missing = 0.2) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  nonfinite <- !is.finite(mat)
  n_neg <- sum(mat == -Inf, na.rm = TRUE)
  n_pos <- sum(mat == Inf, na.rm = TRUE)
  n_na <- sum(is.na(mat))
  drop_missing <- rowMeans(nonfinite) > max_missing
  rng <- if (any(!nonfinite)) range(mat[!nonfinite]) else c(0, 0)
  mat[which(mat == -Inf)] <- rng[1]
  mat[which(mat == Inf)] <- rng[2]
  if (anyNA(mat)) {
    med <- apply(mat, 1, stats::median, na.rm = TRUE)
    idx <- which(is.na(mat), arr.ind = TRUE)
    mat[idx] <- med[idx[, 1]]
  }
  v <- apply(mat, 1, stats::var)
  drop_const <- is.na(v) | v == 0
  keep <- !(drop_missing | drop_const)
  list(mat = mat[keep, , drop = FALSE],
       report = c(neg_inf_capped = n_neg, pos_inf_capped = n_pos,
                  missing_imputed = n_na,
                  genes_removed_missing = sum(drop_missing),
                  genes_removed_constant = sum(drop_const & !drop_missing)))
}
