#' Entropy Weight Calculation (C++)
#'
#' Internal function called by NECorr for entropy-based weighting of matrices.
#' @param mat A numeric matrix.
#' @return A numeric vector of weights.
#' @keywords internal
#' @export
entropy_weight_cpp <- function(mat) {
  .Call(`_NECorr_entropy_weight_cpp`, mat)
}

#' Multi-Correlation Calculation (C++)
#'
#' Internal function called by NECorr to compute multiple correlation metrics on network edges.
#' @param expression Numeric expression matrix.
#' @param ranks Numeric matrix of ranks.
#' @param src Integer vector of source indices.
#' @param tgt Integer vector of target indices.
#' @param bootstrapIterations Integer, number of bootstrap iterations.
#' @param usebestgcc Logical, whether to use the best correlation method.
#' @param asymmetricgcc Logical, whether to use asymmetric GCC.
#' @param rownames_in Logical, whether row names are provided.
#' @return A data frame with correlation values.
#' @export
#' @importFrom Rcpp evalCpp
multi_corr_necorr <- function(expression, ranks, src, tgt, bootstrapIterations, usebestgcc, asymmetricgcc, rownames_in) {
  .Call(`_NECorr_multi_corr_necorr`, expression, ranks, src, tgt, bootstrapIterations, usebestgcc, asymmetricgcc, rownames_in)
}

#'compute_TSI_TSE_parallel (C++)
#'
#' Internal function called by NECorr to compute Tissue Specificity Index (TSI) and Tissue Specific Expression (TSE) in parallel.
#' @param expression Numeric expression matrix.
#' @return A list containing TSI and TSE values.
#' @keywords internal
#' @export
#' @importFrom Rcpp evalCpp
#' @importFrom parallel detectCores

compute_TSI_TSE_parallel <- function(expression) {
  .Call(`_NECorr_compute_TSI_TSE_parallel`, expression)
}
