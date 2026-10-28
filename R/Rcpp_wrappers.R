#' Entropy weight calculation (C++)
#'
#' @param mat A numeric matrix.
#' @return A numeric vector of weights.
#' @keywords internal
#' @name entropy_weight_cpp
NULL

#' Edge correlations (C++)
#'
#' Computes GCC, PCC, SCC and KCC for each network edge. Use
#' \code{\link{multiCorr}} instead of calling this directly.
#' @param expression Numeric expression matrix (genes x samples).
#' @param ranks Numeric matrix of average ranks, same shape as \code{expression}.
#' @param src,tgt Zero-based row indices of the two genes of each edge.
#' @param bootstrapIterations Number of permutations for the GCC p-value.
#' @param useBestGCC Not used.
#' @param asymmetricGCC Logical, use the asymmetric Gini correlation.
#' @param rownames_in Optional row names for the output.
#' @param seed Seed of the permutation test.
#' @param which Logical vector of length 4: compute GCC, PCC, SCC, KCC.
#' @return A data frame of correlations and p-values.
#' @keywords internal
#' @name multi_corr_necorr
NULL

#' Tissue specificity and exclusion indices (C++)
#'
#' @param expr Numeric matrix (genes x conditions) on the linear scale.
#' @return A list with TSI and TSE.
#' @keywords internal
#' @name compute_TSI_TSE_parallel
NULL
