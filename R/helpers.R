# =========================
# NECorr Helper Functions
# =========================

#' DiscNULLrow
#' @description Remove rows from expression matrix where all values are zero.
#' @param eset Numeric matrix or data frame (genes × samples)
#' @return Filtered matrix/data frame with no all-zero rows
#' @export
DiscNULLrow <- function(eset) {
  eset[rowSums(eset) != 0, , drop = FALSE]
}

#' factOrg
#' @description Reorder expression columns by sample type and generate factor list.
#' Handles pseudo-replicates for samples with only one replicate.
#' @param eset Numeric matrix or data frame (genes × samples)
#' @param pDataTissue Character/factor vector of sample names (length = ncol(eset))
#' @return List with reordered expression matrix and factor vector
#' @export
# if there is only one replicate for a sample, we create two pseudo-replicates by duplicating it
# if there is 2 replicate create a third one that is the average of the two other replicates
factOrg <- function(eset, pDataTissue) {
  # Validate inputs
  if (!is.matrix(eset) && !is.data.frame(eset)) {
    stop("'eset' must be a matrix or data frame")
  }
  if (length(pDataTissue) != ncol(eset)) {
    stop("'pDataTissue' length must match the number of columns in 'eset'")
  }

  # Convert tissue labels to character for consistency
  pDataTissue <- as.character(pDataTissue)
  sample.names <- unique(pDataTissue)

  # Precompute counts for each tissue
  tissue_counts <- table(pDataTissue)

  # Warn if any tissue has fewer than 2 replicates
  low_rep <- names(tissue_counts)[tissue_counts < 2]
  if (length(low_rep) > 0) {
    warning(
      "The following tissue(s) have fewer than 2 replicates: ",
      paste(low_rep, collapse = ", "),
      "\nNECorr will duplicate columns to reach 3 replicates, ",
      "but statistical power will be limited."
    )
  }

  # Build index list: ensure each tissue has at least 3 entries
  idx_list <- unlist(lapply(sample.names, function(sn) {
    idx <- which(pDataTissue == sn)
    n <- length(idx)
    if (n == 1) {
      rep(idx, 3)                # singletons replicated 3 times
    } else if (n == 2) {
      c(idx, idx[1])              # pairs get an extra copy of the first
    } else {
      idx                         # keep as is for 3+ replicates
    }
  }), use.names = FALSE)

  # Build factor for grouping
  replicate_counts <- ifelse(tissue_counts == 1, 3,
                             ifelse(tissue_counts == 2, 3, tissue_counts))
  orgFactors <- factor(rep(names(tissue_counts), times = replicate_counts),
                       levels = sample.names)

  # Return reordered expression matrix + grouping factor
  list(
    eset    = as.matrix(eset)[, idx_list, drop = FALSE],
    orgFact = orgFactors
  )
}
