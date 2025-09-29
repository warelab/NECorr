#' DE.ranking: Differential expression ranking for NECorr
#' @description Uses limma to rank genes in the network by differential expression.
#' @param exps Expression matrix (genes × samples)
#' @param GeneList Vector of gene IDs to include
#' @param factortab Factor vector or data frame of sample groups
#' @param sample.l Condition of interest
#' @param sample.names Vector of all sample group names
#' @return Named vector of scaled DE ranks
#' @export
DE.ranking <- function(exps, GeneList, factortab, sample.l, sample.names) {
  #requireNamespace("limma")
  #requireNamespace("Biobase")

  # Select relevant genes (matrix subsetting is fast)
  #sel.Genes <- intersect(rownames(exps), GeneList)
  sel.idx <- match(GeneList, rownames(exps), nomatch = 0)
  sel.Genes <- rownames(exps)[sel.idx[sel.idx > 0]]
  expr <- exps[sel.Genes, , drop = FALSE]

  # Build design matrix directly
  f <- factor(factortab, levels = sample.names)
  design <- model.matrix(~ 0 + f)
  colnames(design) <- sample.names

  # Fit linear model directly on the matrix (no ExpressionSet needed)
  fit <- limma::lmFit(expr, design)
  #eset <- ExpressionSet(assayData = expr)
  #fit <- lmFit(eset, design)

  # Build contrast matrix for sample.l vs all other samples
  contrast <- paste0(sample.l, " - ((", paste(sample.names[sample.names != sample.l], collapse = "+"),
                     ")/", length(sample.names) - 1, ")")
  cont.mat <- limma::makeContrasts(contrasts = contrast, levels = design)

  # Apply contrasts and empirical Bayes smoothing
  fit2 <- limma::contrasts.fit(fit, cont.mat)
  fit2 <- suppressWarnings(limma::eBayes(fit2))

  # Extract results
  de <- limma::topTable(fit2, coef = 1, adjust.method = "fdr", number = length(sel.Genes))

  # Rank
  ranks <- ScalN(setNames(de$B, rownames(de)))
  colnames(ranks)[2] <- "DE_rank"
  return(ranks)
}
