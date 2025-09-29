#' fast_tsi_tse
#'
#' @description Fast computation of Tissue Specificity Index (TSI) and Tissue Specific Expression (TSE)
#' for genes across tissues, handling replicates and normalizing results. Collapse replicates,
#' run parallel TSI/TSE, and normalise from 0 to 1
#' @param expr Numeric expression matrix (genes x samples).
#' @param tissues Character vector of tissue types for each sample (column).
#' @return Data frame with columns: gene, TSI, TSE (normalized from 0 to 1).
#' @export
#' @importFrom parallel detectCores
#'
fast_tsi_tse <- function(expr, tissues) {
  expr <- as.matrix(expr)

  # Auto-detect and fix orientation if necessary
  if (length(tissues) == nrow(expr) && length(tissues) != ncol(expr)) {
    message("fast_tsi_tse: Detected genes in columns, transposing matrix...")
    expr <- t(expr)
  }

  stopifnot(length(tissues) == ncol(expr))

  gene_ids <- rownames(expr) # store before collapsing
  conds <- unique(tissues)

  # Collapse replicates: mean per tissue

  collapsed <- sapply(conds, function(cond) {
    mat <- expr[, tissues == cond, drop = FALSE]
    if (is.null(dim(mat))) {
      # only one replicate, return as-is
      as.numeric(mat)
    } else {
      rowMeans(mat)
    }
  })
  # Ensure correct orientation (genes X conditions)
  if (nrow(collapsed) != length(gene_ids)) {
    collapsed <- t(collapsed)
  }

  rownames(collapsed) <- gene_ids

  # Compute TSI/TSE in parallel
  res <- compute_TSI_TSE_parallel(collapsed)

  # Normalise between 0 and 1
  norm01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(rep(0, length(x)))
    (x - rng[1]) / diff(rng)
  }

  data.frame(
    gene = gene_ids,
    TSI = norm01(res$TSI),
    TSE = norm01(res$TSE),
    row.names = NULL
  )
}


# ========================
# Tissue Specificity Index
# ========================
#' TSI
#' @description Tissue Specificity Index (Yanai et al., 2005) for upregulated genes.
#' This measures how much a gene is expressed in one tissue compared to others.
#' @param x Numeric vector of expression values for a gene across tissues.
#' @return Numeric TSI value (higher = more tissue-specific)
#' @export
TSI <- function(x) {
  sum(1 - x / max(x)) / (length(x) - 1)
}

# ========================
# Tissue Specificity Exclusion Index
# ========================
#' TSE
#' @description Tissue Specificity Exclusion Index (Liseron et al., 2013) for downregulated genes.
#' This measures how much a gene is excluded in one tissue compared to others.
#' @param x Numeric vector of expression values for a gene across tissues.
#' @return Numeric TSE value
#' @export
TSE <- function(x) {
  sum(1 - (min(x) / max(x))) / (length(x) - 1)
}

# ========================
# IUTtsi: Filter tissue-specific genes
# ========================
#' IUTtsi
#' @description Define tissue-specific genes above an expression threshold.
#' Filters out very lowly expressed genes that could be artifacts.
#' @param data Expression matrix (genes x tissues).
#' @param geneTS Logical vector marking tissue-specific genes.
#' @param targetSample Name of the target tissue/condition.
#' @param expThreshold Minimum expression threshold for targetSample.
#' @param tissues Optional vector of tissue types for columns (for handling replicates).
#' @return Data frame of filtered genes with `tsi.order` as last column.
#' @export
IUTtsi <- function(data, geneTS, targetSample, expThreshold, tissues = NULL) {
  tryCatch({
    idx <- which(geneTS)
    Glen <- length(idx)

    if (Glen == 0) {
      return(data.frame(matrix(ncol = ncol(data) + 1, nrow = 0,
                               dimnames = list(NULL, c(colnames(data), "tsi.order")))))
    }

    tsi <- if (Glen > 1) {
      apply(data[idx, , drop = FALSE], 1, TSI)
    } else {
      setNames(TSI(data[idx, , drop = FALSE]), rownames(data)[idx])
    }

    tsi <- sort(tsi, decreasing = TRUE)
    y <- data[names(tsi), , drop = FALSE]

    # Handle replicate groups
    if (!is.null(tissues)) {
      group_cols <- which(tissues == targetSample)
      if (length(group_cols) == 0) {
        stop("IUTtsi: targetSample '", targetSample, "' not found in tissues vector")
      }
      keep <- rowMeans(y[, group_cols, drop = FALSE]) > expThreshold
    } else {
      if (!(targetSample %in% colnames(y))) {
        stop("IUTtsi: targetSample '", targetSample, "' not found in column names of expression data")
      }
      keep <- y[, targetSample] > expThreshold
    }

    y <- y[keep, , drop = FALSE]
    df <- cbind(y, tsi.order = tsi[rownames(y)])
    as.data.frame(df)
  }, error = function(e) {
    message("Error in IUTtsi: ", e$message)
    return(data.frame(matrix(ncol = ncol(data) + 1, nrow = 0,
                             dimnames = list(NULL, c(colnames(data), "tsi.order")))))
  })
}

# ========================
# IUTtse: Filter tissue-excluded genes
# ========================
#' IUTtse
#' @description Define tissue-excluded genes from the tissue/condition of interest.
#' @param data Expression matrix (genes x tissues).
#' @param geneTS Logical vector marking tissue-excluded genes.
#' @param targetSample Name of the target tissue/condition.
#' @param maxExpThreshold Maximum expression threshold for targetSample.
#' @param tissues Optional vector of tissue types for columns (for handling replicates).
#' @return Data frame with `tse.order` as last column.
#' @export
IUTtse <- function(data, geneTS, targetSample = NULL, maxExpThreshold = NULL, tissues = NULL) {
  tryCatch({
    idx <- which(geneTS)
    Glen <- length(idx)

    if (Glen == 0) {
      return(data.frame(matrix(ncol = ncol(data) + 1, nrow = 0,
                               dimnames = list(NULL, c(colnames(data), "tse.order")))))
    }

    tse <- if (Glen > 1) {
      apply(data[idx, , drop = FALSE], 1, TSE)
    } else {
      setNames(TSE(data[idx, , drop = FALSE]), rownames(data)[idx])
    }

    tse <- sort(tse, decreasing = TRUE)
    y <- data[names(tse), , drop = FALSE]

    # Handle replicate groups if filtering by target low-expression
    if (!is.null(maxExpThreshold) && !is.null(targetSample)) {
      if (!is.null(tissues)) {
        group_cols <- which(tissues == targetSample)
        if (length(group_cols) == 0) {
          stop("IUTtse: targetSample '", targetSample, "' not found in tissues vector")
        }
        keep <- rowMeans(y[, group_cols, drop = FALSE]) < maxExpThreshold
      } else {
        if (!(targetSample %in% colnames(y))) {
          stop("IUTtse: targetSample '", targetSample, "' not found in column names of expression data")
        }
        keep <- y[, targetSample] < maxExpThreshold
      }
      y <- y[keep, , drop = FALSE]
    }

    df <- cbind(y, tse.order = tse[rownames(y)])
    as.data.frame(df)
  }, error = function(e) {
    message("Error in IUTtse: ", e$message)
    return(data.frame(matrix(ncol = ncol(data) + 1, nrow = 0,
                             dimnames = list(NULL, c(colnames(data), "tse.order")))))
  })
}

# ========================
# IUTest: Intersection-Union Test
# Optimised IUTest() with optional pre-expanded matrix
# ========================
#' IUTest
#' @description Perform the Intersection-Union Test (IUT) for tissue specificity.
#' @param m.eset Expression matrix (genes x samples).
#' @param sFactors Vector of sample replicate counts.
#' @param sIndex Index of the target sample in the factor list.
#' @param alpha Significance level.
#' @param preExpanded Logical indicating if m.eset is already expanded for replicates.
#' @param adj_sFactors Adjusted sample factors (for pre-expanded data).
#' @return List with TS (tissue-specific) and TE (tissue-excluded) logical vectors.
#' @export
IUTest <- function(m.eset, sFactors, sIndex, alpha,
                   preExpanded = FALSE, adj_sFactors = NULL) {
  tryCatch({
    # Validate inputs
    m.eset <- as.matrix(m.eset)
    gene <- nrow(m.eset)
    nrarrays <- ncol(m.eset)
    # Check if pre-expanded and adj_sFactors provided
    if (!preExpanded) {
      nbTreat <- length(sFactors)
      if (nbTreat < 2) {
        warning("IUTest requires at least 2 conditions")
        return(list(TS = rep(FALSE, gene), TE = rep(FALSE, gene)))
      }
      if (sum(sFactors) != nrarrays) {
        stop("Sum of sFactors does not match ncol(m.eset)")
      }
      if (sIndex < 1 || sIndex > nbTreat) {
        stop("Invalid sIndex: must be between 1 and ", nbTreat)
      }
      # Expand m.eset for replicates
      idx_list <- unlist(lapply(seq_len(nbTreat), function(i) {
        start_col <- sum(sFactors[seq_len(i - 1)]) + 1
        end_col <- start_col + sFactors[i] - 1
        idx <- start_col:end_col
        if (sFactors[i] == 1) rep(idx, 3)
        else if (sFactors[i] == 2) c(idx, idx[1])
        else idx
      }), use.names = FALSE)
      adj_sFactors <- ifelse(sFactors < 3, 3, sFactors)
      m.eset <- m.eset[, idx_list, drop = FALSE]
      nrarrays <- ncol(m.eset)
    }
    # Now m.eset is expanded, and adj_sFactors is set
    nbTreat <- length(adj_sFactors)
    # Calculate means and sum of squares
    meanTab <- matrix(0, nrow = gene, ncol = nbTreat)
    SSdevs <- numeric(gene)
    fL <- 1
    for (i in seq_len(nbTreat)) {
      upper <- fL + adj_sFactors[i] - 1
      m.eset.i <- m.eset[, fL:upper, drop = FALSE]
      meanTab[, i] <- rowMeans(m.eset.i)
      SSdevs <- SSdevs + rowSums((m.eset.i - meanTab[, i])^2)
      fL <- upper + 1
    }
    # Perform IUT
    Nin <- adj_sFactors[sIndex]
    Noth <- adj_sFactors[-sIndex]
    meanSamp <- as.matrix(meanTab[, sIndex])
    meanAll <- as.matrix(meanTab[, -sIndex, drop = FALSE])
    # Compute standard deviations and t-statistics
    SD <- SSdevs / (nrarrays - nbTreat)
    # Avoid zero SD
    SD[SD == 0] <- 0.05
    # Matrix of standard errors
    SE <- sapply(Noth, function(noth) sqrt(SD * (1 / Nin + 1 / noth)))
    meanSamp_rep <- matrix(meanSamp, nrow = nrow(meanSamp), ncol = length(Noth))
    tstattab <- (meanSamp_rep - meanAll) / SE
    # Determine critical t-value
    # thresh <- 1 - (1 - alpha)^(1 / gene)
    #crit.value <- qt(1 - thresh, nrarrays - adj_sFactors[sIndex])
    # Degrees of freedom
    df <- nrarrays - adj_sFactors[sIndex]
    # Calculate p-values for all comparisons
    pvals <- 2 * pt(-abs(tstattab), df = df)
    # Adjust using BH FDR
    pvals_adj <- matrix(p.adjust(as.vector(pvals), method = "BH"),
                        nrow = nrow(pvals), ncol = ncol(pvals))
    # Proportion threshold (e.g., 0.8 means 80% of comparisons must be sig.)
    prop_required <- 0.8
    # TS: significantly higher in target tissue
    resultpos <- rowSums(pvals_adj < alpha & tstattab > 0) >= prop_required * (nbTreat - 1)
    # TE: significantly lower in target tissue
    resultmin <- rowSums(pvals_adj < alpha & tstattab < 0) >= prop_required * (nbTreat - 1)
    # DEBUG
    # cat("Max t-stat:", max(tstattab, na.rm = TRUE), "\n")
    # cat("TS hits:", sum(resultpos), "\n")
    # cat("TE hits:", sum(resultmin), "\n")
    list(TS = resultpos, TE = resultmin)
  }, error = function(e) {
    message("Error in IUTest: ", e$message)
    list(TS = rep(FALSE, nrow(m.eset)), TE = rep(FALSE, nrow(m.eset)))
  })
}

# ========================
# ts.IUT: Adaptive thresholding for NECorr
# ========================
#' ts.IUT
#' @description Combine TSI/TSE tests and IUT to determine tissue-specific and tissue-excluded genes.
#' @param name Analysis name.
#' @param eset Expression matrix.
#' @param tissues Vector of tissue types for columns.
#' @param target Target tissue.
#' @param threshold p-value threshold for IUT.
#' @param filter Minimum expression threshold.
#' @param minGenes Minimum number of genes to find before stopping.
#' @param maxIter Maximum number of iterations to relax criteria.
#' @param loosenFactor Factor to increase p-value threshold each iteration.
#' @param filterDecrease Factor to decrease expression threshold each iteration.
#' @details This function iteratively applies the Intersection-Union Test (IUT)
#' to identify tissue-specific (TS) and tissue-excluded (TE) genes.
#' It starts with stringent criteria and relaxes them if not enough genes are found.
#' The process continues until at least `minGenes` are found or `maxIter` is reached.
#' The function returns filtered TS and TE gene lists.
#' @return List with filtTSI and filtTSE matrices.
#' @export
# Integrated ts.IUT with adaptive threshold/filter and fast TSI/TSE
ts.IUT <- function(name,
                   eset,
                   tissues,
                   target,
                   threshold = 0.05,
                   filter = 1,
                   minGenes = 1,
                   maxIter = 5,
                   loosenFactor = 2,
                   filterDecrease = 0.5) {
  eset <- as.matrix(eset)
  if (length(tissues) != ncol(eset)) {
    stop("'tissues' length must match ncol(eset)")
  }

  sample.names <- unique(tissues)
  sFactors <- table(tissues)
  sIndex <- match(target, sample.names)
  if (is.na(sIndex)) {
    stop("Target condition '", target, "' not found in tissues vector")
  }

  nbTreat <- length(sFactors)
  idx_list <- unlist(lapply(seq_len(nbTreat), function(i) {
    start_col <- sum(sFactors[seq_len(i - 1)]) + 1
    end_col <- start_col + sFactors[i] - 1
    idx <- start_col:end_col
    if (sFactors[i] == 1) rep(idx, 3)
    else if (sFactors[i] == 2) c(idx, idx[1])
    else idx
  }), use.names = FALSE)
  adj_sFactors <- ifelse(sFactors < 3, 3, sFactors)
  eset_expanded <- eset[, idx_list, drop = FALSE]

  iter <- 1
  found <- FALSE
  finalThreshold <- threshold
  finalFilter <- filter
  filtTSI <- NULL
  filtTSE <- NULL

  while (iter <= maxIter && !found) {
    #message(sprintf("ts.IUT: Iteration %d - threshold=%.5f, filter=%.3f",
    #                iter, finalThreshold, finalFilter))

    t.iut <- IUTest(eset_expanded,
                    sFactors = NULL,
                    sIndex = sIndex,
                    alpha = finalThreshold,
                    preExpanded = TRUE,
                    adj_sFactors = adj_sFactors) ######

    ts_scores <- fast_tsi_tse(eset, tissues)

    group_cols <- which(tissues == target)
    ts_keep <- rowMeans(eset[, group_cols, drop = FALSE]) > finalFilter &
      t.iut$TS
    filtTSI <- ts_scores[ts_keep, c("gene", "TSI"), drop = FALSE]

    te_keep <- rowMeans(eset[, group_cols, drop = FALSE]) < finalFilter &
      t.iut$TE
    filtTSE <- ts_scores[te_keep, c("gene", "TSE"), drop = FALSE]

    # message("ts.IUT: Found ", nrow(filtTSI), " TS and ", nrow(filtTSE), " TE after filtering.")

    if (nrow(filtTSI) >= minGenes || nrow(filtTSE) >= minGenes) {
      found <- TRUE
    } else {
      finalThreshold <- min(1, finalThreshold * loosenFactor)
      finalFilter <- max(0, finalFilter * filterDecrease)
      iter <- iter + 1
    }
  }

  if (!found) {
    # message("ts.IUT: No TS/TE found after ", maxIter, " iterations.")
  } else {
    # message(sprintf("ts.IUT: Final parameters - threshold=%.5f, filter=%.3f",
    #                finalThreshold, finalFilter))
  }

  return(list(
    filtTSI = filtTSI,
    filtTSE = filtTSE,
    finalThreshold = finalThreshold,
    finalFilter = finalFilter
  ))
}
