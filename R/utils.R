# =========================
# Edge ranking and regulator ranking
# =========================

#' Rank network edges by the NECorr scores of their two genes
#'
#' @description Each edge receives the mean NECorr score of its two genes and
#' its percentile among all edges (1 = highest-scoring edge). This is a
#' ranking, not a statistical test.
#' @param network.int Data frame of edges (source in column 1, target in column 2).
#' @param gene.rank.vec Named numeric vector (or list) of gene scores.
#' @return Data frame with columns sourceIDs, targetIDs, edge_score, percentile.
#' @export
hub_edge_percentile <- function(network.int, gene.rank.vec) {
  network.int <- as.data.frame(network.int)
  stopifnot(!is.null(names(gene.rank.vec)))
  gene.rank.vec <- setNames(as.numeric(unlist(gene.rank.vec)), names(gene.rank.vec))
  src <- as.character(network.int[, 1]); tgt <- as.character(network.int[, 2])
  if (!all(c(src, tgt) %in% names(gene.rank.vec))) {
    stop("All genes in network.int must be present in names of gene.rank.vec")
  }
  edge_score <- unname((gene.rank.vec[src] + gene.rank.vec[tgt]) / 2)
  data.frame(sourceIDs = src, targetIDs = tgt, edge_score = edge_score,
             percentile = rank(edge_score, ties.method = "max") / length(edge_score),
             stringsAsFactors = FALSE)
}

#' @rdname hub_edge_percentile
#' @description \code{hub_edge_significance()} is the former name and is kept
#' for backward compatibility; the former column \code{p2} was a percentile and
#' is no longer reported as a p-value.
#' @export
hub_edge_significance <- function(network.int, gene.rank.vec) {
  .Deprecated("hub_edge_percentile")
  hub_edge_percentile(network.int, gene.rank.vec)
}

#' Decide whether an edge list is directed
#'
#' @description Heuristic used when \code{directed = NULL}. A regulatory
#' network (transcription factor -> target) has few genes acting as sources and
#' few reciprocal edges; an interaction network listed in arbitrary order does
#' not. Set \code{directed} explicitly whenever the network type is known.
#' @param network Data frame of edges (source, target).
#' @return A list with \code{directed} (logical), \code{source_fraction} and
#' \code{reciprocity}.
#' @export
detect_direction <- function(network) {
  network <- as.data.frame(network)
  src <- as.character(network[, 1]); tgt <- as.character(network[, 2])
  source_fraction <- length(unique(src)) / length(unique(c(src, tgt)))
  reciprocity <- mean(paste(tgt, src) %in% paste(src, tgt))
  list(directed = source_fraction <= 0.5 && reciprocity < 0.2,
       source_fraction = source_fraction, reciprocity = reciprocity)
}

#' Rank candidate regulators of the top-scoring genes
#'
#' @description For every gene with network partners, tests whether its
#' partners are enriched in top-scoring (hub) genes with a hypergeometric test,
#' which accounts for the number of partners. In a \strong{directed} network
#' partners are the gene's targets (source -> target), so enriched genes are
#' candidate upstream regulators of the hubs. In an \strong{undirected} network
#' (protein-protein interactions, yeast two-hybrid) direction is unknown:
#' partners are all neighbours and enriched genes are reported as
#' "associated", never as regulators.
#' @param network.int Data frame of edges (source, target).
#' @param scores Named numeric vector of NECorr gene scores.
#' @param Desc Optional data frame of gene descriptions with a \code{gene_id} column.
#' @param directed TRUE, FALSE, or NULL to decide with \code{detect_direction()}.
#' @param hub_fraction Fraction of top-scoring genes defining the hub set.
#' @param min_partners Minimum number of partners for a gene to be tested.
#' @param order_by "score" (default) orders the candidates by their own NECorr
#' score, which the benchmark on CRISPR dependencies favours; "enrichment"
#' orders them by the hypergeometric p-value of their partners.
#' @return Data frame with columns Gene, role, n_partners,
#' n_hub_partners, fold_enrichment, pvalue, FDR, mean_partner_score, Score,
#' followed by the description columns.
#' @export
rank_regulators <- function(network.int, scores, Desc = NULL, directed = NULL,
                            hub_fraction = 0.05, min_partners = 3,
                            order_by = c("score", "enrichment")) {
  order_by <- match.arg(order_by)
  ord <- function(d) if (order_by == "score") order(-d$Score, d$pvalue) else order(d$pvalue, -d$fold_enrichment)
  network.int <- unique(as.data.frame(network.int)[, 1:2])
  src <- as.character(network.int[, 1]); tgt <- as.character(network.int[, 2])
  keep <- src != tgt
  src <- src[keep]; tgt <- tgt[keep]
  if (is.null(directed)) directed <- detect_direction(network.int)$directed
  scores <- setNames(as.numeric(unlist(scores)), names(scores))
  genes <- names(scores)
  n_hub <- max(min(10, length(genes)), ceiling(hub_fraction * length(genes)))
  hubs <- names(sort(scores, decreasing = TRUE))[seq_len(min(n_hub, length(genes)))]
  if (directed) {
    partners <- split(tgt, src)
  } else {
    partners <- split(c(tgt, src), c(src, tgt))
  }
  partners <- lapply(partners, unique)
  partners <- partners[lengths(partners) >= min_partners]
  if (length(partners) == 0) return(NULL)
  N <- length(genes)
  res <- data.frame(
    Gene = names(partners),
    role = if (directed) "regulator" else "associated",
    n_partners = lengths(partners),
    n_hub_partners = vapply(partners, function(p) sum(p %in% hubs), numeric(1)),
    stringsAsFactors = FALSE, row.names = NULL)
  # universe excludes the tested gene itself
  K <- length(hubs) - as.integer(res$Gene %in% hubs)
  res$fold_enrichment <- (res$n_hub_partners / res$n_partners) / (K / (N - 1))
  res$pvalue <- stats::phyper(res$n_hub_partners - 1, K, (N - 1) - K, res$n_partners, lower.tail = FALSE)
  res$FDR <- stats::p.adjust(res$pvalue, method = "BH")
  res$mean_partner_score <- vapply(partners, function(p) mean(scores[p], na.rm = TRUE), numeric(1))
  res$Score <- unname(scores[res$Gene])
  res <- res[ord(res), , drop = FALSE]
  if (!is.null(Desc)) {
    Desc <- as.data.frame(Desc)
    if ("gene_id" %in% names(Desc)) {
      res <- merge(res, Desc, by.x = "Gene", by.y = "gene_id", all.x = TRUE, sort = FALSE)
      res <- res[ord(res), , drop = FALSE]
    }
  }
  rownames(res) <- NULL
  res
}

#' @rdname rank_regulators
#' @description \code{activator_significant()} is the former name, kept for
#' backward compatibility.
#' @param hub.int.significant Data frame of top edges from \code{hub_edge_percentile()}.
#' @export
activator_significant <- function(hub.int.significant, network.int, Desc = NULL,
                                  scores = NULL, directed = NULL) {
  .Deprecated("rank_regulators")
  if (is.null(scores)) {
    genes <- unique(c(as.character(as.data.frame(network.int)[, 1]), as.character(as.data.frame(network.int)[, 2])))
    scores <- setNames(rep(0, length(genes)), genes)
    top <- unique(c(hub.int.significant$sourceIDs, hub.int.significant$targetIDs))
    scores[top] <- 1
  }
  rank_regulators(network.int, scores, Desc, directed)
}

# =========================
# NECorr Utility Functions
# =========================

#' ScalN
#' @description Scales a numeric vector to a range from 0 to 1 using min-max scaling.
#' This is useful for normalising data so that all values fall within the same range.
#' @param x Numeric vector to be scaled.
#' @return A numeric vector of the same length as `x`, with values rescaled to fall between 0 and 1.
#' If all values in `x` are the same, returns a vector of zeros.
#' @export
ScalN <- function(x) {
  # Handle empty or NULL input
  if (is.null(x) || length(x) == 0) return(NULL)
  # Force numeric conversion
  x <- setNames(as.numeric(x), names(x))
  # Scale between 0 and 1
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) {
    scaled <- rep(0, length(x))  # avoid division by zero
  } else {
    scaled <- (x - rng[1]) / diff(rng)
  }
  # Preserve names if present
  if (!is.null(names(x))) {
    df <- data.frame(Gene = names(x), value = scaled, stringsAsFactors = FALSE)
  } else {
    stop("ScalN requires named vector so it can merge by Gene")
  }
  return(df)
}


#' scaling_param
#' @description Scale named numeric vector from 0 to 100
#' @param x Named numeric vector
#' @return Data frame with scaled values
#' @export
scaling_param <- function(x) {
  rng <- range(x[!is.infinite(x)], na.rm = TRUE)
  x[is.infinite(x) & sign(x) < 0] <- rng[1]
  x[is.infinite(x)] <- rng[2]
  data.frame(val = scales::rescale(x, to = c(0, 100)))
}

#' fishersMethod
#' @description Combine p-values using Fisher's method.
#' @param x Numeric vector of p-values
#' @return Combined p-value
#' @export
# Fisher's method to combine p-values
fishersMethod <- function(x) {
  x <- x[!is.na(x) & x > 0]  # remove NAs and zeros to avoid log(0)
  if (length(x) == 0) return(NA_real_)
  pchisq(-2 * sum(log(x)), df = 2 * length(x), lower.tail = FALSE)
}

#' cauchyCombine
#' @description Cauchy (ACAT) combination of p-values, valid under arbitrary dependence
#' (Liu & Xie 2020, JASA).
#' @param p Numeric vector of p-values
#' @return Combined p-value
#' @export
cauchyCombine <- function(p) {
  p <- p[!is.na(p)]
  if (length(p) == 0) return(NA_real_)
  p <- pmin(pmax(p, 1e-15), 1 - 1e-15)
  stat <- mean(tan((0.5 - p) * pi))
  0.5 - atan(stat) / pi
}

#' indexing.network
#' @description Replace gene names in network with their row indices in expression matrix
#' @param tab Expression matrix (rownames = genes)
#' @param network Two-column data frame of edges (source, target)
#' @return Matrix of same size as network with integer indices
#' @export
indexing.network <- function(tab, network) {
  cbind(match(network[, 1], rownames(tab)),
        match(network[, 2], rownames(tab)))
}

#' mergeMetricsByGene
#' @description Merge multiple gene metric data frames by 'Gene' column.
#' Handles NULL inputs and missing genes.
#' @param ... Data frames with 'Gene' column and metric columns
#' @return Merged data frame with all genes and metrics
#' @export
mergeMetricsByGene <- function(...) {
  dfs <- list(...)
  # Remove NULLs
  dfs <- Filter(Negate(is.null), dfs)
  if (length(dfs) == 0) {
    return(data.frame(Gene = character(0)))
  }
  # Check all have Gene column
  for (df in dfs) {
    if (!"Gene" %in% colnames(df)) {
      stop("All metric data frames must have a 'Gene' column")
    }
  }
  merged <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), dfs)
  return(merged)
}

#' safe_as_data_table
#' @description Safely convert a list to a data.table, padding shorter elements with NA if needed.
#' @param x List to convert
#' @param pad_with_NA Logical, if TRUE pads shorter elements with NA to match longest length
#' @return data.table
#' @export
#' @examples
#' #' # Example with equal lengths
#' lst_equal <- list(a = 1:3, b = letters[1:3])
#' safe_as_data_table(lst_equal)
#' #' # Example with differing lengths, padding with NA
#' lst_diff <- list(a = 1:3, b = letters[1:2])

safe_as_data_table <- function(x, pad_with_NA = TRUE) {
  # Check that input is a list
  if (!is.list(x)) {
    stop("Input must be a list.")
  }
  # Get lengths of each item
  lens <- lengths(x)
  # If all lengths are equal, convert directly
  if (length(unique(lens)) == 1) {
    return(as.data.table(x))
  }
  # If lengths differ
  message("Column lengths differ: ", paste(lens, collapse = ", "))
  if (!pad_with_NA) {
    stop("Lengths differ. Set pad_with_NA = TRUE to pad shorter columns with NA.")
  }
  # Pad shorter elements with NA to match the max length
  max_len <- max(lens)
  x <- lapply(x, function(col) {
    length(col) <- max_len
    col
  })
  return(as.data.table(x))
}
