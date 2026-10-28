# =========================
# Criterion weighting
# =========================

#' Default AHP pairwise-comparison matrix for the NECorr criteria
#'
#' @description Saaty pairwise judgements (1 = equal importance, 3 = moderately
#' more important, 5 = strongly more important; intermediate values allowed).
#' Each criterion has an importance level and the judgement between two
#' criteria is their level difference plus one.
#' The default levels (\code{preset = "validated"}) were optimised on human
#' cancer cell lines, where genome-wide CRISPR screens give a complete ground
#' truth, by leave-one-cancer-type-out cross-validation: degree is most
#' important (weight 0.35), condition specificity and differential expression
#' follow (0.21 each), and the two co-expression criteria and bridging carry
#' 0.08 each. \code{preset = "wiring"} gives the judgements of NECorr 3.2,
#' selected on simulated data, in which degree carries no weight.
#' Entry \code{[i, j]} is the importance of criterion i relative to criterion j.
#' @param criteria Character vector of criteria to include (subset of the
#' default set). Criteria absent from a run are dropped and the weights
#' recomputed from the remaining judgements.
#' @param preset "validated" (default) or "wiring".
#' @return A reciprocal pairwise-comparison matrix.
#' @export
#' @examples
#' round(ahp_weights(necorr_ahp_matrix())$weights, 2)
#' ahp_weights(necorr_ahp_matrix(c("Coexpr", "TS", "CondCoexpr"), preset = "wiring"))$weights
necorr_ahp_matrix <- function(criteria = NULL, preset = c("validated", "wiring")) {
  preset <- match.arg(preset)
  crit <- c("Coexpr", "CondCoexpr", "TS", "DE", "Bridging", "Degree")
  # importance level per criterion (higher = more important)
  # level difference d between two criteria -> Saaty judgement d + 1
  level <- switch(preset,
    validated = c(Coexpr = 0, CondCoexpr = 0, TS = 2, DE = 2, Bridging = 0, Degree = 3),
    wiring = c(Coexpr = 2, CondCoexpr = 1, TS = 2, DE = 1, Bridging = 0, Degree = 0))
  A <- outer(level[crit], level[crit], function(a, b) {
    d <- a - b
    ifelse(d >= 0, d + 1, 1 / (-d + 1))
  })
  dimnames(A) <- list(crit, crit)
  if (!is.null(criteria)) {
    miss <- setdiff(criteria, crit)
    if (length(miss)) stop("Unknown criteria: ", paste(miss, collapse = ", "))
    A <- A[criteria, criteria, drop = FALSE]
  }
  A
}

#' Analytic Hierarchy Process weights
#'
#' @description Computes criterion weights as the principal eigenvector of a
#' Saaty pairwise-comparison matrix and reports the consistency ratio (CR).
#' A CR above 0.1 indicates contradictory judgements.
#' @param A Square reciprocal pairwise-comparison matrix with dimnames.
#' @return A list with \code{weights} (named, sum to 1), \code{lambda_max},
#' \code{CI} (consistency index) and \code{CR} (consistency ratio).
#' @export
#' @examples
#' A <- matrix(c(1, 3, 5,  1/3, 1, 2,  1/5, 1/2, 1), 3, byrow = TRUE,
#'             dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
#' ahp_weights(A)
ahp_weights <- function(A) {
  A <- as.matrix(A)
  n <- nrow(A)
  if (n != ncol(A)) stop("Pairwise matrix must be square")
  if (any(A <= 0)) stop("Pairwise judgements must be positive")
  if (max(abs(A * t(A) - 1)) > 1e-6) stop("Pairwise matrix must be reciprocal: A[j, i] == 1 / A[i, j]")
  if (n == 1) return(list(weights = setNames(1, rownames(A)), lambda_max = 1, CI = 0, CR = 0))
  e <- eigen(A)
  k <- which.max(Re(e$values))
  w <- abs(Re(e$vectors[, k]))
  w <- w / sum(w)
  names(w) <- rownames(A)
  lambda <- Re(e$values[k])
  CI <- (lambda - n) / (n - 1)
  # Saaty random consistency index
  RI <- c(0, 0, 0.58, 0.90, 1.12, 1.24, 1.32, 1.41, 1.45, 1.49)
  ri <- if (n <= length(RI)) RI[n] else 1.98 * (n - 2) / n
  CR <- if (ri > 0) CI / ri else 0
  list(weights = w, lambda_max = lambda, CI = CI, CR = CR)
}

#' CRITIC weights
#'
#' @description CRiteria Importance Through Intercriteria Correlation
#' (Diakoulaki et al. 1995). A criterion receives more weight when it varies
#' across genes and is not redundant with the other criteria, so several
#' correlated criteria (for example centralities) share weight instead of
#' accumulating it.
#' @param P Numeric matrix or data frame (genes x criteria), ideally scaled 0-1.
#' @return Named numeric vector of weights summing to 1.
#' @export
critic_weight <- function(P) {
  P <- as.matrix(P)
  s <- apply(P, 2, stats::sd)
  R <- suppressWarnings(stats::cor(P))
  R[is.na(R)] <- 0
  diag(R) <- 1
  info <- s * rowSums(1 - R)
  if (sum(info) == 0) return(setNames(rep(1 / ncol(P), ncol(P)), colnames(P)))
  setNames(info / sum(info), colnames(P))
}

# Resolve the `weights` argument of NECorr() into a named weight vector.
# Returns list(weights, scheme, CR).
resolve_weights <- function(weights, P) {
  crit <- colnames(P)
  CR <- NA_real_
  if (is.matrix(weights)) {
    if (is.null(rownames(weights)) || !all(crit %in% rownames(weights))) {
      stop("Pairwise matrix must have dimnames covering: ", paste(crit, collapse = ", "))
    }
    a <- ahp_weights(weights[crit, crit, drop = FALSE])
    w <- a$weights; CR <- a$CR; scheme <- "ahp (user matrix)"
  } else if (is.numeric(weights)) {
    if (is.null(names(weights)) || !all(crit %in% names(weights))) {
      stop("Numeric weights must be named and cover: ", paste(crit, collapse = ", "))
    }
    w <- weights[crit] / sum(weights[crit]); scheme <- "custom"
  } else {
    scheme <- match.arg(tolower(weights), c("ahp", "critic", "equal", "entropy"))
    w <- switch(scheme,
      ahp = { a <- ahp_weights(necorr_ahp_matrix(crit)); CR <- a$CR; a$weights },
      critic = critic_weight(P),
      equal = setNames(rep(1 / length(crit), length(crit)), crit),
      entropy = setNames(as.numeric(entropy_weight(P)), crit))
  }
  if (!is.na(CR) && CR > 0.1) {
    warning(sprintf("AHP consistency ratio is %.2f (> 0.10): the pairwise judgements contradict each other.", CR))
  }
  list(weights = w[crit], scheme = scheme, CR = CR)
}
