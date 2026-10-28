#' NECorr: rank genes of a molecular network for one condition
#'
#' @author Christophe Liseron-Monfils, Andrew Olson
#' @description NECorr ranks the genes of a molecular network for their
#' importance in one condition (tissue, treatment, stage) from criteria
#' computed on expression data and on the network:
#' \describe{
#'   \item{DE}{differential expression of the target condition against the
#'     mean of the others (limma moderated t)}
#'   \item{Degree}{number of network partners}
#'   \item{TS}{selective expression in the target condition
#'     (see \code{\link{condition_specificity}})}
#'   \item{Coexpr}{how much more the gene is co-expressed with its network
#'     partners than with random genes (see \code{\link{gene_coexpression}})}
#'   \item{CondCoexpr}{the same for the part of the co-expression that comes
#'     from the target-condition samples}
#'   \item{Bridging}{betweenness centrality beyond what the gene's degree
#'     predicts}
#' }
#' The default (NECorr 3.4) is the rank product of two criteria: up-regulation
#' in the target condition (DE with \code{direction = "up"}) and Degree, the
#' number of targets of the gene in a directed network. It is
#' the scheme that recovered most lineage-selective CRISPR dependencies in the
#' human benchmark; a supervised model given all criteria did not improve on
#' it. The six-criteria scheme of NECorr 3.3 is obtained with
#' \code{criteria = c("Coexpr", "CondCoexpr", "TS", "DE", "Bridging", "Degree")},
#' \code{weights = "ahp"}, \code{aggregation = "sum"}, \code{direction = "both"}
#' and \code{scaling = "softcap"}.
#' @param networkFile File path or data.frame/data.table of edges: source in the
#' first column, target in the second (a three-column file is read as
#' source, interaction type, target).
#' @param expression File path or data.frame/data.table: gene IDs in the first
#' column, one column per sample. Log2 scale expected (see \code{log_scale}).
#' @param description.file File path or data.frame/data.table of gene
#' descriptions, gene IDs in the first column.
#' @param condition Name of the condition to rank genes for.
#' @param meta_data File path or data.frame/data.table giving the condition of
#' each sample: sample IDs in the first column (or as row names) and the
#' condition in the next column.
#' @param permutation Number of permutations for the GCC p-value (0 = none).
#' @param method Correlation method(s): any of "GCC", "PCC", "SCC", "KCC".
#' @param visualize TRUE/FALSE, build the summary plots.
#' @param output_dir Directory for results and plots (required if \code{save_results}).
#' @param interactive_net TRUE/FALSE, return an interactive network (needs visNetwork).
#' @param top_n Number of top genes shown in plots.
#' @param save_results TRUE/FALSE, write result tables as CSV.
#' @param useBestGCC Kept for backward compatibility; not used.
#' @param asymmetricGCC TRUE/FALSE, use the asymmetric Gini correlation.
#' @param weights Criterion weighting: "equal" (default), "ahp", "critic",
#' "entropy", a named numeric vector, or a Saaty pairwise-comparison matrix
#' with criteria as dimnames.
#' @param directed Is the network directed (source regulates target)? TRUE,
#' FALSE, or NULL to decide with \code{\link{detect_direction}}. Direction is
#' used to rank regulators; hub criteria use the undirected network.
#' @param criteria Criteria to combine, among Coexpr, CondCoexpr, TS, DE,
#' Bridging and Degree. Default: \code{c("DE", "Degree")}, the pair validated
#' on human CRISPR dependencies (NECorr 3.4). "Degree" makes the score increase
#' with connectivity; leave it out for a ranking that does not depend on how
#' many partners a gene has.
#' @param direction "up" (default): DE and TS count only genes expressed more
#' in the target condition than elsewhere (positive moderated t), because a
#' gene must be expressed in a condition to be required there. "both": absolute
#' statistics, so that genes excluded from the condition also score (NECorr 3.3).
#' @param scaling "rank" (default): each criterion is replaced by rank / n,
#' which with \code{aggregation = "product"} gives a rank product. "softcap":
#' 0-1 scaling with a bounded soft cap above the 99th percentile (NECorr 3.3).
#' @param degree_type "targets" (default): in a directed network the Degree
#' criterion is the number of targets of the gene (out-degree), so that genes
#' without targets share the lowest degree rank. "total": the number of
#' partners, regulators and targets together; this is always used for
#' undirected networks.
#' @param aggregation "product" (default): weighted geometric mean, so a gene
#' must do well on all criteria together. "sum": weighted sum, where a high value on one
#' criterion compensates for a low value on another.
#' "waspas": the mean of the two.
#' @param epsilon Small constant added to the 0-1 criteria before the product.
#' @param coexpr_stat,tail_q How co-expression with partners is summarised, see
#' \code{\link{gene_coexpression}}.
#' @param n_background Number of random non-partner genes used as each gene's
#' co-expression baseline.
#' @param full_netstat TRUE to also compute betweenness for the \code{netstat}
#' table when the Bridging criterion is not used (slow on large networks).
#' @param edge_top Fraction of top-ranked edges returned in \code{necorrEdges}.
#' @param n_null Number of null permutations for empirical gene p-values (0 = none).
#' Expression profiles are shuffled among the network genes, which keeps the
#' network and the expression data intact but breaks their correspondence.
#' p-values are computed within groups of genes of similar degree.
#' @param log_scale TRUE if expression is log2-transformed.
#' @param ts_alpha FDR threshold used to report the number of selective and
#' excluded genes. It does not affect the ranking.
#' @param ts_relax Kept for backward compatibility; not used.
#' @param seed Integer seed for all permutation steps.
#' @param verbose TRUE/FALSE, show progress.
#' @param ... Ignored.
#' @return A list: \code{necorrHub_nodes} (ranked genes), \code{necorrReg}
#' (candidate regulators, see \code{\link{rank_regulators}}), \code{necorrEdges}
#' (top-ranked edges), \code{edge_ranking} (all edges), \code{netstat},
#' \code{coexpres} (edge correlations), \code{hub.m.param} (criteria matrix),
#' \code{criterion_weights}, \code{settings}, and \code{plots} if requested.
#' @examples
#' data(necorr_example)
#' res <- NECorr(networkFile = necorr_example$networkFile,
#'               expression = necorr_example$expression,
#'               description.file = necorr_example$description.file,
#'               condition = "flower", meta_data = necorr_example$meta_data,
#'               method = "PCC", permutation = 0,
#'               visualize = FALSE, save_results = FALSE, verbose = FALSE)
#' head(res$necorrHub_nodes)
#' res$criterion_weights
#' @export
NECorr <- function(networkFile, expression, description.file,
                   condition, meta_data, permutation = 1000,
                   method = "GCC", visualize = TRUE,
                   output_dir = NULL, interactive_net = FALSE,
                   top_n = 20, save_results = TRUE,
                   useBestGCC = FALSE, asymmetricGCC = FALSE,
                   weights = "equal", criteria = c("DE", "Degree"),
                   aggregation = c("product", "sum", "waspas"), epsilon = 0.01,
                   direction = c("up", "both"), scaling = c("rank", "softcap"),
                   degree_type = c("targets", "total"),
                   coexpr_stat = c("tail", "mean"), tail_q = 0.95,
                   directed = NULL, n_background = 200, full_netstat = FALSE,
                   edge_top = 0.05,
                   n_null = 0, log_scale = TRUE,
                   ts_alpha = 0.05, ts_relax = FALSE,
                   seed = 12345, verbose = TRUE, ...) {

  n_threads <- max(1, detectCores(logical = TRUE) - 2)
  steps <- c("Loading network", "Loading expression", "Loading description",
             "Loading metadata", "Preparing expression", "Co-expression",
             "Condition specificity", "Differential expression", "Network topology",
             "Scoring", "Null model", "Edges and regulators", "Visualization", "Saving results")
  update_progress <- function(step) {
    if (!verbose) return(invisible())
    bar <- paste0(strrep("#", round(30 * step / length(steps))),
                  strrep("-", 30 - round(30 * step / length(steps))))
    cat(sprintf("\r[%s] %-26s", bar, steps[step]))
    if (step == length(steps)) cat("\n")
    flush.console()
  }
  if (save_results && is.null(output_dir)) stop("Please provide an 'output_dir' to save results")
  aggregation <- match.arg(aggregation); coexpr_stat <- match.arg(coexpr_stat)
  direction <- match.arg(direction); scaling <- match.arg(scaling); degree_type <- match.arg(degree_type)
  all_criteria <- c("Coexpr", "CondCoexpr", "TS", "DE", "Bridging", "Degree")
  if (!all(criteria %in% all_criteria)) stop("criteria must be among: ", paste(all_criteria, collapse = ", "))
  cond_coexpr <- "CondCoexpr" %in% criteria
  combine <- function(P, w) {
    eps <- if (scaling == "rank") 0 else epsilon   # ranks are never zero
    wsum <- drop(P %*% w); wprod <- drop(exp(log(P + eps) %*% w))
    switch(aggregation, sum = wsum, product = wprod, waspas = 0.5 * wsum + 0.5 * wprod)
  }

  #### 1. Load network ####
  # --- Load network ---
  if (is.character(networkFile) && file.exists(networkFile)) {
    network.int <- fread(networkFile, header = FALSE,
                         encoding = "Latin-1",
                         colClasses = c("character", "character"),
                         nThread = n_threads, showProgress = FALSE)
    if (ncol(network.int) == 2) {
      setnames(network.int, c("source", "target"))
    } else if (ncol(network.int) >= 3) {
      # Take only columns 1 and 3
      network.int <- network.int[, .(source = V1, target = V3)]
    } else {
      stop("networkFile must have at least 2 columns.")
    }
  } else if (is.data.frame(networkFile) || is.data.table(networkFile)) {
    network.int <- safe_as_data_table(networkFile)
    # Ensure correct colnames
    if (!all(c("source", "target") %in% names(network.int))) {
      setnames(network.int, 1:2, c("source", "target"))
    }
  } else {
    stop("networkFile must be a file path or a data.frame/data.table")
  }
  network.int <- network.int[, .(source = as.character(source), target = as.character(target))]
  network.int <- network.int[!is.na(source) & !is.na(target) & source != "" & target != ""]
  update_progress(1)
  #### 2. Load expression ####
  if (is.character(expression) && file.exists(expression)) {
    eset <- fread(expression, header = TRUE,
                  nThread = n_threads, showProgress = FALSE)
    gene_ids <- eset[[1]]
    eset[, 1 := NULL]
    eset[, gene_id := gene_ids]
    setkey(eset, gene_id)
  } else {
    eset <- data.table::as.data.table(expression)  # copy
    if (!"gene_id" %in% names(eset)) {
      if (!is.character(eset[[1]]) && !is.factor(eset[[1]])) stop("expression: first column must hold gene IDs (or name it 'gene_id')")
      setnames(eset, 1, "gene_id")
    }
    eset[, gene_id := as.character(gene_id)]
    setkey(eset, gene_id)
  }
  update_progress(2)
  #### 3. Load description ####
  if (is.character(description.file) && file.exists(description.file)) {
    description_df <- fread(description.file, header = TRUE,
                            nThread = n_threads, showProgress = FALSE)
    colnames(description_df)[1] <- "gene_id"
  } else {
    description_df <- data.table::as.data.table(description.file)  # copy
    if (!"gene_id" %in% colnames(description_df)) setnames(description_df, 1, "gene_id")
  }
  description_df[, gene_id := as.character(gene_id)]
  setkey(description_df, gene_id)
  update_progress(3)
  #### 4. Load metadata ####
  # Accepted layouts (file, data.frame or data.table):
  #  (a) column 1 = sample IDs matching expression columns, column 2 = condition
  #  (b) sample IDs as row names, column 1 = condition
  expr_sample_cols <- setdiff(colnames(eset), "gene_id")
  if (is.character(meta_data) && length(meta_data) == 1 && file.exists(meta_data)) {
    md <- as.data.frame(fread(meta_data, header = TRUE, nThread = n_threads, showProgress = FALSE))
  } else if (is.data.frame(meta_data)) {
    rn <- rownames(meta_data)  # read before coercion: data.table row names are lost by as.data.frame()
    md <- as.data.frame(meta_data, stringsAsFactors = FALSE)  # copy; never modify the user's object
    rownames(md) <- rn
  } else {
    stop("meta_data must be a file path, data.frame or data.table")
  }
  if (ncol(md) >= 2 && all(as.character(md[[1]]) %in% expr_sample_cols)) {
    ids <- as.character(md[[1]]); cond <- as.character(md[[2]])
  } else if (!is.null(rownames(md)) && all(rownames(md) %in% expr_sample_cols)) {
    ids <- rownames(md); cond <- as.character(md[[1]])
  } else {
    stop("meta_data: sample IDs (column 1 or row names) must match the expression column names")
  }
  factortab <- data.table(sample_id = ids, condition = cond)
  if (!setequal(factortab$sample_id, expr_sample_cols)) {
    stop("Mismatch between metadata sample IDs and expression data column names.")
  }
  update_progress(4)

  #### 5. Prepare expression ####
  methods <- toupper(method)
  if (!all(methods %in% c("GCC", "PCC", "SCC", "KCC"))) stop("method must be among GCC, PCC, SCC, KCC")
  network.int <- unique(network.int[source != target, .(source, target)])
  Genelist <- unique(c(network.int$source, network.int$target))

  sample.names <- unique(factortab$condition)
  if (!condition %in% sample.names) stop("Condition '", condition, "' not found in meta_data")
  if (length(sample.names) < 2) stop("At least two conditions are needed")
  # columns grouped by condition, no duplication of samples
  cond_by_col <- factortab$condition[match(expr_sample_cols, factortab$sample_id)]
  idx_list <- unlist(lapply(sample.names, function(sn) which(cond_by_col == sn)))
  tissues <- cond_by_col[idx_list]
  treatment.f <- factor(tissues, levels = sample.names)
  if (any(table(tissues) < 2)) {
    warning("Some conditions have a single sample: differential expression and tissue specificity rely on pooled variance only.")
  }

  eset <- eset[gene_id %in% Genelist]
  if (anyDuplicated(eset$gene_id)) {
    warning("Duplicated gene IDs in expression data: the first occurrence is kept.")
    eset <- eset[!duplicated(gene_id)]
  }
  m.raw <- as.matrix(eset[, expr_sample_cols, with = FALSE])[, idx_list, drop = FALSE]
  storage.mode(m.raw) <- "double"
  rownames(m.raw) <- eset$gene_id
  clean <- sanitize_expression(m.raw)
  m.eset <- clean$mat
  if (sum(clean$report) > 0 && verbose) {
    message("\nExpression cleaning: ", paste(names(clean$report), clean$report, sep = "=", collapse = ", "))
  }
  if (nrow(m.eset) < 3) stop("Fewer than 3 network genes have usable expression data")
  if (ncol(m.eset) < 8) {
    warning("Only ", ncol(m.eset), " samples: edge correlations are imprecise below about 8-10 samples.")
  }
  target_cols <- tissues == condition
  update_progress(5)

  # local random stream: results are reproducible and the caller's RNG state is left untouched
  if (exists(".Random.seed", envir = globalenv())) {
    old_seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
  set.seed(seed + 104729L)

  #### 6. Co-expression on the network edges ####
  dir_info <- detect_direction(network.int)
  if (is.null(directed)) {
    directed <- dir_info$directed
    if (verbose) message("\nNetwork treated as ", if (directed) "directed" else "undirected",
                         " (", round(100 * dir_info$source_fraction), "% of genes are sources, ",
                         round(100 * dir_info$reciprocity), "% reciprocal edges). Set 'directed' to override.")
  }
  tc <- if (cond_coexpr) target_cols else NULL
  coex <- gene_coexpression(m.eset, network.int, methods, permutation, tc,
                            asymmetricGCC = asymmetricGCC, seed = seed, n_background = n_background,
                            stat = coexpr_stat, tail_q = tail_q, directed = directed)
  int.sig.full <- coex$edges
  update_progress(6)

  #### 7-8. Condition specificity and differential expression ####
  spec <- condition_specificity(m.eset, tissues, condition, alpha = ts_alpha, log_scale = log_scale)
  TS <- setNames(rep(0, length(Genelist)), Genelist); TS[spec$Gene] <- if (direction == "up") pmax(spec$t_selective, 0) else spec$TS
  update_progress(7)
  DE <- setNames(rep(0, length(Genelist)), Genelist); DE[spec$Gene] <- if (direction == "up") pmax(spec$DE_t, 0) else abs(spec$DE_t)
  update_progress(8)

  #### 9. Network topology ####
  g <- igraph::simplify(graph_from_data_frame(network.int, directed = FALSE))
  netstat <- data.table(Gene = V(g)$name, Degree = degree(g),
                        ClusCoef = transitivity(g, type = "local", isolates = "zero"),
                        PageRank = page_rank(g)$vector)
  # Betweenness is slow on large networks: computed only when it is used or asked for
  if ("Bridging" %in% criteria || full_netstat) {
    netstat[, Betweenness := betweenness(g)]
    # Bridging: betweenness beyond what degree predicts (residual on the log scale)
    lb <- log1p(netstat$Betweenness); ld <- log1p(netstat$Degree)
    netstat[, Bridging := if (stats::var(ld) > 0 && stats::var(lb) > 0) pmax(stats::residuals(stats::lm(lb ~ ld)), 0) else 0]
  } else netstat[, Bridging := 0]
  if (directed) {
    gd <- graph_from_data_frame(network.int, directed = TRUE)
    netstat[, OutDegree := degree(gd, mode = "out")[Gene]]
    netstat[, InDegree := degree(gd, mode = "in")[Gene]]
  }
  deg <- setNames(netstat$Degree, netstat$Gene)
  if (degree_type == "targets") {
    if (directed) {
      deg <- setNames(netstat$OutDegree, netstat$Gene)
    } else {
      degree_type <- "total"   # an undirected network has no targets: all partners are counted
      if (verbose) message("Undirected network: the Degree criterion counts all partners (degree_type = \"total\").")
    }
  }
  bridging <- setNames(netstat$Bridging, netstat$Gene)
  update_progress(9)

  #### 10. Criteria matrix and score ####
  sc <- if (scaling == "rank") rank01 else scale01
  build_criteria <- function(cx, TS, DE) {
    P <- cbind(Coexpr = sc(cx$Coexpr[Genelist]),
               CondCoexpr = if (!is.null(cx$CondCoexpr)) sc(cx$CondCoexpr[Genelist]),
               TS = sc(TS[Genelist]), DE = sc(DE[Genelist]),
               Bridging = sc(bridging[Genelist]),
               Degree = sc(log1p(deg[Genelist])))
    rownames(P) <- Genelist
    P[is.na(P)] <- 0
    P[, all_criteria[all_criteria %in% intersect(criteria, colnames(P))], drop = FALSE]
  }
  P <- build_criteria(coex$gene, TS, DE)
  constant <- apply(P, 2, function(x) max(x) == min(x))
  if (any(constant)) {
    if (verbose) message("\nCriteria without variation dropped: ", paste(colnames(P)[constant], collapse = ", "))
    P <- P[, !constant, drop = FALSE]
  }
  w <- resolve_weights(weights, P)
  score <- combine(P, w$weights)
  gene_rank_df <- data.frame(Gene = names(score), Score = unname(score), stringsAsFactors = FALSE)
  update_progress(10)

  #### 11. Empirical null: shuffle expression profiles among network genes ####
  # The network and the expression data are kept; their correspondence is broken.
  # Each gene is compared with null scores of genes of similar degree.
  if (n_null > 0) {
    rn <- rownames(m.eset)
    bins <- cut(rank(deg[Genelist], ties.method = "average"),
                breaks = min(10, length(unique(deg[Genelist]))), labels = FALSE)
    names(bins) <- Genelist
    null_scores <- matrix(NA_real_, length(Genelist), n_null, dimnames = list(Genelist, NULL))
    for (b in seq_len(n_null)) {
      new <- sample(rn)
      mp <- m.eset; rownames(mp) <- new
      bgp <- data.table::copy(coex$background); bgp[, Gene := new[match(Gene, rn)]]
      cxp <- gene_coexpression(mp, network.int, methods, 0, tc, asymmetricGCC = asymmetricGCC,
                               seed = seed + b, background = bgp, stat = coexpr_stat, tail_q = tail_q,
                               directed = directed)$gene
      TSp <- TS; TSp[new] <- TS[rn]
      DEp <- DE; DEp[new] <- DE[rn]
      Pp <- build_criteria(cxp, TSp, DEp)[, colnames(P), drop = FALSE]
      null_scores[, b] <- combine(Pp, w$weights)
    }
    pvals <- vapply(Genelist, function(gn) {
      nb <- null_scores[bins == bins[gn], , drop = FALSE]
      (sum(nb >= score[gn]) + 1) / (length(nb) + 1)
    }, numeric(1))
    gene_rank_df$p_empirical <- pvals[gene_rank_df$Gene]
    gene_rank_df$FDR <- p.adjust(gene_rank_df$p_empirical, method = "BH")
  }
  gene_rank_df <- gene_rank_df[order(-gene_rank_df$Score), , drop = FALSE]
  rownames(gene_rank_df) <- NULL
  update_progress(11)

  #### 12. Edge ranking and candidate regulators ####
  edge_ranking <- hub_edge_percentile(network.int, score)
  edge_ranking <- edge_ranking[order(-edge_ranking$edge_score), , drop = FALSE]
  top_edges <- edge_ranking[edge_ranking$percentile > 1 - edge_top, , drop = FALSE]
  actres <- rank_regulators(network.int, score, description_df, directed = directed)
  update_progress(12)

  #### 13. Results ####
  results <- list(
    necorrHub_nodes = gene_rank_df,
    necorrReg = actres,
    necorrEdges = top_edges,
    edge_ranking = edge_ranking,
    netstat = as.data.frame(netstat),
    coexpres = as.data.frame(int.sig.full),
    hub.m.param = as.data.frame(P),
    coexpression_stats = coex$stats,
    specificity = spec,
    criterion_weights = w$weights,
    settings = list(condition = condition, methods = methods, permutation = permutation,
                    weights = w$scheme, ahp_consistency_ratio = w$CR,
                    aggregation = aggregation, epsilon = epsilon, direction = direction, scaling = scaling, degree_type = degree_type,
                    coexpr_stat = coexpr_stat, tail_q = tail_q, directed = directed, edge_top = edge_top, n_null = n_null, seed = seed,
                    n_background = n_background, ts_alpha = ts_alpha,
                    n_selective = sum(spec$direction == "selective" & spec$FDR < ts_alpha),
                    n_excluded = sum(spec$direction == "excluded" & spec$FDR < ts_alpha),
                    expression_cleaning = clean$report,
                    n_genes = length(Genelist), n_edges = nrow(network.int), n_samples = ncol(m.eset))
  )
  if (visualize) {
    results$plots <- if (nrow(top_edges) == 0) NULL else
      visualize_necorr_results(results, top_n = top_n, interactive_net = interactive_net,
                               highlight_regulators = TRUE, output_dir = output_dir)
  }
  update_progress(13)

  if (save_results) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    fwrite(results$necorrHub_nodes, file.path(output_dir, "necorrHub_nodes.csv"))
    if (!is.null(results$necorrReg)) fwrite(results$necorrReg, file.path(output_dir, "necorrReg.csv"))
    fwrite(results$necorrEdges, file.path(output_dir, "necorrEdges.csv"))
    fwrite(results$coexpres, file.path(output_dir, "coexpres.csv"))
    fwrite(data.table(Gene = rownames(results$hub.m.param), results$hub.m.param), file.path(output_dir, "hub.m.param.csv"))
    fwrite(results$netstat, file.path(output_dir, "netstat.csv"))
    fwrite(data.table(criterion = names(w$weights), weight = w$weights), file.path(output_dir, "criterion_weights.csv"))
  }
  update_progress(14)
  results
}

# 0-1 scaling that keeps names. Values above the 99th percentile of the positive
# values are compressed into (q, 1.25 q] by a bounded increasing function
# (soft cap), so that one extreme gene cannot compress all others while the
# order of the top genes is preserved.
# A hard cap would tie the top 1% of genes and make the head of the ranking,
# the part that is followed up experimentally, depend on input order.
# Rank scaling: rank / n in (0, 1], ties averaged, missing values ranked lowest.
# With the product aggregation this gives a rank product of the criteria.
rank01 <- function(x) {
  x[!is.finite(x)] <- -Inf
  setNames(rank(x, ties.method = "average") / length(x), names(x))
}

scale01 <- function(x, cap = 0.99) {
  x[!is.finite(x)] <- NA
  pos <- x[!is.na(x) & x > 0]
  if (length(pos) >= 20 && cap < 1) {
    q <- stats::quantile(pos, cap, names = FALSE)
    over <- !is.na(x) & x > q
    x[over] <- q * (1 + 0.25 * (1 - q / x[over]))
  }
  r <- range(x, na.rm = TRUE)
  if (!all(is.finite(r)) || r[2] == r[1]) return(setNames(rep(0, length(x)), names(x)))
  (x - r[1]) / (r[2] - r[1])
}

# Lower confidence bound of a standardised excess: (mean of k partner values -
# background mean) / background sd, minus its one-sided 95% margin. Zero when
# the excess is not supported by the number of partners.
excess_lcb <- function(mean_partner, k, bg_mean, bg_sd, n_bg, conf = 0.95) {
  d <- (mean_partner - bg_mean) / bg_sd
  lcb <- d - stats::qnorm(conf) * sqrt(1 / k + 1 / n_bg)
  lcb[!is.finite(lcb)] <- 0
  pmax(lcb, 0)
}

# Contribution of the target samples to the Pearson correlation of each pair
pair_contribution <- function(Z, i, j, target_cols, chunk = 200000L) {
  out <- numeric(length(i)); n <- ncol(Z)
  for (s in seq(1L, length(i), by = chunk)) {
    e <- s:min(length(i), s + chunk - 1L)
    out[e] <- rowSums(Z[i[e], target_cols, drop = FALSE] * Z[j[e], target_cols, drop = FALSE]) / (n - 1)
  }
  out
}

#' Per-gene background of co-expression with non-partner genes
#'
#' @description For every gene, the correlation with \code{n_background}
#' random genes that are not its network partners. Random partners are drawn
#' as in a degree-preserving randomisation of the network (configuration
#' model): with probability proportional to their degree, and in a directed
#' network according to the role of the gene (random targets for a regulator,
#' random regulators for a target). The background is therefore what the gene
#' would show in a network with the same degrees and no biological wiring. It
#' removes two effects from the network co-expression: genes whose profile
#' follows the main tissue differences correlate with many genes, and
#' high-degree genes are over-represented among anyone's partners.
#' @param m Expression matrix (genes x samples).
#' @param edges data.table with columns source and target.
#' @param methods Correlation methods.
#' @param target_cols Logical vector marking the target-condition samples, or NULL.
#' @param n_background Number of random non-partner genes per gene.
#' @param asymmetricGCC Logical.
#' @param tail_q Quantile of the background that a partner must exceed to count
#' as strongly co-expressed.
#' @param directed TRUE if edges go from regulator (source) to target.
#' @return data.table with one row per gene: bg_mean, bg_sd, bg_q (absolute
#' correlation) and bgc_mean, bgc_sd, bgc_q (target-condition contribution).
#' @export
coexpression_background <- function(m, edges, methods = "PCC", target_cols = NULL,
                                    n_background = 200, asymmetricGCC = FALSE, tail_q = 0.95,
                                    directed = FALSE) {
  genes <- rownames(m); G <- length(genes)
  nb <- min(n_background, G - 1)
  ei <- match(edges$source, genes); ej <- match(edges$target, genes)
  ok <- !is.na(ei) & !is.na(ej); ei <- ei[ok]; ej <- ej[ok]
  a <- rep(seq_len(G), each = nb)
  # Degree-preserving null (configuration model): a random partner is drawn with
  # probability proportional to its degree, so that well-connected genes are as
  # frequent among the random partners as among the real ones.
  if (directed) {
    kout <- tabulate(ei, G); kin <- tabulate(ej, G)
    as_source <- stats::runif(length(a)) < (kout / pmax(kout + kin, 1))[a]   # role of the gene on this draw
    b <- integer(length(a))
    if (any(as_source)) b[as_source] <- sample.int(G, sum(as_source), replace = TRUE, prob = kin + 1e-9)
    if (any(!as_source)) b[!as_source] <- sample.int(G, sum(!as_source), replace = TRUE, prob = kout + 1e-9)
  } else {
    b <- sample.int(G, length(a), replace = TRUE, prob = tabulate(c(ei, ej), G) + 1e-9)
  }
  key <- function(x, y) (pmin(x, y) - 1) * G + pmax(x, y)
  ok <- a != b & !(key(a, b) %in% key(ei, ej))
  a <- a[ok]; b <- b[ok]
  bg <- data.table(source = genes[a], target = genes[b])
  r <- multiCorr(m, bg, methods = methods, pernum = 0, asymmetricGCC = asymmetricGCC)
  bg[, absr := abs(rowMeans(as.matrix(r[, methods, with = FALSE]), na.rm = TRUE))]
  if (!is.null(target_cols)) {
    Z <- t(scale(t(m)))
    bg[, absc := abs(pair_contribution(Z, a, b, target_cols))]
  } else bg[, absc := NA_real_]
  bg[, .(bg_mean = mean(absr), bg_sd = stats::sd(absr), bg_q = stats::quantile(absr, tail_q, names = FALSE),
         bgc_mean = mean(absc), bgc_sd = stats::sd(absc),
         bgc_q = if (all(is.na(absc))) NA_real_ else stats::quantile(absc, tail_q, names = FALSE, na.rm = TRUE),
         n_bg = .N), by = .(Gene = source)]
}

# z-score of the number of partners above the gene's own background quantile
tail_excess <- function(count, k, tail_q) {
  e <- 1 - tail_q
  z <- (count - k * e) / sqrt(k * e * (1 - e))
  z[!is.finite(z)] <- 0
  pmax(z, 0)
}

#' Gene-level co-expression along network edges
#'
#' @description Computes the edge correlations, then asks, for each gene,
#' whether it is more co-expressed with its network partners than with random
#' non-partner genes.
#' \describe{
#'   \item{Coexpr}{with \code{stat = "tail"} (default): z-score of the number
#'     of partners whose absolute correlation with the gene exceeds the
#'     \code{tail_q} quantile of the gene's own background, against the number
#'     expected by chance. Large networks from binding assays contain many
#'     edges without regulatory effect; counting the strongly co-expressed
#'     partners is not diluted by them. With \code{stat = "mean"}: excess of the
#'     mean absolute correlation over the background, as the one-sided 95\%
#'     lower confidence bound of the standardised difference. Neither grows
#'     with degree when there is no signal.}
#'   \item{CondCoexpr}{the same for the contribution of the target-condition
#'     samples to the correlation: high when the gene and its partners depart
#'     from their means together in the target condition.}
#' }
#' @param m Expression matrix (genes x samples).
#' @param edges data.table with columns source and target.
#' @param methods Correlation methods.
#' @param pernum Number of permutations for the GCC p-value.
#' @param target_cols Logical vector marking the samples of the target condition.
#' @param asymmetricGCC Logical.
#' @param seed Integer seed.
#' @param background Output of \code{\link{coexpression_background}}; computed if NULL.
#' @param n_background Number of random non-partner genes per gene.
#' @param stat "tail" or "mean", see above.
#' @param tail_q Background quantile used by \code{stat = "tail"}.
#' @param directed TRUE if edges go from regulator (source) to target; used for the background.
#' @return List with \code{edges} (edge table with correlations, combined
#' \code{Correlation}, \code{pvalue}, \code{cond_contribution}), \code{gene}
#' (named vectors Coexpr, CondCoexpr), \code{stats} (per-gene details) and
#' \code{background}.
#' @export
gene_coexpression <- function(m, edges, methods = "PCC", pernum = 0, target_cols = NULL,
                              asymmetricGCC = FALSE, seed = 12345, background = NULL,
                              n_background = 200, stat = c("tail", "mean"), tail_q = 0.95,
                              directed = FALSE) {
  stat <- match.arg(stat)
  edges <- data.table::as.data.table(edges)
  full <- multiCorr(m, edges, methods = methods, pernum = pernum, asymmetricGCC = asymmetricGCC, seed = seed)
  full[, Correlation := rowMeans(as.matrix(.SD), na.rm = TRUE), .SDcols = methods]
  pm <- as.matrix(full[, paste0(methods, "_pvalue"), with = FALSE])
  # Cauchy combination: valid for the dependent p-values of the different methods
  full[, pvalue := if (all(is.na(pm))) NA_real_ else apply(pm, 1, cauchyCombine)]
  genes <- rownames(m)
  i <- match(full$source, genes); j <- match(full$target, genes)
  if (!is.null(target_cols)) {
    Z <- t(scale(t(m)))
    cc <- rep(NA_real_, nrow(full)); ok <- !is.na(i) & !is.na(j)
    cc[ok] <- pair_contribution(Z, i[ok], j[ok], target_cols)
    full[, cond_contribution := cc]
  } else full[, cond_contribution := NA_real_]
  if (is.null(background)) {
    background <- coexpression_background(m, edges, methods, target_cols, n_background, asymmetricGCC, tail_q, directed)
  }
  long <- rbind(full[, .(Gene = source, r = abs(Correlation), cc = abs(cond_contribution))],
                full[, .(Gene = target, r = abs(Correlation), cc = abs(cond_contribution))])
  long <- long[is.finite(r)]
  long <- merge(long, background[, .(Gene, bg_q, bgc_q)], by = "Gene", all.x = TRUE)
  agg <- long[, .(k = .N, mean_abs_r = mean(r), mean_abs_c = mean(cc),
                  n_strong = sum(r > bg_q, na.rm = TRUE), n_strong_c = sum(cc > bgc_q, na.rm = TRUE)), by = Gene]
  agg <- merge(agg, background, by = "Gene", all.x = TRUE)
  if (stat == "tail") {
    agg[, Coexpr := tail_excess(n_strong, k, tail_q)]
    agg[, CondCoexpr := if (is.null(target_cols)) NA_real_ else tail_excess(n_strong_c, k, tail_q)]
  } else {
    agg[, Coexpr := excess_lcb(mean_abs_r, k, bg_mean, bg_sd, n_bg)]
    agg[, CondCoexpr := if (is.null(target_cols)) NA_real_ else excess_lcb(mean_abs_c, k, bgc_mean, bgc_sd, n_bg)]
  }
  list(edges = full,
       gene = list(Coexpr = setNames(agg$Coexpr, agg$Gene),
                   CondCoexpr = if (is.null(target_cols)) NULL else setNames(agg$CondCoexpr, agg$Gene)),
       stats = as.data.frame(agg), background = background)
}

#' Condition specificity of each gene
#'
#' @description A continuous version of the intersection-union test. The target
#' condition is compared with every other condition (limma moderated t). A gene
#' selectively expressed in the target must be higher than the other conditions;
#' a gene excluded from it must be lower. The score is the t-statistic that
#' \code{prop_required} of the comparisons reach (with the default 0.8, one
#' similar condition in five is tolerated), so every gene receives a value and
#' no significance threshold decides which genes enter the ranking.
#' @param m Expression matrix (genes x samples), log scale.
#' @param tissues Condition label of each column.
#' @param target Target condition.
#' @param prop_required Fraction of the pairwise comparisons that must agree.
#' @param alpha FDR threshold used only to report how many genes are selective or excluded.
#' @param log_scale TRUE if \code{m} is log2-transformed (used for the tau index).
#' @return Data frame: Gene, TS (score, zero when neither selective nor
#' excluded), direction, t_selective, t_excluded, FDR, tau, DE_t (target
#' against the mean of the other conditions).
#' @export
condition_specificity <- function(m, tissues, target, prop_required = 0.8, alpha = 0.05, log_scale = TRUE) {
  conds <- unique(tissues)
  if (!target %in% conds) stop("Target condition '", target, "' not found")
  others <- setdiff(conds, target)
  f <- factor(tissues, levels = conds)
  design <- stats::model.matrix(~ 0 + f); colnames(design) <- paste0("c", seq_along(conds))
  tn <- paste0("c", match(target, conds)); on <- paste0("c", match(others, conds))
  fit <- limma::lmFit(m, design)
  cm <- limma::makeContrasts(contrasts = c(paste(tn, on, sep = " - "),
                                           paste0(tn, " - (", paste(on, collapse = " + "), ")/", length(on))), levels = design)
  fit2 <- suppressWarnings(limma::eBayes(limma::contrasts.fit(fit, cm)))
  tt <- fit2$t[, seq_along(others), drop = FALSE]
  # value reached by prop_required of the comparisons: order statistic from the bottom
  kth <- max(1, length(others) - ceiling(prop_required * length(others)) + 1)
  pick <- function(x) apply(x, 1, function(v) sort(v)[kth])
  up <- pick(tt); down <- pick(-tt)
  df <- fit2$df.total[1]
  stat <- pmax(up, down)
  p <- stats::pt(stat, df = df, lower.tail = FALSE)
  lin <- if (log_scale) 2^m else m
  means <- vapply(conds, function(cn) rowMeans(lin[, tissues == cn, drop = FALSE]), numeric(nrow(m)))
  tau <- rowSums(1 - means / pmax(apply(means, 1, max), 1e-12)) / (length(conds) - 1)
  data.frame(Gene = rownames(m), TS = pmax(stat, 0),
             direction = ifelse(stat <= 0, "none", ifelse(up >= down, "selective", "excluded")),
             t_selective = up, t_excluded = down, FDR = stats::p.adjust(p, "BH"), tau = tau,
             DE_t = fit2$t[, length(others) + 1], stringsAsFactors = FALSE, row.names = NULL)
}
