test_that("the bundled example runs with every weighting scheme", {
  data(necorr_example)
  for (w in c("ahp", "critic", "equal", "entropy")) {
    r <- NECorr(necorr_example$networkFile, necorr_example$expression, necorr_example$description.file,
                "flower", necorr_example$meta_data, permutation = 0, method = "PCC", weights = w,
                visualize = FALSE, save_results = FALSE, verbose = FALSE)
    expect_equal(sum(r$criterion_weights), 1)
    expect_false(anyNA(r$necorrHub_nodes$Score))
  }
})

test_that("the default is the rank product of up-regulation and the number of targets", {
  d <- sim_data(1)
  r <- NECorr(d$net, d$expr, d$desc, d$target, d$meta, permutation = 0, method = "PCC",
              visualize = FALSE, save_results = FALSE, verbose = FALSE)
  expect_setequal(colnames(r$hub.m.param), c("DE", "Degree"))
  expect_equal(r$settings$aggregation, "product"); expect_equal(r$settings$direction, "up")
  expect_equal(r$settings$scaling, "rank"); expect_equal(r$settings$weights, "equal")
  expect_equal(r$settings$degree_type, if (isTRUE(r$settings$directed)) "targets" else "total")
  P <- as.matrix(r$hub.m.param); n <- nrow(P)
  # criteria are ranks / n and the score is their geometric mean
  sp <- r$specificity; up <- setNames(rep(0, n), rownames(P)); up[sp$Gene] <- pmax(sp$DE_t, 0)
  expect_equal(unname(P[, "DE"]), unname(rank(up[rownames(P)]) / n))
  sc <- setNames(r$necorrHub_nodes$Score, r$necorrHub_nodes$Gene)
  expect_equal(unname(sc[rownames(P)]), unname(sqrt(P[, "DE"] * P[, "Degree"])), tolerance = 1e-10)
  # genes expressed less in the target condition get the lowest DE rank
  down <- sp$Gene[sp$DE_t < 0]; expect_true(all(P[down, "DE"] <= min(P[setdiff(rownames(P), down), "DE"])))
  expect_lte(match(d$hub, r$necorrHub_nodes$Gene), 10)
  rb <- NECorr(d$net, d$expr, d$desc, d$target, d$meta, permutation = 0, method = "PCC", direction = "both",
               visualize = FALSE, save_results = FALSE, verbose = FALSE)
  expect_false(isTRUE(all.equal(rb$hub.m.param$DE, r$hub.m.param$DE)))
})

test_that("a planted condition-specific hub is ranked at the top", {
  d <- sim_data(1)
  r <- run_necorr(d)
  expect_lte(match(d$hub, r$necorrHub_nodes$Gene), 3)
  expect_equal(r$settings$weights, "ahp")
  expect_lt(r$settings$ahp_consistency_ratio, 0.1)
  expect_equal(r$settings$aggregation, "sum")
  expect_setequal(colnames(r$hub.m.param), c("Coexpr", "CondCoexpr", "TS", "DE", "Bridging", "Degree"))
  # the wiring-only ranking of NECorr 3.2 is still available
  c3 <- c("Coexpr", "TS", "CondCoexpr")
  r3 <- run_necorr(d, criteria = c3, aggregation = "product", weights = necorr_ahp_matrix(c3, preset = "wiring"))
  expect_lte(match(d$hub, r3$necorrHub_nodes$Gene), 3)
  expect_equal(unname(round(r3$criterion_weights[c3], 2)), c(0.4, 0.4, 0.2))
})

test_that("condition-specific genes receive a specificity score without any gate", {
  d <- sim_data(1, effect = 1.5)                 # weak contrast: the old gated criterion was empty here
  r <- run_necorr(d)
  expect_true("TS" %in% colnames(r$hub.m.param))
  expect_gt(r$hub.m.param[d$hub, "TS"], stats::median(r$hub.m.param$TS))
  expect_gt(mean(r$hub.m.param$TS > 0), 0.3)
  expect_equal(nrow(r$specificity), nrow(r$hub.m.param))
})

test_that("NECorr leaves the caller's random number stream untouched", {
  d <- sim_data(3, G = 60)
  set.seed(42); a <- runif(1); set.seed(42); invisible(run_necorr(d)); b <- runif(1)
  expect_equal(a, b)
})

test_that("results do not depend on condition names or their order", {
  d1 <- sim_data(2, labels = c("zeta", "alpha", "mid", "beta"))
  d2 <- sim_data(2, labels = c("a", "b", "c", "d"))
  r1 <- run_necorr(d1); r2 <- run_necorr(d2)
  expect_equal(r1$necorrHub_nodes$Score, r2$necorrHub_nodes$Score)
  expect_equal(r1$necorrHub_nodes$Gene, r2$necorrHub_nodes$Gene)
})

test_that("metadata and expression are accepted in every layout", {
  d <- sim_data(3, G = 80)
  ref <- run_necorr(d)$necorrHub_nodes
  meta_rn <- data.frame(Condition = d$meta$Condition, row.names = d$meta$SampleID)
  meta_before <- d$meta
  variants <- list(
    run_necorr(d, meta_data = meta_rn),
    run_necorr(d, meta_data = data.table::as.data.table(d$meta)),
    run_necorr(d, expression = data.table::as.data.table(d$expr), networkFile = data.table::as.data.table(d$net)),
    run_necorr(d, meta_data = d$meta[sample(nrow(d$meta)), ]))   # row order must not matter
  for (v in variants) expect_equal(v$necorrHub_nodes, ref)
  expect_identical(d$meta, meta_before)   # user objects are never modified
  tmp <- tempfile(); dir.create(tmp)
  f <- function(x, n, ...) { p <- file.path(tmp, n); data.table::fwrite(x, p, sep = "\t", ...); p }
  rf <- run_necorr(d, networkFile = f(d$net, "net.txt", col.names = FALSE), expression = f(d$expr, "expr.txt"),
                   description.file = f(d$desc, "desc.txt"), meta_data = f(d$meta, "meta.txt"))
  expect_equal(rf$necorrHub_nodes$Gene, ref$Gene)
})

test_that("infinite and missing expression values do not break the ranking", {
  d <- sim_data(4, G = 80)
  d$expr[5, 3] <- -Inf; d$expr[6, 4] <- NA
  r <- run_necorr(d)
  expect_false(anyNA(r$necorrHub_nodes$Score))
  expect_equal(unname(r$settings$expression_cleaning["neg_inf_capped"]), 1)
})

test_that("regulators are ranked using direction when the network has one", {
  d <- sim_data(1, directed = TRUE)
  r <- run_necorr(d, directed = TRUE)
  expect_equal(unique(r$necorrReg$role), "regulator")
  expect_lte(match(d$hub, r$necorrReg$Gene), 3)
  u <- run_necorr(sim_data(1, directed = FALSE), directed = FALSE)
  expect_equal(unique(u$necorrReg$role), "associated")
})

test_that("top edges are a percentile selection", {
  r <- run_necorr(sim_data(5), edge_top = 0.1)
  expect_true(all(r$necorrEdges$percentile > 0.9))
  expect_false("p2" %in% names(r$necorrEdges))
  expect_equal(nrow(r$edge_ranking), r$settings$n_edges)
})

test_that("empirical p-values are returned and favour the planted hub", {
  d <- sim_data(6, G = 150)
  r <- run_necorr(d, n_null = 30)
  h <- r$necorrHub_nodes
  expect_true(all(c("p_empirical", "FDR") %in% names(h)))
  expect_lt(h$p_empirical[h$Gene == d$hub], 0.05)
  r2 <- run_necorr(d, n_null = 30)
  expect_identical(r$necorrHub_nodes, r2$necorrHub_nodes)   # seeded
})

test_that("criteria, aggregation, custom weights and pairwise matrices are honoured", {
  d <- sim_data(7, G = 80)
  r0 <- run_necorr(d, criteria = c("Coexpr", "TS", "CondCoexpr", "DE", "Bridging", "Degree"), aggregation = "sum")
  crit <- colnames(r0$hub.m.param)
  expect_setequal(crit, c("Coexpr", "CondCoexpr", "TS", "DE", "Bridging", "Degree"))
  expect_equal(r0$necorrHub_nodes$Score[1], max(drop(as.matrix(r0$hub.m.param) %*% r0$criterion_weights)))
  rp <- run_necorr(d, aggregation = "product"); rw <- run_necorr(d, aggregation = "waspas")
  expect_false(isTRUE(all.equal(rp$necorrHub_nodes$Score, rw$necorrHub_nodes$Score)))
  w <- setNames(rep(1, length(crit)), crit); w["Degree"] <- 0
  r <- run_necorr(d, criteria = crit, weights = w)
  expect_equal(unname(r$criterion_weights["Degree"]), 0)
  c3 <- c("Coexpr", "TS", "CondCoexpr")
  A <- necorr_ahp_matrix(c3); A[1, 2] <- 9; A[2, 1] <- 1 / 9; A[2, 3] <- 9; A[3, 2] <- 1 / 9; A[1, 3] <- 1 / 9; A[3, 1] <- 9
  expect_warning(run_necorr(d, criteria = c3, weights = A), "consistency ratio")
  expect_error(run_necorr(d, weights = c(Coexpr = 1)), "named and cover")
  expect_error(run_necorr(d, criteria = c("Coexpr", "Nope")), "criteria must be")
})

test_that("degree: the wiring-only score ignores it, the default score includes it, p-values are calibrated", {
  set.seed(11); G <- 600; genes <- paste0("G", seq_len(G))
  g <- igraph::sample_pa(G, m = 2, directed = TRUE); el <- igraph::as_edgelist(g)
  d <- list(net = data.frame(source = genes[el[, 2]], target = genes[el[, 1]]),
            expr = data.frame(GeneID = genes, matrix(rnorm(G * 16, 7), G, 16, dimnames = list(NULL, paste0("S", 1:16)))),
            desc = data.frame(gene_id = genes), meta = data.frame(SampleID = paste0("S", 1:16), Condition = rep(paste0("c", 1:4), each = 4)), target = "c1")
  # wiring-only ranking (NECorr 3.2 setting): the score is independent of degree
  c3 <- c("Coexpr", "TS", "CondCoexpr")
  r3 <- run_necorr(d, directed = TRUE, criteria = c3, aggregation = "product", weights = necorr_ahp_matrix(c3, preset = "wiring"))
  h3 <- r3$necorrHub_nodes
  expect_lt(abs(cor(h3$Score, igraph::degree(g)[match(h3$Gene, genes)], method = "spearman")), 0.25)
  # validated default: the score includes degree, the empirical p-values stay calibrated and independent of degree
  r <- run_necorr(d, directed = TRUE, n_null = 40)
  h <- r$necorrHub_nodes; deg <- igraph::degree(g)[match(h$Gene, genes)]
  expect_gt(cor(h$Score, deg, method = "spearman"), 0.3)
  expect_lt(abs(cor(-log10(h$p_empirical), deg, method = "spearman")), 0.15)
  expect_lt(mean(h$p_empirical < 0.05), 0.09)
  expect_gt(mean(h$p_empirical < 0.05), 0.02)
})

test_that("invalid input is reported clearly", {
  d <- sim_data(8, G = 60)
  expect_error(run_necorr(d, condition = "nope"), "not found")
  bad <- d$meta; bad$SampleID <- paste0("X", bad$SampleID)
  expect_error(run_necorr(d, meta_data = bad), "sample IDs")
  expect_error(run_necorr(d, method = "XYZ"), "method must be")
})

test_that("degree counts targets in a directed network and all partners otherwise", {
  d <- sim_data(1)
  r <- NECorr(d$net, d$expr, d$desc, d$target, d$meta, permutation = 0, method = "PCC", directed = TRUE,
              visualize = FALSE, save_results = FALSE, verbose = FALSE)
  expect_equal(r$settings$degree_type, "targets")
  P <- as.matrix(r$hub.m.param); n <- nrow(P); net <- unique(as.data.frame(d$net)[, 1:2]); net <- net[net[[1]] != net[[2]], ]
  out <- setNames(rep(0, n), rownames(P)); tb <- table(net[[1]]); out[names(tb)] <- as.numeric(tb)
  expect_equal(unname(P[, "Degree"]), unname(rank(out) / n))
  # genes without targets share the lowest degree rank
  expect_true(all(P[out == 0, "Degree"] < min(P[out > 0, "Degree"])))
  # total degree on request, and always for an undirected network
  rt <- NECorr(d$net, d$expr, d$desc, d$target, d$meta, permutation = 0, method = "PCC", directed = TRUE, degree_type = "total",
               visualize = FALSE, save_results = FALSE, verbose = FALSE)
  g <- igraph::simplify(igraph::graph_from_data_frame(net, directed = FALSE)); dg <- igraph::degree(g)[rownames(rt$hub.m.param)]
  expect_equal(unname(rt$hub.m.param[, "Degree"]), unname(rank(dg) / length(dg)))
  ru <- NECorr(d$net, d$expr, d$desc, d$target, d$meta, permutation = 0, method = "PCC", directed = FALSE,
               visualize = FALSE, save_results = FALSE, verbose = FALSE)
  expect_equal(ru$settings$degree_type, "total")
  expect_equal(unname(ru$hub.m.param[rownames(rt$hub.m.param), "Degree"]), unname(rt$hub.m.param[, "Degree"]))
})
