test_that("tissue indices separate specificity from exclusion", {
  x <- rbind(excluded_in_one = c(8, 8, 8, 8, 1), specific_to_one = c(1, 8, 1, 1, 1))
  r <- compute_TSI_TSE_parallel(x)
  expect_equal(r$TSE[1], 1)
  expect_lt(r$TSE[2], 0.5)
  expect_gt(r$TSI[2], r$TSI[1])
})

test_that("tissue indices are computed on the linear scale and stay in range", {
  m <- matrix(c(-2, -1, 0.5, 1, -3,  2, 2, 2, 9, 2), 2, byrow = TRUE,
              dimnames = list(c("a", "b"), paste0("s", 1:5)))
  r <- fast_tsi_tse(m, tissues = paste0("t", 1:5), log_scale = TRUE)
  expect_true(all(is.finite(r$TSI) & r$TSI >= 0 & r$TSI <= 1))
})

test_that("ts.IUT finds planted selective and excluded genes without relaxing", {
  set.seed(1)
  tissues <- rep(c("root", "leaf", "stem", "seed"), each = 4)
  m <- matrix(rnorm(300 * 16, 8, 0.3), 300, 16, dimnames = list(paste0("g", 1:300), NULL))
  m[1:10, tissues == "leaf"] <- m[1:10, tissues == "leaf"] + 4
  m[11:20, tissues == "leaf"] <- m[11:20, tissues == "leaf"] - 4
  r <- ts.IUT(eset = m, tissues = tissues, target = "leaf")
  expect_false(r$relaxed)
  expect_equal(r$finalThreshold, 0.05)
  expect_gte(sum(paste0("g", 1:10) %in% r$filtTSI$gene), 9)
  expect_gte(sum(paste0("g", 11:20) %in% r$filtTSE$gene), 9)
  expect_lte(nrow(r$filtTSI) + nrow(r$filtTSE), 25)
})

test_that("sanitize_expression caps infinities at the observed range", {
  m <- matrix(c(1, 2, -Inf, 4, 5, 6,   3, Inf, 5, 6, 7, 8,   2, 2, 2, 2, 2, 2,   NA, NA, NA, 1, 2, 3,   1, NA, 3, 4, 5, 6),
              5, byrow = TRUE, dimnames = list(paste0("g", 1:5), NULL))
  s <- sanitize_expression(m, max_missing = 0.2)
  expect_true(all(is.finite(s$mat)))
  expect_equal(max(s$mat), 8)   # +Inf -> finite maximum
  expect_equal(min(s$mat), 1)   # -Inf -> finite minimum
  expect_setequal(rownames(s$mat), c("g1", "g2", "g5"))
  expect_equal(unname(s$mat["g5", 2]), 4)  # NA -> gene median
  expect_equal(unname(s$report["genes_removed_constant"]), 1)
})

test_that("AHP weights and consistency ratio are correct", {
  a <- ahp_weights(necorr_ahp_matrix())
  expect_equal(sum(a$weights), 1)
  expect_lt(a$CR, 0.1)
  # validated default: degree first, specificity and differential expression next, the rest equal
  expect_equal(unname(round(a$weights[c("Degree", "TS", "DE", "Coexpr", "CondCoexpr", "Bridging")], 2)), c(0.35, 0.21, 0.21, 0.08, 0.08, 0.08))
  # NECorr 3.2 judgements remain available
  b <- ahp_weights(necorr_ahp_matrix(preset = "wiring"))
  expect_equal(unname(b$weights["Coexpr"]), unname(b$weights["TS"]))
  expect_gt(b$weights["CondCoexpr"], b$weights["Bridging"])
  d <- ahp_weights(necorr_ahp_matrix(c("Coexpr", "TS", "CondCoexpr"), preset = "wiring"))$weights
  expect_equal(unname(round(d, 2)), c(0.4, 0.4, 0.2))
  # perfectly consistent matrix built from known weights
  w <- c(a = 0.6, b = 0.3, c = 0.1); A <- outer(w, w, "/")
  b <- ahp_weights(A)
  expect_equal(b$weights, w, tolerance = 1e-8)
  expect_equal(b$CR, 0, tolerance = 1e-8)
  bad <- A; bad[1, 2] <- 9
  expect_error(ahp_weights(bad), "reciprocal")
})

test_that("CRITIC shares weight between redundant criteria", {
  set.seed(1)
  x <- runif(200); y <- runif(200)
  P <- cbind(a = x, a_copy = x + rnorm(200, 0, 0.01), b = y)
  w <- critic_weight(P)
  expect_equal(sum(w), 1)
  expect_gt(w["b"], w["a"])
})

test_that("edge ranking reports a percentile, not a p-value", {
  net <- data.frame(source = c("a", "a", "b"), target = c("b", "c", "c"))
  e <- hub_edge_percentile(net, c(a = 1, b = 0.5, c = 0))
  expect_named(e, c("sourceIDs", "targetIDs", "edge_score", "percentile"))
  expect_equal(e$edge_score, c(0.75, 0.5, 0.25))
  expect_equal(e$percentile, c(1, 2 / 3, 1 / 3))
})

test_that("direction is detected for regulatory and interaction networks", {
  expect_true(detect_direction(sim_data(1, directed = TRUE)$net)$directed)
  expect_false(detect_direction(sim_data(1, directed = FALSE)$net)$directed)
})

test_that("condition specificity is continuous and finds selective and excluded genes", {
  set.seed(1)
  tissues <- rep(c("root", "leaf", "stem", "seed", "flower"), each = 3)
  m <- matrix(rnorm(300 * 15, 8, 0.5), 300, 15, dimnames = list(paste0("g", 1:300), NULL))
  m[1:10, tissues == "leaf"] <- m[1:10, tissues == "leaf"] + 2
  m[11:20, tissues == "leaf"] <- m[11:20, tissues == "leaf"] - 2
  m[21:30, tissues %in% c("leaf", "stem", "seed")] <- m[21:30, tissues %in% c("leaf", "stem", "seed")] + 2   # high in three: not specific
  s <- condition_specificity(m, tissues, "leaf")
  expect_true(all(s$direction[1:10] == "selective"))
  expect_true(all(s$direction[11:20] == "excluded"))
  expect_gt(min(s$TS[1:20]), stats::quantile(s$TS[31:300], 0.95))
  expect_lt(mean(s$TS[21:30]), mean(s$TS[1:10]) / 2)
  expect_gt(mean(s$TS[31:300] > 0), 0.3)     # every gene gets a value: no gate
  expect_true(all(s$tau >= 0 & s$tau <= 1))
})

test_that("co-expression with partners is separated from degree and from global correlation", {
  set.seed(2); G <- 400; n <- 20
  genes <- paste0("g", 1:G)
  E <- matrix(rnorm(G * n), G, n, dimnames = list(genes, NULL))
  f <- rnorm(n); E[1:11, ] <- outer(rep(1, 11), f) + matrix(rnorm(11 * n, 0, 0.4), 11)       # module: g1 with partners g2..g11
  h <- rnorm(n); E[201:300, ] <- outer(rep(1, 100), h) + matrix(rnorm(100 * n, 0, 0.6), 100)  # broad cluster unrelated to wiring
  net <- data.frame(source = c(rep("g1", 10), rep("g201", 10), rep("g50", 60), sample(genes, 600, TRUE)),
                    target = c(genes[2:11], sample(genes[12:200], 10), sample(genes[12:200], 60), sample(genes, 600, TRUE)))
  net <- unique(net[net$source != net$target, ])
  cx <- gene_coexpression(E, net, "PCC", directed = TRUE)
  co <- cx$gene$Coexpr
  expect_gt(co["g1"], stats::quantile(co, 0.99))
  expect_lt(co["g201"], co["g1"] / 3)   # correlated with a broad cluster, but not specifically with its partners
  expect_lt(co["g50"], co["g1"] / 3)    # many partners, no co-expression
})

test_that("criterion scaling keeps the order of the top genes (no ties above the 99th percentile)", {
  set.seed(1); x <- setNames(c(rexp(500), 0, 0), paste0("g", 1:502))
  s <- NECorr:::scale01(x)
  expect_equal(range(s), c(0, 1))
  expect_equal(order(-s), order(-x))
  expect_equal(sum(s == 1), 1L)
  # an extreme value does not compress the rest
  y <- x; y[1] <- 1e6
  expect_gt(stats::median(NECorr:::scale01(y)[-1]), 0.05)
})
