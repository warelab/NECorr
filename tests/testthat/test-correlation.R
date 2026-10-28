null_edges <- function(n, K, seed) {
  set.seed(seed)
  X <- matrix(rnorm(2 * K * n), 2 * K, n, dimnames = list(paste0("g", seq_len(2 * K)), NULL))
  list(X = X, ed = data.frame(source = paste0("g", seq(1, 2 * K, 2)), target = paste0("g", seq(2, 2 * K, 2))))
}

test_that("GCC permutation p-values are calibrated under the null", {
  d <- null_edges(10, 1500, 1)
  r <- multiCorr(d$X, d$ed, methods = c("GCC", "PCC"), pernum = 500)
  expect_gt(mean(r$GCC_pvalue < 0.05), 0.03)
  expect_lt(mean(r$GCC_pvalue < 0.05), 0.08)
  expect_gt(mean(r$PCC_pvalue < 0.05), 0.03)
  expect_lt(mean(r$PCC_pvalue < 0.05), 0.08)
})

test_that("GCC detects strong correlation and is reproducible", {
  set.seed(2); n <- 20; K <- 200
  a <- matrix(rnorm(K * n), K, n); b <- 0.9 * a + sqrt(0.19) * matrix(rnorm(K * n), K, n)
  X <- rbind(a, b); rownames(X) <- paste0("g", seq_len(2 * K))
  ed <- data.frame(source = paste0("g", 1:K), target = paste0("g", (K + 1):(2 * K)))
  r1 <- multiCorr(X, ed, methods = "GCC", pernum = 500, seed = 7)
  r2 <- multiCorr(X, ed, methods = "GCC", pernum = 500, seed = 7)
  expect_gt(mean(r1$GCC_pvalue < 0.05), 0.95)
  expect_identical(r1$GCC_pvalue, r2$GCC_pvalue)
})

test_that("rank correlations handle ties like stats::cor", {
  set.seed(3)
  X <- matrix(rpois(40 * 12, 1.5), 40, 12, dimnames = list(paste0("g", 1:40), NULL))  # many ties
  X <- X[apply(X, 1, var) > 0, ]
  ed <- data.frame(source = rownames(X)[-nrow(X)], target = rownames(X)[-1])
  r <- multiCorr(X, ed, methods = c("PCC", "SCC", "KCC"))
  ref <- function(m) mapply(function(a, b) cor(X[a, ], X[b, ], method = m), ed$source, ed$target)
  expect_equal(r$PCC, unname(ref("pearson")), tolerance = 1e-8)
  expect_equal(r$SCC, unname(ref("spearman")), tolerance = 1e-8)
  expect_equal(r$KCC, unname(ref("kendall")), tolerance = 1e-8)
})

test_that("edges with a gene absent from the expression matrix return NA", {
  d <- null_edges(8, 5, 4)
  ed <- rbind(d$ed, data.frame(source = "g1", target = "missing"))
  r <- multiCorr(d$X, ed, methods = "PCC")
  expect_equal(nrow(r), 6)
  expect_true(is.na(r$PCC[6]))
})

test_that("Cauchy combination is calibrated for dependent p-values", {
  d <- null_edges(10, 3000, 5)
  r <- multiCorr(d$X, d$ed, methods = c("PCC", "SCC", "KCC"))
  p <- apply(as.matrix(r[, c("PCC_pvalue", "SCC_pvalue", "KCC_pvalue")]), 1, cauchyCombine)
  expect_lt(mean(p < 0.05), 0.07)
})
