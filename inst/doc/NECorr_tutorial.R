## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 4.5)
library(NECorr)

## -----------------------------------------------------------------------------
data(necorr_example)
head(necorr_example$networkFile, 3)
necorr_example$expression[1:3, 1:5]
head(necorr_example$meta_data, 3)
head(necorr_example$description.file, 3)

## -----------------------------------------------------------------------------
list.files(system.file("extdata", package = "NECorr"), pattern = "^example_")

## -----------------------------------------------------------------------------
res <- NECorr(
  networkFile      = necorr_example$networkFile,
  expression       = necorr_example$expression,
  description.file = necorr_example$description.file,
  meta_data        = necorr_example$meta_data,
  condition        = "flower",
  method           = "GCC",
  permutation      = 200,
  directed         = TRUE,
  n_null           = 50,
  visualize        = FALSE,
  save_results     = FALSE,
  verbose          = FALSE
)

## -----------------------------------------------------------------------------
head(res$necorrHub_nodes, 8)

## -----------------------------------------------------------------------------
round(head(res$hub.m.param[res$necorrHub_nodes$Gene, ], 5), 2)
round(res$criterion_weights, 3)

## -----------------------------------------------------------------------------
head(res$necorrReg[, c("Gene", "role", "n_partners", "n_hub_partners", "fold_enrichment", "FDR")], 5)

## -----------------------------------------------------------------------------
head(res$necorrEdges, 5)

## -----------------------------------------------------------------------------
str(res$settings[c("weights", "aggregation", "direction", "scaling", "degree_type", "directed")])

## -----------------------------------------------------------------------------
c3 <- c("Coexpr", "TS", "CondCoexpr")
A <- necorr_ahp_matrix(c3, preset = "wiring")
round(A, 2)
ahp <- ahp_weights(A)
round(ahp$weights, 2)
ahp$CR   # consistency ratio, should be below 0.1
res_wiring <- NECorr(
  necorr_example$networkFile, necorr_example$expression, necorr_example$description.file,
  "flower", necorr_example$meta_data, method = "PCC", permutation = 0,
  criteria = c3, aggregation = "product", scaling = "softcap", weights = A,
  directed = TRUE, visualize = FALSE, save_results = FALSE, verbose = FALSE)
round(res_wiring$criterion_weights, 2)
head(res_wiring$coexpression_stats[order(-res_wiring$coexpression_stats$Coexpr), c("Gene", "k", "n_strong", "Coexpr", "CondCoexpr")], 3)

## -----------------------------------------------------------------------------
res_alt <- NECorr(
  necorr_example$networkFile, necorr_example$expression, necorr_example$description.file,
  "flower", necorr_example$meta_data, method = "PCC", permutation = 0,
  criteria = c("Coexpr", "TS", "CondCoexpr", "DE"), weights = "critic", aggregation = "sum",
  directed = TRUE, visualize = FALSE, save_results = FALSE, verbose = FALSE)
round(res_alt$criterion_weights, 3)
head(res_alt$necorrHub_nodes$Gene, 8)

## ----fig.height = 4-----------------------------------------------------------
plots <- visualize_necorr_results(res, top_n = 15)
plots$top_hubs_plot

## ----fig.height = 5-----------------------------------------------------------
plots$network_plot

