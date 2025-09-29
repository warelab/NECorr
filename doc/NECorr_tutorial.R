## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)
library(NECorr)

## ----eval=FALSE---------------------------------------------------------------
# networkFile <- system.file("extdata", "example_network.tsv", package = "NECorr")
# expressionFile <- system.file("extdata", "example_expression.tsv", package = "NECorr")
# descriptionFile <- system.file("extdata", "example_desc.tsv", package = "NECorr")
# metadataFile <- system.file("extdata", "example_metadata.tsv", package = "NECorr")
# condition <- "Leaf"
# 
# head(net)
# head(expr)
# head(desc)
# head(meta)
# 

## ----eval=FALSE---------------------------------------------------------------
# results <- NECorr(
#   networkFile = networkFile,
#   expression = expressionFile,
#   description.file = descriptionFile,
#   condition = condition,
#   metadata = metadataFile,
#   methods = c("PCC", "SCC", "KCC"),
#   visualize = TRUE,
#   save_results = TRUE,
#   output_dir = "NECorr_results",
#   top_n = 15,
#   interactive_net = FALSE
# )

## ----eval=FALSE---------------------------------------------------------------
# # View top hub genes
# head(results$necorrHub_nodes, 10)
# 
# # View top regulators
# results$necorrReg
# 
# # View significant edges
# results$necorrEdges
# 

## ----eval=FALSE---------------------------------------------------------------
# results$plots$top_hubs_plot

## ----eval=FALSE---------------------------------------------------------------
# results$plots$network_plot

## ----eval=FALSE---------------------------------------------------------------
# if (!is.null(results$plots$coexpression_plot)) results$plots$coexpression_plot

## ----eval=FALSE---------------------------------------------------------------
# results$plots$degree_distribution

## ----eval=FALSE---------------------------------------------------------------
# results <- NECorr(
#   networkFile, expressionFile, descriptionFile, condition, metadataFile,
#   visualize = TRUE, save_results = TRUE, output_dir = "NECorr_results",
#   interactive_net = TRUE
# )
# results$plots$network_plot
# 

