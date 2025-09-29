library(data.table)
library(usethis)

# Create folder for raw data if it doesn't exist
dir.create("data-raw", showWarnings = FALSE)

set.seed(123)  # For reproducibility

# ---- 1. Expression data ----
genes <- paste0("Gene", 1:50)
samples <- paste0("S", 1:9)

# Create expression table as data.table
expression <- data.table(
  gene_id = genes,
  matrix(
    round(runif(50 * 9, min = 2, max = 10), 2),
    nrow = 50,
    dimnames = list(NULL, samples)
  )
)
setnames(expression, old = names(expression)[-1], new = samples)  # keep sample names

# ---- Metadata ----
meta_data <- data.table(
  sample_id = samples,
  Condition = rep(c("f1", "f2", "f3"), each = 3)
)
sample_ids <- meta_data[[1]]
# Remove first column from data.table
meta_data[, 1 := NULL]
rownames(meta_data) <- sample_ids

# ---- Description file ----
categories <- c("TF", "Kinase", "Enzyme", "Receptor", "Other")
description.file <- data.table(
  gene_id = genes,
  Category = sample(categories, 50, replace = TRUE),
  Description = paste("Description for", genes)
)

# ---- Network ----
num_edges <- 80
edges <- data.table(
  source = sample(genes, num_edges, replace = TRUE),
  target = sample(genes, num_edges, replace = TRUE)
)
edges <- edges[source != target] # remove self-loops
edges <- unique(edges)           # remove duplicates

# ---- Save as internal list object ----
#' Example NECorr dataset
#'
#' A toy dataset containing a small gene expression matrix, sample description,
#' gene metadata, and a gene network, for demonstrating NECorr functionality.
#'
#' @format A list with 4 elements:
#' \describe{
#'   \item{expression}{A data.table (50 genes x 10 columns: gene_id + 9 samples) with random expression values.}
#'   \item{meta_data}{A data.table (9 rows) mapping sample IDs to experimental conditions (`f1`, `f2`, `f3`).}
#'   \item{description.file}{A data.table (50 rows) with gene categories and descriptions.}
#'   \item{networkFile}{A data.table of gene-gene edges (2 columns: source, target).}
#' }
#'
#' @details
#' This dataset is useful for:
#' \itemize{
#'   \item Running NECorr without external files
#'   \item Testing different conditions (`f1`, `f2`, `f3`)
#'   \item Demonstrating input flexibility (data.table or file paths)
#' }
#'
#' @examples
#' data(necorr_example)
#' str(necorr_example)
#' head(necorr_example$expression)
#'
necorr_example <- list(
  expression       = expression,
  meta_data        = meta_data,
  description.file = description.file,
  networkFile      = edges
)

# Save to data/ folder as .rda (compressed R object)
usethis::use_data(necorr_example, overwrite = TRUE)

cat("✅ Internal NECorr example dataset created as data.tables: data/necorr_example.rda\n")
