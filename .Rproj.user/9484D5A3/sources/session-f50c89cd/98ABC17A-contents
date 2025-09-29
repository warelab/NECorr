# Load necessary libraries
library(Rcpp)
library(data.table)
setwd("~/Documents/1_projects/NECorr/")
# Set seed for reproducibility
set.seed(42)

# Compile and load the C++ code
Rcpp::sourceCpp("src/multi_corr_necorr.cpp")

# Generate and prepare test data
set.seed(42)
num_genes <- 50
num_samples <- 10
num_edges <- 35

expression_data <- matrix(rnorm(num_genes * num_samples), nrow = num_genes)
rownames(expression_data) <- paste0("Gene", 1:num_genes)
colnames(expression_data) <- paste0("Sample", 1:num_samples)

edges <- data.frame(
  source = paste0("Gene", sample(1:num_genes, num_edges, replace = TRUE)),
  target = paste0("Gene", sample(1:num_genes, num_edges, replace = TRUE))
)

edges <- edges[edges$source != edges$target, ]

expression_data[1, ] <- expression_data[2, ] + rnorm(num_samples, sd = 0.1)
expression_data[3, ] <- expression_data[4, ] * 0.1 + rnorm(num_samples, sd = 1)
expression_data[5, ] <- -expression_data[6, ] + rnorm(num_samples, sd = 0.1)

ranks <- matrix(
  apply(expression_data, 1, function(row) rank(row, ties.method = "first") - 1),
  nrow = nrow(expression_data), byrow = TRUE
)
# Convert gene names to indices
gene_index <- setNames(seq_len(num_genes) - 1, rownames(expression_data))
src_indices <- as.integer(gene_index[edges$source])
tgt_indices <- as.integer(gene_index[edges$target])

# Call the multi_corr_necorr function
results <- multi_corr_necorr(
  expression = expression_data,
  ranks = ranks,
  src = src_indices,
  tgt = tgt_indices,
  bootstrapIterations = 1000,
  useBestGCC = FALSE,
  asymmetricGCC = TRUE
)

# Print the results
print(results)
