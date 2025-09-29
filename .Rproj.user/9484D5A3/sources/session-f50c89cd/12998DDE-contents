# ============================
# NECorr Package Test Script
# ============================
remove.packages("NECorr")
rm(list=ls())
# rm -rf /Users/liseronc/Library/R/arm64/4.3/library/NECorr
setwd("~/Documents/1_projects/NECorr/")
devtools::clean_dll()
# !! restart session within the RStudio IDE (Session > Restart R)
# install.packages("devtools")
devtools::document()
devtools::build_vignettes()
#devtools::compileAttributes()
devtools::check()
devtools::install()

# 1. Load your package
library(NECorr)
#run the tutorial
#vignette("NECorr_tutorial")

# 2. Load example data
source("data-raw/make_necorr_example.R")
data(necorr_example)

# 3. Run NECorr
results <- NECorr(
  networkFile     = necorr_example$network,
  expression      = necorr_example$expression,
  description.file= necorr_example$description,
  condition       = "f1",
  meta_data       = necorr_example$meta_data,
  permutation     = 100,
  save_results = FALSE,
  interactive_net = FALSE
)
print(results)

# test with mutiple correlation score
res <- NECorr(
  networkFile     = necorr_example$network,
  expression      = necorr_example$expression,
  description.file= necorr_example$description,
  condition       = "f1",
  meta_data       = necorr_example$meta_data,
  method = c("PCC", "SCC", "KCC", "GCC"),
  permutation     = 100,
  save_results = FALSE,
  interactive_net = FALSE
)
print(res)
# 4. Visualize interactively
plots <- NECorr::visualize_necorr_results(res, top_n = 20, interactive_net = TRUE, output_dir = "plots")
# Explore interactively
plots$network_plot

# 5. Visualize the results
plots <- visualize_necorr_results(res, top_n = 20, interactive_net = FALSE, output_dir = "NECorr_plots")
# Show top hubs
plots$top_hubs_plot
# Show network
plots$network_plot
# Co-expression scatter
plots$coexpression_plot
# Degree distribution
plots$degree_distribution

