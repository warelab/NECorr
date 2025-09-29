#' Example NECorr dataset
#'
#' A toy dataset containing a small gene expression matrix, sample description,
#' gene metadata, and a gene network, for demonstrating NECorr functionality.
#'
#' @format A list with 4 elements:
#' \describe{
#'   \item{expression}{A data frame (50 genes x 10 columns: GeneID + 9 samples) with random expression values.}
#'   \item{description}{A data frame (9 rows) mapping sample IDs to experimental conditions (`f1`, `f2`, `f3`).}
#'   \item{metadata}{A data frame (50 rows) with gene categories and descriptions.}
#'   \item{network}{A data frame of gene-gene edges (2 columns: Gene1, Gene2).}
#' }
#'
#' @details
#' This dataset is useful for:
#' \itemize{
#'   \item Running NECorr without external files
#'   \item Testing different conditions (`f1`, `f2`, `f3`)
#'   \item Demonstrating input flexibility (data frames or file paths)
#' }
#'
#' @examples
#' data(necorr_example)
#' str(necorr_example)
#' head(necorr_example$expression)
#'
"necorr_example"
