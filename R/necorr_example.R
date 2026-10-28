#' Example NECorr data set
#'
#' A small simulated data set: a regulatory network of 12 transcription factors
#' and their targets, and log2 expression of 120 genes in four conditions
#' (root, leaf, flower, seed; four samples each). \code{Gene001} and its
#' targets form a co-expressed module that is induced in flower. The same data
#' are available as text files in \code{system.file("extdata", package = "NECorr")}
#' (\code{example_network.txt}, \code{example_expression.txt},
#' \code{example_gene_description.txt}, \code{example_sample_metadata.txt}).
#'
#' @format A list with 4 data frames, named after the arguments of \code{\link{NECorr}}:
#' \describe{
#'   \item{networkFile}{edges: \code{source} (transcription factor) and \code{target}}
#'   \item{expression}{\code{GeneID} followed by one column per sample}
#'   \item{description.file}{\code{gene_id}, \code{Category}, \code{Description}}
#'   \item{meta_data}{\code{SampleID} and \code{Condition}}
#' }
#' @examples
#' data(necorr_example)
#' str(necorr_example, max.level = 1)
"necorr_example"
