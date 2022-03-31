
#' multiCorr
#' @description the function is calculating co-expression
#' @param x the table of expression
#' @param net the network to test for co-expression
#' @param verbose enumerate all steps
#' @param pernum number of permutation for significance test
#' @return corMAT
#' @export
multiCorr <- function(x="" ,net="", verbose=TRUE,  pernum=0){
  tryCatch(
    expr = {
      corMAT <- c()
      # RCPP Gini Correlation calculations
      corMAT <- giniR(edges=net, expression=x, bootstrapIterations=pernum, statCutoff=0.6)
      return(corMAT)
    },
    error = function(e){
      message("Error during multiCorr running")
      message(e)
    },
    warning = function(w){
      message("Warning during multiCorr running")
      message(w)
    },
    finally = {
    }
  )
}

#' GiniR
#' @description conversion of Rcpp into R function
#' @param edges the network edges
#' @param expression the expression tables
#' @param bootstrapIterations the parameter of bootstrap for randomization
#' @param statCutoff significance cut off
#' @useDynLib NECorr
#' @export
giniR <- function(edges, expression, bootstrapIterations, statCutoff) {
  res <- gini(edges, expression, bootstrapIterations, statCutoff)
  return(res)
}
