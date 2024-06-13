#' multiCorr
#' @description the function is calculating co-expression using Gini Correlation
#' @param x the table of expression
#' @param net the network to test for co-expression
#' @param sigmethod significance of the p-value for the chosen method
#' @param verbose enumerate all steps
#' @param pernum number of permutation for significance test
#' @return corMAT a correlation matrix
#' @export
multiCorr <- function(x ,net= NA,
                      nblocks = 10, verbose = TRUE, pernum = 0, sigmethod = 0.6){
  tryCatch(
    expr = {
      corMAT <- c()
      runtime_giniR <- system.time(corMAT <- giniR(edges=net, expression=x, bootstrapIterations=pernum, statCutoff=sigmethod))
      cat("Runtime giniR (seconds):", runtime_giniR["elapsed"], "\n")
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

#' expand.grid.unique
#' @description expand the data
#' @param x x parameter
#' @param y y parameter
#' @param include.equals including equals
#' @return res
#' @export
expand.grid.unique <- function(x, y, include.equals=FALSE){
  x <- unique(x)
  y <- unique(y)
  g <- function(i){
    z <- setdiff(y, x[seq_len(i-include.equals)])
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  res <- do.call(rbind, lapply(seq_along(x), g))
  return(res)
}


#' GiniR
#' @description convertion of Rcpp into R function
#' @param edges the network edges
#' @param expression the expression tables
#' @param bootstrapIterations the parameter of boothtrap for randomization
#' @param statCutoff significance cut off
#' @useDynLib NECorr
#' @export
giniR <- function(edges, expression, bootstrapIterations, statCutoff) {
  res <- gini(edges, expression, bootstrapIterations, statCutoff)
  return(res)
}


