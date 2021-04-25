#' cor.matrix.NECorr
#' @param GEMatrix expression matrix
#' @param cpus number of cpu used
#' @param cormethod correlation method:"GCC", "PCC", "SCC", "KCC"
#' @param style type of pairs: "pairs.between","pairs.only"
#' @param var1.id variable vector source
#' @param var2.id varaible vector target
#' @param pernum permutation number
#' @param sigmethod significance
#' @param output matrix or paired
#' @return resList
#' @export
cor.matrix.NECorr <- function (GEMatrix, cpus = cpus,
                               cormethod = c("GCC", "PCC", "SCC", "KCC"),
                               style = c("pairs.between","pairs.only"),
                               var1.id = NA, var2.id = NA,
                               pernum = 0,
                               sigmethod = c("two.sided","one.sided"),
                               output = c("matrix","paired")){
  tryCatch(
    expr = {
      if (cpus > 1) {
        #   library(snowfall)
      }
      if (pernum == 0){
        sigmethod <- "two.sided"
      }
      if (!is.matrix(GEMatrix) | !is.numeric(GEMatrix)){
        stop("Error: GEMatrix in cor.matrix function is not matrix or Error: GEMatrix is not numeric")
      }
      if (length(rownames(GEMatrix)) == 0) {
        rownames(GEMatrix) <- seq(1, dim(GEMatrix)[1], by = 1)
      }
      VariableNum <- nrow(GEMatrix)
      SampleSize <- ncol(GEMatrix)
      if (VariableNum <= 1 || SampleSize <= 1) {
        stop("Error:the number of variable is less than 2, or the number of observation is less than 2")
      }
      if (style == "pairs.between" || style =="pairs.only") {
        if (length(which(is.na(var1.id) == TRUE)) > 0 | length(which(is.na(var1.id) ==  TRUE)) > 0) {
          stop("Error: no variable IDs are given")
        }
        if (length(which(is.numeric(var1.id) == FALSE)) > 0 | length(which(is.numeric(var2.id) == FALSE)) > 0) {
          stop("Error:var1.id and var2.id should be numeric vector")
        }
      }
      if(style == "pairs.only"){
        taskmatrix <- cbind(var1.id, var2.id)
      }
      if(style == "pairs.between") {
        df <- expand.grid.unique(var1.id, var2.id, include.equals = TRUE)
        df2 = t(apply(df, 1, sort))
        taskmatrix <- df2[!duplicated(df2),]
      }
      pernum = as.numeric(pernum)
      results <- apply(taskmatrix, 1, cor.pair, GEMatrix = GEMatrix,
                       rowORcol = "row", cormethod = cormethod, pernum = pernum,
                       sigmethod = sigmethod)
      if(output == "paired"){
        kk <- 0
        corpvalueMatrix <- matrix(NA, nrow = dim(taskmatrix)[1], ncol = 4)
        for (i in 1:dim(taskmatrix)[1]) {
          if (taskmatrix[i, 1] == taskmatrix[i, 2]) {
            next
          }
          kk <- kk + 1
          corpvalueMatrix[kk, 1:2] <- rownames(GEMatrix)[taskmatrix[i, ]]
          if (cormethod == "GCC") {
            fGCC <- gcc.corfinal(results[i][[1]])
            corpvalueMatrix[kk, 3] <- fGCC$gcc.fcor
            corpvalueMatrix[kk, 4] <- fGCC$gcc.fpvalue
          }
          else {
            corpvalueMatrix[kk, 3] <- results[i][[1]]$cor
            corpvalueMatrix[kk, 4] <- results[i][[1]]$pvalue
          }
        }
        return(corpvalueMatrix[1:kk, ])
      }else{
        UniqueRow <- sort(unique(taskmatrix[, 1]))
        UniqueCol <- sort(unique(taskmatrix[, 2]))
        corMatrix <- matrix(0, nrow = length(UniqueRow), ncol = length(UniqueCol))
        rownames(corMatrix) <- rownames(GEMatrix)[UniqueRow]
        colnames(corMatrix) <- rownames(GEMatrix)[UniqueCol]
        pvalueMatrix <- corMatrix
        pvalueMatrix[] <- NA
        for (i in 1:dim(taskmatrix)[1]) {
          rowidx <- which(UniqueRow == taskmatrix[i, 1])
          colidx <- which(UniqueCol == taskmatrix[i, 2])
          if (cormethod == "GCC") {
            fGCC <- gcc.corfinal(results[i][[1]])
            corMatrix[rowidx, colidx] <- fGCC$gcc.fcor
            pvalueMatrix[rowidx, colidx] <- fGCC$gcc.fpvalue
          }
          else {
            corMatrix[rowidx, colidx] <- results[i][[1]]$cor
            pvalueMatrix[rowidx, colidx] <- results[i][[1]]$pvalue
          }
          # fill up the inverse part of the table
          corMatrix[colidx, rowidx] <- corMatrix[rowidx, colidx]
          pvalueMatrix[colidx, rowidx] <- pvalueMatrix[rowidx, colidx]
        }
        resList <- list(corMatrix = corMatrix, pvalueMatrix = pvalueMatrix)
        return(resList)
      }
    },
    error = function(e){ 
      message("Error in cor.matrix.NECorr")
      message(e)
    },
    warning = function(w){
      message("Warning in cor.matrix.NECorr ")
      message(w)
    },
    finally = {
    }
  )
  
 
}

#' multiCorr
#' @description the function is calculating co-expression
#' @param x the table of expression
#' @param net the network to test for co-expression
#' @param nsockets number of parallel sockets
#' @param methods correlation method
#' @param sigmethod significance of the p-value for the chosen method
#' @param nblocks number of blocks to be cut
#' @param verbose enumerateall steps
#' @param cpus number of cpu used
#' @param pernum number of permutation for singificance test
#' @param ... other parameters
#' @return corMAT a correlation matrix
#' @export
multiCorr <- function(x ,net= NA, nsockets= 4, methods = c("GCC","PCC","SCC","KCC"),
                      sigmethod = c("two.sided", "one.sided"),
                      nblocks = 10, verbose = TRUE, cpus = 1, pernum = 0, ...){
  tryCatch(
    expr = {
      corMAT <- c()
      if (methods == "GCC"){
        #### NEW GINI CORRELATION CALCULATION
        corMAT <- giniR(edges=net, expression=x, bootstrapIterations=pernum, statCutoff=0.6)
      }else{
        #suppressWarnings(suppressPackageStartupMessages(require(foreach)))
        #suppressWarnings(suppressPackageStartupMessages(require(doSNOW)))
        nsockets <- as.numeric(nsockets)
        cl <- makeCluster(nsockets, type="SOCK")
        registerDoSNOW(cl)
        
        Nrow = nrow(net)
        fc <- gl(nblocks, ceiling(Nrow/nblocks), length = Nrow)
        matnet <- split(net,fc)
        j <- 1:nblocks
        #corMAT<-foreach(j=1:nblocks, .combine='rbind',
        corMAT<-foreach(j, .combine='rbind',
                        .export=c('indexing.network','cor.matrix.NECorr'))%dopar%{
                          netindexed <- indexing.network(as.matrix(x) ,matnet[[j]])
                          G1 <- as.numeric(as.vector(netindexed[,1]))
                          G2 <- as.numeric(as.vector(netindexed[,3]))
                          cor.matrix.NECorr(x, var1.id=G1, var2.id=G2, sigmethod=sigmethod, cormethod=methods,
                                            pernum=pernum, output="paired", cpus=cpus, style="pairs.only")
                        }
        stopCluster(cl)
      }
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

#' cor.pair
#' @param idxvec vector of the index
#' @param GEMatrix Gene expression network
#' @param rowORcol row or columns
#' @param cormethod correlation type "GCC","PCC", "SCC" or "KCC"
#' @param pernum permutation number
#' @param sigmethod p-value threshold
#' @return res paired correlations
cor.pair <- function (idxvec, GEMatrix, rowORcol = c("row", "col"),
                      cormethod = c("GCC","PCC", "SCC", "KCC"),
                      pernum = 0, sigmethod = c("two.sided", "one.sided")){
  if (!is.vector(idxvec) | length(idxvec) != 2) {
    stop("Error: idxvec must be a vector with two elements indicating the indexs(rows) in GEMatrix")
  }
  if (class(GEMatrix) != "matrix") {
    stop("Error: GEMatrix should be a numeric data matrix")
  }
  if (rowORcol == "row") {
    x1 <- GEMatrix[idxvec[1], ]
    y1 <- GEMatrix[idxvec[2], ]
  }
  else if (rowORcol == "col") {
    x1 <- GEMatrix[, idxvec[1]]
    y1 <- GEMatrix[, idxvec[2]]
  }
  else {
    stop("Error: rowORcol must be \"row\" or \"col\"")
  }
  if (length(cormethod) > 1) {
    stop("Error: length of cormethod must be of length 1")
  }
  if (!is.vector(x1) || !is.vector(y1)) {
    stop("Error: input two vectors for gcc.corpair")
  }
  if (!is.numeric(x1) || !is.numeric(y1)) {
    stop("Error: x should be numeric")
  }
  if (length(which(is.na(x1) == TRUE)) > 0 || length(which(is.na(y1) ==
                                                           TRUE)) > 0) {
    stop("Error: There are Na(s) in x")
  }
  if (length(x1) != length(y1)) {
    stop("Error: the lengths of each row in x are different.\n")
  }
  realcor <- getcor(x1, y1, cormethod)
  if (pernum <= 0) {
    if (cormethod == "GCC") {
      res <- list(gcc.rankx = realcor$gcc.rankx, gcc.ranky = realcor$gcc.ranky,
                  gcc.rankx.pvalue = NA, gcc.ranky.pvalue = NA)
    }
    else {
      res <- list(cor = realcor, pvalue = NA)
    }
  }
  else {
    pGCCMatrix <- matrix(0, nrow = pernum, ncol = 2)
    colnames(pGCCMatrix) <- c("gcc.rankx", "gcc.ranky")
    rownames(pGCCMatrix) <- paste("permut", seq(1, pernum,
                                                by = 1), sep = "")
    GenePairXY <- t(matrix(c(x1, y1), ncol = 2))
    pGenePairXY <- GenePairXY
    Length <- length(x1)
    for (i in 1:pernum) {
      curtime <- format(Sys.time(), "%H:%M:%OS4")
      XXX <- unlist(strsplit(curtime, ":"))
      curtimeidx <- (as.numeric(XXX[1]) * 3600 + as.numeric(XXX[2]) *
                       60 + as.numeric(XXX[3])) * 10000
      set.seed(curtimeidx)
      TT = sort(runif(Length), index.return = TRUE)$ix
      pGenePairXY[1, ] <- GenePairXY[1, TT]
      if (cormethod == "GCC") {
        cortmp <- getcor(pGenePairXY[1, ], pGenePairXY[2, ],
                         cormethod)
        pGCCMatrix[i, 1] <- cortmp$gcc.rankx
        pGCCMatrix[i, 2] <- cortmp$gcc.ranky
      }
      else {
        pGCCMatrix[i, 1] <- getcor(pGenePairXY[1, ],
                                   pGenePairXY[2, ], cormethod)
      }
    }
    if (cormethod == "GCC") {
      res <- list(gcc.rankx = realcor$gcc.rankx, gcc.ranky = realcor$gcc.ranky,
                  gcc.rankx.pvalue = getpvalue(pGCCMatrix[, 1],
                                               pernum, realcor$gcc.rankx, sigmethod),
                  gcc.ranky.pvalue = getpvalue(pGCCMatrix[,2],
                                               pernum, realcor$gcc.ranky, sigmethod))
    }
    else {
      res <- list(cor = realcor, pvalue = getpvalue(pGCCMatrix[,1],
                                                    pernum, realcor, sigmethod))
    }
  }
  return(res)
}

#' gcc.corfinal
#' @param gcccor table of the rank gcc
#' @return res table of gcc score and gcc pvalue
gcc.corfinal <- function (gcccor){
  if (is.numeric(gcccor$gcc.rankx.pvalue) & is.numeric(gcccor$gcc.ranky.pvalue)){
    fpvalue <- gcccor$gcc.rankx.pvalue
    fgcc <- gcccor$gcc.rankx
    if (gcccor$gcc.ranky.pvalue < fpvalue) {
      fpvalue <- gcccor$gcc.ranky.pvalue
      fgcc <- gcccor$gcc.ranky
    }
  }
  else {
    fpvalue <- NA
    x <- c(gcccor$gcc.rankx, gcccor$gcc.ranky)
    fgcc <- x[which(abs(x) == max(abs(x)))]
    if (length(fgcc) > 1) {
      fgcc <- fgcc[1]
    }
  }
  res <- list(gcc.fcor = fgcc, gcc.fpvalue = fpvalue)
  return(res)
}

#' getcor
#' @description from package XXX - reference
#' @param g1 gcc first
#' @param g2 gcc second
#' @param cormethod method
#' @return res
getcor <- function(g1, g2, cormethod){
  if (cormethod == "PCC") {
    res <- cor.test(g1, g2, method = "pearson")$estimate
  }
  else if (cormethod == "SCC") {
    res <- cor.test(g1, g2, method = "spearman")$estimate
  }
  else if (cormethod == "KCC") {
    res <- cor.test(g1, g2, method = "kendall")$estimate
  }
  else if (cormethod == "GCC") {
    res <- list(gcc.rankx = onegcc(g1, g2),
                gcc.ranky = onegcc(g2,g1))
  }
  return(res)
}

#' getpvalue
#' @description from package XXX - reference
#' @param percorvec percentage corelation vector
#' @param pernum permutaion number
#' @param realcor real correlation
#' @param sigmethod significance method
#' @return pvalue
getpvalue <- function(percorvec, pernum, realcor, sigmethod){
  pvalue <- length(which(percorvec >= realcor))/pernum
  if (pvalue == 0)
    pvalue <- 1/pernum
  if (pvalue > 0.5)
    pvalue <- 1 - pvalue
  if (sigmethod == "two.sided") {
    pvalue <- 2 * pvalue
  }
  return(pvalue)
}


#' gcc
#' @description fom package XXX - reference
#' @param weightvec weight
#' @param vectsort vector sorted
#' @param vectselfsort selection of sorted
#' @return gcccor
gcc <- function(weightvec, vectsort, vectselfsort){
  Sum1 <- sum(weightvec * vectsort)
  Sum2 <- sum(weightvec * vectselfsort)
  if (Sum2 == 0) {
    gcccor <- 0
  }
  else {
    gcccor <- Sum1/Sum2
  }
  return(gcccor)
}

#' onegcc
#' @description from package XXXX  reference
#' @param x first vector
#' @param y second vector
#' @return gcc.rankx
onegcc <- function(x, y){
  Length <- length(x)
  Wt <- t(2 * seq(1, Length, by = 1) - Length - 1)
  GenePairXY <- t(matrix(c(x, y), ncol = 2))
  SortYlist <- getrank(GenePairXY)
  gcc.rankx <- gcc(Wt, SortYlist$Sort2By1, SortYlist$Sort2By2)
  return(gcc.rankx)
}

#' getrank
#' @param datamatrix matrix with the data
#' @return res ranked data
getrank <- function(datamatrix){
  if (dim(datamatrix)[1] != 2) {
    stop("Error: the row num of datamatrix must be 2")
  }
  OrderIndex <- order(datamatrix[1, ], decreasing = FALSE)
  SortGenePair <- datamatrix[, OrderIndex]
  Sort2By1 <- SortGenePair[2, ]
  res <- list(Sort2By1 = Sort2By1, Sort2By2 = sort(datamatrix[2,
  ], decreasing = FALSE))
  return(res)
}
