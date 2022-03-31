#' DiscNULLrow
#' @description Discard rows of expression with only NULL expression
#' @param eset entry table
#' @return res table with null column
#' @export
DiscNULLrow <- function(eset){
  notNULLrows <-  apply(eset, 1, function(row) sum(row)!=0)
  res <- eset[notNULLrows, ]
  return(res)
}


#' factOrg
#' @description reorganize the column from the list of Tissue
# from phenoData !!!! organized in the same order
#' @param eset expression data
#' @param pDataTissue vector containing the sample names
#' @return out from list the factor list of the sample type
#' @export
factOrg  <- function(eset, pDataTissue){
  # @retrun a list with the nset reorganized by sample type
  # defined by pDataTissue)
  # @return from list the factor list of the sample type
  tryCatch(
    expr = {
      sample.names <- unique(pDataTissue)
      xcol <-c()
      orgFactors <- factor()
      for (i in 1:length(sample.names)){
        nrep <- length(which(pDataTissue==sample.names[i]))
        if(nrep == 1){
          # if a condition has no replicate generate false replicates.
          # this will need to be simulated usign the variance from all samples to be
          # accurate
          xcol <- c(xcol,rep(i,3))
          orgFactors <-
            factor(c(as.character(orgFactors),as.character(rep(sample.names[i],3))))
        } else if(nrep > 1){
          xcol <- c(xcol,which(pDataTissue==sample.names[i]))
          orgFactors <-
            factor(c(as.character(orgFactors),as.character(rep(sample.names[i],nrep))))
        }
      }
      orgFactors <- factor(orgFactors, levels = sample.names)
      m.eset = as.matrix(eset)
      m.eset <- m.eset[,xcol]
      out <- list(eset=m.eset, orgFact=orgFactors)
      return(out)
    },
    error = function(e){
      message("Error in the expression column reorganization")
      message(e)
    },
    warning = function(w){
      message("Warning in the expression column reorganization")
      message(w)
    },
    finally = {
    }
  )
  sample.names <- unique(pDataTissue)
  xcol <-c()
  orgFactors <- factor()
  for (i in 1:length(sample.names)){
    nrep <- length(which(pDataTissue==sample.names[i]))
    if(nrep == 1){
      # if a condition has no replicate generate false replicates.
      # this will need to be simulated usign the variance from all samples to be
      # accurate
      xcol <- c(xcol,rep(i,3))
      orgFactors <-
        factor(c(as.character(orgFactors),as.character(rep(sample.names[i],3))))
    } else if(nrep > 1){
      xcol <- c(xcol,which(pDataTissue==sample.names[i]))
      orgFactors <-
        factor(c(as.character(orgFactors),as.character(rep(sample.names[i],nrep))))
    }
  }
  orgFactors <- factor(orgFactors, levels = sample.names)
  m.eset = as.matrix(eset)
  m.eset <- m.eset[,xcol]
  out <- list(eset=m.eset, orgFact=orgFactors)
  return(out)
}


#' TSI
#' @description tissue-selective upregulated genes from Tissue Specificty index (Yanai et al, 2005)
#' @param x vector of value needed to be "TSIed"
#' @return res vector scaled to the TSI
#' @export
TSI <- function(x){
  tryCatch(
    expr = {
      res <- sum(1 - x/max(x))/(length(x) - 1)
      return(res)
    },
    error = function(e){
      message("Error in the Tissue Specific Index calculation")
      message(e)
    },
    warning = function(w){
      message("Warning in the Tissue Specific Index calculation")
      message(w)
    },
    finally = {
    }
  )

}


#' TSE
#' @description tissue-selective downregulated genes (Liseron et al., 2013)
#' @param x vector of value needed to be "TSEed"
#' @return res
#' @export
TSE <- function(x){
  tryCatch(
    expr = {
      res <- sum( 1 - (min(x)/max(x)))/(length(x) - 1)
      return(res)
    },
    error = function(e){
      message("Error in the Tissue Specific Exclusion")
      message(e)
    },
    warning = function(w){
      message("Warning in Tissue Specific Exclusion")
      message(w)
    },
    finally = {
    }
  )
}


#' IUTtsi
#' @description define tissue specific expression threshold
#' to apply to detect tissue-specificity
#' This could avoid very lowly expressed genes that could be artifact
#' @param data expression data
#' @param geneTS genes with tissue specific proprieties
#' @param targetSample target sample name
#' @param expThreshold threshold of the expression for RNAseq around 10 for microarray
#' could be different according to the platform
#' @return y.sort
#' @export
IUTtsi <- function(data, geneTS, targetSample, expThreshold){
  tryCatch(
    expr = {
      tsi.order <- c()
      Glen <- length(geneTS[which(geneTS == TRUE)])
      if(Glen > 0){
        geneTS <- geneTS[which(geneTS == TRUE)]
        if(Glen > 1){
          tsi <- apply(data[names(geneTS), ], 1, TSI)
          tsi.order <- tsi[order(tsi, decreasing = TRUE)]
          y <- data[names(tsi.order),]
          y.sort <- as.matrix(y[which(y[, targetSample] > expThreshold), ])
          #
          y.sort <- as.data.frame(y.sort)
          y.sort$RowNames <- row.names(y.sort)
          #
          tsi.order <- as.data.frame(tsi.order);
          tsi.order$RowNames <- row.names(tsi.order)
          #
          y.sort <- left_join(y.sort, tsi.order, by = 'RowNames')
          y.sort <- as.data.frame(y.sort)
          row.names(y.sort) <- y.sort$RowNames
          y.sort <- y.sort[,-which(colnames(y.sort) == "RowNames")]
          # y.sort <- as.matrix(y.sort)
        }else if(Glen == 1){ # if there is only one tissue-selective genes
          tsi.order <- TSI(data[which(geneTS==TRUE), ])
          names(tsi.order) <- names(geneTS)
          y.sort <- t(as.matrix(data[names(tsi.order), ]))
          #
          rownames(y.sort) <- names(geneTS)
          y.sort <- as.data.frame(y.sort)
          y.sort$RowNames <- row.names(y.sort)
          #
          tsi.order <- as.data.frame(tsi.order);
          tsi.order$RowNames <- row.names(tsi.order)
          #
          y.sort <- left_join(y.sort, tsi.order, by = 'RowNames')
          y.sort <- as.data.frame(y.sort)
          row.names(y.sort) <- y.sort$RowNames
          y.sort <- y.sort[,-which(colnames(y.sort) == "RowNames")]
          y.sort <- as.matrix(y.sort)
        }

      }else{
        y.sort <- matrix(nrow = 0, ncol = length(colnames(data))+1,
                         dimnames = list(c(), c(colnames(data), "tsi.order")))
      }
      return(y.sort)
    },
    error = function(e){
      message("Error in IUTtsi")
      message(e)
    },
    warning = function(w){
      message("Warning in IUTtsi")
      message(w)
    },
    finally = {
    }
  )
}


#' IUTtse
#' @description define the tissue excluded genes from the tissue/condition of interest
#' @param data expression data
#' @param geneTS genes with tissue specific proprieties
#' @param targetSample target sample name
#' @param expThreshold threshold of the expression for RNAseq around 10 for microarray
#' could be different according to the platform
#' @return y.sort
#' @export
IUTtse <- function(data, geneTS, targetSample, expThreshold){
  tryCatch(
    expr = {
      tse.order <- c()
      Glen <- length(geneTS[which(geneTS == TRUE)])
      if(Glen > 0){
        geneTS <- geneTS[which(geneTS == TRUE)]
        if(Glen > 1){
          tse <- apply(data[names(geneTS), ], 1, TSE)
          tse.order <- tse[order(tse, decreasing = TRUE)]
          y <- data[names(tse.order),]
          y.sort <- as.matrix(y[which(y[, targetSample] > expThreshold), ])
          #
          y.sort <- as.data.frame(y.sort)
          y.sort$RowNames <- row.names(y.sort)
          #
          tse.order <- as.data.frame(tse.order);
          tse.order$RowNames <- row.names(tse.order)
          #
          y.sort <- left_join(y.sort, tse.order, by = 'RowNames')
          y.sort <- as.data.frame(y.sort)
          row.names(y.sort) <- y.sort$RowNames
          y.sort <- y.sort[,-which(colnames(y.sort) == "RowNames")]
        }else if(Glen == 1){ # if there is only one tissue-selective genes
          tse.order <- TSE(data[which(geneTS==TRUE), ])
          names(tse.order) <- names(geneTS)
          y.sort <- t(as.matrix(data[names(tse.order), ]))
          rownames(y.sort) <- names(geneTS)
          #
          y.sort <- as.data.frame(y.sort)
          y.sort$RowNames <- row.names(y.sort)
          #
          tse.order <- as.data.frame(tse.order)
          tse.order$RowNames <- row.names(tse.order)
          #
          y.sort <- left_join(y.sort, tse.order, by = 'RowNames')
          y.sort <- as.data.frame(y.sort)
          row.names(y.sort) <- y.sort$RowNames
          y.sort <- y.sort[,-which(colnames(y.sort) == "RowNames")]
          # y.sort <- as.matrix(y.sort)
        }
      }else{
        y.sort <- matrix(nrow = 0, ncol = length(colnames(data))+1,
                         dimnames = list(c(), c(colnames(data),"tse.order")))
      }
      return(y.sort)
    },
    error = function(e){
      message("Error in IUTtse")
      message(e)
    },
    warning = function(w){
      message("Warning in  IUTtse")
      message(w)
    },
    finally = {
    }
  )

}


#' IUTest
#' @description perform the complete IUT for a gene expression table
#' @param m.eset expression matrix
#' @param sFactors number of sample replicates
#' @param sIndex index of the sample in the table
#' @param alpha significance level for the IUT
#' @return res
#' @export
IUTest <- function(m.eset, sFactors, sIndex, alpha){
  tryCatch(
    expr = {
      gene <- nrow(m.eset)
      nrarrays <- ncol(m.eset)
      nbTreat <- length(sFactors)
      # variance calculation for each gene
      fL <- 1
      meanTab <- matrix(rep(0, gene*nbTreat), nrow=gene)
      SSdevs <- rep(0, gene)
      for (i in (1:nbTreat)){
        upper <- fL + sFactors[i]-1
        m.eset.i <- m.eset[, fL:upper]
        # mean calculation
        meanTab[, i] <- apply(m.eset.i, 1, mean)
        # variance calculation
        SSdev <- apply(m.eset.i, 1, function(x) var(x)*(length(x)-1))
        SSdevs <- as.matrix(SSdevs + SSdev)
        fL <- upper + 1
      }

      # t-test between sample and each other sample
      Nin  <- sFactors[sIndex]
      Noth <- sFactors[-sIndex]
      meanSamp <- meanTab[, sIndex]
      meanAll <- meanTab[, -sIndex]
      SD <- SSdevs/(nrarrays-nbTreat)
      SD[SD == 0] <- 0.05 # arbitrary SD if equal to 0
      SE <- (SD %*% (1/Nin + 1/Noth))^(1/2)
      tstattab <- (meanSamp - meanAll)/SE
      # critical value for the two side t-test to compare with each t-test
      #alpha <- 0.05
      thresh <- 1-(1 - alpha)^(1/gene)
      crit.value <- qt(1-thresh, nrarrays-sFactors[sIndex])
      # find genes with all significant t-test for the chosen sample compared each
      # other ones.
      if(isTRUE(nrow(tstattab) > 0)){
        StudentSumP <- rowSums(tstattab > crit.value)
        resultpos <- StudentSumP == (nbTreat-1)
        StudentSumN <- rowSums(tstattab < - crit.value)
        resultmin <- StudentSumN == (nbTreat-1)
      }else{
        resultpos <- NULL
        resultmin <- NULL
      }
      res <- list(TS = resultpos, TE = resultmin)
      return(res)
    },
    error = function(e){
      message("Error in IUTest")
      message(e)
    },
    warning = function(w){
      message("Warning in IUTest")
      message(w)
    },
    finally = {
    }
  )

}

#' ts.IUT
#' @description combine the TSI/TSE test and IUT
#' to determine the tissue specific genes and the tissue excluded genes.
#' list of list with 4 methods raw and filtered results matrix
#' of tissue specific genes for each method; upload all libraries gplots
#' TCC
#' upload function below
#' TissueTransform, DiscNULLrow, factOrg, IUTest, numeratorTSI, denominator,
#' TSI, tsi.matrix,
#' plotlOnegene, plotLfewGenes, heatmap.ts, tindex, tindex.matrix,
#' @param name name of the data
#' @param eset matrix of expression data in matrix format in TPM for example
#' but any other expression matrix
#' @param tissues character vector defining the the tissue type of each
#' column present in the expression table
#' @param target the target tissues
#' @param threshold represent the p-value for the IUT test after Sidak
#' correction; the outlier percentage for the ROKU test
#' the percentage of (1 -threshold for the Tindex and the TSI tests)
#' @param filter the minimal expression level for the tissue/condition of
#' interest; the default is set at 100; expression of 100 TPM
#' @return list.all csv.table output for each type of analysis will be written in the
#' table folder
#' @export
ts.IUT <- function(name = "exp", eset, tissues, target,
                   threshold = 0.05, filter = 50){
  tryCatch(
    expr = {
      # Reorganize Data
      eset <- DiscNULLrow(eset) ##if eset in log 2 scale retransform it 2**eset to
      ##have the right scale
      # Discard  NULL expression
      # Reorganize expression data for the tissue specific tests
      org.nset <- factOrg(eset, tissues)
      m.eset <- org.nset$eset
      # Create means of expression for each tissue to get a compact
      # representation
      m.eset <- as.data.frame(m.eset)
      meansFactor <- do.call("cbind", tapply(1:ncol(m.eset), org.nset$orgFact,
                                             function(x) rowMeans(m.eset[x])))
      m.eset <- as.matrix(m.eset)
      # Intersection union test - using Intersection-Union Test (Berger et al) -
      # TSI
      rankedTSE <- c()
      rankedTSI <- c()
      factors <- table(org.nset$orgFact)
      samples <- unique(org.nset$orgFact)
      t.iut <- IUTest(m.eset, factors, as.numeric(which(samples==target)),  ############
                      threshold)
      tm.iut <- meansFactor[names(t.iut$TS)[which(t.iut$TS == TRUE)], ]
      m.iut.tsi.filt <- IUTtsi(meansFactor, t.iut$TS, target, filter)
      # Rank the tissue selective genes using the Tissue Selective Index
      # for each gene TSI for activation and modify for tissue repression
      if((length(rownames(m.iut.tsi.filt)) > 1) == TRUE){
        m.iut.tsi.filt <- m.iut.tsi.filt[which(m.iut.tsi.filt[,target] > filter), ]
        others1 <- meansFactor[setdiff(rownames(meansFactor), rownames(m.iut.tsi.filt)), ]
        others1 <- cbind(others1, tsi.order = rep(0, nrow(others1)))
        rankedTSI <- rbind(m.iut.tsi.filt, others1)
      }else{
          rankedTSI <- cbind(meansFactor, tsi.order = rep(0, nrow(meansFactor)))
      }
      m.iut.e <- meansFactor[names(t.iut$TE)[which(t.iut$TE == TRUE)], ]
      m.iut.tse.filt <- IUTtse(meansFactor, t.iut$TE, target, filter)
      # Rank the tissue selective genes using the Tissue Selective Index
      # for each gene TSI for activation and modify for tissue repression
      if((length(rownames(m.iut.tse.filt)) > 1) == TRUE){
        m.iut.tse.filt <- m.iut.tse.filt[which(m.iut.tse.filt[,target] > filter), ]
        others2 <- meansFactor[setdiff(rownames(meansFactor), rownames(m.iut.tse.filt)), ]
        others2 <- cbind(others2, tse.order = rep(0, nrow(others2)))
        rankedTSE <- rbind(m.iut.tse.filt, others2)
      }else{
          rankedTSE <- cbind(meansFactor, tse.order = rep(0, nrow(meansFactor)))
      }
      # Results
      listall <- list(filtTSI = rankedTSI, filtTSE = rankedTSE)
      return(listall)
    },
    error = function(e){
      message("Error in ts.IUT")
      message(e)
    },
    warning = function(w){
      message("Warning in ts.IUT")
      message(w)
    },
    finally = {
    }
  )
}

