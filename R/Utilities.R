#' indexing.network
#' @param tab expression table
#' @param network network needed to be indexed
#' @return netindex
#' @export
indexing.network <- function (tab,network) {
  # create table with all the gene in network
  t <- as.data.frame(unique(apply(network,1, function(x) c(x[1],x[2]))))
  tt1 <- t(t)
  #replace the gene name by their index in expression table
  indexing <- function (x,y) {
    if ( length(which(rownames(y) == x[1])) != 0 && length(which(rownames(y) == x[2])) != 0) {
      c(which(rownames(y) == x[1]),which(rownames(y) == x[2]))
    }
  }
  netindex <-  do.call(rbind,apply(tt1, 1, indexing,tab))
  return(netindex)
}


#' DE.ranking
#' @description differential expression ranking using the network gene
#'   1-expression file with NECorr format
#'   2-list of network genes
#'   3-expression factor table generate automatically
#'   4-name of the studied sample
#'   read expressiontable from microarrays or RNAseq
#' @param exps expression table
#' @param GeneList gene list present in the network
#' @param factortab table with factors
#' @param sample.l the sample chosen for tissue specificity
#' @param sample.names all the samples in the metadata
#' @param exps.file the expression file
#' @import Biobase
#' @import stats
#' @import limma
#' @return ranks rank of the gene in DEG analysis
#' @export
DE.ranking <- function(exps, GeneList, factortab, sample.l,
                       sample.names, exps.file = FALSE){
  if (exps.file == TRUE){
    df <- as.matrix(read.table(exps, header=TRUE,sep="\t",row.names = 1))
    expression <- df[,colSums(is.na(df)) != nrow(df)]
    targets <- read.table(factortab,header=TRUE,sep="\t",row.names=1)
    sample.names <- unique(targets[,1])
    f <- factor(targets[,1], levels=sample.names)
  }
  else {
    expression <- as.matrix(exps)
    f <- factor(factortab)
  }
  # parse expression data to only have the gene that are
  # in the network - load the genes that are in the networks
  sel.Genes <- intersect(rownames(expression),GeneList)
  no.duplicate.col <- as.character(colnames(expression))
  colnames(expression) <- make.unique(no.duplicate.col)
  expressio <- expression[sel.Genes,]
  design <- model.matrix(~ 0 + f)
  # DE (differential expression)
  eset <- ExpressionSet(assayData=expressio)
  colnames(design) <- sample.names # Assigns column names.
  fit <- lmFit(eset, design)
  other.names = sample.names[which(sample.names != sample.l)]
  equation <- noquote(paste0(noquote(sample.l), " - ((",
           paste(other.names, collapse = '+'), ") /",
           length(sample.names)-1,")"))
  assign(".equation", equation, envir = .GlobalEnv)
  #contrast.matrix <- makeContrasts(equation , levels=design)
  contrast.matrix <- makeContrasts(.equation, levels=design)
  remove(".equation", envir= .GlobalEnv)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  # Computes estimated coefficients and standard errors for a given set of contrasts.
  fit2 <- eBayes(fit2)
  #write.table(topTable(fit2, coef=1, adjust="fdr", sort.by="B",
  # number=50000), file="limma_complete.xls", row.names=F, sep="\t")
  # Exports complete limma statistics table for first comparison group ('coef=1')
  # to tab delimited text file.
  de <- topTable(fit2, coef=1, adjust.method="fdr", sort.by="P", number=length(sel.Genes))
  ranks <- ScalN(de$B)
  names(ranks) = rownames(de)
  return(ranks)
}


#' NetTopology
#' @description network topology statistics using igraph methods;
#' script partly adapted from https://gist.github.com/mhawksey/16-2306
#' @author Christophe Liseron-Monfils, Martin Hawksey
#' @param network the network should be a list of interactions space by a line
#' @return datagrid
#' @export
NetTopology <-function(network){
  #suppressWarnings(suppressPackageStartupMessages(require(igraph)))
  # pass to igraph the network that would transform in a graph object
  g <- graph.data.frame(network, directed = T)
  ## calculate the topology stats necessary for network statistics
  betweenness_centrality <- betweenness(g,v=V(g),directed = F, normalized = T)
  #betweenness_centrality <- betweenness.estimate(g,vids=V(g),directed = F, normalized = T)
  eigenvector_centrality<-evcent(g, scale = TRUE)
  # algorithm used by google to find webpages; this is an evolution of eigen-vector
  # not sure if applicable to biological data yet
  pagerank<-page.rank(g)$vector
  # degree centrality
  degree<-degree(g, v=V(g), mode = "total",normalized = T)
  # in and out degree
  degree_in<-degree(g, v=V(g), mode = "in")
  degree_out<-degree(g, v=V(g), mode = "out")
  # transitivity of gene in the network, also called coefficient cluster
  transitivity_centrality <- transitivity( g, vids=V(g),type ="local", isolates= "zero")
  # preparation to transfrom the network statistics in a table
  gene_name<-V(g)$name
  # bind into matrice
  datagrid <- data.frame(I(gene_name),degree,degree_in,degree_out,
                         betweenness_centrality,eigenvector_centrality[[1]],
                         transitivity_centrality,pagerank)
  cc <- c("gene_name", "EdgeCount","degree_in","degree_out","BetweennessCentrality",
          "eigenvector_centrality","ClusteringCoefficient","pagerank")
  colnames(datagrid) <- cc
  rownames(datagrid) <- gene_name
  #eliminate the first column being a duplicate of rownames now
  datagrid <- datagrid[,-1]
  return(datagrid)
}


#' ScalN
#' @description rescaling data in a range between 0 to 1
#' @param x vector that need to be rescaled
#' @param ... other paramter that can be passed on
#' @return res
#' @export
ScalN <- function(x, ...) {
  res <- (x - min(x, ...))/(max(x, ...) - min(x, ...))
  return(res)
}


#' fishersMethod
#' @description combine p-values Fisher method to have overall importance of
#' this node in the studied tissue (using code from Michael Love code)
#' (http://mikelove.wordpress.com/2012/03/12/combining-p-values-fishers-method-sum-of-p-values-binomial/)
#' @param x vector of p-values
#' @return res
#' @export
fishersMethod <- function(x) {
  res <- pchisq(-2 * sum(log(x)),df=2*length(x),lower.tail=FALSE)
  return(res)
}
