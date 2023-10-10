#' effector_significance
#'
#' @param eff.m.param table with scale parameters for gene ranking
#' @param Desc description table of the genes 
#' @param j.nB predictive model based on Naives Bayes
#' @param sample.l sample name
#' @nGenes number of genes in hub gene ranking
#' @description define the ranking of the effector genes 
#' @return gene.rank.e.description table of ranking and descriptions
effector_significance <- function(eff.m.param, Desc, j.nB, sample.l){
  tryCatch(
    expr = {
      prob.pred <- suppressWarnings(predict(j.nB, type="prob",
                                            newdata = eff.m.param,
                                            threshold = 0.01))
      gene.rank.eff <- (prob.pred$posterior)[,2]
      #
      gene.rank.eff = gene.rank.eff[order(gene.rank.eff, decreasing=TRUE)]
      gene.rank.eff = as.data.frame(gene.rank.eff)
      # merge description and ranking
      gene.rank.e.description <- left_join(rownames_to_column(gene.rank.eff),
                                           rownames_to_column(Desc),
                                           by = ("rowname" = "rowname"))
      
      gene.rank.e.description <- as.data.frame(gene.rank.e.description)
      colnames(gene.rank.eff)[1] <- sample.l
      return(gene.rank.e.description)
    },
    error = function(e){ 
      message("NECorr error in the significance")
      message(e)
    },
    warning = function(w){
      message("NECorr warning in the significance ranking of the effectors")
      message(w)
    },
    finally = {
      
    }
  )
  
}


#' hub_edge_significance
#'
#' @param network.int network file
#' @param gene.rank.hash gene ranking hash
#' @description define the hub edge ranking
#' @return res 
hub_edge_significance <- function(network.int=network.int, gene.rank.hash=gene.rank.hash){
  tryCatch(
    expr = {
      sourceIDs <- as.vector(network.int[,1])
      targetIDs <- as.vector(network.int[,2])
      ranks.sum <- rep(1,length(targetIDs))
      for(i in 1:nrow(network.int)){
        ranks.sum[i] <- sum(gene.rank.hash[[sourceIDs[i]]], gene.rank.hash[[targetIDs[i]]])
      }
      hub.int.ranks <- data.frame(sourceIDs,targetIDs,ranks.sum)
      break.points <- c(-Inf, unique(sort(as.numeric(hub.int.ranks[,3]))), Inf)
      p2 <- cut( as.numeric(hub.int.ranks[,3]), breaks=break.points, labels=FALSE ) 
      p2 <- 1 - p2/length(break.points)
      hub.int.ranks <- as.data.frame(cbind(hub.int.ranks,p2))
      hub.int.ranks$p2 <- as.numeric(as.character(hub.int.ranks$p2))
      res <- hub.int.ranks
      return(res)
    },
    error = function(e){ 
      message("Error in the hub significance calculation:")
      message(e)
    },
    warning = function(w){
      message("Warning in hub significance calculation:")
      message(w)
    },
    finally = {}
  )
}


#' effector_edge_significance
#'
#' @param network.int network data
#' @param gene.rank.eff.hash gene ranking for the effector in a hash
#' @description define the edge ranking around the effectors
#' @return res
effector_edge_significance <- function(network.int=network.int, 
                                       gene.rank.eff=gene.rank.eff, nGenes=nGenes){
  tryCatch(
    expr = {
      targetIDs <- as.vector(network.int[,2])
      eff_ranks.sum <- rep(1,length(targetIDs))
      
      gene.rank.eff.hash <- hash()
      gene.rank.eff = as.data.frame(gene.rank.eff)
      geneIDs <- row.names(gene.rank.eff)
      for (i in 1:nGenes){
        geneRank <- as.numeric(gene.rank.eff[i,1])
        gene.rank.eff.hash[[geneIDs[i]]] <- geneRank
      }
      sourceIDs <- as.vector(network.int[,1])
      targetIDs <- as.vector(network.int[,2])
      eff_ranks.sum <- rep(1,length(targetIDs))
      for(i in 1:nrow(network.int)){
        eff_ranks.sum[i] <- sum(gene.rank.eff.hash[[sourceIDs[i]]], gene.rank.eff.hash[[targetIDs[i]]])
      }
      eff.int.ranks <- data.frame(sourceIDs,targetIDs,eff_ranks.sum)
      break.points <- c(-Inf, unique(sort(as.numeric(eff.int.ranks[,3]))), Inf)
      p2 <- cut( as.numeric(eff.int.ranks[,3]), breaks=break.points, labels=FALSE )
      p2 <- 1 - p2/length(break.points)
      eff.int.ranks <- as.data.frame(cbind(eff.int.ranks,p2))
      eff.int.ranks$p2 <- as.numeric(as.character(eff.int.ranks$p2))
      res <- eff.int.ranks
      return(res)
    },
    error = function(e){ 
      message("Error in effector edge significance:")
      message(e)
    },
    warning = function(w){
      message("Warnings in effector edge significance:")
      message(w)
    },
    finally = {
    }
  )
}

#' activator_significant
#'
#' @param hub.int.significant significance of the hub genes
#' @param network.int network edges
#' @Desc description file genes and gene names
#' @description define the activator significance
#' @return res 
activator_significant <- function(hub.int.significant = hub.int.significant, 
                                  network.int = network.int,
                                  Desc = Desc){
  tryCatch(
    expr = {
      # find genes that are significant in the hub subnetwork in the complete network
      sig.hub <- unique(c(as.character(hub.int.significant$sourceIDs), 
                          as.character(hub.int.significant$targetIDs)))
      sc.sig.hub <- subset(network.int, network.int[,1] %in% sig.hub)
      tg.sig.hub <- subset(network.int, network.int[,2] %in% sig.hub)
      net.extension.sig.hub <- rbind(sc.sig.hub,tg.sig.hub)
      g <- graph.data.frame(net.extension.sig.hub, directed = T)
      ## calculate the pagerank known as regulator
      pagerank<-page.rank(g)$vector
      # degree centrality
      degree<-degree(g, v=V(g), mode = "total",normalized = T)
      # out degree
      degree_out<-degree(g, v=V(g), mode = "out")
      # preparation to transform the network statistics in a table
      gene_name<-V(g)$name
      # bind into matrices
      datagrid <- data.frame(I(gene_name), pagerank, degree, degree_out)
      cc <- c("gene_name", "pagerank", "EdgeCount", "degree_out")
      colnames(datagrid) <- cc
      rownames(datagrid) <- gene_name
      datagrid <- datagrid[order(datagrid$pagerank, decreasing = T),]
      datagrid$pagerank <- ScalN(datagrid$pagerank)
      gene.rank.act.description <- left_join(rownames_to_column(datagrid), 
                                             rownames_to_column(Desc), 
                                             by = ("rowname" = "rowname"))
      rownames(gene.rank.act.description) <- gene.rank.act.description$rowname
      gene.rank.act.description <- gene.rank.act.description[, -c(1,2)]
      res <- gene.rank.act.description
      return(res)
    },
    error = function(e){ 
      message("Error in activator significant search:")
      message(e)
    },
    warning = function(w){
      message("Warning in activator significant search:")
      message(w)
    },
    finally = {
    }
  )
}

#' linked_act_hub_net
#'
#' @param hub.int.significant significant network of the hub genes
#' @param gene.rank.act.significant the ranking of the activator genes
#' @param network.int full gene network
#' @return act.net
linked_act_hub_net <- function(hub.int.significant=hub.int.significant,
                               gene.rank.act.significant=gene.rank.act.significant,
                               network.int=network.int){
  tryCatch(
    expr = {
      sig.hub <- unique(c(as.character(hub.int.significant[,1]), 
                          as.character(hub.int.significant[,2])))
      sc.sig.hub <- subset(network.int, network.int[,1] %in% sig.hub)
      tg.sig.hub <- subset(network.int, network.int[,2] %in% sig.hub)
      net.extension.sig.hub <- rbind(sc.sig.hub,tg.sig.hub)
      act.net.1 <- subset(
        net.extension.sig.hub,
        net.extension.sig.hub[,1] %in% as.vector(as.character(gene.rank.act.significant[,1])))
      act.net.2 <- subset(
        net.extension.sig.hub,
        net.extension.sig.hub[,2] %in% as.vector(as.character(gene.rank.act.significant[,1])))
      act.net.pre <- rbind(act.net.1,act.net.2)
      meanSig <- mean(as.vector(as.numeric(as.character(hub.int.significant[,3]))))
      act.net <- cbind(act.net.pre,
                       rep(meanSig,nrow(act.net.pre)),
                       rep("act",nrow(act.net.pre)))
      colnames(act.net) <- c("source","target","score","node.type")
      return(act.net)
    },
    error = function(e){ 
      message("Error in linking hub genes significant and activator significant genes:")
      message(e)
    },
    warning = function(w){
      message("Warning in linking hub genes significant and activator significant genes:")
      message(w)
    },
    finally = {
    }
  )
}

#' linked_eff_hub_net
#'
#' @param hub.int.significant significant network of the hub genes
#' @param eff.int.significant significant network of the effector genes
#'
#' @return eff.net
linked_eff_hub_net <- function(hub.int.significant=hub.int.significant,
                               eff.int.significant=eff.int.significant){
  tryCatch(
    expr = {
      sig.hub <- unique(c(as.character(hub.int.significant[,1]), 
                          as.character(hub.int.significant[,2])))
      sc.sig.eff <- subset(eff.int.significant, eff.int.significant[,1] %in% sig.hub)
      tg.sig.eff <- subset(eff.int.significant, eff.int.significant[,2] %in% sig.hub)
      eff.net.pre <- rbind(sc.sig.eff,tg.sig.eff)
      eff.net <- cbind(eff.net.pre[,1],eff.net.pre[,2],eff.net.pre[,3], rep(nrow(eff.net.pre)))
      return(eff.net)
    },
    error = function(e){ 
      message("Error in linking hub genes significant and effector significant genes:")
      message(e)
    },
    warning = function(w){
      message("Warning in linking hub genes significant and effector significant genes:")
      message(w)
    },
    finally = {
    }
  )
}

#' indexing.network
#' @param tab expression table
#' @param network network needed to be indexed
#' @return netindex
#' @export
indexing.network <- function (tab,network) {
  tryCatch(
    expr = {
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
    },
    error = function(e){ 
      message("Error in the Network indexing")
      message(e)
    },
    warning = function(w){
      message("Warning in Network indexing")
      message(w)
    },
    finally = {
    }
  )

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
  tryCatch(
    expr = {
      if (exps.file == TRUE){
        df <- as.matrix(read.table(exps, header=TRUE,sep="\t",row.names = 1))
        expression <- df[,colSums(is.na(df)) != nrow(df)]
        targets <- read.table(factortab,header=TRUE,sep="\t",row.names=1)
        sample.names <- unique(targets[,1])
        f <- factor(targets[,1], levels=sample.names)
      }else{
        expression <- as.matrix(exps)
        f <- factor(factortab)
      }
      # parse expression data to only have the gene that are
      # in the network - load the genes that are in the networks
      sel.Genes <- intersect(rownames(expression),GeneList)
      no.duplicate.col <- as.character(colnames(expression))
      colnames(expression) <- make.unique(no.duplicate.col)
      expressio <- expression[sel.Genes, ]
      design <- model.matrix(~ 0 + f)
      # DE (differential expression)
      eset <- ExpressionSet(assayData=expressio)
      colnames(design) <- gsub("^f", "", colnames(design)) # Assigns column names.
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
      fit2 <- suppressWarnings(eBayes(fit2))
      #write.table(topTable(fit2, coef=1, adjust="fdr", sort.by="B",
      # number=50000), file="limma_complete.xls", row.names=F, sep="\t")
      # Exports complete limma statistics table for first comparison group ('coef=1')
      # to tab delimited text file.
      de <- topTable(fit2, coef=1, adjust.method="fdr", sort.by="P", number=length(sel.Genes))
      ranks <- ScalN(de$B)
      names(ranks) = rownames(de)
      #print(head(ranks))
      return(ranks)
    },
    error = function(e){ 
      message("Error in the Differential Expression ranking")
      message(e)
    },
    warning = function(w){
      message("Warning in the Differential Expression ranking")
      message(w)
    },
    finally = { 
    }
  )
}


#' NetTopology
#' @description network topology statistics using igraph methods;
#' script partly adapted from https://gist.github.com/mhawksey/16-2306
#' @author Christophe Liseron-Monfils, Martin Hawksey
#' @param network the network should be a list of interactions space by a line
#' @return datagrid
#' @export
NetTopology <- function(network){
  tryCatch(
    expr = {
      # pass to igraph the network that would transform in a graph object
      g <- graph.data.frame(network, directed = T)
      ## calculate the topology stats necessary for network statistics
      betweenness_centrality <- betweenness(g,v=V(g), directed = F, normalized = T)
      #betweenness_centrality <- betweenness.estimate(g,vids=V(g), directed = F, normalized = T)
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
    },
    error = function(e){ 
      message("Error in the network topology parameters")
      message(e)
    },
    warning = function(w){
      message("Warning in the network topology parameters")
      message(w)
    },
    finally = {
    }
  )

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

#' scaling_param
#' @description
#' rescaling the data from the parameter
#' @param x vector with names of the genes
#' @return res
#' @export
scaling_param <- function(x){
  naming <- deparse(substitute(x))
  x1 <- range(x[!is.infinite(x)])
  x[is.infinite(x) & sign(x) < 0] <- x1[1]
  x[is.infinite(x)] <- x1[2]
  x <- as.data.frame(x)
  x$x <- rescale(x$x, to=c(0,100))
  colnames(x) <- naming
  return(x)
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


