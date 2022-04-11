#' effector_significance
#'
#' @param eff.m.param table with scale parameters for gene ranking
#' @param Desc description table of the genes
#' @param j.nB predictive model based on Naives Bayes
#' @param sample.l sample name
#' @nGenes number of genes in hub gene ranking
#' @description define the ranking of the "effector" genes, downstream to the hub genes
#' @return gene.rank.e.description table of ranking and descriptions
effector_significance <- function(eff.m.param, Desc, j.nB, sample.l){
  tryCatch(
    expr = {
      # Calculate the probability that the gene are involved in the studied process
      # using a Naive Bayes ML frame trained on Arabidopsis data
      prob.pred <- suppressWarnings(predict(j.nB, type="prob", newdata = eff.m.param, threshold = 0.01))
      gene.rank.eff <- (prob.pred$posterior)[,2]
      #
      gene.rank.eff = gene.rank.eff[order(gene.rank.eff, decreasing=TRUE)]
      gene.rank.eff = as.data.frame(gene.rank.eff)
      # merge description and ranking
      gene.rank.e.description <- left_join(rownames_to_column(gene.rank.eff),
                                           rownames_to_column(Desc),
                                           by = ("rowname" = "rowname"))

      gene.rank.e.description <- as.data.frame(gene.rank.e.description)
      colnames(gene.rank.e.description)[c(1,2)] <- c("Gene", "Effector.prob")
      #colnames(gene.rank.eff)[1] <- sample.l
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
#' @description define the hub edge ranking, interaction ranking
#' @return res
hub_edge_significance <- function(network.int=network.int, gene.rank.hash=gene.rank.hash){
  tryCatch(
    expr = {
      sourceIDs <- as.vector(network.int[,1])
      targetIDs <- as.vector(network.int[,2])
      score <- rep(1,length(targetIDs))
      # create the sum of the significance from each node from the linked by a edge
      for(i in 1:nrow(network.int)){
        score[i] <- sum(gene.rank.hash[[sourceIDs[i]]], gene.rank.hash[[targetIDs[i]]])
      }
      hub.int.ranks <- data.frame(sourceIDs,targetIDs,score)
      # generate a ranking "tree"
      break.points <- c(-Inf, unique(sort(as.numeric(hub.int.ranks[,3]))), Inf)
      p2 <- cut( as.numeric(hub.int.ranks[,3]), breaks=break.points, labels=FALSE )
      # use the ranking to generate a p-value of the hub
      p2 <- 1 - p2/length(break.points)
      hub.int.ranks <- as.data.frame(cbind(hub.int.ranks,p2))
      # add the p-value to final table
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
#' @description define the edge ranking around the downstream effectors of the hub genes
#' @return res
effector_edge_significance <- function(network.int=network.int,
                                       gene.rank.eff=gene.rank.eff,
                                       nGenes=nGenes){
  tryCatch(
    expr = {
      targetIDs <- as.vector(network.int[,2])
      eff_ranks.sum <- rep(1,length(targetIDs))

      gene.rank.eff.hash <- hash()
      gene.rank.eff = as.data.frame(gene.rank.eff[,c(1,2)])
      geneIDs <- gene.rank.eff[,1]
      for (i in 1:nGenes){
        geneRank <- as.numeric(gene.rank.eff[i,2])
        gene.rank.eff.hash[[geneIDs[i]]] <- geneRank
      }
      sourceIDs <- as.vector(network.int[,1])
      targetIDs <- as.vector(network.int[,2])
      eff_ranks.sum <- rep(1,length(targetIDs))
      for(i in 1:nrow(network.int)){
        eff_ranks.sum[i] <- sum(gene.rank.eff.hash[[sourceIDs[i]]], gene.rank.eff.hash[[targetIDs[i]]])
      }
      eff.int.ranks <- data.frame(sourceIDs,targetIDs,eff_ranks.sum)
      # calculation of the p-values
      break.points <- c(-Inf, unique(sort(as.numeric(eff.int.ranks[,3]))), Inf)
      p2 <- cut( as.numeric(eff.int.ranks[,3]), breaks=break.points, labels=FALSE )
      p2 <- 1 - p2/length(break.points)
      eff.int.ranks <- as.data.frame(cbind(eff.int.ranks,p2))
      eff.int.ranks$p2 <- as.numeric(as.character(eff.int.ranks$p2))
      ###
      #eff_density <- density(eff.int.ranks$eff_ranks.sum)
      #plot(eff_density)
      #plot(eff.int.ranks$eff_ranks.sum, -log10(eff.int.ranks$p2))
      ####
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
#' @description define the activator significance, upstream to the hub genes
#' @return res
activator_significant <- function(hub.int.significant=hub.int.significant,
                                  network.int=network.int,
                                  Desc=Desc){
  tryCatch(
    expr = {
      # find genes that are significant in the hub sub-network in the complete network
      sig.hub <- unique(c(as.character(hub.int.significant$sourceIDs),
                          as.character(hub.int.significant$targetIDs)))
      sc.sig.hub <- subset(network.int, network.int[,1] %in% sig.hub)
      tg.sig.hub <- subset(network.int, network.int[,2] %in% sig.hub)
      net.extension.sig.hub <- rbind(sc.sig.hub,tg.sig.hub)

      ## change the names of the hub genes in the extended hub network using source genes
      change.sig.hubNames <- net.extension.sig.hub[,1] %in% sig.hub
      tmp <- as.character(net.extension.sig.hub[,1])
      tmp[change.sig.hubNames] <- "sig"
      net.extension.sig.hub[,1] <- tmp

      ## change the names of the hub genes in the extended hub network using target genes
      change.sig.hubNames <- net.extension.sig.hub[,2] %in% sig.hub
      tmp <- as.character(net.extension.sig.hub[,2])
      tmp[change.sig.hubNames] <- "sig"
      net.extension.sig.hub[,2] <- tmp

      ## network extended to putative regulators using source genes
      act.m.param <- net.extension.sig.hub
      sc.count <- rle(sort( act.m.param[,1] ))
      act.m.param$Count <- sc.count[ match( act.m.param[,1] , sc.count )]

      # gene ranking of linked to hub nodes
      gene.rank.act <- cbind(sc.count$values, sc.count$lengths)
      gene.rank.act <- gene.rank.act[order(as.numeric(gene.rank.act[,2]), decreasing=TRUE),]
      gene.rank.act <- gene.rank.act[ - which(gene.rank.act[,1] == "sig"),]

      # add the genes that have 20% of the genes linked to hub genes
      gene.rank.act.significant <- gene.rank.act[
        which(gene.rank.act[,2] >= (length(sig.hub)*0.20)),]
      gene.rank.act.significant <- as.data.frame(gene.rank.act.significant)
      # write the putative activator genes
      gene.rank.act.b <- as.data.frame(gene.rank.act)
      colnames(gene.rank.act.b)[c(1,2)] <- c("Gene", "Nb.Connections.with.Hub")
      Desc <- as.data.frame(rownames_to_column(Desc))
      colnames(Desc)[1] <- "Gene"
      gene.rank.act.description <- left_join(gene.rank.act.b, Desc,
                                             by = "Gene")

      gene.rank.act.description <- as.data.frame(gene.rank.act.description)
      res <- list(rank=gene.rank.act.significant, description=gene.rank.act.description)
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
#' @description join the hub significant interactions and the putative upstream activators
#' @return act.net
linked_act_hub_net <- function(hub.int.significant=hub.int.significant,
                               gene.rank.act.significant=gene.rank.act.significant,
                               network.int=network.int){
  tryCatch(
    expr = {
      # isolation of genes present in the hub significant list
      sig.hub <- unique(c(as.character(hub.int.significant[,1]),
                          as.character(hub.int.significant[,2])))
      # take all the network interactions that contains the genes in the significant hubs
      sc.sig.hub <- subset(network.int, network.int[,1] %in% sig.hub)
      tg.sig.hub <- subset(network.int, network.int[,2] %in% sig.hub)
      net.extension.sig.hub <- rbind(sc.sig.hub,tg.sig.hub)
      # find the significant of all interaction containing this gene
      # eliminate the duplicate rows
      net.extension.sig.hub <- net.extension.sig.hub[!duplicated(net.extension.sig.hub), ]
      net.extension.sig.hub$score <- 0
      matching.hubgene.sc <- as.character(hub.int.significant[,1])
      matching.hubgene.tg <- as.character(hub.int.significant[,2])
      for(i in 1:nrow(net.extension.sig.hub)){
        #matching the hub interaction that significant with the extended network
        search.gene <- as.character(net.extension.sig.hub[i,1])
        search1 <- which(matching.hubgene.sc %in% search.gene)
        search2 <- which(matching.hubgene.tg %in% search.gene)
        row.pre.mean <- unique(c(search1, search2))
        # extract the mean for the interactions part of the hub gene significant
        if(length(row.pre.mean) > 0){
          int.premean <- hub.int.significant[row.pre.mean,]
          meanI <- mean(as.numeric(as.character(int.premean[,3])))
          if(net.extension.sig.hub[i,"score"] == 0){
            net.extension.sig.hub[i,"score"] <- meanI
            }else{
            net.extension.sig.hub[i,"score"] <- (net.extension.sig.hub[i,"score"] + meanI)/2
            }
          }
        }
        #
      act.net.1 <- subset(
        net.extension.sig.hub,
        net.extension.sig.hub[,1] %in% as.vector(as.character(gene.rank.act.significant[,1])))
      act.net.2 <- subset(
        net.extension.sig.hub,
        net.extension.sig.hub[,2] %in% as.vector(as.character(gene.rank.act.significant[,1])))
      act.net.pre <- rbind(act.net.1,act.net.2)
      # eliminate putative duplicates
      act.net.pre <- act.net.pre[!duplicated(act.net.pre), ]
      act.net.pre <- act.net.pre[order(act.net.pre$score, decreasing = TRUE),]
      act.net <- cbind(act.net.pre,
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
#' @description find the downstream effector genes that are linked to the significant hub genes
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
      # create table with all the genes in network
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
      fit2 <- suppressWarnings(eBayes(fit2))
      # Exports complete limma statistics table for first comparison group ('coef=1')
      # to tab delimited text file.
      de <- topTable(fit2, coef=1, adjust.method="fdr", sort.by="P", number=length(sel.Genes))
      ranks <- ScalN(de$B)
      names(ranks) = rownames(de)
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
NetTopology <-function(network){
  tryCatch(
    expr = {
      # pass to igraph the network that would transform in a graph object
      g <- graph.data.frame(network, directed = T)
      ## calculate the topology stats necessary for network statistics
      betweenness_centrality <- betweenness(g,v=V(g),directed = F, normalized = T)
      eigenvector_centrality<-evcent(g, scale = TRUE)
      # former Google algorithm for web-page search; this is an evolution of eigen-vector
      pagerank<-page.rank(g)$vector
      # degree centrality
      degree<-degree(g, v=V(g), mode = "total",normalized = T)
      # in and out degree
      degree_in<-degree(g, v=V(g), mode = "in")
      degree_out<-degree(g, v=V(g), mode = "out")
      # transitivity of gene in the network, also called coefficient cluster
      transitivity_centrality <- transitivity( g, vids=V(g),type ="local", isolates= "zero")
      # preparation to transform the network statistics in a table
      gene_name<-V(g)$name
      # bind into matrices
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


