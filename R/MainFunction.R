#' Necorr
#' @author Christophe Liseron-Monfils, Andrew Olson
#' @param expression Expression file in log2 (ratio expression) with row: gene,
#' first column: type of sample,second column: sample names
#' @param networkFile Molecular network file with source in the first column, targets in
#'  the second column
#' @param description.file genome description
#' @param condition Condition from expression to study the network co-expression
#' correlation
#' @param metadata dataframe with the metadata
#' @param permutation permutation number used for all significance calculation
#' #param lmiR List of miRNAs
#' @param method used for co-expression correlation: GCC, MINE, PCC, SCC or KCC
#' @param sigcorr significance of the correlation
#' @param NSockets number of sockets
#' @param verbose the number for the warn message
#' @description NECorr helps discover candidate genes that could be
#' important for specific conditions.
#' The principal inputs are the expression data and the network file.
#' The expression data should start with 3 header columns.
#' The first column describes the conditions. Each condition will be
#' treated separately for the co-expression analysis
#' The output of the program will be generated in a result folder generated
#' in the working path
#' Create the output directory if not existing; generate "./results" dir and
#' "./results/tmp"
#' C.Liseron-Monfils - Ware lab Sept2013 - CSHL
#' partly based on rsgcc package for the GCC, PCC,KCC and SPP Ma et al, 2012, plant Physiology
#' @return res
#' @export
Necorr <- function(networkFile = "",
                   expression = "",
                   description.file = "",
                   condition = "",
                   metadata = "", 
                   method = "GCC", permutation = 1000, 
                   sigcorr = 0.01,
                   NSockets = 2,
                   verbose= 0){
  verbose=0
  #options(warn = warn_verbose)
  # load file and options
  # read network File
  network.int <- read.delim(networkFile, sep ="\t", header = T, fileEncoding="latin1")
  network.int <- network.int[!is.na(network.int[1,]),]
  network.int <- network.int[!is.na(network.int$target),]
  
  # take only in consideration the real interactions
  if(ncol(network.int) == 2){
    network.int <- network.int[,c(1,2)]
  } else if(ncol(network.int) > 2){
    network.int <- network.int[, c(1,3)]
  } 
  # remove NA from the network data
  network.int <- network.int[!is.na(network.int[,1]),]
  network.int <- network.int[!is.na(network.int[,2]),]
  # Condition experiment and/or tissue
  factortab <- read.table(metadata, header = T, sep = "\t", row.names = 1) #factor.file #condition <- "Radial"
  # Description file name with gene name and annotations
  Desc <-  read.csv(description.file, header=T, row.name=1) #description.file <- "1.Ath.GeneDesc.csv"
  # Read expression file
  eset <- read.table(expression, header = T, row.names=1) #!! change made Here
  eset.rows <- rownames(eset)
  eset <- do.call(data.frame,lapply(eset, function(x) replace(x, is.infinite(x), 10^300))) # max number that R can handle
  rownames(eset) <- eset.rows
  # need to see if all these tables are still useful??
  #___________________ MAIN SCRIPT ___________________________________
  
  # calculate the weight for each parameter
  # (use same order for columns in gene stats table)
  #suppressWarnings(suppressPackageStartupMessages(require(pmr)))
  #preftable <-  matrix(data = rep(1,(5*5)), nrow = 5, ncol = 5, byrow = TRUE)

  #preftable[1,2] = 1; preftable[2,1] = 1 # interaction pval and expression
  #preftable[1,3] = 2; preftable[3,1] = 1/2 # interaction pval and betweenness
  #preftable[1,4] = 2; preftable[4,1] = 1/2 # interaction pval and connectivity
  #preftable[1,5] = 2; preftable[5,1] = 1/2 # interaction pval and transitivy

  #preftable[2,3] = 2; preftable[3,2] = 1/2 # expression and betweenness
  #preftable[2,4] = 2; preftable[4,2] = 1/2 # expression and connectivity
  #preftable[2,5] = 2; preftable[5,2] = 1/2 # expression and transitivy

  #preftable[3,4] = 2; preftable[4,3] = 1/2 # Betweenness and connectivity
  #preftable[3,5] = 2; preftable[5,3] = 1/2 # Betweenness and transitivy

  #preftable[4,5] = 2; preftable[5,4] = 1/2 # connectivity and transitivy
  
  # weight of parameter(wpar) give the criteria ranking
  #wpar <- ahp(preftable)$weighting
  #print (wpar)
  wpar <- c(0.2825910, 0.2825910, 0.1864405, 0.1412955, 0.1070820)
  
  # The next section is dealing with the fact
  # that some sample types can have no replicate
  # so if this is the case, pseudo-replicates will be generated
  # to be able to calculate a statistic for
  # the tissue/stress-selectivity of each sample type
  # and order these tissue/stress-selective genes
  sample.names <- unique(factortab[,1])
  #print(sample.names)
  xcol <- c()
  treatment.f <- factor()
  for (i in 1:length(sample.names)){
    #print(paste0("iteration i:", i))
    #print(paste0("iteration j:", j))
    nrep <- length(which(factortab[,1] == sample.names[i]))
    if(i==1){
      j <- i
    } 
    if(nrep == 1){
      # create 3 artificial columns if the sample has only one replicate
      xcol <- c(xcol,rep(j,3))
      j <- j+3
      treatment.f <- factor(c(as.character(treatment.f), as.character(rep(sample.names[i],3))))
    } else if(nrep > 1){
      xcol <- c(xcol,seq(j,j+(nrep-1)))
      treatment.f <- factor(c(as.character(treatment.f), as.character(rep(sample.names[i],nrep))))
      j <- j+(nrep)
    }
    #print(treatment.f)
    #print(xcol)
  }
  #print(xcol)
  treatment.f <- factor(treatment.f, levels = sample.names)
  #-----------------------------------------------------------------------------------------------
  
  #start.time <- Sys.time()
  # Generate all the topology statistics using the function NetTopology
  netname <- basename(networkFile)
  netstat <- NetTopology(network.int)
  Genelist <- rownames(netstat)
  BetwC <- structure(netstat$BetweennessCentrality, names=Genelist)  #3 betweenness
  Conn <- structure(netstat$EdgeCount, names=Genelist) # 4 connectivity
  ClusCoef <- structure(netstat$ClusteringCoefficient, names=Genelist) #5 transitivity
  EigenC <- structure(netstat$eigenvector_centrality, names=Genelist) #6 eigen-vector
  PageRank <- structure(netstat$pagerank, names=Genelist) #7 page rank
  
  # Generate an empty vector to serve as table with all the ranks for conditions
  total.rank <- matrix(data = NA, nrow = length(Genelist))
  rownames(total.rank) <- Genelist
  total.rank[,1] = Genelist
  # R network topology done 

  # Take only the genes that are part of the molecular network
  eset <- eset[intersect(Genelist,rownames(eset)),] #!! change made Here
  #eset <- eset[-grep("NA", rownames(m.eset)),]
  m.eset <- as.matrix(eset)
  m.eset <- m.eset[, xcol]
  #print(m.eset)
  # Loop to measure the importance of gene expression
  conditionList <-factortab[,2]
  factorList <- factortab[,1]
  df <- cbind(as.character(conditionList),as.character(factorList))
  df <- df[!duplicated(df),]
  #fileexp <- as.character(basename(expression))
  #end.time <- Sys.time()
  #time.taken <- end.time - start.time
  #print(time.taken)
  
  #-----------------------------------------------------------------------------------------------
  message("Analysis of the co-expression file and p-value sums")
  #start.time <- Sys.time()
  
  int.sig <- multiCorr(m.eset, net = network.int, nsockets = NSockets,
                       methods = method, output = "paired", sigmethod = "two.sided",
                       pernum = permutation , verbose = T, cpus = NSockets)
  #int.sig.file <- int.sig
  #colnames(int.sig.file) <- c("Source","Target","Correlation","p-value")
  # This adds a column of color values based on the y values
  # remove the line with NA due to not found gene due to no express/ion
  # or error of annotation in the initial network file
  int.sig <- int.sig[complete.cases(int.sig),]
  
  # Transform the p-value in the maximal p-value in function of number of permutations
  # to avoid log-scale = infinity
  # Replace the p-value equal to 0 by the maximal p-value e.g.
  # 1/10000 knowing that we have 10000 randomizations in initial calculations
  max.p.val <- 1/as.numeric(permutation)
  int.sig[which(int.sig[,4] == 0),4] <- max.p.val
  int.sig <-as.data.frame(int.sig)
  p.int.sig  <- as.numeric(as.character(int.sig[,4]))
  pval <- -log(p.int.sig,10)
  #-----------------------------------------------------------------------------------------------
  message("Define gene tissue specificity index")
  #start.time <- Sys.time()
  # Add the weight for the studied tissue
  # tissue selectivity or tissue exclusion from the tissue should be considered
  # as both can be important at a genetic level, repression or activation, of a gene
  # are part of the tissue specificity
  #nreplics <- table(as.character(treatment.f))
  ####  Need to find a way to order factor by order of appearance in expression file
  # e.g nreplics <- c(4,4,4)
  sample.l <- condition
  #m.eset <- as.matrix(m.eset)
  # Rank the tissue selective genes using the Tissue Selective Index
  # for each gene TSI for activation and modify for tissue repression
  m.eset <- as.data.frame(m.eset)
  prets <- ts.IUT(name = paste0(sample.l, netname,"_", method, "_", permutation),
                  eset=m.eset, tissues=sample.names, target=sample.l,
                  threshold = 0.05, filter = 10)
  
  # Vector of tissue-selectivity
  ts <- rbind(prets$filtTSI , prets$filtTSE)
  ts <- data.frame(names = row.names(ts), ts)
  #print(head(ts[order(ts$tsi.order, decreasing=T),])) #######################
  # use data table to select only the best tsi or tse for each genes
  ts <- setDT(ts)[, .SD[which.max(tsi.order)], by=names]
  ts <- as.data.frame(ts)
  rownames(ts) <- ts$names; ts <- ts[,-1] 
  tslist = ts
  #end.time <- Sys.time()
  #time.taken <- end.time - start.time
  #print(time.taken)
  
  #-----------------------------------------------------------------------------------------------
  message("Determine the overall interaction importance for each node gene")
  # for the interactions coming from the same node
  # create hash with all the genes in the network as keys
  coexprs.pvals <- structure(rep(1, length(Genelist)), names=Genelist)   #1 interaction p-values
  #request of the hash package
  h <- hash()
  for (i in 1:length(int.sig[,1])){
    g1 = as.character(int.sig[i,1])
    g2 = as.character(int.sig[i,2])
    twistedcor = 1 - abs(int.sig[i,3]);
    if (twistedcor == 0) {
      twistedcor = max.p.val
    }
    if (all(has.key(g1, h))) {
      h[[g1]] <- append(h[[g1]], twistedcor)
    }else {
      h[[g1]] <- twistedcor
    }
    if (all(has.key(g2, h))) {
      h[[g2]] <- append(h[[g2]], twistedcor)
    }else {
      h[[g2]] <- twistedcor
    }
  }
  # add a vector with all the p-value attached to a gene
  # in the co-expression analysis
  for (i in 1:length(Genelist)){
    gene <- Genelist[i] # get the name of the gene
    trans.cumpvals <- fishersMethod(as.numeric(as.vector(h[[gene]])))
    cumpvals <- -log(trans.cumpvals,10)  # transform pval in factor and Take -log10 of the results
    if (is.na(cumpvals)){
      coexprs.pvals[gene] <- 0;
    }else {
      coexprs.pvals[gene] <- cumpvals
    }
  }
  
  # Differential Expression and ranking for the network genes for the factor")
  DE.ranks <- DE.ranking(m.eset, Genelist, treatment.f, sample.l, sample.names)
  DE.ranks <- as.data.frame(DE.ranks)
  #print(head(DE.ranks))  ####
  #the results are already scaled
  
  # scaling the network parameters 
  coexprs.pvals <- ScalN(coexprs.pvals)  #1 interaction p-values
  ts <- structure(ts$tsi.order, names=rownames(ts))
  tsi.order <- ScalN(ts ) #2 expression specificity
  BetwC <- ScalN(BetwC) #3 betweenness
  Conn <- ScalN(Conn) #4 connectivity
  ClusCoef <- ScalN(ClusCoef)  #5 transitivity
  PageRank <- ScalN(PageRank) #6 PageRank centrality
  
  coexprs.pvals = as.data.frame(coexprs.pvals)
  tsi.order <- as.data.frame(tsi.order)
  BetwC <- as.data.frame(BetwC)
  Conn <- as.data.frame(Conn)
  ClusCoef <- as.data.frame(ClusCoef) 
  EigenC <- as.data.frame(EigenC)
  PageRank <- as.data.frame(PageRank)
  
  ###-------------------------------------------------------------------------------------------------
  message("Merging and scaling the network parameters")
  # Hub NECorr merge all the sub-network statistics
  m.param <- list(coexprs.pvals, tsi.order, BetwC, Conn, ClusCoef, EigenC, PageRank, DE.ranks) %>%  
    map(~ .x %>% as.data.frame %>% rownames_to_column('rn')) %>% 
    reduce(full_join, by = 'rn') %>%
    column_to_rownames('rn')
  m.param[is.na(m.param)] <- 0
  m.param <- as.data.frame(m.param)
  m.param$BetwC <- ifelse(m.param$DE.ranks >= 0.4, m.param$BetwC, 0)
  m.param$Conn <- ifelse(m.param$DE.ranks >= 0.4, m.param$Conn, 0)
  m.param$ClusCoef <- ifelse(m.param$DE.ranks >= 0.4, m.param$ClusCoef, 0)
  #colnames(m.param) <- c("interaction.pvals",
  #                      "tissue.treatment.specificity",
  #                      "Betweenness.Centrality",
  #                      "Connectivity", "Transitivity",
  #                      ,"DE")
  #--------------------------------------------------------------------------------------------------
  # calculate the weight for: (use ahp to generate weight)")
  #done before the loop wpar variable
  
  # Hub NECorr - pick the top genes for validations - best alternatives")
  # calculate the weighted ranking for each gene = alternative ranking
  hub.m.param <- m.param[,c("coexprs.pvals", "tsi.order", "BetwC", "Conn", "ClusCoef")]
  colnames(hub.m.param)  <- c("interaction.pvals", "tissue.treatment.specificity",
                              "Betweenness.Centrality", "Connectivity", "Transitivity")
  for (i in 1:length(colnames(hub.m.param))){
    hub.m.param[,i] <- hub.m.param[,i]*wpar[i]
  }
   #################
  # gene ranking
  gene.rank.h <- rowSums(hub.m.param)*100
  gene.rank.h <- as.data.frame(gene.rank.h[order(gene.rank.h, decreasing=TRUE)])
  colnames(gene.rank.h)[1] <- "gene.rank.h"
  nGenes <- nrow(gene.rank.h)
  #test2 <- left_join(rownames_to_column(gene.rank.h), 
  #          rownames_to_column(hub.m.param),
  #          by = ("rowname" = "rowname"))
  #print("#### test2 ####")
  #print(head(test2[order(test2$Betweenness.Centrality, decreasing = T ),]))
  #################
  #building hub rank hash
  #load the ranks into a hash to make faster with the hash package
  gene.rank.hash <- hash()
  geneIDs <- row.names(gene.rank.h)
  for (i in 1:nGenes) {
    geneRank <- as.numeric(gene.rank.h[i,1])
    gene.rank.hash[[geneIDs[i]]] <- geneRank
  }
  gene.rank.h.description <- left_join(rownames_to_column(gene.rank.h), 
                                       rownames_to_column(Desc),
                                       by = ("rowname" = "rowname"))
  gene.rank.h.description <- as.data.frame(gene.rank.h.description)
  colnames(gene.rank.h)[1] <- sample.l
  total.rank <- as.data.frame(total.rank)
  total.rank.h <- left_join(rownames_to_column(total.rank),
                            rownames_to_column(gene.rank.h),
                            by = ("rowname" = "rowname"))
  total.rank.h <- as.data.frame(total.rank.h)
  
  #end.time <- Sys.time()
  #time.taken <- end.time - start.time
  #print(time.taken)
  #------------------------------------------------------------------------------------------------
  message("Effector NECorr - pick the top genes for validations - best alternatives")
  #start.time <- Sys.time()
  # calculate the probability = alternative ranking
  eff.m.param <-  m.param[ ,c("coexprs.pvals", "tsi.order", "BetwC", 
                              "Conn", "EigenC", "PageRank", "DE.ranks")] 
  
  colnames(eff.m.param) <- c("interaction.pvals", "tissue.treatment.specificity",
                             "Betweenness.Centrality", "Connectivity", "Eigenvector",
                             "PageRank","DE")
  
  ### MACHINE LEARNING with parameter define from validations
  ### Predict probability to be important gene
  
  #!!!! if Naive Bayes
  #suppressWarnings(suppressPackageStartupMessages(library(klaR)))
  #j.nB <- NaiveBayes(phenotype ~ . , data = eff.m.param, kernel = "rectangular", n = 148)
  load(system.file("extdata", "ML_model.RData", package="NECorr", mustWork = TRUE))
  j.nB <- model
  #
  gene.rank.e.description <- effector_significance(eff.m.param=eff.m.param, Desc=Desc, j.nB=j.nB, 
                                  sample.l=sample.l)
  message("effector calculation done")
  gene.rank.eff <- gene.rank.e.description$gene.rank.eff
  names(gene.rank.eff) <- rownames(gene.rank.e.description)
  
  
  #end.time <- Sys.time()
  #time.taken <- end.time - start.time
  #print(time.taken)
  #------------------------------------------------------------------------------------------------
  message("Hub interaction ranking")
 
  # Hub interactions (edges) ranking
  hub.int.ranks <- hub_edge_significance(network.int=network.int, 
                                         gene.rank.hash=gene.rank.hash)
  hub.int.significant <- subset(hub.int.ranks, p2 < 0.05)
  #------------------------------------------------------------------------------------------------
  # Cleaning results for display
  MCDAnode <- gene.rank.h.description[order(gene.rank.h.description$gene.rank.h, decreasing = T),] 
  MLconfirmed <- as.vector(
    gene.rank.e.description$rowname[which(gene.rank.e.description$gene.rank.eff > 0.6)])
  
  MCDAnode_ML <- gene.rank.h.description[gene.rank.h.description$rowname %in% MLconfirmed, ]
  MCDAnode_ML  <- MCDAnode_ML[order(MCDAnode_ML$gene.rank.h, decreasing = T),] 
  
  MCDAedge <- hub.int.significant[unique(which(hub.int.significant$targetIDs %in% MLconfirmed),
                                         which(hub.int.significant$sourceIDs %in% MLconfirmed)), ]
  
  #MCDAedge.significant <- subset(MCDAedge, p2 < 0.005)
  MCDAedge <- MCDAedge[order(MCDAedge$ranks.sum, decreasing = T),]
  
  #print("eff.int.significant")
  prereg.net <- hub.int.significant[ ,c(1,2)]
  
  # Activator NECorr - based on genes in the complete network edges linking the hub genes")
  actres <- activator_significant(hub.int.significant = hub.int.significant, 
                        network.int=network.int, Desc = Desc)
  
  #print("actres")
  #print(head(actres))
  #gene.rank.act.significant <- actres$rank
  #gene.rank.act.description <-  actres$description
  #end.time <- Sys.time()
  #time.taken <- end.time - start.time
  #print(time.taken)
   
  #hub.net <- as.data.frame(
  #  cbind(as.character(hub.int.significant$sourceIDs),
  #        as.character(hub.int.significant$targetIDs),
  #        as.numeric(as.character(hub.int.significant$ranks.sum)),
  #        rep("hub", nrow(hub.int.significant))))
  #colnames(hub.net) <- c("source","target","score","node.type")
  #print("hub.net")
  #print(head(hub.net))
  # significant activator sub-network linked to hub network
  #act.net <- linked_act_hub_net(hub.int.significant, gene.rank.act.significant,
  #                              network.int)
  
  # significant effector sub-network linked to hub network
  #eff.net <- linked_eff_hub_net(hub.int.significant, eff.int.significant)
  
  # merge important networks
  #hub.act.net <- rbind(act.net,hub.net)
  #print("hub.act.net")
  #print(head(hub.act.net))
  
  res <- list(
    necorr_hub_nodesRank = MCDAnode, # hub gene ranking using MCDA linear model
    necorr_hub_edgesRank = MCDAedge, # edge including hubs ranked using MCDA linear model
    necorr_reg = actres, # regulator based on PageRank from extended network around hub genes
    necorr_hub_nodesML = MCDAnode_ML ## hub gene ranking present in ML
    #netstat = netstat, # network statistic
    #coexpres = int.sig, # co-expression ranking
    #deg_rank = DE.ranks, # Differentially expressed ranking
    #NBeff_rank = gene.rank.e.description, # Naives Bayes effector gene ranking
    #edge_rank = coexprs.pvals, # interaction importance
    #tsi_rank = tslist # tissue specificity
    #hub_edge_rank = hub.int.ranks, # hub interaction gene ranking
    #gene.rank.act.significant = gene.rank.act.significant,
    #activator_rank = gene.rank.act.description, # activator gene ranking
    #significant_network = hub.act.net,
    )
  return(res)
}
 
  
#' necorr_graph
#'
#' @param hubnet significant network continuing hub, activators
#' @param hub.int.significant hub interaction significance
#' @param gene.rank.act.significant  activator ranking
#' 
#' @return netgraph
#' @export
necorr_graph <- function(hubnet, hub.int.significant, 
                         gene.rank.act.significant){
  # link activator and hub significant sub-networks
  sig.hub <- unique(c(as.character(hub.int.significant$V1), 
                     as.character(hub.int.significant$V2)))
  if(( nrow(hubnet)>0) == TRUE){
    #print(hub.act.net)  
    g <- graph.data.frame(hubnet, directed = T)
    
    vcolors <- rep("cyan",length(V(g)$name))
    vcolors[which(V(g)$name %in% sig.hub)] <- "red"
    
    vsize.hub <- as.numeric(
      gene.rank.h[V(g)$name[which(V(g)$name %in% sig.hub)],1])
    vsize.hub <- 2^(ScalN(vsize.hub) + 2.21)
    
    temp <- as.data.frame(gene.rank.act.significant)
    rownames(temp) <- temp[,1]
    vsize.act  <- as.numeric(temp[V(g)$name[which(V(g)$name %in% temp[,1])],2])
    vsize.act <- 2^(ScalN(vsize.act) + 2.21)
    vsize <- c(vsize.act,vsize.hub)
    vlabel <- as.character(Desc[V(g)$name,"Associated.Gene.Name"])
    if (is.finite(vsize) & vsize>0){
      # 	    mark.groups <- vcolors
      # 	    mark.col <- visColoralpha(vcolors, alpha=0.2)
      # 	    mark.border <- visColoralpha(vcolors, alpha=0.2)
      #
      # 	    mark.groups= mark.groups,mark.col= mark.col, mark.border=mark.border,
      # pdf(file = title, width = 10, height = 10)
      netgraph <- dnet::visNet(g, glayout=layout.fruchterman.reingold(g) ,
                               vertex.shape="sphere",
                               vertex.label = vlabel, edge.color = "grey",
                               edge.arrow.size = 0.3, vertex.color = vcolors,
                               vertex.frame.color = vcolors, newpage = F)
      return(netgraph)
    }
  }
  else{
    message("no network can be drawn")
    }
  }
