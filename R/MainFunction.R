#' Necorr
#' @author Christophe Liseron-Monfils, Andrew Olson
#' @param expression Expression file in log2 (ratio expression) with row: gene,
#' first column: type of sample,second column: sample names
#' @param networkFile Molecular network file with source in the first column, 
#' targets in the second column
#' @param pretopology optional argument FALSE by default if not
#' a precalculation of the network topology is expected as R RDA files
#' with the same name than the network file and the same localization
#' @param description.file genome description
#' @param condition Condition from expression to study the network co-expression
#' correlation
#' @param metadata dataframe with the metadata
#' @param permutation permutation number used for all significance calculation
#' #param lmiR List of miRNAs
#' @param sigcorr significance of the correlation
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
#' partly based on rsgcc package for the GCC, PCC,KCC and SPP 
#' Ma et al, 2012, plant Physiology
#' @return res
#' @export
Necorr <- function(networkFile = "",
                   expression = "",
                   description.file = "",
                   condition = "",
                   metadata = "", 
                   permutation = 1000, 
                   sigcorr = 0.01,
                   pretopology = FALSE,
                   verbose= 0){
  verbose = 0
  #options(warn = warn_verbose)
  # load file and options
  # Read network File
  network.int <- read.delim(networkFile, sep ="\t", header = T, 
                            fileEncoding="latin1")
  network.int <- network.int[!is.na(network.int[1, ]), ]
  network.int <- network.int[!is.na(network.int$target), ]
  # Take only in consideration the real interactions
  if (ncol(network.int) == 2) {
    network.int <- network.int[,c(1, 2)]
  } else if (ncol(network.int) > 2) {
    network.int <- network.int[, c(1, 3)]
  } 
  # Remove NA from the network data
  network.int <- network.int[!is.na(network.int[, 1]), ]
  network.int <- network.int[!is.na(network.int[, 2]), ]
  # Get the list of genes from the initial network
  Genelist <- unique(c(as.character(network.int[, 1]),
                     as.character(network.int[, 2])))
  # Condition experiment and/or tissue
  factortab <- read.table(metadata, header = T, sep = "\t", row.names = 1)
  # Description file name with gene name and annotations
  Desc <-  read.csv(description.file, header = T, row.name = 1) 
  # Read expression file
  eset <- read.table(expression, header = T, row.names = 1) 
  eset[sapply(eset, is.infinite)] <- 10^300
  #_________________________ MAIN SCRIPT _______________________________________
  # The next section is dealing with the fact
  # that some sample types can have no replicate
  # so if this is the case, pseudo-replicates will be generated
  # to be able to calculate a statistic for
  # the tissue/stress-selectivity of each sample type
  # and order these tissue/stress-selective genes
  message("Progressing on the transcriptomic parameters")
  sample.names <- unique(factortab[,1])
  xcol <- c()
  treatment.f <- factor()
  for (i in 1:length(sample.names)) {
    nrep <- length(which(factortab[, 1] == sample.names[i]))
    if (i == 1) {
      j <- i
    } 
    if (nrep == 1) {
      # Create 3 artificial columns if the sample has only one replicate
      xcol <- c(xcol, rep(j, 3))
      j <- j + 3
      treatment.f <- factor(c(as.character(treatment.f),
                              as.character(rep(sample.names[i], 3))))
    } else if (nrep > 1) {
      xcol <- c(xcol, seq(j, j + (nrep - 1)))
      treatment.f <- factor(c(as.character(treatment.f),
                              as.character(rep(sample.names[i], nrep))))
      j <- j + (nrep)
    }
  }
  treatment.f <- factor(treatment.f, levels = sample.names)
  #-----------------------------------------------------------------------------
  # Take only the genes that are part of the molecular network
  eset <- eset[intersect(Genelist, rownames(eset)), ]
  m.eset <- as.matrix(eset)
  m.eset <- m.eset[, xcol]
  # Loop to measure the importance of gene expression
  conditionList <- factortab[, 2]
  factorList <- factortab[, 1]
  df <- cbind(as.character(conditionList), as.character(factorList))
  df <- df[!duplicated(df), ]
  #-----------------------------------------------------------------------------
  #message("Analysis of the co-expression file and p-value sums")
  int.sig <- multiCorr(m.eset, net = network.int,
                       pernum = permutation,
                       sigmethod = 0.6,
                       verbose = TRUE)
  # This adds a column of color values based on the y values
  # remove the line with NA due to not found gene due to no express/ion
  # or error of annotation in the initial network file
  int.sig <- int.sig[complete.cases(int.sig), ]
  # Transform the p-value in the maximal p-value in function
  # of number of permutations to avoid log-scale = infinity
  # Replace the p-value equal to 0 by the maximal p-value e.g.
  # 1/10000 knowing that we have 10000 randomization in initial calculations
  max.p.val <- 1/as.numeric(permutation)
  int.sig[which(int.sig[,4] == 0), 4] <- max.p.val
  int.sig <-as.data.frame(int.sig)
  #-----------------------------------------------------------------------------
  #message("Define gene tissue specificity index")
  # Add the weight for the studied tissue
  # tissue selectivity or tissue exclusion from the tissue should be considered
  # as both can be important at a genetic level, 
  # repression or activation, of a gene are part of the tissue specificity
  #nreplics <- table(as.character(treatment.f))
  #  Need to order factor by order of appearance in the expression file
  # e.g nreplics <- c(4,4,4)
  sample.l <- condition
  #m.eset <- as.matrix(m.eset)
  # Rank the tissue selective genes using the Tissue Selective Index
  # for each gene TSI for activation and modify for tissue repression
  m.eset <- as.data.frame(m.eset)
  prets <- ts.IUT(name = paste0(sample.l, netname, "_", 
                                "GCC_", permutation),
                  eset = m.eset, tissues = sample.names, target = sample.l,
                  threshold = 0.05, filter = 10)
  # Vector of tissue-selectivity
  ts <- rbind(prets$filtTSI , prets$filtTSE)
  ts <- data.frame(names = row.names(ts), ts)
  ts <- ts[order(ts$tsi.order, decreasing = T), ]
  #saveRDS(ts, "ts.rds") #######################
  # use data table to select only the best tsi or tse for each genes
  ts <- setDT(ts)[, .SD[which.max(tsi.order)], by=names]
  ts <- as.data.frame(ts)
  #print("ts")
  #print(head(ts))
  rownames(ts) <- ts$names; ts <- ts[,-1] 
  tslist = ts
  #-----------------------------------------------------------------------------
  #message("Determine the overall interaction importance for each node gene")
  # for the interactions coming from the same node
  # create hash with all the genes in the network as keys
  # 1 interaction p-values
  coexprs.pvals <- structure(rep(1, length(Genelist)), names=Genelist)   
  coexpressionres <- paste0(expression, "co-express.rda")
  if (file.exists(coexpressionres)) {
    coexprs.pvals <- readRDS(coexpressionres)
  } else {
    # request of the hash package
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
      # transform pval in factor and Take -log10 of the results
      cumpvals <- -log(trans.cumpvals,10)  
      if (is.na(cumpvals)){
        coexprs.pvals[gene] <- 0;
      } else {
        coexprs.pvals[gene] <- cumpvals
      }
    }
    saveRDS(coexprs.pvals, file = coexpressionres)
  }
  
  #-----------------------------------------------------------------------------
  #message("Differential expression filtering")
  # Differential Expression and ranking for the network genes for the factor")
  # the results are already scaled
  DE.ranks <- DE.ranking(m.eset, Genelist, treatment.f, sample.l, sample.names)
  DE.ranks <- as.data.frame(DE.ranks, stringsAsFactors = FALSE)
  #-----------------------------------------------------------------------------
  # scaling the network parameters 
  #1 interaction p-values
  coexprs.pvals <- scaling_param(coexprs.pvals)  
  #2 expression specificity
  tsi.order <- structure(ts$tsi.order, names=rownames(ts))
  tsi.order <- scaling_param(tsi.order)
  #-----------------------------------------------------------------------------
  # Generate all the topology statistics using the function NetTopology
  message("Progressing on the molecular network parameters")
  netname <- basename(networkFile)
  expressionname <- basename(expression)
  de.genes <- rownames(subset(DE.ranks, DE.ranks >= 0.5))
  coexp.genes <- rownames(subset(coexprs.pvals, coexprs.pvals >= 0.5))
  sel_for_net <- unique(c(de.genes, coexp.genes))
  # Filter the network for genes expressed in the experiment 
  # or with significant co-expressions
  network.int.filt <- network.int
  colnames(network.int.filt) <- c("source", "target")
  network.int.filt <- subset(
    network.int.filt,
    source %in% sel_for_net | target %in% sel_for_net)
  net_in_mem <- paste0(networkFile, "-", expressionname, "-", condition, "_netstat.rda")
  if (file.exists(net_in_mem)) {
    netstat <- readRDS(net_in_mem)
  } else {
    netstat <- NetTopology(network.int.filt)
    saveRDS(netstat, file = net_in_mem)
  }
  # 3 betweenness
  BetwC <- structure(netstat$BetweennessCentrality, names = rownames(netstat))  
  BetwC <- scaling_param(BetwC)
  # 4 connectivity
  Conn <- structure(netstat$EdgeCount, names = rownames(netstat)) 
  Conn <- scaling_param(Conn)
  # 5 transitivity
  ClusCoef <- structure(netstat$ClusteringCoefficient, names = rownames(netstat)) 
  ClusCoef <- scaling_param(ClusCoef)
  # 6 eigen-vector
  EigenC <- structure(netstat$eigenvector_centrality, names = rownames(netstat)) 
  EigenC <- scaling_param(EigenC)
  # 7 page rank
  PageRank <- structure(netstat$pagerank, names = rownames(netstat)) 
  PageRank <- scaling_param(PageRank)
  # Generate an empty vector to serve as table with all the ranks for conditions
  total.rank <- matrix(data = NA, nrow = length(Genelist))
  rownames(total.rank) <- Genelist
  total.rank[, 1] <- Genelist
  #-----------------------------------------------------------------------------
  message("NECorrHub: Merging and scaling the network parameters")
  # Hub NECorr merge all the sub-network statistics
  #print(head(coexprs.pvals))
  #print(head(BetwC))
  m.param <- left_join(rownames_to_column(coexprs.pvals), 
                       rownames_to_column(tsi.order), by=c("rowname"))
  m.param <- purrr::reduce(list(
    m.param,
    rownames_to_column(BetwC),
    rownames_to_column(Conn),
    rownames_to_column(ClusCoef),
    rownames_to_column(EigenC),
    rownames_to_column(PageRank),
    rownames_to_column(DE.ranks)), dplyr::left_join, by = 'rowname') %>% column_to_rownames('rowname')
  m.param[is.na(m.param)] <- 0
  m.param <- as.data.frame(m.param)
  #-----------------------------------------------------------------------------
  # message("Gene ranking")
  # Calculate the weight for: (use ahp to generate weight)") done before 
  # the loop wpar variable
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
  #print(wpar)
  wpar <- c(0.2825910, 0.2825910, 0.1864405, 0.1412955, 0.1070820)
  # Hub NECorr - pick the top genes for validations - best alternatives")
  # calculate the weighted ranking for each gene = alternative ranking
  if (nrow(network.int) < 2000 ) {
    # small network
    wpar <- c(0.2825910, 0.2825910, 0.1864405, 0.1412955, 0.1070820)
    hub.m.param <- m.param[, c("coexprs.pvals", "tsi.order",
                               "BetwC", "Conn", "ClusCoef")]
    colnames(hub.m.param) <- c("interaction.pvals", 
                                "tissue.treatment.specificity",
                                "Betweenness.Centrality",
                                "Connectivity", "Transitivity")
  } else {
    # big network option
    wpar <- c(0.310830903472379, 0.0705565032549202, 0.450522631406784, 0.168089961865917)
    hub.m.param <- m.param[, c("tsi.order", "BetwC", "PageRank", "DE.ranks")]
    colnames(hub.m.param) <- c("tissue.treatment.specificity", "Betweenness.Centrality", "PageRank", "DE")
  }
  
  #print(str(hub.m.param))
  #print(data.frame(min=sapply(hub.m.param, min),max=sapply(hub.m.param, max)))
  hub.m.param <- t(t(hub.m.param) * wpar)
  #################
  # gene ranking
  gene.rank.h <- rowSums(hub.m.param)
  gene.rank.h <- as.data.frame(gene.rank.h[order(gene.rank.h, decreasing=TRUE)])
  colnames(gene.rank.h)[1] <- "gene.rank.h"
  nGenes <- nrow(gene.rank.h)
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
  #-----------------------------------------------------------------------------
  #message("NECorrHub_ML")
  #start.time <- Sys.time()
  # calculate the probability = alternative ranking
  eff.m.param <-  m.param[ ,c("coexprs.pvals", "tsi.order", "BetwC", 
                              "Conn", "EigenC", "PageRank", "DE.ranks")] 
  
  colnames(eff.m.param) <- c("interaction.pvals",
                             "tissue.treatment.specificity",
                             "Betweenness.Centrality", "Connectivity",
                             "Eigenvector",
                             "PageRank", "DE")
  
  ### MACHINE LEARNING with parameter define from validations
  ### Predict probability to be important gene
  
  #!!!! if Naive Bayes
  #suppressWarnings(suppressPackageStartupMessages(library(klaR)))
  #j.nB <- NaiveBayes(phenotype ~ . , data = eff.m.param, kernel = "rectangular", n = 148)
  load(system.file("extdata", "ML_model.RData",
                   package="NECorr", mustWork = TRUE))
  j.nB <- model
  #
  gene.rank.e.description <- effector_significance(eff.m.param = eff.m.param, 
                                                   Desc = Desc, j.nB =j.nB,
                                                   sample.l = sample.l)
  gene.rank.eff <- gene.rank.e.description$gene.rank.eff
  names(gene.rank.eff) <- rownames(gene.rank.e.description)
  #-----------------------------------------------------------------------------
  message("NeCorrEdge: interaction ranking")
  # Hub interactions (edges) ranking
  hub.int.ranks <- hub_edge_significance(network.int = network.int, 
                                         gene.rank.hash=gene.rank.hash)
  hub.int.significant <- subset(hub.int.ranks, p2 < 0.05)
  #-----------------------------------------------------------------------------
  # Cleaning results for display
  MCDAnode <- gene.rank.h.description[order(gene.rank.h.description$gene.rank.h,
                                            decreasing = T), ] 
  rownames(MCDAnode) <- MCDAnode$rowname
  MCDAnode <- MCDAnode[, -1]
  # Get table for the hub genes back up by Machine Learning
  MLconfirmed <- as.vector(
    gene.rank.e.description$rowname[
      which(gene.rank.e.description$gene.rank.eff > 0.6)])
  # 
  MCDAnode_ML <- gene.rank.h.description[
    gene.rank.h.description$rowname %in% MLconfirmed, ]
  MCDAnode_ML  <- MCDAnode_ML[order(MCDAnode_ML$gene.rank.h, decreasing = T), ] 
  rownames(MCDAnode_ML) <- MCDAnode_ML$rowname
  MCDAnode_ML <- MCDAnode_ML[, -1]
  #
  MCDAedge <- hub.int.significant
  MCDAedge_ML <- hub.int.significant[
    unique(which(hub.int.significant$targetIDs %in% MLconfirmed),
           which(hub.int.significant$sourceIDs %in% MLconfirmed)), ]
  MCDAedge <- MCDAedge[order(MCDAedge$ranks.sum, decreasing = T),]
  prereg.net <- hub.int.significant[ ,c(1,2)]
  # Activator NECorr - based on genes in the complete network 
  # edges linking the hub genes")
  actres <- activator_significant(hub.int.significant = hub.int.significant, 
                        network.int = network.int, Desc = Desc)
  colnames(MCDAnode) <- c("Gene_rank_h", "Associated_Gene_Name",
                           "Chromosome_Name", "Strand",	
                           "Gene_Biotype",	"TF_family")
  colnames(actres) <- c("PageRank", "EdgeCount",	"Degree_Out",
                       "Associated_Gene_Name","Chromosome_Name", "Strand",
                       "Gene_Biotype",	"TF_family")
  colnames(MCDAnode_ML) <- c("Gene_rank_h", "Associated_Gene_Name",
                              "Chromosome_Name", "Strand",	
                              "Gene_Biotype",	"TF_family")
  res <- list(
    # hub gene ranking using MCDA linear model
    necorrHub_nodes = MCDAnode, 
    # regulator based on PageRank from extended network around hub genes
    necorrReg = actres, 
    # edge including hubs ranked using MCDA linear model
    necorrEdges = MCDAedge,
    ## hub gene ranking present in ML
    necorrHub_nodesML = MCDAnode_ML,
    # edge including hubs ranked using MCDA linear model in ML
    necorrEdgesML = MCDAedge_ML 
    #netstat = netstat, # network statistic
    #coexpres = int.sig, # co-expression ranking
    #NBeff_rank = gene.rank.e.description, # Naives Bayes effector gene ranking
    #edge_rank = coexprs.pvals, # interaction importance
    #tsi_rank = tslist, # tissue specificity
    #hub_edge_rank = hub.int.ranks, # hub interaction gene ranking
    #hub.m.param = m.param
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
