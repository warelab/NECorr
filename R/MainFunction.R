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
#' @param name the name of the
#' @param Filelist condition list see if still necessary with metadata
#' @param type Omics comparative expression type: protein or gene
#' @param permutation permutation number used for all significance calculation
#' #param lmiR List of miRNAs
#' @param method used for co-expression correlation: GCC, MINE, PCC, SCC or KCC
#' @param sigcorr significance of the correlation
#' @param NSockets number of sockets
#' @param fadjacency correlation with all combination (all) or network
#' combination only (only)
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
#' @export
Necorr <- function(networkFile = "",
                   expression = "",
                   description.file = "",
                   condition = "",
                   metadata = "", name = "",
                   Filelist = "",
                   method = "GCC", permutation = 1000, sigcorr = 0.01,
                   fadjacency = "only", type = "gene",
                   NSockets = 2){
  
  # load file and options
  # read networkFile
  network.int <- read.delim(networkFile, sep ="\t", header = T, fileEncoding="latin1")
  # take only in consideration the real interactions
  if(ncol(network.int) == 3){
    network.int <- network.int[,c(1,3)]
  } else if(ncol(network.int) == 5){
    network.int <- network.int[, c(1,3)]
  } else if(ncol(network.int) == 2){
    network.int <- network.int[, c(1,2)]
  }
  # Condition experiment and/or tissue
  factortab <- read.table(metadata, header = T) #factor.file #condition <- "Radial"
  # Description file name with gene name and annotations
  Desc <-  read.csv(description.file, header=T, row.name=1) #description.file <- "1.Ath.GeneDesc.csv"
  # Create the subdirectory for the fianl results
  # need to see if all these tables are still useful??
  #___________________ MAIN SCRIPT ___________________________________
  
  # calculate the weight for each parameter
  # (use same order for columns in gene stats table)
  # suppressWarnings(suppressPackageStartupMessages(require(pmr)))
  # preftable <-  matrix(data = rep(1,(5*5)), nrow = 5, ncol = 5, byrow = TRUE)
  
  # preftable[1,2] = 1; preftable[2,1] = 1 # interaction pval and expression
  # preftable[1,3] = 2; preftable[3,1] = 1/2 # interaction pval and betweeness
  # preftable[1,4] = 4; preftable[4,1] = 1/4 # interaction pval and connectivity
  # preftable[1,5] = 3; preftable[5,1] = 1/3 # interaction pval and transitivy
  
  # preftable[2,3] = 2; preftable[3,2] = 1/2 # expression and betweeness
  # preftable[2,4] = 4; preftable[4,2] = 1/4 # expression and connectivity
  # preftable[2,5] = 3; preftable[5,2] = 1/3 # expression and transitivy
  
  # preftable[3,4] = 3; preftable[4,3] = 1/3 # Betweeness and connectivity
  # preftable[3,5] = 3; preftable[5,3] = 1/3 # Betweeness and transitivy
  
  # preftable[4,5] = 1/2; preftable[5,4] = 2 # connectivity and transitivy
  
  
  # # weight of parameter(wpar) give the criteria ranking
  # wpar <- ahp(preftable)$weighting
  # #print (wpar)
  wpar <- c(0.31744701,0.31744701,0.19772678,0.06747995,0.09989925)
  
  message("Start of R analysis for the Main script part")
  # The next section is dealing with the fact
  # that some sample types can have no replicate
  # so if this is the case, pseudo-replicates will be generated
  # to be able to calculate a statistic for
  # the tissue/stress-selectivity of each sample type
  # and order these tissue/stress-selective genes
  sample.names <- unique(factortab[,1])
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
  #---
  message("I - Loading of the network statistics and the gene list")
  #---
  #start.time <- Sys.time()
  # Generate all the topology statistics using the function NetTopology
  netname <- basename(networkFile)
  netstat <- NetTopology(network.int)
  Genelist <- rownames(netstat)
  BetwC <- structure(netstat$BetweennessCentrality, names=Genelist)  #3 betweeness
  Conn <- structure(netstat$EdgeCount, names=Genelist) # 4 connectivity
  ClusCoef <- structure(netstat$ClusteringCoefficient, names=Genelist) #5 transitivity
  EigenC <- structure(netstat$eigenvector_centrality, names=Genelist) #6 eigenvector
  PageRank <- structure(netstat$pagerank, names=Genelist) #7 pagerank
  
  # Generate an empty vector to serve as table with all the ranks for conditions
  total.rank <- matrix(data = NA, nrow = length(Genelist))
  rownames(total.rank) <- Genelist
  total.rank[,1] = Genelist
  # R network topology done 

  # Read expression file
  eset <- read.table(expression,header = T,row.names=1) #!! change made Here
  # Take only the genes that are part of the molecular network
  eset <- eset[Genelist[Genelist %in% rownames(eset)],] #!! change made Here
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
  
  #---
  message("II - Analysis of the co-expression file and p-value sums")
  #start.time <- Sys.time()
  #---
  int.sig <- multiCorr(m.eset, net = network.int, nsockets = NSockets,
                       methods = method, output = "paired", sigmethod = "two.sided",
                       pernum = permutation , verbose = T, cpus = NSockets)
  #int.sig.file <- int.sig
  #colnames(int.sig.file) <- c("Source","Target","Correlation","p-value")
  # This adds a column of color values based on the y values
  # remove the line with NA due to not found gene due to no expresssion
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
  #---
  message("III - Define gene tissue specificity index (TSI) (Yanai 2011, Bioinformatics)")
  #start.time <- Sys.time()
  #---
  
  # Add the weight for the studied tissue
  # tissue selectivity or tissue exclusion from the tissue should be considered
  # as both can be important at a genetic level, repression or activation, of a gene
  # are part of the tissue specificity
  #nreplics <- table(as.character(treatment.f))
  ####  Need to find a way to order factor by order of appearance in expression file
  # e.g nreplics <- c(4,4,4)
  sample.l <- condition
  m.eset <- as.matrix(m.eset)
  # Rank the tissue selective genes using the Tissue Selective Index
  # for each gene TSI for activation and modify for tissue repression
  m.eset <- as.data.frame(m.eset)
  prets <- ts.IUT(name = paste0(sample.l, netname,"_", method, "_", permutation),
                  eset=m.eset, tissues=sample.names, target=sample.l,
                  threshold = 0.05, filter = 20)
  
  # Vector of tissue-selectivity
  ts <- rbind(prets$filtTSI , prets$filtTSE)
  ts <- data.frame(names = row.names(ts), ts)
  # use data table to select only the best tsi or tse for each genes
  ts <- setDT(ts)[, .SD[which.max(tsi.order)], by=names]
  ts <- as.data.frame(ts)
  rownames(ts) <- ts$names; ts <- ts[,-1] 
  #end.time <- Sys.time()
  #time.taken <- end.time - start.time
  #print(time.taken)
  
  #-----------------------------------------------------------------------------------------------
  message("IV - Determine the overall interaction importance for each node (gene)")
  #-----------------------------------------------------------------------------------------------
  
  #start.time <- Sys.time()
  # for the interactions coming from the same node
  # create hash with all the gene in the network as keys
  int.pvals <- structure(rep(1, length(Genelist)), names=Genelist)   #1 interaction p-values
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
  #start.time <- Sys.time()
  for (i in 1:length(Genelist)){
    gene <- Genelist[i] # get the name of the gene
    trans.cumpvals <- fishersMethod(as.numeric(as.vector(h[[gene]])))
    cumpvals <- -log(trans.cumpvals,10)  # transform pval in factor and Take -log10 of the results
    if (is.na(cumpvals)){
      int.pvals[gene] <- 0;
    }else {
      int.pvals[gene] <- cumpvals
    }
  }
  
  #---------------------------------------------------------------------------------------
  message("V - Differential Expression and ranking for the network genes for the factor")
  #start.time <- Sys.time()
  #----------------------------------------------------------------------------------------
  DE.ranks <- DE.ranking(m.eset, Genelist, treatment.f, sample.l, sample.names)
  DE.ranks <- as.data.frame(DE.ranks)
  #the results are already scaled
  
  #end.time <- Sys.time()
  #time.taken <- end.time - start.time
  #print(time.taken)
  ###-----------------------------------------------------------------------------------------
  message("VI - Normalize each column between 0 and 1 using the some")
  #start.time <- Sys.time()
  ###----------------------------------------------------------------------------------------
  int.pvals <- ScalN(int.pvals)  #1 interaction p-values
  ts <- ScalN(ts) #2 expression specificity
  BetwC <- ScalN(BetwC) #3 betweeness
  Conn <- ScalN(Conn) #4 connectivity
  ClusCoef <- ScalN(ClusCoef)  #5 transitivity
  PageRank <- ScalN(PageRank)
  
  int.pvals = as.data.frame(int.pvals)  #1 interaction p-values
  ts <- as.data.frame(ts); #2 expression specificity
  BetwC <- as.data.frame(BetwC) #3 betweeness
  Conn <- as.data.frame(Conn) #4 connectivity
  ClusCoef <- as.data.frame(ClusCoef)  #5 transitivity
  EigenC <- as.data.frame(EigenC)
  PageRank <- as.data.frame(PageRank)
  
  #end.time <- Sys.time()
  #time.taken <- end.time - start.time
  #print(time.taken)
  ###-------------------------------------------------------------------------------------------------
  message("VII - Hub NECorr merge all the subnetwork statistics")
  #start.time <- Sys.time()
  ###-------------------------------------------------------------------------------------------------
  m.param <- list(int.pvals, ts, BetwC, Conn, ClusCoef, EigenC, PageRank, DE.ranks) %>% 
    map(~ .x %>% as.data.frame %>% rownames_to_column('rn')) %>% 
    reduce(full_join, by = 'rn') %>%
    column_to_rownames('rn')
  print(m.param)
  m.param[is.na(m.param)] <- 0
  #par.tab <- m.param
  m.param <- as.data.frame(m.param)
  colnames(m.param) <- c("interaction.pvals",
                        "tissue.treatment.specificity",
                        "Betweenness.Centrality",
                        "Connectivity","Transitivity",
                        "Eigenvector","PageRank","DE")
  #tail(m.param)
  #end.time <- Sys.time()
  #time.taken <- end.time - start.time
  #print(time.taken)
  #--------------------------------------------------------------------------------------------------
  message("VIII - calculate the weight for: (use ahp to generate weight)")
  #start.time <- Sys.time()
  #--------------------------------------------------------------------------------------------------
  
  #done before the loop wpar variable
  #end.time <- Sys.time()
  #time.taken <- end.time - start.time
  #print(time.taken)
  
  #--------------------------------------------------------------------------------------------------
  message("IX - Hub NECorr - pick the top genes for validations - best alternatives")
  #--------------------------------------------------------------------------------------------------
  
  # calculate the weighted ranking for each gene = alternative ranking
  hub.m.param <-  m.param[,c( "interaction.pvals",
                              "tissue.treatment.specificity",
                              "Betweenness.Centrality",
                              "Connectivity","Transitivity")]
  for (i in 1:length(colnames(hub.m.param))){
    hub.m.param[,i] <- hub.m.param[,i]*wpar[i]
  }
  
  # gene ranking
  gene.rank.h <- apply(hub.m.param,1,sum)*100
  gene.rank.h <- gene.rank.h[order(gene.rank.h, decreasing=TRUE)]
  nGenes <- length(gene.rank.h)
  #print("building hub rank hash")
  gene.rank.h <- as.data.frame(gene.rank.h)
  #load the ranks into a hash to make step XI faster
  #suppressWarnings(suppressPackageStartupMessages(require(hash)))
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
  #print(gene.rank.h)
  total.rank <- as.data.frame(total.rank)
  total.rank.h <- left_join(rownames_to_column(total.rank),
                            rownames_to_column(gene.rank.h),
                            by = ("rowname" = "rowname"))
  total.rank.h <- as.data.frame(total.rank.h)
  
  #end.time <- Sys.time()
  #time.taken <- end.time - start.time
  #print(time.taken)
  #------------------------------------------------------------------------------------------------
  message("X - Effector NECorr - pick the top genes for validations - best alternatives")
  #start.time <- Sys.time()
  #------------------------------------------------------------------------------------------------
  # calculate the probability = alternative ranking
  eff.m.param <-  m.param[ ,c("interaction.pvals",
                              "tissue.treatment.specificity",
                              "Betweenness.Centrality",
                              "Connectivity",
                              "Eigenvector",
                              "PageRank","DE")]
  
  ### MACHINE LEARNING with parameter define from validations
  ### Predict probability to be important gene
  
  #!!!! if Naive Bayes
  #suppressWarnings(suppressPackageStartupMessages(library(klaR)))
  #j.nB <- NaiveBayes(phenotype ~ . , data = eff.m.param, kernel = "rectangular", n = 148)
  load(system.file("extdata", "ML_model.RData", package="NECorr", mustWork = TRUE))
  j.nB <- model
  #data("./R/model.rds")
  #j.nB <- readRDS("./R/model.rds")
  prob.pred <- predict(j.nB, type="prob", newdata = eff.m.param )
  
  gene.rank.eff <- (prob.pred$posterior)[,2]
  gene.rank.eff = gene.rank.eff[order(gene.rank.eff, decreasing=TRUE)]
  gene.rank.eff = as.data.frame(gene.rank.eff)
  # merge description and ranking
  gene.rank.e.description <- left_join(rownames_to_column(gene.rank.eff),
                                       rownames_to_column(Desc),
                                       by = ("rowname" = "rowname"))
  gene.rank.e.description <- as.data.frame(gene.rank.e.description)
  colnames(gene.rank.eff)[1] <- sample.l
  total.rank.eff <- left_join(rownames_to_column(total.rank),
                              rownames_to_column(gene.rank.eff),
                              by = ("rowname" = "rowname"))
  total.rank.eff <- as.data.frame(total.rank.eff)
  #
  gene.rank.eff.hash <- hash()
  #print("building eff rank hash")
  geneIDs <- row.names(gene.rank.eff)
  for (i in 1:nGenes){
    geneRank <- as.numeric(gene.rank.eff[i,1])
    gene.rank.eff.hash[[geneIDs[i]]] <- geneRank
  }
  
  #end.time <- Sys.time()
  #time.taken <- end.time - start.time
  #print(time.taken)
  #------------------------------------------------------------------------------------------------
  message("XI - Hub interaction ranking")
  #start.time <- Sys.time()
  #------------------------------------------------------------------------------------------------
  # hub sub-network
  # cl2 <- makeCluster(nsockets, type="SOCK")
  # registerDoSNOW(cl2)
  # hub.int.ranks <- foreach(i=1:nrow(network.int),.combine='rbind' )%dopar%{
  sourceIDs <- as.vector(network.int[,1])
  targetIDs <- as.vector(network.int[,2])
  ranks.sum <- rep(1,length(targetIDs))
  for(i in 1:nrow(network.int)){
    ranks.sum[i] <- sum(gene.rank.hash[[sourceIDs[i]]], gene.rank.hash[[targetIDs[i]]])
  }
  hub.int.ranks <- data.frame(sourceIDs,targetIDs,ranks.sum)
  # stopCluster(cl2)
  break.points <- c(-Inf, unique(sort(as.numeric(hub.int.ranks[,3]))), Inf)
  p2 <- cut( as.numeric(hub.int.ranks[,3]), breaks=break.points, labels=FALSE )  ####### !!!!!
  p2 <- 1 - p2/length(break.points)
  hub.int.ranks <- as.data.frame(cbind(hub.int.ranks,p2))
  hub.int.ranks$p2 <- as.numeric(as.character(hub.int.ranks$p2))
  hub.int.significant <- subset(hub.int.ranks, p2 < 0.005)
  # ranking of the effector per edges
  eff_ranks.sum <- rep(1,length(targetIDs))
  for(i in 1:nrow(network.int)) {
    SourceID <- as.character(network.int[i,1])
    TargetID <- as.character(network.int[i,2])
    eff_ranks.sum <- sum(gene.rank.eff.hash[[sourceIDs[i]]], gene.rank.eff.hash[[targetIDs[i]]])
  }
  eff.int.ranks <- data.frame(sourceIDs, targetIDs, eff_ranks.sum)
  # stopCluster(cl3)
  
  break.points <- c(-Inf, unique(sort(as.numeric(eff.int.ranks[,3]))), Inf)
  p2 <- cut( as.numeric(eff.int.ranks[,3]), breaks=break.points, labels=FALSE )
  p2 <- 1 - p2/length(break.points)
  eff.int.ranks <- as.data.frame(cbind(eff.int.ranks,p2))
  eff.int.ranks$p2 <- as.numeric(as.character(eff.int.ranks$p2))
  eff.int.significant = subset(eff.int.ranks, p2 < 0.005)
  
  #end.time <- Sys.time()
  #time.taken <- end.time - start.time
  #print(time.taken)
  ###----------------------------------------------------------------------------------------------
  message("XII - Activator NECorr - based on genes in the complete network edges linking the hub genes")
  #start.time <- Sys.time()
  ###----------------------------------------------------------------------------------------------
  # find genes that are significant in the hub subnetwork in the complete network
  sig.hub <- unique(c(as.character(hub.int.significant$V1), as.character(hub.int.significant$V2)))
  sc.sig.hub <- subset(network.int, network.int[,1] %in% sig.hub)
  tg.sig.hub <- subset(network.int, network.int[,2] %in% sig.hub)
  net.extension.sig.hub <- rbind(sc.sig.hub,tg.sig.hub)
  
  ## change the names of the hub gene in the extented hub network using source genes
  change.sig.hubNames <- net.extension.sig.hub[,1] %in% sig.hub
  tmp<-as.character(net.extension.sig.hub[,1])
  tmp[change.sig.hubNames]<-"sig"
  net.extension.sig.hub[,1] <- tmp
  
  ## change the names of the hub gene in the extented hub network using target genes
  change.sig.hubNames <- net.extension.sig.hub[,2] %in% sig.hub
  tmp <- as.character(net.extension.sig.hub[,2])
  tmp[change.sig.hubNames] <- "sig"
  net.extension.sig.hub[,2] <- tmp
  
  ## network extended to putative regulator using source genes
  act.m.param <- net.extension.sig.hub
  sc.count <- rle( sort( act.m.param[,1] ))
  act.m.param$Count <- sc.count[ match( act.m.param[,1] , sc.count ) ]
  
  # gene ranking of linked to hub nodes
  gene.rank.act <- cbind(sc.count$values, sc.count$lengths)
  gene.rank.act <- gene.rank.act[order(as.numeric(gene.rank.act[,2]), decreasing=TRUE),]
  gene.rank.act <- gene.rank.act[ - which(gene.rank.act[,1] == "sig"),]
  
  # add the genes that have 25% of the genes linked to hub genes
  gene.rank.act.significant <- gene.rank.act[which(gene.rank.act[,2] >= (length(sig.hub)*0.20)),]
  
  # write the putative activator genes
  gene.rank.act.b <- as.data.frame(gene.rank.act)
  gene.rank.act.description <- left_join(rownames_to_column(gene.rank.act.b), 
                                         rownames_to_column(Desc), 
                                         by = ("rowname" = "rowname"))
  gene.rank.act.description <- as.data.frame(gene.rank.act.description)
  
  #end.time <- Sys.time()
  #time.taken <- end.time - start.time
  #print(time.taken)
  ###-------------------------------------------------------------------------------------------------
  message("XIII  - Vizualization of network and list of interactions that could affected and color")
  ###-------------------------------------------------------------------------------------------------
  #start.time <- Sys.time()
  #suppressWarnings(suppressPackageStartupMessages(library(igraph)))
  #suppressWarnings(suppressPackageStartupMessages(library(dnet)))
  hub.net <- as.data.frame(
    cbind(as.character(hub.int.significant$V1),
          as.character(hub.int.significant$V2),
          as.numeric(as.character(hub.int.significant$V3)),
          rep("hub", nrow(hub.int.significant))))
  colnames(hub.net) <- c("source","target","score","node.type")
  # significant activator sub-netowrk linked to hub network
  net.extension.sig.hub <- rbind(sc.sig.hub,tg.sig.hub)
  act.net.1 <- subset(
    net.extension.sig.hub,
    net.extension.sig.hub[,1] %in% as.vector(as.character(gene.rank.act.significant[,1])))
  act.net.2 <- subset(
    net.extension.sig.hub,
    net.extension.sig.hub[,2] %in% as.vector(as.character(gene.rank.act.significant[,1])))
  act.net.pre <- rbind(act.net.1,act.net.2)
  meanSig <- mean(as.vector(as.numeric(as.character(hub.int.significant$V3))))
  act.net <- cbind(act.net.pre,
                   rep(meanSig,nrow(act.net.pre)),
                   rep("act",nrow(act.net.pre)))
  rm(act.net.1,act.net.2,act.net.pre,meanSig)
  colnames(act.net) <- c("source","target","score","node.type")
  # significant effector sub-netowrk linked to hub network
  sc.sig.eff <- subset(eff.int.significant, eff.int.significant$V1 %in% sig.hub)
  tg.sig.eff <- subset(eff.int.significant, eff.int.significant$V2 %in% sig.hub)
  eff.net.pre <- rbind(sc.sig.eff,tg.sig.eff)
  eff.net <- cbind(eff.net.pre[,1],eff.net.pre[,2],eff.net.pre[,3], rep(nrow(eff.net.pre)))
  rm(sc.sig.eff,tg.sig.eff,eff.net.pre)
  
  # link activator and hub significant sub-networks
  if(!is.null(act.net) & !is.null(hub.net)){
    hub.act.net <- rbind(act.net,hub.net)
  }else if(!is.null(hub.net)){
    hub.act.net <- hub.net
  }
  if(( nrow(hub.act.net)>0) == TRUE){
    #print(hub.act.net)  
    g <- graph.data.frame(hub.act.net, directed = T)
    
    vcolors <- rep("cyan",length(V(g)$name))
    vcolors[which(V(g)$name %in% sig.hub)] <- "red"
    
    vsize.hub <- as.numeric(gene.rank.h[V(g)$name[which(V(g)$name %in% sig.hub)],1])
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
      #dev.off()
    }else{
      g <- graph.data.frame(hub.act.net, directed = T)
      vcolors <- rep("cyan",length(V(g)$name))
      temp <- as.data.frame(gene.rank.act.significant)
      rownames(temp) <- temp[,1]
      vlabel <- as.character(Desc[V(g)$name,"Associated.Gene.Name"])
      if (is.finite(vsize) & vsize>0){
        # 	    mark.groups <- vcolors
        # 	    mark.col <- visColoralpha(vcolors, alpha=0.2)
        # 	    mark.border <- visColoralpha(vcolors, alpha=0.2)
        #
        # 	    mark.groups= mark.groups,mark.col= mark.col, mark.border=mark.border,
        # pdf(file = title, width = 10, height = 10)
        netgraph <- dnet::visNet(g, glayout=layout.fruchterman.reingold(g) ,
                     vertex.shape="sphere", vertex.size = vsize,
                     vertex.label = vlabel, edge.color = "grey",
                     edge.arrow.size = 0.3, vertex.color = vcolors,
                     vertex.frame.color = vcolors, newpage = F)
        #dev.off()
        
      }
    }
    res <- list(
      netstat = netstat, # network statistic
      coexpres = int.sig, # co-expresion ranking
      deg_rank = DE.ranks, # Diffeentially expressed ranking
      edge_rank <- int.pvals,# interaction importance
      tsi_rank = tslist, # tissue specificity
      hub_rank = gene.rank.h.description, # hub gene ranking
      activator_rank = gene.rank.act.description, # activator gene ranking
      effector_rank = gene.rank.e.description, # effector gene ranking
      hub_edge_rank = hub.int.ranks # interaction gene rankin
      # graph preparation
    )
    }
    return(res)
  }
 
  
