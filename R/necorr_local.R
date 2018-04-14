# core functions necessary to apply the NECorr pipeline
# variable $1 : name of the file with variables
# variable $2 : network file -n (Ath_Y1H-RR-interactions.txt)
# variable $3 : condition file -c (Colmutant)
# variable $4 : expression file -e (Mutant-ratios.txt)
# variable $5 : miR list file -l (Ath_miRNA_TAIR_Seq.txt to transform to Ath_affinity_gene.txt)
# variable $6 : expression type -t (gene or protein); 
# variable $7 : permutation number -p for ACscore and correlations
# variable $8 : method for correlation -m (GCC, MINE, PCC, SPP, KCC)
# variable $9 : type of combination all (possible within gene of network) or only (network) -f
# variable $10: significance of p-value correlation obtained thanks to the permutation  -s (0.01)

# variable $11: Description file  "Ath.gene.desc.csv"
# variable $12: Expression Factor file "ColdFactor.txt"

# variable $13: file names of co-expression
# variable $14: nb
# variable $15: script

Necorr <- function(network.file, description.file, factor.file, 
                   metadata, name, 
                   permutation = 1000,
                   CoorMetric= 'GCC'
                   Filelist #condition list see if still necessary with metadata
                   ){
  
  # AT THE END ADD THE OTHER FUNCTIONS FROM SOURCE
  #source("./src/NECorr.functions.v3.R")
  #load("./src/ML_model.RData")
  suppressWarnings(suppressPackageStartupMessages(require(RColorBrewer)))
  suppressWarnings(suppressPackageStartupMessages(require(gplots)))
  library(Rcpp)
  sourceCpp("gini.cpp")
  double_me3(5)
  # load file and options
  # read network file
  network.int <- read.delim(network.file, sep ="\t", header = T,fileEncoding="latin1")
  #
  factortab <- read.table(metadata,header = T) #factor.file
  # Description file name with gene name and annotations
  Desc <-  read.csv(description.file,header=T, row.name=1) #description.file <- "1.Ath.GeneDesc.csv"
  # Condition experiment and/or tissue
  factortab <-  read.table(factor.file,header = T) #condition <- "Radial"
  # Create the subdirectory for the fianl results
  # need to see if all these tables are still useful??
  subDirTS = paste0("results/",condition,"/4_TS_file/")
  subDirGraph = paste0("results/",condition,"/6_TS_graph/")
  subDirFile = paste0("results/",condition,"/7_gene_ranking_per_condition/")
  dir.create(file.path(mainDir, subDirTS), showWarnings = FALSE)
  dir.create(file.path(mainDir, subDirGraph), showWarnings = FALSE)
  dir.create(file.path(mainDir, subDirFile), showWarnings = FALSE)
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
  
  print("Start of R analysis for the Main script part")
  # The next section is dealing with the fact 
  # that some sample types can have no replicate 
  # so if this is the case, pseudo-replicates will be generated 
  # to be able to calculate a statistic for 
  # the tissue/stress-selectivity of each sample type
  # and order these tissue/stress-selective genes
  
  sample.names <- unique(factortab$Treatment)
  coexpression.cond <- unique(factortab$Condition)
  xcol <-c()
  treatment.f <- factor()
  for (i in 1:length(sample.names)){ 
    nrep <- length(which(factortab$Treatment== sample.names[i]))
    if(nrep == 1){
      xcol <- c(xcol,rep(i,3))
      treatment.f <- factor(c(as.character(treatment.f),as.character(rep(sample.names[i],3))))
    } else if(nrep > 1){
      xcol <- c(xcol,seq(i,i+(nrep-1)))
      treatment.f <- factor(c(as.character(treatment.f),as.character(rep(sample.names[i],nrep))))
    }
    #print(treatment.f)
  }
  #print(xcol)
  treatment.f <- factor(treatment.f, levels = sample.names)
  ###------------------------------------------------------------------------------------------------------------------------
  print("### I - Load the network statistics and the gene list")
  ###------------------------------------------------------------------------------------------------------------------------
  #start.time <- Sys.time()
  # Generate all the topology statistics using the function NetTopology
  netstat <- NetTopology(network.int)
  write.csv(netstat,paste0(mainDir,"/results/",condition,"/5_",netname,"_NetStat.csv"))
  Genelist <- rownames(netstat)
  BetwC <- structure(netstat$BetweennessCentrality, names=Genelist)  #3 betweeness
  Conn <- structure(netstat$EdgeCount, names=Genelist) # 4 connectivity
  ClusCoef <- structure(netstat$ClusteringCoefficient, names=Genelist) #5 transitivity
  EigenC <- structure(netstat$eigenvector_centrality, names=Genelist) #6 eigenvector
  PageRank <- structure(netstat$pagerank, names=Genelist) #7 pagerank
  
  # Generate an empty vector to serve as table with all the ranks for conditions 
  total.rank = matrix(data = NA, nrow = length(Genelist))
  rownames(total.rank) <- Genelist
  total.rank[,1] = Genelist
  print("R network topology done")
 
  # Read expression file
  eset <- read.csv(expression,header = T) #!! change made Here
  # Take only the genes that are part of the molecular network
  eset <- eset[Genelist,] #!! change made Here
  m.eset <- as.matrix(eset)  
  m.eset <- m.eset[-grep("NA",rownames(m.eset)),xcol]
  # Loop to measure the importance of gene expression 
  conditionList <-factortab$Condition
  factorList <- factortab$Treatment
  df <- cbind(as.character(conditionList),as.character(factorList))
  df <- df[!duplicated(df),]
  
  
  suppressWarnings(suppressPackageStartupMessages(require(foreach)))
  suppressWarnings(suppressPackageStartupMessages(require(doSNOW)))
  files.list <- paste0(mainDir,"/",Filelist)
  print(files.list)
  files <- read.table(files.list)
  condnb <- nrow(files)
  
  ## Loop to do the anlysis per condition defined in the first row of the expression data
  #CondNB <-1  ### test of the code without loop
  for (CondNB in 1:condnb ){
    fileexp <- as.character(files[CondNB,1])
    x.exp <- as.matrix(read.table(paste0(mainDir,"/",fileexp)))
    
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(time.taken)
    
    ###------------------------------------------------------------------------------------------------------------------------
    print("### II - Analysis of the co-expression file and p-value sums")
    #start.time <- Sys.time()
    ###------------------------------------------------------------------------------------------------------------------------
    int.sig <- bigcorGCC(x.exp, net = network.int, nsockets = NSockets,
                         methods= CoorMetric, output = "paired", sigmethod = "two.sided", 
                         pernum = permutation , verbose = FALSE, cpus = NSockets, type = incrtype)
    
  }
}
#################### TO merge with preceding function
necorr <- function(expsion , network, condition,
                   lmiR = "1.Ath.affinity.gene.txt",
                   method = "GCC", permutation = 1000, sigcorr = 0.01,
                   fadjacency = "only",
                   gdesc,
                   type = "gene",
                   dirtmp="./results/tmp", dirout = './results'){
#' @author Christophe Liseron-Monfils
#' @param expsion Expression file in log2 (ratio expression) with row: gene,
#' first column: type of sample,second column: sample names
#' @param network Molecular network file with source in the first column, targets in
#'  the second column
#' @param condition Condition from expression to study the network co-expression
#' correlation
#' @param type Omics comparative expression type: protein or gene
#' @param permutation permutation number used for all significance calculation
#' @param lmiR List of miRNAs
#' @param method Method used for co-expression correlation: GCC, MINE, PCC,
#' SCC or KCC
#' @param gdesc genome description
#' @param dirtmp directory for the temporary results
#' @param dirout directory for the results
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

# Create a random number used for the rest of the analysis
# this will distingish analysis done with the same dataset.
nb <- RANDOM

# Create the output directory if not existing; generate "./results" dir and
# "./results/tmp"
if(!dir.exists(dirout)){
  dir.create(dirout)
}
if(!dir.exists(dirtmp)){
  dir.create(dirtmp, recursive = T)
}

# Generate the expression factor file ( see ./src/a.1.factorfile.v2.pl)
factorstab <- factorfile(exps = expsion, nb = nb,  condition = condition, dirout = dirtmp)

stdout <- mainNecorr( )

}
########################
factorfile <- function(exps, nb, condition , dirout = "./results/tmp/"){
#' @author Christophe Liseron-Monfils
#' @param exps Expression file in log2(ratio expression) with row:gene,
#' first column: type of sample, second colum: sample names
#' @param condition Different sample present in the expression data: stress
#' condition, tissue types and/or developmental stages
#' @param nb
#' @param dirout results directory
#' @description Generate the factor file for the co-expression analysis

# Keep only the core filename by erasing the prefix dir and the suffix
# (file extension)
name <- extractname(exps)
out <- paste0(dirout, "/", condition, ".", name, ".", nb, ".Factor.txt")
tab <- read.delim2(exps)
Conditions <- as.character(tab[1, ])
Factors <- as.character(tab[2, ])
Samples <- as.character(tab[3, ])

# Print the vector file
tabFact <- do.call("rbind", list(Samples, Factors, Conditions))
rownames(tabFact) <- tabFact$Samples
tabFact <- tabFact[, -1]
return(tabFact)

}

extractname <- function(filestring){
  #' @author Christophe Liseron-Monfils
  #' @description Extract the core name of a file by erasing
  #' the directory and the file extension. This allows to use the file
  #' in another the naming of other files.
  #' @param filestring directory and filename (../tmp/file.txt)
  pfile <- strsplit(filestring, "[.]")
  namefile <- strsplit(pfile[[1]][length(pfile[[1]]) -1], "[/]")
  namefile <- namefile[[1]][length(namefile)]
  return(namefile)
}

# C.Liseron-Monfils - Ware lab Sept2013
# CSHL
# partly based on rsgcc package for the GCC, PCC,KCC and SPP Ma et al, 2012, plant Physiology
# 

cor.matrix.NECorr <- function (GEMatrix, cpus = cpus, 
                               cormethod = c("GCC", "PCC", "SCC", "KCC", "BiWt"), 
                               style = c("pairs.between","pairs.only"), 
                               var1.id = NA, var2.id = NA, 
                               pernum = 0, 
                               sigmethod = c("two.sided","one.sided"), 
                               output = c("matrix","paired")) {
  if (cpus > 1) {
    #   library(snowfall)
  }
  if (pernum == 0) {
    sigmethod <- "two.sided"
  }
  if (!is.matrix(GEMatrix) | !is.numeric(GEMatrix)) {
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
  if( style =="pairs.only"){
    taskmatrix <- cbind(var1.id, var2.id)
  }
  # if (style == "all.pairs") {
  # var1.id <- seq(1, dim(GEMatrix)[1], by = 1)
  # var2.id <- var1.id
  # }
  #if(style == "pairs.between" || style == "all.pairs") {
  if(style == "pairs.between") {
    df <- expand.grid.unique(var1.id, var2.id, include.equals = TRUE)
    df2 = t(apply(df, 1, sort))
    taskmatrix <- df2[!duplicated(df2),]
  }
  #if (cpus == 1 | cormethod == "BiWt") {
  #  suppressWarnings(suppressPackageStartupMessages(require(biwt)))
  pernum = as.numeric(pernum)
  results <- apply(taskmatrix, 1, cor.pair, GEMatrix = GEMatrix,
                   rowORcol = "row", cormethod = cormethod, pernum = pernum, 
                   sigmethod = sigmethod)
  #}else(){
  #sfnit(parallel = TRUE, cpus = cpus)
  #print(sprintf("%s cpus to be used", sfCpus()))
  #     results <- sfApply(taskmatrix, 1, cor.pair, GEMatrix = GEMatrix, 
  #                        rowORcol = "row", cormethod = cormethod, pernum = pernum, 
  #                        sigmethod = sigmethod)
  #     sfStop()
  #   }
  if (output == "paired") {
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
  } else {
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
    return(list(corMatrix = corMatrix, pvalueMatrix = pvalueMatrix))
  }
}

bigcorGCC <- function(x ,net= NA, nsockets= 4, methods = c("GCC","PCC","SCC","KCC"),
                      sigmethod = c("two.sided", "one.sided"),
                      nblocks = 10, verbose = TRUE, cpus = 1, pernum = 0, ...){
  suppressWarnings(suppressPackageStartupMessages(require(foreach)))
  suppressWarnings(suppressPackageStartupMessages(require(doSNOW)))
  nsockets <- as.numeric(nsockets)
  cl <- makeCluster(nsockets, type="SOCK")
  registerDoSNOW(cl)
  corMAT <- c()
  Nrow = nrow(net)
  fc <- gl(nblocks, ceiling(Nrow/nblocks), length = Nrow)
  matnet <- split(net,fc)
  #corMAT <-foreach(i=1:nblocks, .combine='rbind') %do% {
  corMAT<-foreach(j=1:nblocks, .combine='rbind', .export=c('indexing.network','cor.matrix.NECorr','cor.pair','gcc.corfinal'))%dopar%{
    netindexed <- indexing.network(as.matrix(x) ,matnet[[j]])
    G1 <- as.numeric(as.vector(netindexed[,1]))
    G2 <- as.numeric(as.vector(netindexed[,2]))
    cor.matrix.NECorr(x, var1.id=G1, var2.id=G2, sigmethod=sigmethod, cormethod=methods, 
                      pernum=pernum, output=output, cpus=cpus, style="pairs.only")
    }
  stopCluster(cl)
  return(corMAT)
  
}


expand.grid.unique <- function(x, y, include.equals=FALSE){
  x <- unique(x)
  y <- unique(y)
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}

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



#differential expression ranking uisng the network gene
DE.ranking <- function(exps,GeneList,factortab,sample.l, exps.file = FALSE){
  suppressWarnings(suppressPackageStartupMessages(require(Biobase)))
  suppressWarnings(suppressPackageStartupMessages(require(limma)))
  # 1-expression file with NECorr format
  # 2-list of netowrk genes
  # 3-expression factor table generate automatically
  # 4-name of the studied sample
  # read expressiontable from microarrays or RNAseq
  if (exps.file == TRUE){
    df<- as.matrix(read.table(exps, header=TRUE,sep="\t",row.names = 1))
    expression <- df[,colSums(is.na(df)) != nrow(df)]
    targets <- read.table(factortab,header=TRUE,sep="\t",row.names=1)
    sample.names <- unique(targets$Treatment)
    f <- factor(targets$Treatment, levels=sample.names)
  }
  else { 
    expression <- as.matrix(exps)
    f <- factor(factortab)
  }
  # parse expression data to only have the gene that are in the network - load the genes that are in the networks
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
  equation <- noquote(paste0(noquote(sample.l)," - ((", paste(other.names, collapse = '+'),") /", length(sample.names)-1,")"))
  assign(".equation", equation, envir=.GlobalEnv)
  #contrast.matrix <- makeContrasts(equation , levels=design)
  contrast.matrix <- makeContrasts(.equation, levels=design)
  remove(".equation", envir=.GlobalEnv)
  fit2 <- contrasts.fit(fit, contrast.matrix) # Computes estimated coefficients and standard errors for a given set of contrasts.
  fit2 <- eBayes(fit2)
  #write.table(topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=50000), file="limma_complete.xls", row.names=F, sep="\t")
  
  # Exports complete limma statistics table for first comparison group ('coef=1') to tab delimited text file.
  de <- topTable(fit2, coef=1, adjust="fdr", sort.by="P", number=length(sel.Genes))
  ranks <- ScalN(de$B)
  names(ranks) = rownames(de)
  return(ranks)
}

# network topology statistics using igraph methods
NetTopology <-function(network)
{
  #the network should be a list of interactions space by a line
  suppressWarnings(suppressPackageStartupMessages(require(igraph)))
  # script partly adapted from https://gist.github.com/mhawksey/1682306
  # author: Martin Hawksey 
  
  # pass to igraph the network that would transform in a graph object
  g <- graph.data.frame(network, directed = T)
  
  ## calculate the topology stats necessary for network statistics
  betweenness_centrality <- betweenness(g,v=V(g),directed = F, normalized = T)
  #betweenness_centrality <- betweenness.estimate(g,vids=V(g),directed = F, normalized = T)
  
  eigenvector_centrality<-evcent(g, scale = TRUE)
  # algorithm used by google to find webpages; this is an evolution of eigen-vector
  # not sure if applicable to bological data yet
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
  cc <- c("gene_name", "EdgeCount","degree_in","degree_out","BetweennessCentrality","eigenvector_centrality","ClusteringCoefficient","pagerank")
  colnames(datagrid) <- cc
  rownames(datagrid) <- gene_name
  #eliminate the first column being a duplicate of rownames now
  datagrid <- datagrid[,-1]
  return(datagrid)
}

# tissue selectivity from Van Deun 2009
IUT<-function(DATA,nrreplics,target,updown,alpha,mc)
{
  nrprobesets<-nrow(DATA)
  nrarrays<-ncol(DATA)
  nrconditions<-length(nrreplics)
  
  #Calculation of condition-specific averages and variances for all probesets
  lower<-1
  AV<-matrix(rep(0,nrprobesets*nrconditions),nrow=nrprobesets)
  SSWITHIN<-rep(0,nrprobesets)
  for (i in (1:nrconditions))
  {
    upper<-lower+nrreplics[i]-1
    v<-rep(1,nrreplics[i])
    DATA_i<-DATA[,lower:upper]
    average<-DATA_i%*%(v%*%((t(v)%*%v)^(-1)))
    AV[,i]<-t(average)
    sswithin<-((DATA_i-average%*%v)^2)%*%v
    SSWITHIN<-SSWITHIN+sswithin
    lower<-upper+1
  }
  #Calculation of t-statistics with s≤ from ANOVA context
  MTARGET<-AV[,target]%*%t(rep(1,nrconditions-1))
  MOTHERS<-AV
  MOTHERS<-MOTHERS[,-target]
  invN_target<-1/nrreplics[target]
  invN_OTHERS<-nrreplics^(-1)
  invN_OTHERS<-invN_OTHERS[-target]
  SE2<-(SSWITHIN/(nrarrays-nrconditions))%*%(invN_target+invN_OTHERS)
  SE<-SE2^(.5)
  TMATRIX<-(MTARGET-MOTHERS)/SE
  
  #Critical t-value using multiple comparison correction for testing if mc=1
  #multiple probes (Sidak ~= Bonferonni)
  if (mc==1)
    alpha<-1-(1-alpha)^(1/nrprobesets)#Sidak
  critT<-qt(1-alpha,nrarrays-nrreplics[target])
  
  #SUMMARY
  SUMSTAT<-rowSums(updown*TMATRIX>critT)
  IUT<-SUMSTAT==nrconditions-1
}

# tissue-selective upregulated genes from Tissue Specificty index (Yanai et al, 2005)
numeratorTSI = function(x) 2**x/max(2**x)
denominator = function(x) length(2**x)-1
TSI = function(x) sum(1-numeratorTSI (x))/denominator(x)

# tissue-selective downregulated genes (Liseron, 2013)
numeratorTSR = function(x) (2**x - 2**min(x))/(2**max(x))
TSR = function(x) sum(numeratorTSR(x))/denominator(x)

# rescaling data in a range between 0 to 1
ScalN = function(x) (x - min(x))/(max(x) - min(x))


# combine p-values Fisher method to have overall importance of 
# this node in the studied tissue (using code from Michael Love code) 
# (http://mikelove.wordpress.com/2012/03/12/combining-p-values-fishers-method-sum-of-p-values-binomial/)
fishersMethod = function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE) 

Aff.table <- function (type, species, file_ddG,file_misM,file_ini) {
  # create 3 table for each of the factor
  ddG <- as.matrix(read.table(file_ddG))
  misM <- as.matrix(read.table(file_misM))
  InibT <- as.matrix(read.table(file_ini))
  
  # use this table to modify the value in the affinity table
  tab <- 100*((((((2^(misM+0.5))*exp(-2))/factorial(misM+0.5))+0.7120482)*.3)+(1/exp((max(ddG)/100)*ddG)*.7))*(InibT/10)
  tab[tab < 0.1] <- 0
  write.table(tab,paste(species,"_affinity_",type,".txt", sep =""), sep = "\t", row.names = TRUE, col.names=NA,quote =FALSE)
}

Aff.table.v3 <- function (type, species, file_ddG,file_misM,file_ini) {
  # create 3 table for each of the factor
  ddG <- as.matrix(read.table(file_ddG))
  misM <- as.matrix(read.table(file_misM))
  InibT <- as.matrix(read.table(file_ini))
  
  # use this table to modify the value in the affinity table
  tab <- 10*((((((2^(misM+0.5))*exp(-2))/factorial(misM+0.5))+0.7120482)*.5) + (1/exp((max(ddG)/100)*ddG)*.5))
  tab[tab ==  min(tab)] <- 0
  write.table(tab,paste(species,"_affinity_",type,".txt", sep =""), sep = "\t", row.names = TRUE, col.names=NA,quote =FALSE)
}


## subdivision of file in subfile to parallelize using qsub
## all combinations
sub.parallel.matrix <- function (file,tmp) {
  tab <- read.table(file)
  mat <-as.matrix(tab)
  d = length(rownames(mat))
  t <- combn(1:d,2)
  tt <- as.data.frame(t(t))
  # number of files increasing with nber of rows and combinations (100 genes to 4000)
  # setting number of files to 2000 for 4000 genes
  # x <- c(9900,15996000); y <- c(10,2000); lm(y~x)
  nb = (4.977e-04)*length(rownames(tt)) + 5.0726459
  nb <- ceiling(nb)
  if (nb >= 100){ nb = 100 }
  #print(nb)
  z <- split (tt,c(rep(1:(nb-1), each = floor(length(rownames(tt))/nb)),
                   rep(nb,length(rownames(tt)) - ((nb-1)*floor(length(rownames(tt))/nb)) )
  ))
  # print each combination in number of data.frames from z
  for (i in 1:nb) {
    #write.table(z[[i]], paste("tmp/sub_tab",i, sep =""), quote = FALSE)
    write.table(z[[i]], paste(tmp,"/",file,".",i, sep =""), quote = FALSE, row.names = FALSE,col.names=FALSE)
  }
  return(nb)     
}

## subdivision of file in subfile to parallelize using qsub
## network combinations only
sub.parallel.network <- function (file,netfile,tmp) {
  tab <- read.table(file)
  network <- read.table(netfile)
  network <- network[-1,]
  # create table with all the gene in network 
  t <- as.data.frame(unique(apply(network,1, function(x) c(x[1],x[3]))))
  tt1 <- t(t)
  #replace the gene name by their index in expression table
  #which(rownames(tab) == tt[2,1])
  indexing <- function (x,y) {
    if ( length(which(rownames(y) == x[1])) != 0 && 
         length(which(rownames(y) == x[2])) != 0) {
      c(which(rownames(y) == x[1]),which(rownames(y) == x[2]))
    }
  }
  tt <- as.data.frame(do.call(rbind,apply(tt1, 1, indexing,tab)))
  nb = (4.977e-04)*length(rownames(tt)) + 5.0726459
  nb <- ceiling(nb)
  #print(nb)
  z <- split (tt,c(rep(1:(nb-1), each = floor(length(rownames(tt))/nb)),
                   rep(nb,length(rownames(tt)) - ((nb-1)*floor(length(rownames(tt))/nb)) )
  ))
  # print each combination in number of data.frames from z
  for (i in 1:nb) {
    #write.table(z[[i]], paste("tmp/sub_tab",i, sep =""), quote = FALSE)
    write.table(z[[i]], paste(tmp,"/",file,".",i, sep =""), quote = FALSE, row.names = FALSE,col.names=FALSE)
  }
  return(nb)     
}

## local run of the network under the form of interactionlist
sub.parallel.network.local <- function (file,netfile,tmp) {
  tab <- read.table(file)
  network <- read.table(netfile)
  network <- network[-1,]
  # create table with all the gene in network 
  t <- as.data.frame(unique(apply(network,1, function(x) c(x[1],x[3]))))
  tt1 <- t(t)
  #replace the gene name by their index in expression table
  #which(rownames(tab) == tt[2,1])
  indexing <- function (x,y) {
    if ( length(which(rownames(y) == x[1])) != 0 && 
         length(which(rownames(y) == x[2])) != 0) {
      c(which(rownames(y) == x[1]),which(rownames(y) == x[2]))
    }
  }
  ns <- apply(tt1, 1, indexing,tab)
  if(typeof(ns) == "integer"){ 
    ns = as.data.frame(ns)
  }
  tt <- as.data.frame(do.call(rbind,ns))
  nb = 6
  #print(nb)
  z <- split (tt,c(rep(1:(nb-1), each = floor(length(rownames(tt))/nb)),
                   rep(nb,length(rownames(tt)) - ((nb-1)*floor(length(rownames(tt))/nb)) )
  ))
  # print each combination in number of data.frames from z
  for (i in 1:nb) {
    #write.table(z[[i]], paste("tmp/sub_tab",i, sep =""), quote = FALSE)
    write.table(z[[i]], paste("./",tmp,"/",file,".",i, sep =""), quote = FALSE, row.names = FALSE,col.names=FALSE)
  }
  return(nb)     
}

## calculation of correlation for matrix creating adjacency table
NetCor <- function(file,method = c("MINE","GCC", "PCC", "SCC", "KCC"),nbperm,cpu) {
  tab <- read.table(file)
  tab <-as.matrix(tab)
  if (method == "MINE") {
    x <- RMINE(tab,nbperm)
    write.table(x$corMatrix,paste(file,"_CorrM.txt", sep =""), sep = "\t", row.names = TRUE, col.names=NA,quote =FALSE)
    write.table(x$pvalueMatrix, paste(file,"_PvalM.txt", sep=""), sep = "\t", row.names =TRUE, col.names=NA,quote =FALSE)
  } 
  else {
    x <- Rscorr(tab,method = method,nbperm,1)
    write.table(x$corMatrix,paste(file,"_CorrM.txt", sep =""), sep = "\t", row.names =TRUE, col.names=NA,quote =FALSE)
    write.table(x$pvalueMatrix, paste(file,"_PvalM.txt", sep=""), sep = "\t", row.names =TRUE, col.names=NA,quote =FALSE)
  }
}

## calculation of correlation for sub table of pairs obtained with sub.parallel.matrix
NetCor.pair <- function(file,method = c("MINE","GCC", "PCC", "SCC", "KCC"),nbperm,subfile) {
  tab <- read.table(file,comment.char = "")   #add comment.char = "" to fasten the process
  mat <-as.matrix(tab)
  zi <- read.table(subfile)
  res <- c()
  if (method == "MINE") {
    suppressWarnings(suppressPackageStartupMessages(require(minerva)))
    for (i in 1:length(rownames(zi))) {
      # MIC score calculation
      cor.value <- mine(mat[zi[i,1],],mat[zi[i,2],], alpha=0.6)$MIC
      # p-value calculation
      rowx <- rownames(mat)[zi[i,1]]
      rowy <- rownames(mat)[zi[i,2]]
      x <- as.numeric(mat[rowx,])
      y <- as.numeric(mat[rowy,])
      perm.values = numeric(nbperm)
      for (k in 1:nbperm) {
        # loop for permutations tests
        rx <- sample(x, length(x),replace=FALSE)
        ry <- sample(y, length(y),replace=FALSE)        
        # record the value of the test
        perm.values[k] = mine(rx,ry)$MIC 
      }   
      perm.values = abs(perm.values)
      vals = c(cor.value,perm.values)
      rank.vals = sort(vals, decreasing = TRUE)
      pval <- mean(which(rank.vals == cor.value))/length(vals)
      content <- c(rowx,rowy,cor.value,pval)
      res <- rbind(res,content)
    }
    write.table(res,paste(subfile,"_CorrM_PvalM", sep =""), sep = "\t", row.names = FALSE, col.names=FALSE,quote =FALSE)
  } 
  else {
    for (i in 1:length(rownames(zi))) {
      content <- corr.gcc(mat,cormethod=method,sigmethod="one.sided", style="one.pair", var1.id = zi[i,1], var2.id = zi[i,2], output="paired", pernum = nbperm)
      content <- as.vector(content)
      res <- rbind(res,content)
    }
    #print(res)
    write.table(res,paste(subfile,"_CorrM_PvalM", sep =""), sep = "\t", row.names =FALSE, col.names=FALSE,quote =FALSE)
  }
}


## correlation rsgcc
Rscorr <- function (mat, method = c("GCC", "PCC", "SCC", "KCC"), nbperm, cpu) {
  x <- corr.gcc( mat, cormethod = method, cpus = 1, style = "all.pairs", sigmethod = "one.sided",pernum = nbperm, output = "matrix")
  return(x)
}

## probalistic correlation
RMINE <- function (mat,nbperm) {
  suppressWarnings(suppressPackageStartupMessages(require(minerva)))
  #nbperm= 10000
  # extract all row with a variance inferior to 1e-05 as mine not working with variance equal to 0
  tmat <- t(mat)
  tmat  <- tmat[apply(tmat,1,var) > 1e-05,apply(tmat,2,var) > 1e-05]
  m <- mine(tmat, alpha = 0.7)$MIC
  rownames(m) <- rownames(mat)
  colnames(m) <- rownames(mat)
  # Create an empty table for p-value with same dimension than m with value = 1
  pm <- matrix(1/nbperm, ncol = length(colnames(m)), nrow = length(rownames(m)))
  rownames(pm) <- rownames(mat)
  colnames(pm) <- rownames(mat)    
  # Combine all vector row in table
  z <- as.data.frame(combn(rownames(mat), 2))
  for ( i in 1:length(colnames(z))) {
    rowx <- as.character(z[1,i])
    rowy <- as.character(z[2,i])
    #print(paste(i,z[1,i],z[2,i],sep =" "))
    x <- as.numeric(mat[rowx,])
    y <- as.numeric(mat[rowy,])
    perm.values = numeric(nbperm)
    for (k in 1:nbperm) {
      # loop for permutations tests
      rx <- sample(x, length(x),replace=FALSE)
      ry <- sample(y, length(y),replace=FALSE)        
      # record the value of the test
      perm.values[k] = mine(cbind(rx,ry))$MIC[1,2] 
    }
    
    real.test <- m[rowx,rowy] 
    #print(paste(i,real.test,sep =" "))   
    perm.values = abs(perm.values)
    vals = c(real.test,perm.values)
    rank.vals = sort(vals, decreasing = TRUE)
    pval <- mean(which(rank.vals == real.test))/length(perm.values) 
    pm[rowx,rowy] <- pval
    pm[rowy,rowx] <- pval
    #if (pval < 0.05){
    #    hist(perm.values,col = "lightblue",border ="lightgrey", breaks = 1000, xlim=c(0,1))
    #    curve(dnorm(x, mean=mean(perm.values), sd=sd(perm.values)), add=TRUE, col="darkblue", lwd=2)
    #    abline(v = real.test, col = "red", lwd = 2)
    #}
  }
  
  return(list(corMatrix = m, pvalueMatrix = pm))
}

#affinity table generation

corr.gcc <- function (GEMatrix, cpus = 1, 
                      cormethod = c("GCC", "PCC", "SCC", "KCC", "BiWt"), 
                      style = c("all.pairs", "pairs.between", 
                                "adjacent.pairs", "one.pair"), 
                      var1.id = NA, var2.id = NA, pernum = 0, 
                      sigmethod = c("two.sided","one.sided"), output = c("matrix", "paired")) 
{
  if (cpus > 1) {
    suppressWarnings(suppressPackageStartupMessages(require(snowfall)))
  }
  if (length(cormethod) > 1) {
    cat("Warning: one correlation method should be specified. Default:GCC")
  }
  if (length(style) > 1) {
    cat("Warning: one style should be specified. Default: all.pairs")
  }
  if (pernum > 0 & length(sigmethod) > 1) {
    sigmethod = "two.sided"
  }
  if (pernum == 0) {
    sigmethod <- "two.sided"
  }
  if (!is.matrix(GEMatrix)) {
    stop("Error: GEMatrix in cor.matrix function is not matrix")
  }
  if (!is.numeric(GEMatrix)) {
    stop("Error: GEMatrix is not numeric")
  }
  if (length(rownames(GEMatrix)) == 0) {
    rownames(GEMatrix) <- seq(1, dim(GEMatrix)[1], by = 1)
  }
  VariableNum <- nrow(GEMatrix)
  SampleSize <- ncol(GEMatrix)
  if (VariableNum <= 1 || SampleSize <= 1) {
    stop("Error:the number of variable is less than 2, or the number of observation is less than 2 ")
  }
  if (style == "one.pair") {
    if (length(var1.id) != 1 || length(var2.id) != 1 || is.na(var1.id) == 
        TRUE || is.na(var2.id) == TRUE) {
      stop("Error: Not define the var1.id or var1.id")
    }
  }
  if (style == "adjacent.pairs") {
    var1.id <- seq(1, (VariableNum - 1), by = 1)
    var2.id <- var1.id + 1
  }
  if (style == "one.pair" || style == "adjacent.pairs") {
    if (length(var1.id) != length(var2.id)) {
      stop("Error: var1.id and var2.id should be vectors with the same length")
    }
    taskmatrix <- matrix(c(var1.id, var2.id), ncol = 2)
  }
  if (style == "pairs.between") {
    if (length(which(is.na(var1.id) == TRUE)) > 0 | length(which(is.na(var1.id) == 
                                                                 TRUE)) > 0) {
      stop("Error: no variable IDs are given")
    }
    if (length(which((var1.id != var2.id) == TRUE)) > 0) {
      stop("Error: var1.id should be the same with var2.id for the pairs.between style")
    }
    if (length(which(is.numeric(var1.id) == FALSE)) > 0 | 
        length(which(is.numeric(var2.id) == FALSE)) > 0) {
      stop("Error:var1.id and var2.id should be numeric vector")
    }
  }
  if (style == "all.pairs") {
    var1.id <- seq(1, dim(GEMatrix)[1], by = 1)
    var2.id <- var1.id
  }
  if (style == "pairs.between" || style == "all.pairs") {
    CurLen <- length(var1.id)
    taskmatrix <- matrix(0, nrow = length(var1.id) * (length(var1.id) + 
                                                        1)/2, ncol = 2)
    kk = 0
    for (i in 1:length(var1.id)) {
      j <- 1
      while (j <= i) {
        kk <- kk + 1
        taskmatrix[kk, ] <- c(i, j)
        j <- j + 1
      }
    }
  }
  if (cpus == 1 | cormethod == "BiWt") {
    results <- apply(taskmatrix, 1, cor.pair, GEMatrix = GEMatrix, 
                     rowORcol = "row", cormethod = cormethod, pernum = pernum, 
                     sigmethod = sigmethod)
  }
  else {
    sfInit(parallel = TRUE, cpus = cpus)
    #print(sprintf("%s cpus to be used", sfCpus()))
    results <- sfApply(taskmatrix, 1, cor.pair, GEMatrix = GEMatrix, 
                       rowORcol = "row", cormethod = cormethod, pernum = pernum, 
                       sigmethod = sigmethod)
    sfStop()
  }
  if (output == "paired") {
    kk <- 0
    corpvalueMatrix <- matrix(NA, nrow = dim(taskmatrix)[1], 
                              ncol = 4)
    for (i in 1:dim(taskmatrix)[1]) {
      if (taskmatrix[i, 1] == taskmatrix[i, 2]) 
        next
      kk <- kk + 1
      corpvalueMatrix[kk, 1:2] <- rownames(GEMatrix)[taskmatrix[i, 
                                                                ]]
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
  }
  else {
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
      if (style == "pairs.between" | style == "all.pairs") {
        corMatrix[colidx, rowidx] <- corMatrix[rowidx, 
                                               colidx]
        pvalueMatrix[colidx, rowidx] <- pvalueMatrix[rowidx, 
                                                     colidx]
      }
    }
    return(list(corMatrix = corMatrix, pvalueMatrix = pvalueMatrix))
  }
}
##
cor.pair <- function (idxvec, GEMatrix, rowORcol = c("row", "col"), 
                      cormethod = c("GCC", "PCC", "SCC", "KCC", "BiWt"), 
                      pernum = 0, 
                      sigmethod = c("two.sided","one.sided")) 
{
  onegcc <- function(x, y) {
    getrank <- function(datamatrix) {
      if (dim(datamatrix)[1] != 2) {
        stop("Error: the row num of datamatrix must be 2")
      }
      OrderIndex <- order(datamatrix[1, ], decreasing = FALSE)
      SortGenePair <- datamatrix[, OrderIndex]
      Sort2By1 <- SortGenePair[2, ]
      return(list(Sort2By1 = Sort2By1, Sort2By2 = sort(datamatrix[2, 
                                                                  ], decreasing = FALSE)))
    }
    gcc <- function(weightvec, vectsort, vectselfsort) {
      Sum1 <- sum(weightvec * vectsort)
      Sum2 <- sum(weightvec * vectselfsort)
      if (Sum2 == 0) {
        cat("\n", x, "\n", y, "\n")
        cat("Warning: the Denominator is ZRRO, the value of one variable is consistent.")
        gcccor <- 0
      }
      else {
        gcccor <- Sum1/Sum2
      }
      return(gcccor)
    }
    Length <- length(x)
    Wt <- t(2 * seq(1, Length, by = 1) - Length - 1)
    GenePairXY <- t(matrix(c(x, y), ncol = 2))
    SortYlist <- getrank(GenePairXY)
    gcc.rankx <- gcc(Wt, SortYlist$Sort2By1, SortYlist$Sort2By2)
    return(gcc.rankx)
  }
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
  getcor <- function(g1, g2, cormethod) {
    if (cormethod == "PCC") {
      return(cor.test(g1, g2, method = "pearson")$estimate)
    }
    else if (cormethod == "SCC") {
      return(cor.test(g1, g2, method = "spearman")$estimate)
    }
    else if (cormethod == "KCC") {
      return(cor.test(g1, g2, method = "kendall")$estimate)
    }
    else if (cormethod == "BiWt") {
      return(biwt.cor(t(matrix(c(g1, g2), ncol = 2)), output = "vector")[1])
    }
    else if (cormethod == "GCC") {
      ## TO MODIFY ####### return(list(gcc.rankx = onegcc(g1, g2), gcc.ranky = onegcc(g2, g1)))
    }
  }
  getpvalue <- function(percorvec, pernum, realcor, sigmethod) {
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
  # if (pernum <= 0) {
  #   if (cormethod == "GCC") {
  #     ## TO MODIFY ####### return(list(gcc.rankx = realcor$gcc.rankx, gcc.ranky = realcor$gcc.ranky, 
  #     ## TO MODIFY #######             gcc.rankx.pvalue = NA, gcc.ranky.pvalue = NA))
  #   }
  #   else {
  #     return(list(cor = realcor, pvalue = NA))
  #   }
  # }
  # else {
  #   pGCCMatrix <- matrix(0, nrow = pernum, ncol = 2)
  #   colnames(pGCCMatrix) <- c("gcc.rankx", "gcc.ranky")
  #   rownames(pGCCMatrix) <- paste("permut", seq(1, pernum, 
  #                                               by = 1), sep = "")
  #   GenePairXY <- t(matrix(c(x1, y1), ncol = 2))
  #   pGenePairXY <- GenePairXY
  #   Length <- length(x1)
  #   for (i in 1:pernum) {
  #     curtime <- format(Sys.time(), "%H:%M:%OS4")
  #     XXX <- unlist(strsplit(curtime, ":"))
  #     curtimeidx <- (as.numeric(XXX[1]) * 3600 + as.numeric(XXX[2]) * 
  #                      60 + as.numeric(XXX[3])) * 10000
  #     set.seed(curtimeidx)
  #     TT = sort(runif(Length), index.return = TRUE)$ix
  #     pGenePairXY[1, ] <- GenePairXY[1, TT]
  #     if (cormethod == "GCC") {
  #       cortmp <- getcor(pGenePairXY[1, ], pGenePairXY[2, 
  #                                                      ], cormethod)
  #       pGCCMatrix[i, 1] <- cortmp$gcc.rankx
  #       pGCCMatrix[i, 2] <- cortmp$gcc.ranky
  #     }
  #     else {
  #       pGCCMatrix[i, 1] <- getcor(pGenePairXY[1, ], 
  #                                  pGenePairXY[2, ], cormethod)
  #     }
  #   }
    if (cormethod == "GCC") {
      ## TO MODIFY #######  return(list(gcc.rankx = realcor$gcc.rankx, gcc.ranky = realcor$gcc.ranky, 
      ## TO MODIFY #######             gcc.rankx.pvalue = getpvalue(pGCCMatrix[, 1], 
      ## TO MODIFY #######                                          pernum, realcor$gcc.rankx, sigmethod), gcc.ranky.pvalue = getpvalue(pGCCMatrix[, 
                                                                                                                              2], pernum, realcor$gcc.ranky, sigmethod)))
    }
    else {
      return(list(cor = realcor, pvalue = getpvalue(pGCCMatrix[, 
                                                               1], pernum, realcor, sigmethod)))
    }
  }
}

