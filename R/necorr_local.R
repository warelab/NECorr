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
#network.file <- "../data/network.txt"
#expression <- "../data/expression.txt"
#description.file <- "../data/description.csv"
#metadata <- "../data/metadata.txt"

  network.file <- "../data/grnmetanet.txt"
  expression <- "../data/gene_expression_matrix.txt"
  description.file <-"../data/1.Ath.GeneDesc.csv"
  metadata <- "../data/grnmeta.metadata.txt"
  

# setwd("/Users/cliseron/Documents/1_Repository/NECorr/R/")

Necorr <- function(network.file, expression, description.file,
                   condition, metadata, name, 
                   Filelist, #condition list see if still necessary with metadata
                   method = "GCC", permutation = 1000, sigcorr = 0.01,
                   fadjacency = "only",type = "gene",
                   dirtmp="./results/tmp", dirout = './results',
                   NSockets = 8){
  #' @author Christophe Liseron-Monfils
  #' @param expression Expression file in log2 (ratio expression) with row: gene,
  #' first column: type of sample,second column: sample names
  #' @param network.file Molecular network file with source in the first column, targets in
  #'  the second column
  #' @param description.file genome description
  #' @param condition Condition from expression to study the network co-expression
  #' correlation
  #' @param type Omics comparative expression type: protein or gene
  #' @param permutation permutation number used for all significance calculation
  #' @param lmiR List of miRNAs
  #' @param method Method used for co-expression correlation: GCC, MINE, PCC,
  #' SCC or KCC
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
  
  # Create the output directory if not existing; generate "./results" dir and
  # "./results/tmp"
  if(!dir.exists(dirout)){
    dir.create(dirout)
  }
  if(!dir.exists(paste0(dirout,"/",condition))){
    dir.create(paste0(dirout,"/",condition), recursive = T)
  }
  if(!dir.exists(dirtmp)){
    dir.create(dirtmp, recursive = T)
  }
  # AT THE END ADD THE OTHER FUNCTIONS FROM SOURCE
  #source("./src/NECorr.functions.v3.R")
  load("./ML_model.RData")
  suppressWarnings(suppressPackageStartupMessages(require(RColorBrewer)))
  suppressWarnings(suppressPackageStartupMessages(require(gplots)))
  library(Rcpp)
  sourceCpp("../src/gini.cpp")
  # load file and options
  # read network file
  network.int <- read.delim(network.file, sep ="\t", header = T,fileEncoding="latin1")
  # Condition experiment and/or tissue
  
  factortab <- read.table(metadata,header = T) #factor.file #condition <- "Radial"
  # Description file name with gene name and annotations
  Desc <-  read.csv(description.file,header=T, row.name=1) #description.file <- "1.Ath.GeneDesc.csv"
  # Create the subdirectory for the fianl results
  # need to see if all these tables are still useful??
  subDirTS = paste0(dirout, "/", condition,"/4_TS_file/")
  subDirGraph = paste0(dirout, "/", condition, "/6_TS_graph/")
  subDirFile = paste0(dirout, "/", condition, "/7_gene_ranking_per_condition/")
  dir.create(file.path(subDirTS), showWarnings = FALSE)
  dir.create(file.path(subDirGraph), showWarnings = FALSE)
  dir.create(file.path(subDirFile), showWarnings = FALSE)
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
  
  sample.names <- unique(factortab[,1]) 
  print(sample.names)
  xcol <-c()
  treatment.f <- factor()
  for (i in 1:length(sample.names)){ 
    nrep <- length(which(factortab[,1] == sample.names[i]))
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
  netname <- basename(network.file)
  netstatFile <- paste0(dirout,"/",condition,"/5_",netname,"_NetStat.csv")
  if(file.exists(netstatFile)){
    print(paste0("reusing netstat output file ",netstatFile))
    netstat <- read.csv(netstatFile,header=T, row.name=1)
  }
  else {
    print("running NetTopology")
    netstat <- NetTopology(network.int)
    write.csv(netstat,netstatFile)
  }
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
  #print(Genelist)
  print("R network topology done")
 
  # Read expression file
  eset <- read.table(expression,header = T,row.names=1) #!! change made Here
  # Take only the genes that are part of the molecular network
  eset <- eset[Genelist,] #!! change made Here
  m.eset <- as.matrix(eset)  
  #print(eset)
  m.eset <- m.eset[-grep("NA", rownames(m.eset)), xcol]
  # Loop to measure the importance of gene expression 
  conditionList <-factortab[,2]
  factorList <- factortab[,1]
  df <- cbind(as.character(conditionList),as.character(factorList))
  df <- df[!duplicated(df),]

  suppressWarnings(suppressPackageStartupMessages(require(foreach)))
  suppressWarnings(suppressPackageStartupMessages(require(doSNOW)))
  #files.list <- paste0(dirout,"/",Filelist)
  #print(files.list)
  #files <- read.table(files.list)
  #condnb <- nrow(files)
  condnb <- length(conditionList)
  ## Loop to do the anlysis per condition defined in the first row of the expression data
  #CondNB <- 1  ### test of the code without loop
  for (CondNB in 1:condnb){
    #fileexp <- as.character(files[CondNB,1])
    fileexp <- as.character(basename(expression))
    #x.exp <- as.matrix(read.table(paste0(dirout,"/",fileexp)))
    #end.time <- Sys.time()
    #time.taken <- end.time - start.time
    #print(time.taken)
    
    ###------------------------------------------------------------------------------------------------------------------------
    print("### II - Analysis of the co-expression file and p-value sums")
    start.time <- Sys.time()
    ###------------------------------------------------------------------------------------------------------------------------
    #int.sig <- bigcorGCC(x.exp, net = network.int, nsockets = NSockets,
    filecoexp <- paste0(as.character(basename(expression)),"_CorrM_PvalM.txt")
    pathcoexp <- paste0(dirout, "/", condition, "/", filecoexp)
    if (file.exists(pathcoexp)) {
      print(paste0("reusing coexp output file ",pathcoexp))
      int.sig <- read.csv(pathcoexp,header=T, row.name=1)
    }
    else {
      print("calculating coexpression")
      int.sig <- bigcorGCC(m.eset, net = network.int, nsockets = NSockets,
                          methods= method, output = "paired", sigmethod = "two.sided", 
                          pernum = permutation , verbose = FALSE, cpus = NSockets)
   
    
    #filecoexp <- paste0(as.character(files[CondNB,1]),"_CorrM_PvalM.txt")
      int.sig.file <- int.sig
      colnames(int.sig.file) <- c("Source","Target","Correlation","p-value")
      write.csv(int.sig.file, pathcoexp)
    }
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(time.taken)
    # create a function to generate a continuous color palette
    rbPal <- colorRampPalette(c('yellow','blue'))
    
    # This adds a column of color values based on the y values
    # remove the line with NA due to not found gene due to no expresssion
    # or error of annotation in the initial network file
    int.sig <- int.sig[complete.cases(int.sig),]
    
    # Transform the p-value in the maximal p-value in function of number of permutations to avoid log-scale = infinity
    # Replace the p-value equal to 0 by the maximal p-value e.g. 
    #1/10000 knowing that we have 10000 randomizations in initial calculations
    max.p.val = 1/as.numeric(permutation)
    int.sig[which(int.sig[,4] == 0),4] = max.p.val
    int.sig <-as.data.frame(int.sig)
    p.int.sig  = as.numeric(as.character(int.sig[,4]))
    pval <- -log(p.int.sig,10)
    int.sig$Col <- rbPal(10)[as.numeric(cut(pval,breaks = 10))]
    title <- paste0(as.character(basename(expression)),"Sig.Interaction.pdf")
    
    pdf(title)
    xlabel = paste0(method," score")
    #plot(as.numeric(as.character(int.sig[,3])),pval,pch = ".", xlab=xlabel,ylab= "Significance (-log scale)",col= int.sig$Col )
    #abline( h=1.3, col ="red",lty='dotted')
    #legend("top", pch= c(21,21,NA), title="Interaction Significance",
    #       bty = "n",
    #       col = c("yellow", "blue", "red"),
    #       pt.bg= c("yellow", "blue", NA),
    #       lty = c(0, 0, 3), lwd = c(0, 0, 1),
    #       legend=c("low","high", "pval = 0.05"),
    #       cex =0.8)
    #dev.off()
    
    # main analysis and weight
    for (i in 1:length(df[,1])){
      nsockets = NSockets
      #i <- 1
      CondName <- df[i,1]
      CondPutFile <- paste0("3_Expression",name,"_",netname,"_",condition, 
                            "_gene_",method,"_",permutation,
                            "_",CondName,"_CorrM_PvalM.txt")
      sample.l <- df[i,2]
      print(paste( "the condition is",CondName,"the factor is:",sample.l,sep=" "))
      print(paste( "correlation",CondPutFile,"files:",filecoexp,sep=" "))
      # interaction file
      
      #end.time <- Sys.time()
      #time.taken <- end.time - start.time
      #print(time.taken)
      ###------------------------------------------------------------------------------------------------------------------------
      print("### III - Define gene tissue specificity index (TSI) (Yanai 2011, Bioinformatics)")
      #start.time <- Sys.time()
      ###------------------------------------------------------------------------------------------------------------------------
      
      # Add the weight for the studied tissue
      # tissue selectivity or tissue exclusion from the tissue should be considered
      # as both can important at genetic level repression or activation of a gene
      # are part of the tissue specificity
      
      nreplics <- table(treatment.f) 
      ####  Need to find a way to order factor by order of appearance in expression file
      # e.g nreplics <- c(4,4,4)
      
      m.eset = as.matrix(m.eset)
      
      # Comparison to find tissue selectivity for the selected condition
      # using Intersection-Union Test (Berger et al)
      act <- IUT(m.eset, nreplics, which(sample.names == sample.l), 1, 0.5, -1)
      rep <- IUT(m.eset, nreplics, which(sample.names == sample.l), -1, 0.5, -1)
      #### do following steps if the expression if act and/or rep are not null
      table(treatment.f)
      
      # Rank the tissue selective genes using the Tissue Selective Index
      # for each gene TSI for activation and modify for tissue repression
      m.eset = as.data.frame(m.eset)
      meansFactor <- do.call("cbind", tapply(1:ncol(m.eset), 
                                             treatment.f, 
                                             function(x) rowMeans(m.eset[x])))
      # Tissue-selective genes(up)
      if(length(act[which(act==TRUE)])>0){
        act <- act[which(names(act) %in% rownames(meansFactor))]
        StA <- length(meansFactor[names(act[which(act==TRUE)]),])
        StB <- length(colnames(meansFactor))
        if( StA > StB){
          tsi = apply(meansFactor[which(act==TRUE),],1,TSI)
          tsi.order <- tsi[order(tsi,decreasing = TRUE)]
        }else if( StA == StB){ # if there is only one tissue-selective genes
          tsi = TSI(meansFactor[which(act==TRUE),])
          names(tsi) <- rownames(meansFactor)[which(act==TRUE)]
          tsi.order =  tsi
        }else if(StA < StB){ # if there is none tissue-selective genes
          tsi = c()
          tsi.order = c()
        }
      }else{
        tsi = c()
        tsi.order = c()       
      }
      # Tissue-selective genes(down)
      if(length(rep[which(rep==TRUE)])>0){
        rep <- rep[which(names(rep) %in% rownames(meansFactor))]
        StA <- length(meansFactor[which(rep==TRUE),])
        StB <- length(colnames(meansFactor))
        if( StA > StB){
          tsr = apply(meansFactor[which(rep==TRUE),],1,TSR)
          tsr.order <- tsr[order(tsr,decreasing = TRUE)]
        }else if( StA == StB){ # if there is only one tissue selective gene
          tsr = TSR(meansFactor[which(rep==TRUE),])
          names(tsr) <- rownames(meansFactor)[which(rep==TRUE)]
          tsr.order =  tsr
        }else if(StA < StB){ # if there is none tissue selective genes
          tsr = c()
          tsr.order = c()
        }
      }else{
        tsr = c()
        tsr.order = c()
      }
      # Vector of tissue-selectivity
      ts <- c(tsr,tsi) #2 expression specificity
      if(!is.null(ts)){
        write.csv(meansFactor[names(c(tsr.order,tsi.order)),],
                  paste0(subDirTS,sample.l,CondName,netname,"_",method,"_",permutation,"TS_ranking.csv"))
        
        title <- paste0(subDirGraph,sample.l,"_",CondName,"_",method,"_",permutation,"TS_ranking.pdf")
        y <- as.matrix(meansFactor[names(c(tsr.order,tsi.order)),])
        if(length(names(c(tsr.order,tsi.order))) > 10){
          decal.y <- -3*(length(names(c(tsr.order,tsi.order)))/length(Genelist))
        }else{
          decal.y <- 0.2
        }
        # if((nrow(meansFactor[names(c(tsr.order,tsi.order)),] > 0) == TRUE) | !is.null(meansFactor[names(c(tsr.order,tsi.order)), ])){
        #   pdf(title, height=10, width=15)
        #   heatmap.2(as.matrix(meansFactor[names(c(tsr.order,tsi.order)),]), Rowv =FALSE,
        #             Colv = FALSE,scale ="row",
        #             main = paste0(sample.l,CondName," Condition selective genes"),
        #             density.info= "none",cexCol= 0.8,labRow='',
        #             trace= "none",dendrogram = "none",
        #             margins = c(15, 4),
        #             lhei = c(1, 8),
        #             keysize= .6,labCol = '',
        #             add.expr = text(x = seq_along(colnames(y)),
        #                             y = decal.y,
        #                             srt = 60,cex=.8,
        #                             labels = colnames(y), xpd = NA,adj=0, pos =2),
        #             col=colorpanel(121,"lightyellow","yellow","darkblue")
        #   )
        #   dev.off()
        #}

      }
      #end.time <- Sys.time()
      #time.taken <- end.time - start.time
      #print(time.taken)
      
      ###------------------------------------------------------------------------------------------------------------------------
      print("### IV - Determine the overall interaction importance for each node (gene)")
      start.time <- Sys.time()
      ###------------------------------------------------------------------------------------------------------------------------
      importanceFile = paste0(subDirTS,sample.l,CondName,netname,"_",method,"_",permutation,"importance.csv")
      if (file.exists(importanceFile)) {
        print(paste0("reusing importance output file ",importanceFile))
        int.pvals <- read.csv(importanceFile,header=T, row.name=1)
      }
      else {
        # for the interactions coming from the same node
        # create hash with all the gene in the network as keys
        int.pvals = structure(rep(1, length(Genelist)), names=Genelist)   #1 interaction p-values
        suppressWarnings(suppressPackageStartupMessages(require(hash)))
        h <- hash()
        for (i in 1:length(int.sig[,1])) {
          g1 = as.character(int.sig[i,1])
          g2 = as.character(int.sig[i,2])
          twistedcor = 1 - abs(int.sig[i,3]);
          if (twistedcor == 0) {
            twistedcor = max.p.val
          }
          if (all(has.key(g1, h))) {
            h[[g1]] <- append(h[[g1]], twistedcor)
          }
          else {
            h[[g1]] <- twistedcor
          }
          if (all(has.key(g2, h))) {
            h[[g2]] <- append(h[[g2]], twistedcor)
          }
          else {
            h[[g2]] <- twistedcor
          }
        }
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        print(time.taken)
        print(" -- populated hash")
        # add a vector with all the p-value attached to a gene
        # in the co-expression analysis
        start.time <- Sys.time()
        for (i in 1:length(Genelist)){
          gene = Genelist[i] # get the name of the gene
          trans.cumpvals = fishersMethod(as.numeric(as.vector(h[[gene]])))
          cumpvals = - log(trans.cumpvals,10)  # transform pval in factor and Take -log10 of the results
          int.pvals[gene]=cumpvals
        }
        print(" -- calculated fishersMethod")
        write.csv(int.pvals,importanceFile)
      }
      end.time <- Sys.time()
      time.taken <- end.time - start.time
      print(time.taken)
      ###------------------------------------------------------------------------------------------------------------------------
      print("### V Differential Expression ranking for the gene in the network for the factor")
      #start.time <- Sys.time()
      ###------------------------------------------------------------------------------------------------------------------------
      
      DE.ranks <- DE.ranking(m.eset,Genelist,treatment.f,sample.l,sample.names)
      #the results are already scaled
      
      #end.time <- Sys.time()
      #time.taken <- end.time - start.time
      #print(time.taken)
      ###------------------------------------------------------------------------------------------------------------------------
      print("### VI Normalize each column between [0,1] using the some")
      #start.time <- Sys.time()
      ###------------------------------------------------------------------------------------------------------------------------
      
      int.pvals = ScalN(int.pvals)  #1 interaction p-values
      ts = ScalN(ts) #2 expression specificity
      BetwC = ScalN(BetwC) #3 betweeness
      Conn = ScalN(Conn) #4 connectivity
      ClusCoef = ScalN(ClusCoef)  #5 transitivity
      PageRank <- ScalN(PageRank)
      
      int.pvals = as.data.frame(int.pvals)  #1 interaction p-values
      ts = as.data.frame(ts); #2 expression specificity
      BetwC = as.data.frame(BetwC) #3 betweeness
      Conn = as.data.frame(Conn) #4 connectivity
      ClusCoef = as.data.frame(ClusCoef)  #5 transitivity
      EigenC = as.data.frame(EigenC)
      PageRank = as.data.frame(PageRank)
      
      #end.time <- Sys.time()
      #time.taken <- end.time - start.time
      #print(time.taken)
      ###------------------------------------------------------------------------------------------------------------------------
      print("### VII - Hub NECorr merge all the subnetwork statistics for 
            the studied tissue or condition in one table column")
      #start.time <- Sys.time()
      ###------------------------------------------------------------------------------------------------------------------------
      
      m.param = merge(int.pvals,ts,by="row.names",all=T)
      rownames(m.param) = m.param$Row.names
      m.param = merge(m.param[,-1],BetwC,by="row.names",all=T)
      rownames(m.param) = m.param$Row.names
      m.param = merge(m.param[,-1],Conn,by="row.names",all=T)
      rownames(m.param) = m.param$Row.names
      m.param = merge(m.param[,-1],ClusCoef,by="row.names",all=T)
      rownames(m.param) = m.param$Row.names
      m.param = merge(m.param[,-1],EigenC,by="row.names",all=T)
      rownames(m.param) = m.param$Row.names
      m.param = merge(m.param[,-1],PageRank,by="row.names",all=T)
      rownames(m.param) = m.param$Row.names
      m.param = merge(m.param[,-1],DE.ranks,by="row.names",all=T)
      rownames(m.param) = m.param$Row.names
      m.param[is.na(m.param)] = 0
      par.tab <- m.param
      colnames(m.param) = c("GeneID","interaction.pvals",
                            "tissue.treatment.specificity","Betweenness.Centrality",
                            "Connectivity","Transitivity","Eigenvector","PageRank","DE")
      #tail(m.param)
      par.tab <- m.param
      write.table(par.tab,
                  paste0("results/",condition,"/8_hub_gene_Std_param.",
                         sample.l,"_",CondName,"_",netname,"_",method,
                         "_",permutation,".txt"),
                  sep = "\t",quote = FALSE,row.names = FALSE)
      rm(par.tab)
      m.param = m.param[,-1]
      #end.time <- Sys.time()
      #time.taken <- end.time - start.time
      #print(time.taken)
      ###------------------------------------------------------------------------------------------------------------------------
      print("### VIII - calculate the weight for: (use ahp to generate weight)")
      #start.time <- Sys.time()
      ###------------------------------------------------------------------------------------------------------------------------
      
      #done before the loop wpar variable
      #end.time <- Sys.time()
      #time.taken <- end.time - start.time
      #print(time.taken)
      
      ###------------------------------------------------------------------------------------------------------------------------
      print("### IX Hub NECorr - pick the top genes for validations (best alternatives)")
      ###------------------------------------------------------------------------------------------------------------------------
      
      # calculate the weighted ranking for each gene = alternative ranking
      hub.m.param <-  m.param[,c( "interaction.pvals",
                                  "tissue.treatment.specificity","Betweenness.Centrality","Connectivity","Transitivity")]
      for (i in 1:length(colnames(hub.m.param))){
        hub.m.param[,i] = hub.m.param[,i]*wpar[i]
      }
      
      # gene ranking
      gene.rank.h = apply(hub.m.param,1,sum)*100
      gene.rank.h = gene.rank.h[order(gene.rank.h, decreasing=TRUE)]
      nGenes = length(gene.rank.h)
      print("building hub rank hash")
      gene.rank.h = as.data.frame(gene.rank.h)
      #load the ranks into a hash to make step XI faster
      suppressWarnings(suppressPackageStartupMessages(require(hash)))
      gene.rank.hash <- hash()
      geneIDs = row.names(gene.rank.h)
      for (i in 1:nGenes) {
        geneRank = as.numeric(gene.rank.h[i,1])
        gene.rank.hash[[geneIDs[i]]] = geneRank
      }
      
      gene.rank.h.description = merge(gene.rank.h,Desc,all.x = T, by = "row.names",sort=FALSE)
      colnames(gene.rank.h.description)[1]<- "GeneID"
      write.csv(as.data.frame(gene.rank.h.description),paste0(subDirFile,sample.l,
                                                              netname,permutation,
                                                              "_Hub_gene_prioritization.csv")
                ,quote=FALSE,row.names=FALSE)
      colnames(gene.rank.h)[1] <- sample.l
      total.rank.h <- merge(total.rank,gene.rank.h,all.x=T, by= "row.names")
      rownames(total.rank.h) <- total.rank.h[,1]
      total.rank.h <- total.rank.h[,-1]
      
      #end.time <- Sys.time()
      #time.taken <- end.time - start.time
      #print(time.taken)
      ###------------------------------------------------------------------------------------------------------------------------
      print("### X - Effector NECorr - pick the top genes for validations (best alternatives)")
      #start.time <- Sys.time()
      ###------------------------------------------------------------------------------------------------------------------------
      
      # calculate the probability = alternative ranking
      eff.m.param <-  m.param[,c( "interaction.pvals","tissue.treatment.specificity","Betweenness.Centrality",
                                  "Connectivity","Eigenvector","PageRank","DE")]
      
      ### MACHINE LEANRING with paramter deine from validation
      ### Predict probality to be important gene
      
      #!!!! if Naive Bayes
      suppressWarnings(suppressPackageStartupMessages(library(klaR)))
      #j.nB <- NaiveBayes(phenotype ~ . , data = eff.m.param, kernel = "rectangular", n = 148)
      j.nB <- model
      prob.pred  <- predict(j.nB, type="prob",newdata = eff.m.param )
      
      gene.rank.eff <- (prob.pred$posterior)[,2]
      gene.rank.eff = gene.rank.eff[order(gene.rank.eff, decreasing=TRUE)]
      gene.rank.eff = as.data.frame(gene.rank.eff)
      
      gene.rank.e.description = merge(gene.rank.eff,Desc,all.x = T, by = "row.names",sort=FALSE)
      colnames(gene.rank.e.description)[1]<- "GeneID"
      write.csv(as.data.frame(gene.rank.e.description),paste0(subDirFile,sample.l,netname,permutation,"_effector_gene_prioritization.csv")
                ,quote=FALSE,row.names=FALSE)
      colnames(gene.rank.eff)[1] <- sample.l
      total.rank.eff <- merge(total.rank,gene.rank.eff,all.x=T, by= "row.names")
      rownames(total.rank.eff) <- total.rank.eff[,1]
      total.rank.eff <- total.rank.eff[,-1]
      
      gene.rank.eff.hash <- hash()
      print("building eff rank hash")
      geneIDs = row.names(gene.rank.eff)
      for (i in 1:nGenes) {
        geneRank = as.numeric(gene.rank.eff[i,1])
        gene.rank.eff.hash[[geneIDs[i]]] = geneRank
      }
      
      #end.time <- Sys.time()
      #time.taken <- end.time - start.time
      #print(time.taken)
      ###------------------------------------------------------------------------------------------------------------------------
      print("### XI - Hub interaction ranking")
      #start.time <- Sys.time()
      ###------------------------------------------------------------------------------------------------------------------------
      # hub sub-network
      # cl2 <- makeCluster(nsockets, type="SOCK")
      # registerDoSNOW(cl2)
      # hub.int.ranks <- foreach(i=1:nrow(network.int),.combine='rbind' )%dopar%{
        sourceIDs <- as.vector(network.int[,1])
        targetIDs <- as.vector(network.int[,3])
        ranks.sum <- rep(1,length(targetIDs))
      for(i in 1:nrow(network.int)) {
        ranks.sum[i] <- sum(gene.rank.hash[[sourceIDs[i]]], gene.rank.hash[[targetIDs[i]]])
      }
      hub.int.ranks <- data.frame(sourceIDs,targetIDs,ranks.sum)
      # stopCluster(cl2)
      break.points <- c(-Inf, unique(sort(as.numeric(hub.int.ranks[,2]))), Inf)
      p2 <- cut( as.numeric(hub.int.ranks[,2]), breaks=break.points, labels=FALSE )
      p2 <- 1 - p2/length(break.points)
      hub.int.ranks <- as.data.frame(cbind(hub.int.ranks,p2))
      hub.int.ranks$p2 <- as.numeric(as.character(hub.int.ranks$p2))
      hub.int.significant = subset(hub.int.ranks, p2 < 0.005)
      write.csv(hub.int.ranks,paste0(subDirFile,sample.l,netname,permutation,"_interactions_prioritization.csv"),quote=FALSE,row.names=FALSE)
      # effector sub-network
      # cl3 <- makeCluster(nsockets, type="SOCK")
      # registerDoSNOW(cl3)
      # eff.int.ranks <- foreach(i=1:nrow(network.int),.combine='rbind' )%dopar%{
      eff_ranks.sum <- rep(1,length(targetIDs))
      for(i in 1:nrow(network.int)) {
        # SourceID <- as.character(network.int[i,1])
        # TargetID <- as.character(network.int[i,2])
        eff_ranks.sum <- sum(gene.rank.eff.hash[[sourceIDs[i]]], gene.rank.eff.hash[[targetIDs[i]]])
      }
      eff.int.ranks <- data.frame(sourceIDs,targetIDs,eff_ranks.sum)
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
      ###------------------------------------------------------------------------------------------------------------------------
      print("### XII - Activator NECorr - based on genes in the complete network that are linked to the hub genes")
      #start.time <- Sys.time()
      ###------------------------------------------------------------------------------------------------------------------------
      # find genes that are significant in the hub subnetwork in the complete network
      sig.hub <- unique(c(as.character(hub.int.significant$V1), as.character(hub.int.significant$V2)))
      sc.sig.hub <- subset(network.int, network.int$source %in% sig.hub)
      tg.sig.hub <- subset(network.int, network.int$target %in% sig.hub)
      net.extension.sig.hub <- rbind(sc.sig.hub,tg.sig.hub)
      
      ## change the names of the hub gene in the extented hub network using source genes
      change.sig.hubNames <- net.extension.sig.hub[,1] %in% sig.hub
      tmp<-as.character(net.extension.sig.hub[,1])
      tmp[change.sig.hubNames]<-"sig"
      net.extension.sig.hub[,1] <- tmp
      
      ## change the names of the hub gene in the extented hub network using target genes
      change.sig.hubNames <- net.extension.sig.hub[,2] %in% sig.hub
      tmp<-as.character(net.extension.sig.hub[,2])
      tmp[change.sig.hubNames]<-"sig"
      net.extension.sig.hub[,2] <- tmp
      
      ## network extended to putative regulator using source genes
      act.m.param <- net.extension.sig.hub
      sc.count <- rle( sort( act.m.param[,1] ))
      act.m.param$Count <- sc.count[ match( act.m.param[,1] , sc.count ) ]
      
      # gene ranking of linked to hub nodes
      gene.rank.act <-  cbind(sc.count$values, sc.count$lengths)
      gene.rank.act <- gene.rank.act[order(as.numeric(gene.rank.act[,2]), decreasing=TRUE),]
      gene.rank.act <- gene.rank.act[ - which(gene.rank.act[,1] == "sig"),]
      
      # add the gene that have 25% of the gene linked to hub genes
      gene.rank.act.significant <- gene.rank.act[which(gene.rank.act[,2] >= (length(sig.hub)*0.20)),]
      
      # write the putative activator genes
      gene.rank.act.b <- as.data.frame(gene.rank.act)
      gene.rank.act.description = merge(gene.rank.act.b,Desc,all.x = T, by.x = "V1", by.y = 0,sort=FALSE)
      colnames(gene.rank.act.description)[1]<- "GeneID"
      write.csv(as.data.frame(gene.rank.act.description),paste0(subDirFile,sample.l,netname,permutation,"_act_gene_prioritization.csv")
                ,quote=FALSE,row.names=FALSE)
      
      #end.time <- Sys.time()
      #time.taken <- end.time - start.time
      #print(time.taken)
      ###------------------------------------------------------------------------------------------------------------------------
      print("### XIII  - Vizualization of network and list of interaction that could affected and color")
      ###------------------------------------------------------------------------------------------------------------------------
      #start.time <- Sys.time()
      suppressWarnings(suppressPackageStartupMessages(library(igraph)))
      suppressWarnings(suppressPackageStartupMessages(library(dnet)))
      hub.net <- as.data.frame(cbind(as.character(hub.int.significant$V1), as.character(hub.int.significant$V2),
                                     as.numeric(as.character(hub.int.significant$V3)), rep("hub", nrow(hub.int.significant))))
      colnames(hub.net) <- c("source","target","score","node.type")
      # significant activator sub-netowrk linked to hub network
      net.extension.sig.hub <- rbind(sc.sig.hub,tg.sig.hub)
      act.net.1 <- subset(net.extension.sig.hub, net.extension.sig.hub[,1] %in% as.vector(as.character(gene.rank.act.significant[,1])))
      act.net.2 <- subset(net.extension.sig.hub, net.extension.sig.hub[,2] %in% as.vector(as.character(gene.rank.act.significant[,1])))
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
        print(hub.act.net)  ############################################
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
          title <- paste0(subDirGraph,sample.l,"_",CondName,"_",
                          method,"_",permutation,"interaction_graph.pdf")
          # 	    mark.groups <- vcolors
          # 	    mark.col <- visColoralpha(vcolors, alpha=0.2)
          # 	    mark.border <- visColoralpha(vcolors, alpha=0.2)
          #
          # 	    mark.groups= mark.groups,mark.col= mark.col, mark.border=mark.border,
          pdf(file = title, width = 10, height = 10)
          dnet::visNet(g, glayout=layout.fruchterman.reingold(g) , 
                       vertex.shape="sphere", vertex.size = vsize, 
                       vertex.label = vlabel, edge.color = "grey",
                       edge.arrow.size = 0.3, vertex.color = vcolors,
                       vertex.frame.color = vcolors, newpage = F)
          dev.off()
        }
      }
    }
  }
}
########################

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
                               cormethod = c("GCC", "PCC", "SCC", "KCC"), 
                               style = c("pairs.between","pairs.only"), 
                               var1.id = NA, var2.id = NA, 
                               pernum = 0, 
                               sigmethod = c("two.sided","one.sided"), 
                               output = c("matrix","paired")){
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
    return(list(corMatrix = corMatrix, pvalueMatrix = pvalueMatrix))
  }
}

bigcorGCC <- function(x ,net= NA, nsockets= 4, methods = c("GCC","PCC","SCC","KCC"),
                      sigmethod = c("two.sided", "one.sided"),
                      nblocks = 10, verbose = TRUE, cpus = 1, pernum = 0, ...){
  corMAT <- c()
  if (methods == "GCC"){
    #### NEW GINI CORRLATION CALCULATION
    corMAT <- gini(edges=net, expression=x, bootstrapIterations=pernum, statCutoff=0.6)
  }else{
    suppressWarnings(suppressPackageStartupMessages(require(foreach)))
    suppressWarnings(suppressPackageStartupMessages(require(doSNOW)))
    nsockets <- as.numeric(nsockets)
    cl <- makeCluster(nsockets, type="SOCK")
    registerDoSNOW(cl)
    
    Nrow = nrow(net)
    fc <- gl(nblocks, ceiling(Nrow/nblocks), length = Nrow)
    matnet <- split(net,fc)
    corMAT<-foreach(j=1:nblocks, .combine='rbind', .export=c('indexing.network','cor.matrix.NECorr','cor.pair'))%dopar%{
      netindexed <- indexing.network(as.matrix(x) ,matnet[[j]])
      G1 <- as.numeric(as.vector(netindexed[,1]))
      G2 <- as.numeric(as.vector(netindexed[,2]))
      cor.matrix.NECorr(x, var1.id=G1, var2.id=G2, sigmethod=sigmethod, cormethod=methods, 
                        pernum=pernum, output="paired", cpus=cpus, style="pairs.only")
    }
    stopCluster(cl)
  }
  return(corMAT)
}

expand.grid.unique <- function(x, y, include.equals=FALSE){
  x <- unique(x)
  y <- unique(y)
  g <- function(i){
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
DE.ranking <- function(exps,GeneList,factortab,sample.l,sample.names, exps.file = FALSE){
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
    sample.names <- unique(targets[,1])
    f <- factor(targets[,1], levels=sample.names)
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
NetTopology <-function(network){
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
IUT<-function(DATA,nrreplics,target,updown,alpha,mc){
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

## calculation of correlation for matrix creating adjacency table
NetCor <- function(file,method = c("GCC", "PCC", "SCC", "KCC"),nbperm,cpu) {
  tab <- read.table(file)
  tab <-as.matrix(tab)
  x <- Rscorr(tab,method = method,nbperm,1)
  write.table(x$corMatrix,paste(file,"_CorrM.txt", sep =""), sep = "\t", row.names =TRUE, col.names=NA,quote =FALSE)
  write.table(x$pvalueMatrix, paste(file,"_PvalM.txt", sep=""), sep = "\t", row.names =TRUE, col.names=NA,quote =FALSE)
}
##
cor.pair <- function (idxvec, GEMatrix, rowORcol = c("row", "col"), 
                      cormethod = c("GCC", "PCC", "SCC", "KCC"), 
                      pernum = 0, 
                      sigmethod = c("two.sided","one.sided")) {
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
  if (pernum <= 0) {
    return(list(cor = realcor, pvalue = NA))
  }else {
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
       pGCCMatrix[i, 1] <- getcor(pGenePairXY[1, ], 
                                  pGenePairXY[2, ], cormethod)
       }
     return(list(cor = realcor, pvalue = getpvalue(pGCCMatrix[,1], pernum, realcor, sigmethod)))
  }
}

Necorr(network.file=network.file, expression=expression, 
       description.file=description.file,
       condition="meristem_young_leaves_4_week_old",metadata=metadata, name="test")
