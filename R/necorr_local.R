# core functions necessary to apply the NECorr pipeline

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

factorfile <- function(exps, nb, condition , dirout = "./results/tmp/"){
#' @author Christophe Liseron-Monfils
#' @param exps Expression file in log2(ratio expression) with row:gene,
#' first column: type of sample, second colum: sample names
#' @param condition Different sample present in the expression data: stress
#' condition, tissue types and/or developmental stages
#' @param nb
#' @param dirout results directory
#' @description Generate the factor file for the co-expression analysis

# Keep only the core filename by reasing the prefix dir and the suffix
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


mainNecorr <- function(lmiR = lmiR,  ){
   #' @author Christophe Liseron-Monfils
   #' @param
   #' @param
   #' @param
   #'

   prepareNECorr(lmiR = lmiR, condition, exps, net, outdir)
   rankingNEcorr()
   cleaningNECorr()

}

prepareNECorr <- function(lmiR, condition, exps, net, type, outdir,
                          method, fadjacency, perm){
  #' @author Christophe Liseron-Monfils
  #' @param
  #'

  # Preparation for the co-expression analysis using
  # the whole adjacency matrix (-all)
  # or the part that are in the network (-only)

  # 1 - Load the expression files

  #    1a) - Collect miR present in the expression file
  BsAfftab <- read.delim(lmiR, sep = "\t", row.names = 1)
  miRnames <- colnames(BsAfftab)

  #    1b) Generate condition list that will be used to do specific co-expression
  # for each condition studied in the analysis
  tab <- read.delim2(exps, row.names = 1)
  Conditions <- as.character(tab[1,])
  rm(tab)
  #    1c) Collect all genes with expression from the network
  tab <- read.delim2(exps, skip = 2, row.names = 1)
  logtab <- log(data.matrix(tab), 2)
  #????? need to log transform if table is already in log2 formats
  # Eliminate possible miRNA expression from the data
  rowLessmiR <- rownames(logtab)[rownames(logtab) %in% miRnames]
  if(length(rowLessmiR) > 0){
   logtab <- logtab[-rowLessmiR, ]
  }
  #       1d) Print expression table with ratio as needed
  #       for MIR tool to work properly
  nameexps <- extractname(exps)
  namenet <- extractname(net)
  myexpmiR <- paste0(outdir, "1_Expression", nameexps, "_", namenet, type,
                     method, fadjacency, method, perm, "nomiRNA.txt")
  # Print the file without miRNA expression
  # execute miRNA expression effects
  ######### END OFCODING SO FAR NEED TO FINISH REST OF THE CODE
  system(.)
}

rankingNEcorr <- function(){

}

cleaningNECorr <- function(){

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

