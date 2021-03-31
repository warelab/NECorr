# clean working place and load required libraries and functions
#rm(list=ls())
#setwd("./")
#setwd("/Users/cliseron/files2share/-a-clm_Code_Repository/NECorr_validations")

#install libraries if not present

library(limma)
library(Biobase)
library(ROCR)
library(plyr)
#source("ML_validation_function.R")
ScalN = function(x) (x - min(x))/(max(x) - min(x))
####
####    A - creation of positive set and negative set for the ROC curve analysis
####
# load positive and negative set (generate a file that can be reused for purpose)
# list of gene in the network here we used the expressiion file for but this can be changed in a vector
Netgene.file <- "SteleNet.txt"
Netgene <- read.table(Netgene.file,header = T, row.names = 1)
#exps = exps[-1,]
Netgene <- data.matrix(Netgene)
# positive set for the condition
positive.table = "Ath_vascular_TAIR.txt"
d.positive.table <- read.table(positive.table,header = F, sep ="\t")
d.positive <- as.vector(d.positive.table$V1)
d.positive <- unique(d.positive)
# negative set: randomly choose genes that are not part of the positive set
dn <- Netgene[-(which(d.positive %in% rownames(Netgene))),]
dn <- dn[sample(nrow(dn), length(d.positive)), ]
d.negative <- rownames(dn)
# extra step to be sure that the positive and negative sets are unique
d.negative <- setdiff(d.negative,d.positive)
d.positive <- setdiff(d.positive,d.negative)
# generating class labels (1 for positive to condition and 0 for negative)
data.lab <- cbind(c(d.positive,d.negative),c(rep(1,length(d.positive)),rep(0,length(d.negative))))
rownames(data.lab) <- data.lab[,1]
data.lab <- data.lab[,-1]
####
####    B - Differential expression
####
# expression file from microarrays in this case
exprs<- as.matrix(read.table("Radial_expression.txt", header=TRUE,sep="\t",row.names=1))
# parse expression data to only have the gene that are in the network - load the genes that are in the networks
exprs <- exprs[rownames(Netgene),]
# DE (differential expression)
eset <- ExpressionSet(assayData=exprs)
design <- model.matrix(~ 0 + factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6)))
colnames(design) <- c("Stele", "Endodermis", "Cortex", "Epidermis","Columella", "QCenter") # Assigns column names.
fit <- lmFit(eset, design)
contrast.matrix <- makeContrasts(Stele - (Endodermis + Cortex + Epidermis + Columella + QCenter)/5, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix) # Computes estimated coefficients and standard errors for a given set of contrasts.
fit2 <- eBayes(fit2)
#write.table(topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=50000), file="limma_complete.xls", row.names=F, sep="\t")
# Exports complete limma statistics table for first comparison group ('coef=1') to tab delimited text file.
de <- topTable(fit2, coef=1, adjust="fdr", sort.by="P", number=length(rownames(Netgene)))
de$B <- ScalN(de$B)
####
####    C - Parameter (factor) use to define the phenotype in the model
####
# load parameter table: tissue specificity, co-expression, topology
tab.Param <- read.table("8_stele_root.tissue_GRNnPnov13_GCC_10000_gene_Std_param.txt",
                        header = T, sep ="\t",row.names = 1)
# merge data label (phenotype) and the parameter table
tab.Param <- merge(data.lab,tab.Param, by =0, all =F)
colnames(tab.Param)[2] <- "phenotype"
rownames(tab.Param) <- tab.Param$Row.names
tab.Param <- tab.Param[,-1]
# add tissue specific and Differential Expression
tab.ParamDE <- merge(tab.Param,de, by=0, all=F)
rownames(tab.ParamDE) <- tab.ParamDE$Row.names
tab.ParamDE <- tab.ParamDE[,-1]
tab.ParamDE <- tab.ParamDE[,c(1:8,14)]
# add the hierachy table calculated by the Gerstein lab script (Yale, New Haven)
# ROC AUC preparation for comparison tissue specificity and DE for estimation of the phenotype
tab.ParamDE <- na.omit(tab.ParamDE)
datam <- tab.ParamDE[,-which(colnames(tab.ParamDE) == "Transitivity")]
colnames(datam)[ncol(datam)] <- "DE"
library(klaR)
Kernel = "rectangular"
N = 148
model <- NaiveBayes(phenotype ~ . , data = datam, kernel = Kernel, n = N)

rm(eset,design,fit2,fit,tab.Param,datam,tab.ParamDE,contrast.matrix,d.positive.table)
rm(exprs, de, dn, Netgene, d.negative, d.positive, Netgene.file,positive.table,ScalN,data.lab)
save.image("/Users/cliseron/files2share/1_NECorr_server_v2/src/ML_model.RData")
