setwd("/Users/olson/src/warelab/NECorr/R/")

network.file <- "../data/network_DAPonly.txt"
expression <- "../data/expression.txt"
description.file <- "../data/description.csv"
metadata <- "../data/metadata.txt"
#network.file <- "../data/grnmetanet.txt"
#expression <- "../data/gene_expression_matrix.txt"
#description.file <-"../data/1.Ath.GeneDesc.csv"
#metadata <- "../data/grnmeta.metadata.txt"

source("./necorr_local.R")
Necorr(network.file=network.file, expression=expression, 
       description.file=description.file,
       condition="flower",metadata=metadata, name="test")

# Necorr(network.file=network.file, expression=expression,
#        description.file=description.file,
#        condition="cond1",metadata=metadata, name="test")
#
