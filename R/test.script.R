setwd("/Users/olson/src/warelab/NECorr/R/")

<<<<<<< HEAD
network.file <- "../data/grnmetanet.txt"
expression <- "../data/gene_expression_matrix.txt"
description.file <-"../data/1.Ath.GeneDesc.csv"
metadata <- "../data/grnmeta.metadata.txt"

Necorr(network.file=network.file, expression=expression, 
       description.file=description.file,
       condition="flower",metadata=metadata, name="test2", NSockets = 4)

#########
network.file <- "../data/network.txt"
expression <- "../data/expression.txt"
description.file <- "../data/description.csv"
metadata <- "../data/metadata.txt"

Necorr(network.file=network.file, expression=expression, 
       description.file=description.file,
       condition="cond1",metadata=metadata, name="test", NSockets = 4)
=======
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
>>>>>>> f6eeba4daa66f6eebaeb83ccf1d199a29b063a51

# Necorr(network.file=network.file, expression=expression,
#        description.file=description.file,
#        condition="cond1",metadata=metadata, name="test")
#
