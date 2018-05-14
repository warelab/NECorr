setwd("/Users/cliseron/Documents/1_Repository/NECorr/R/")

network.file <- "../data/network.txt"
expression <- "../data/expression.txt"
description.file <- "../data/description.csv"
metadata <- "../data/metadata.txt"
#network.file <- "../data/grnmetanet.txt"
#expression <- "../data/gene_expression_matrix.txt"
#description.file <-"../data/1.Ath.GeneDesc.csv"
#metadata <- "../data/grnmeta.metadata.txt"

Necorr(network.file=network.file, expression=expression, 
       description.file=description.file,
       condition="meristem_young_leaves_4_week_old",metadata=metadata, name="test")

Necorr(network.file=network.file, expression=expression, 
       description.file=description.file,
       condition="cond1",metadata=metadata, name="test")

