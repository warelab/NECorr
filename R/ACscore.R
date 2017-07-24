rm(list=ls())
setwd("/Users/cliseron/Documents/1_Repository/NECorr_server_v2/")


# 1 - ranked expression order of the miRNA targets for each ratio
# matrix multiply by number of ratio vector
# Eprime = > ( E_1 to E_n)
expr.ranking <- function(Etab, BSaff.tab){
  Esrt<-c()
  Bmir <- c()
  for(i in seq_len(ncol(Etab))){
    colI <- colnames(Etab)[i]
    VectI <- Etab[, i]
    # Ranked Expression
    OrderI <- rev(grr::order2(VectI))
    Esrt[[colI]] <- as.data.frame(Etab[ OrderI, i])
    # Restitute the names to the files
    rownames(Esrt[[colI]]) <- oldrow <- rownames(Etab)[OrderI]
    colnames(Esrt[[colI]]) <- colI
    # Order the miRNA affinity score
    # remove the genes with no affinity for any of the miRNA
    a <- matches(rownames(Esrt[[colI]]), rownames(BSaff.tab), all.y=F, all.x=F)
    Bmir[[colI]] <- BSaff.tab[a$y[order2(a$x)], ]
    e <- Esrt[[colI]]
    #Esrt[[colI]] <- as.data.frame(e[a$x[order2(a$x)], ]) #v0.100
    Esrt[[colI]] <- as.data.table(e[a$x[order2(a$x)], ])
    rownames(Esrt[[colI]]) <- oldrow[a$x[order2(a$x)]]
  }
  out <- list(expr=Esrt, aff=Bmir)
  return(out)
}


#2 - binding score matrix of miRNA targets
# and ranked the binding affinity in function of the expression order
# Bprime => (Bi_1 to Bi_n)

prescore <- function(e.prime, b.prime){
  # find Arg max for a function in R
  # function f
  e.prime <- rep(seq(e.prime), ncol(b.prime))
  product.e.b <- b.prime * rep(e.prime, ncol(b.prime))
  #fi <- cumsum(product.e.b)/sum(product.e.b)
  sum.e.b = apply(product.e.b, 2, sum)
  #sum.e.b.tab <- abs(sum.e.b, rep(nrow(product.e.b),length(sum.e.b)))
  sum.e.b.tab <- rep(sum.e.b, nrow(product.e.b))
  #cumulative sums
  fi <- cumsum(abs(product.e.b)) / sum.e.b.tab
  # function g
  gi <- cumsum(e.prime) / sum(e.prime)
  # pre score
  # ps(i) <- max(gi) - max(fi)
  ps.i <- apply(gi-fi, 2, which.max)
  return(ps.i)
}
#random prescore with the bprime and the other outcome from prescore
rand.ps <- function(pset, bprime, e){
  b.prime <- bprime[pset,]
  #print("check3")
  product.e.b <- b.prime * rep(e, ncol(b.prime))
  #fi <- cumsum(product.e.b)/sum(product.e.b)
  sum.e.b = apply(product.e.b, 2, sum)
  #print(head(product.e.b))
  #sum.e.b.tab <- abs(sum.e.b, rep(nrow(product.e.b),length(sum.e.b)))
  sum.e.b.tab <- rep(sum.e.b, nrow(product.e.b))
  #print(head(sum.e.b.tab))
  #cumulative sum
  fi <- cumsum(abs(product.e.b)) / sum.e.b.tab
  # function g
  gi <- cumsum(e) / sum(e)
  # pre score
  # ps(i) <- max(gi) - max(fi)
  ps.i <- apply(gi-fi, 2, which.max)
  return(ps.i)
}

ps.test <- function(b.prime, e, score,nperm = 5, cores = 2) {
  library(doParallel)
  library(permute)
  ## check x and group are of same length
  #stopifnot(all.equal(length(x), length(e)))
  ## number of observations
  #print("process")
  N <- nobs(b.prime)
  ## generate the required set of permutations
  pset <- shuffleSet(N, nset = nperm)
  condCores <- cores > 1 #take out condition from the loop
  if (condCores) {
    ## initiate a cluster
    cl <- makeCluster(cores)
    #on.exit(stopCluster(cl = cl))
    ## iterate over the set of permutations applying meanDif
    getDoParWorkers()
    D <- parApply(cl, pset, 1, rand.ps, bprime = b.prime, e = e)
    stopCluster(cl)
  } else {
    D <- apply(pset, 1, rand.ps, bprime = b.prime, e = e)
  }
  ## add on the observed mean difference
  #print(dim(D))
  p.Means <- rowMeans(D)  #### meanRows()
  p.sd <- apply(D, 1, sd)  ### sd(Rows)
  D <- cbind(score, D) ### concatenate with eahc row (merge fast)
  ## compute & return the p-value
  Ds <- apply(D, 1, function(x) sum(x[2:length(x)] >= x[1])) # how many >= to the observed diff? sum(D >= D[1])
  pval <- Ds / (nperm + 1)# what proportion of perms is this (the pval)?
  pval[pval == 0] = 1/nperm
  AC <- (score - p.Means) / p.sd
  return(list(AC, pval))
}
ACsctab <- function(Etab, BSaff.tab, nperm){
  ranked <- expr.ranking(Etab, BSaff.tab)
  Esrt <- ranked$expr
  Bmir <- ranked$aff
  out <- vector("list", ncol(Etab)) 
  for(i in seq_len(ncol(Etab))){
    e.prime.tab <- Esrt[[i]][,1]
    e.prime.i <- as.numeric(e.prime.tab)
    b.prime.i <- Bmir[[i]]
    ps.i <- prescore(e.prime.i, b.prime.i)
    # t-test like of the score
    res <- ps.test(b.prime.i, e.prime.i, ps.i,nperm = nperm, cores = 4)
    # permutation of Bprime number of time required
    # FIND efficient R permutation algorithm
    out[[i]] <- res
    # significance
  }
  return(res)
}


#---------------------------
# main

BSaff.tab <- read.table("1.Ath.affinity.gene.txt")
Etab <- read.table("Tem_Tiss_Roots.txt", skip = 2)
library(grr)
library(data.table)
####
# number of permutation
np <- 2
system.time(res <- ACsctab(Etab, BSaff.tab, np))
#
#library(microbenchmark)
#microbenchmark(res <- ACsctab(Etab, BSaff.tab, np))






