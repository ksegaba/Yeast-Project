#!/usr/bin/env Rscript
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --constraint=intel16|intel18

rm(list = ls())
library(tidyverse)
library(BGLR)
library(SFSI)
library(corrplot)

Y <- read.csv('pheno.csv')
G <- read.csv("G.csv")
G <- G %>% column_to_rownames(var = "X")
G <- as.matrix(G)
Y <- Y[,-1]

h2 <- c()

for (i in 1:ncol(Y)) {
  fm <- fitBLUP(y=Y[,i], K = G)
  h2[i] <- fm$h2
}


names(h2) <- colnames(Y)

ETA <- list(list(K=G,model="RKHS"))

h2_bglr <- c()

for (i in 1:ncol(Y)) {
fm= BGLR(y=Y[,i],ETA=ETA,nIter=4000,burnIn=400)
varU=scan('ETA_1_varU.dat')
varE=scan('varE.dat')
h2_ma = varU/(varU+varE)
h2_bglr[i] <- mean(h2_ma)
}

names(h2_bglr) <- colnames(Y)

save(h2,
     file=paste0("h2_SFSI.RData"))
     
save(h2_bglr,
     file=paste0("h2_BGLR.RData"))

