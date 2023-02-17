# Marco's code to perform LD pruning (In Gustavo's lab)
# Gustavo said this code was under development

rm(list=ls())
setwd("/mnt/research/quantgen/projects/G2F/tools/LD_prune")
source("pruning_functions.R")
load("sample_data.RData")

dim(geno)
geno[1:10,1:6]
table(MAP$chr)

thr <- 0.90

#================================================
# EXAMPLE 1
# Pruning using a maximum R2
out1 <- LD_prune(geno, MAP, threshold=thr, mc.cores=1)

# Check chromosome-wise (e.g., chr 10)
indexCHR <- which(MAP$chr=='10')
X <- geno[,indexCHR]
R2 <- cor(X)^2

# Before pruning (some SNPs are 'linked' with other SNPs)
A <- (R2>thr);  diag(A) <- FALSE
range(colSums(A))

# After pruning (no SNPs are connected with others)
tmp <- colnames(R2) %in% out1$pruneIn[[1]]$NAME
A <- (R2[tmp,tmp]>thr); diag(A) <- FALSE
range(colSums(A))

#================================================
# EXAMPLE 2.
# Pruning using a a maximum R2 and a maximum distance (e.g., 100 Kb)
dMax <- 1E5
out2 <- LD_prune(geno, MAP, threshold=thr, d.max=dMax, mc.cores=1)

D <- as.matrix(dist(MAP[indexCHR,"pos",drop=F], method='manhattan'))

# Check again (within chromosome)
tmp <- colnames(R2) %in% out2$pruneIn[[1]]$NAME
A <- (R2[tmp,tmp]>thr) & (D[tmp,tmp]<=dMax); diag(A) <- FALSE
range(colSums(A))
