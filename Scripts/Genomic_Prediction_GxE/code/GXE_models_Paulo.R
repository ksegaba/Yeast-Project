#!/usr/bin/env Rscript
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --constraint=intel16|intel18
#SBATCH --array=1-1000

library(tidyverse)
library(BGLR)
library(SFSI)
library(corrplot)

rm(list = ls())

setwd('/mnt/gs21/scratch/izquier7/yeast/Data_for_Paulo')

Y <- read_csv('pheno.csv')
#geno <- fread('../geno.csv')

##### Add ID col as rownames
#M <- geno %>% remove_rownames() %>% column_to_rownames(var="ID")
#M <- as.matrix(M)

##### Genomic relationship
#G <- tcrossprod(scale(M))/ncol(M) # additice relationship
#D<-as.matrix(dist(M,method="euclidean"))^2 # euclidian distance
#D<-D/mean(D) # Gaussian kernel

# write.csv(G, "../G.csv")
# write.csv(D, "../D.csv")

G <- read_csv("G.csv")
G <- G %>% column_to_rownames(var = "...1")
D <- read_csv("D.csv")
D <- D %>% column_to_rownames(var = "...1")

G <- as.matrix(G)
D <- as.matrix(D)

Y <-  column_to_rownames(Y, var="ID")
cor_phe <- cor(Y)

pvalues_cor <- cor.mtest(Y)
pvalues_cor <- pvalues_cor$p # p.values - pearson

##### There are some strong correlation, but also negative correlation across enviroments. 

##* 1. Try using all the env (maybe not to promising)
##* 2. Using only correlated environments (> 0.4?)

job <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
JOBS <- expand.grid(rep=1:50, trait=1:20)

trait <- as.vector(JOBS[job,"trait"])
rep <- as.vector(JOBS[job,"rep"])

trait_name <- colnames(Y[trait])
# index <- which(pvalues_cor[,trait]< 0.05)

index <- which(cor_phe[,trait]>= 0.4)
Ycor <- Y[,index]


#_______________________________
# Training and testing sets 

set.seed(rep)

YNA <- Ycor

n <- nrow(Ycor)
percTST <- 0.3
nTST <- round(percTST*n)
indexTST <- sample(1:n,size=nTST,replace=FALSE)
indexTRN <- seq_along(1:n)[-indexTST]
YNA[indexTST,trait_name] <- NA


#---------------------------------------
# Design matrix ana Interaction termns
#---------------------------------------

# Design matrix for individuals. It connects individuals with environments
GID <- factor(rep(rownames(YNA),ncol(YNA)),
              levels=rownames(YNA))

Zg <- model.matrix(~GID-1) # leave out the intercept

# Design matrix for environments. Used in the multi-environment R-Norm model
envID <- factor(rep(colnames(YNA),
                    each=nrow(YNA)),
                levels=colnames(YNA))
ZE <- model.matrix(~envID-1) 

#  Covariance structure for effects
ZgGZgt <- Zg%*%G%*%t(Zg)    # Genetic effect  
ZEZEt <- tcrossprod(ZE)     # Environmental effect
GE <- ZgGZgt*ZEZEt          # GxE interaction term (R-Norm model)

# Eigen decomposition (to speed computational time)
eigen_G <- eigen(G)
eigen_G0 <- eigen(ZgGZgt)
eigen_GE <- eigen(GE)

# Interaction terms (MxE model)
MxE_eigen <- vector("list",ncol(YNA))

for(env in 1:ncol(YNA)){ 
  tmp <- rep(0,ncol(YNA)) ; 
  tmp[env] <- 1; 
  G1 <- kronecker(diag(tmp),G)
  MxE_eigen[[env]] <- eigen(G1)
}

#---------------------------------------
# Single environment
#---------------------------------------
nIter=10000
burnIn=1000

#### GBLUP


ETA <- list(G=list(V=eigen_G$vectors,
                   d=eigen_G$values,model='RKHS'))

fm_se <-BGLR(y=YNA[,trait_name],ETA=ETA,
          nIter=nIter,burnIn=burnIn,
          saveAt = paste0(trait_name,"_GBLUP",rep,"_"))

indexTST <- which(is.na(YNA[,trait_name]))
acc_GBLUP = cor(Y[indexTST,trait_name], fm_se$yHat[indexTST])


#### Gaussian kernel
h<- c(.02,1,5)
KList<-list()

for(i in 1:length(h)){
  KList[[i]]<-list(K=exp(-h[i]*D),model='RKHS')
}

ETA <- list(list(K=D,model="RKHS"))

fmKA <-BGLR(y=YNA[,trait_name],ETA=KList,
           nIter=nIter,burnIn=burnIn,
           saveAt = paste0(trait_name,"_KA",rep,"_"))

acc_KA = cor(Y[indexTST,trait_name], fmKA$yHat[indexTST])


#### SSI
fm_ssi = SSI.CV(YNA[,trait_name],K=G,trn=indexTRN,
             nCV=10, nfolds = 5)

lambda = summary(fm_ssi)$optCOR["lambda"]

# Fit the index with the obtained lambda
fm_ssi2 = SSI(Y[,trait_name],K=G,trn=indexTRN, 
          tst=indexTST,lambda=lambda)

acc_SSI = summary(fm_ssi2)$optCOR


#---------------------------------------
# Multi-environment models
#---------------------------------------
yNA <- as.vector(as.matrix(YNA))

# Across-environments model =============
ETA <- list(list(~envID-1,model="FIXED"))
ETA[[2]] <- list(V=eigen_G0$vectors,d=eigen_G0$values,
                 model='RKHS')

fm_acEnv <- BGLR(y=yNA,ETA=ETA,nIter=nIter,burnIn=burnIn,
                 saveAt = paste0(trait_name,"_accenv",rep,"_"))

acenv <- rep(NA,ncol(YNA))

YHat_acenv <- matrix(fm_acEnv$yHat,ncol=ncol(YNA))
colnames(YHat_acenv) <- colnames(YNA)

acc_acEnv <- cor(Y[indexTST,trait_name],
                 YHat_acenv[indexTST,trait_name])

#MxE interaction model

ETA <- list(list(~envID-1,model="FIXED"))
ETA[[2]] <- list(V=eigen_G0$vectors,d=eigen_G0$values,
                 model='RKHS')

# Interaction terms
for(env in 1:ncol(YNA)){
  eigen_G1 <- MxE_eigen[[env]]
  ETA[[(env+2)]] <- list(V=eigen_G1$vectors,
                         d=eigen_G1$values,model='RKHS')
}

# Model 
fm_me <- BGLR(y=yNA,ETA=ETA,nIter=nIter,burnIn=burnIn,
              saveAt = paste0(trait_name,"_me",rep,"_"))

YHat_me <- matrix(fm_me$yHat,ncol=ncol(YNA))
colnames(YHat_me) <- colnames(YNA)

acc_me <- cor(Y[indexTST,trait_name],
                      YHat_me[indexTST,trait_name])

# Reaction-Norm model
ETA <- list(list(~envID-1,model="FIXED"))
ETA[[2]] <- list(V=eigen_G0$vectors,
                 d=eigen_G0$values,model='RKHS')
ETA[[3]] <- list(V=eigen_GE$vectors,d=eigen_GE$values,
                 model="RKHS")

# Model 
fm_rnorm <- BGLR(y=yNA,ETA=ETA,nIter=nIter,burnIn=burnIn,
                 saveAt = paste0(trait_name,"_rn",rep,"_"))

YHat_rnorm <- matrix(fm_rnorm$yHat,ncol=ncol(YNA))

colnames(YHat_rnorm) <- colnames(YNA)

acc_rNorm <- cor(Y[indexTST,trait_name],
                 YHat_rnorm[indexTST,trait_name])

#---------------------------------------
# Multi-trait
#---------------------------------------


# GBLUP
ETA <- list(list(K=G,model="RKHS"))

YNA0 <- as.matrix(YNA)

fm_mt_gblup <- Multitrait(y=YNA0, 
                     ETA=ETA, 
                     resCov=list(type="DIAG"), 
                     nIter=nIter, burnIn=burnIn,
                     saveAt = paste0(trait_name,"_mt",rep,"_"))

acc_MT = cor(Y[indexTST,trait_name],
             fm_mt_gblup$ETAHat[indexTST,trait_name])


#---------------------------------------
# Results
#---------------------------------------
results <- c(acc_GBLUP, acc_KA, acc_SSI[c(1,3)],
                   acc_acEnv, acc_me, acc_rNorm, 
                   acc_MT)

names(results) <- c('GBLUP','KA','SSI', "df",
                   'Across-env', 'MxE', 
                   'R-Norm', "MT")

save(results,
     file=paste0("results_",trait_name,"_rep_",rep,".RData"))
