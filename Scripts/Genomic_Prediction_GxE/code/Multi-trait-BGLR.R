################################################################################
# Fit a multi-trait Bayesian kernel BLUP with a Gaussian response variable
# with effects of G + E + GxE in the predictor

################################################################################

library(data.table)
library(tidyr)
library(BGLR)
source('Kernel_Construction.R')

# Test instances
Test <- scan('/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt',
            what='character')

# Genotype data
X <- fread('/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv')
X <- as.matrix(X, rownames=1, colnames=1)

# Phenotype data
y <- read.csv('/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv')
y_melt <- melt(y, id.vars='Unnamed..0') # pivot longer
colnames(y_melt) <- c('ID', 'Env', 'y')

# train-test split
`%notin%` <- Negate(`%in%`)
X_train <- X[rownames(X) %notin% Test,]
X_test <- X[Test,]
y_train <- y_melt[y_melt$ID %notin% Test,] 
y_test <- y_melt[y_melt$ID %in% Test,]

# Genomic relationship matrix
G <- getG(X_train, center=T, scaleG=T, scale=T)

# Matrix design for markers
K.gauss <- Kernel_computation(X=X_train, name='Gaussian', degree=NULL, nL=NULL)

# Environment design matrix
Z_E <- model.matrix(~0+Env, data=y_train)
K_E <- Z_E %*% t(Z_E)

ETA <- list(list(model='FIXED', X=

y_NA <- data.frame(y_melt) # test set mask (don't know if this is helpful)
y_NA$y[y_NA$ID %in% Test] <- NA 