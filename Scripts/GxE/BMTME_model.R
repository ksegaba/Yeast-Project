# Bayesian multiple-environment model for GxE anlysis

library(data.table)
library(BMTME)

### Read in data ###
# Test instances
Test <- scan('/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt',
            what='character')

# Genotype data
X <- fread('/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv')
X <- as.matrix(X, rownames=1, colnames=1)

# Phenotype data
y <- read.csv('/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv')
#y_melt <- melt(y, id.vars='Unnamed..0') # pivot longer
#colnames(y_melt) <- c('ID', 'Env', 'y')

# train-test split
`%notin%` <- Negate(`%in%`)
X_train <- X[rownames(X) %notin% Test,]
X_test <- X[Test,]
#y_train <- y_melt[y_melt$ID %notin% Test,] 
#y_test <- y_melt[y_melt$ID %in% Test,]
y_train <- y[y$Unnamed..0 %notin% Test,]
y_test <- y[y$Unnamed..0 %in% Test,]

########### Calculate the Genomic Relationship Matrix ###########
# Method 1: tcrossprod(X)/# of markers
# G <- tcrossprod(X)/ncol(X) 
G_train <- tcrossprod(X_train)/ncol(X_train) 

# Method 2: tcrossprod(centered X)/(2*sum(pq))
# phat <- colMeans(X_train)/2 # minor allele frequency
# X_train2 <- scale(X_train, center=TRUE, scale=FALSE)
# k <- 2*sum(phat*(1-phat)) # denominator
# G_train2 <- tcrossprod(X2)/k

# # Method 3: tcrossprod(Z)/# of markers
# Z_train <- scale(X_train, center=TRUE, scale=TRUE)
# G_train3 <- tcrossprod(Z_train)/ncol(Z_train)

########### Fit Bayesian Multi-Environment Model ###########
LG <- cholesky(G_train) #?
ZG <- model.matrix(~0 + as.factor(y_train$Unnamed..0)) #?
Z.G <- ZG %*% LG #?

fm <- BME(Y=as.matrix(y_train[-1]), # pheno data in each env 
    Z1=Z.G, # genetic effects
    nIter=10000, # iterations
    burnIn=5000,
    thin=2, 
    bs=625)