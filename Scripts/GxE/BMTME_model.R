# Bayesian multiple-trait multiple-environment model for GxE anlysis

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
y_melt <- melt(y, id.vars='Unnamed..0') # pivot longer
colnames(y_melt) <- c('ID', 'Env', 'y')

# train-test split
`%notin%` <- Negate(`%in%`)
X_train <- X[rownames(X) %notin% Test,]
X_test <- X[Test,]
y_train <- y_melt[y_melt$ID %notin% Test,] 
y_test <- y_melt[y_melt$ID %in% Test,]


### Bayesian Multi-Environment Model ###
G <- X%*%t(X)# VanRaden Genomic Relationship Matrix

LG <- cholesky(X)

BME(Y=y_train, Z1=, nIter=10000, burnIn=1000,)