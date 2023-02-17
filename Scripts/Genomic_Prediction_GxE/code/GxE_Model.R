# Multi-Environment Models
# Y = Xb + Zu
# where Y_ij is the fitness value of the ith genotype 
# grown in the jth environment. mu is the grand mean.
# G_i is the random effect of the ith genotype. E_j is 
# the fixed effect of the jth environment. GE_ij is a 
# two-way interaction term. e_ij is the error term.

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(sommer))

r2_score <- function(preds, actual) {
	# This function is comparable to sklearn's r2_score function.
	# Computes the coefficient of determination (R^2)
	rss <- sum((preds - actual) ^ 2) # residual sum of squares
	tss <- sum((actual - mean(actual)) ^ 2) # total sum of squares
	return(1 - (rss/tss)) # return R^2 value
}

### Read in data ###
# Test instances
Test <- scan('/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt',
            what='character')

# Genotype data
X <- fread('/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv')
X <- as.matrix(X, rownames=1, colnames=1)

# Phenotype data
y <- read.csv('/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv')
y <- y[,c(1,3,7,19,35)] # to test out functions
y_melt <- melt(y, id.vars='Unnamed..0') # pivot longer
colnames(y_melt) <- c('ID', 'Env', 'y')

# train-test split
`%notin%` <- Negate(`%in%`)
X_train <- X[rownames(X) %notin% Test,]
X_test <- X[Test,]
y_train <- y_melt[y_melt$ID %notin% Test,] 
y_test <- y_melt[y_melt$ID %in% Test,]
# y_train <- y[y$Unnamed..0 %notin% Test,]
# y_test <- y[y$Unnamed..0 %in% Test,]

########### Calculate the Genomic Relationship Matrix (VanRaden) ###########
# Method 1: tcrossprod(X)/# of markers
# G <- tcrossprod(X)/ncol(X) 
#G_train <- tcrossprod(X_train)/ncol(X_train) 

# Method 2: tcrossprod(centered X)/(2*sum(pq))
# phat <- colMeans(X_train)/2 # minor allele frequency
# X_train2 <- scale(X_train, center=TRUE, scale=FALSE)
# k <- 2*sum(phat*(1-phat)) # denominator
# G_train2 <- tcrossprod(X_train2)/k

# # Method 3: tcrossprod(Z)/# of markers
# Z_train <- scale(X_train, center=TRUE, scale=TRUE)
# G_train3 <- tcrossprod(Z_train)/ncol(Z_train)

# Endelman Method
#G <- A.mat(X)
G_train4 <- A.mat(X_train)
# G_test4 <- A.mat(X_test)

########### Main effect model: assume GxE doesn't exist ###########
# Main genotype random effects + env fixed effects
fit <- mmer(fixed=y~Env, random=~vsr(ID, Gu=G_train4),
            rcov=~units, data=y_train, verbose=F)
summary(fit)
saveRDS(fit, file="GxE_MET_4envs.rds")

f <- fitted.mmer(fit)
y_pred <- f$dataWithFitted$y.fitted # predicted training set values

r2_score(y_pred[1:625], y_train$y[1:625]) # YPDCAFEIN40
r2_score(y_pred[626:1250], y_train$y[626:1250]) # YPDCAFEIN50
r2_score(y_pred[1251:1875], y_train$y[1251:1875]) # YPDCUSO410MM
r2_score(y_pred[1876:2500], y_train$y[1876:2500]) # YPDBENOMYL500

# Evaluate on test set

#p <- predict.mmer(fit, classify="Env")


########### Unstructured model: GxE exists ###########
# Env-specific variances and covariances for each env pair exist
# y_ie = unstructured cov struc of Envs

y <- read.csv('/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv')
y <- y[,c(1,3,35)] # to test out functions
y_melt <- melt(y, id.vars='Unnamed..0') # pivot longer
colnames(y_melt) <- c('ID', 'Env', 'y')
y_train <- y_melt[y_melt$ID %notin% Test,] 
y_test <- y_melt[y_melt$ID %in% Test,]

fitU <- mmer(y~Env, random=~ vsr(usr(Env), ID, Gu=G_train4),
             rcov=~units, data=y_train, verbose=F) #  tolParInv=0.01, for caf 40 & 50
summary(fitU)
#saveRDS(fitU, file="GxE_UMET_4envs.rds")
#saveRDS(fitU, file="GxE_UMET_3envs.rds") # no CuSO4
#saveRDS(fitU, file="GxE_UMET_2envs.rds") # no CuSO4 or Caffeine50
#saveRDS(fitU, file="GxE_UMET_2envsb.rds") # added vsr(ID, Gu=G_train4) +... term

fU <- fitted.mmer(fitU)
y_predU <- fU$dataWithFitted$y.fitted # predicted training set values

r2_score(y_predU[1:625], y_train$y[1:625]) # YPDCAFEIN40
#r2_score(y_predU[751:1500], y$YPDCAFEIN50)
#r2_score(y_predU[1501:2250], y$YPDCUS0410MM)
r2_score(y_predU[626:1250], y_train$y[626:1250]) # YPDBENOMYL500

# Evaluate on test set


#pU <- predict.mmer(fitU, classify="Env")