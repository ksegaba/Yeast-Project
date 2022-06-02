# Genomic BLUP using BGLR 
# Note: This is to compare to the GBLUP using rrblup package

library(BGLR)
library(data.table)
library(parallel)

# Genotype data (-1, 0, 1 encoding)
X <- fread("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv") # genotype
X <- as.matrix(X, rownames=1, colnames=1)
X[1:5,1:5]
Xs <- scale(X)/sqrt(ncol(X)) # centered
Xs[1:5,1:5]

test <- read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt", header=F)
X_train <- X[!(rownames(X) %in% test$V1),]
Xs_train <- scale(X_train, center=TRUE, scale=FALSE) # training set
X_test <- X[test$V1,] # test set
Xs_test <- scale(X_test, center=TRUE, scale=FALSE)

# Genotype data (0, 1, 2 encoding)
W <- X 
W[W==1] <- 2
W[W==0] <- 1
W[W==-1] <- 0
W[1:5,1:5]
Ws <- scale(W)/sqrt(ncol(W))
W_train <- W[!(rownames(W) %in% test$V1),]
Ws_train <- scale(W_train, center=TRUE, scale=FALSE) # training set
W_test <- W[test$V1,] # test set
Ws_test <- scale(W_test, center=TRUE, scale=FALSE)

# Phenotype data
y <- read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv", row.names=1)
y_train <- y[!(rownames(y) %in% test$V1),]
y_test <- y[test$V1,]

# Collect results
file <- "RESULTS_BGLR_GBLUP.csv"
if (!file.exists(file)) {
    cat("Date", "Trait", "yHat_test_SD", "PCC_test_mean", "PCC_test_SD", "R2_test_mean", "R2_test_SD\n", file=file, append=FALSE, sep=",")
} else {message("RESULTS_BLUPs.csv exists") }

# Genomic BLUP
nIter=12000
burnIn=2000
ETA=list(list(X=Xs_train, model="BRR"))
gblup <- function(Y){
    print(Y)
    fit <- BGLR(y=y_train[,Y], ETA=ETA, nIter=nIter, 
                burnIn=burnIn, saveAt=paste(Y, "training_brr_", sep="_"), verbose=FALSE)
    varE=scan(paste(Y, "brr_varE.dat", sep="_"))
    varU=scan(paste(Y, "brr_ETA_mrk_varB.dat", sep="_"))
    h2_1=varU/(varU+varE) # snp heritability
    
    # Append results to file
    cat(Sys.time(), file=file, append=TRUE, sep="")
    cat(trait, ",", file=file, append=TRUE, sep="")
    cat(varE, ",", file=file, append=TRUE, sep="")
    cat(varU, ",", file=file, append=TRUE, sep="")
    cat(h2_1, ",", file=file, append=TRUE, sep="")
    cat(mean(r2), ",", file=file, append=TRUE, sep="")
    cat(r2_sd, ",", file=file, append=TRUE, sep="")
    
    list(noquote(Y), varE, varU, h2_1) # output

}
setwd("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_GBLUP/BGLR_GBLUP")
system.time(tmp <- mclapply(colnames(y), gblup, mc.cores=35))

# write output to file
out <- do.call(rbind, tmp)
colnames(out) <- c("Y", "varE", "varB", "h2")
write.csv(out, "GBLUP_brr_training_results.csv", quote=FALSE, row.names=FALSE)