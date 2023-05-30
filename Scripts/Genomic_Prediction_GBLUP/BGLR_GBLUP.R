################################################################################
# Description: Genomic BLUP using BGLR with cross-validation
#
# Arguments:
#   [1] X_file:     genotype matrix
#   [2] Y_file:     phenotype matrix
#   [3] feat_file   selected features or "all" for all features in X_file
#   [4] test_file   file with samples in test set
#   [5] trait:      name of trait (column in phenotype matrix) or "all" for all traits in Y_file
#   [6] cvf_file:   cross-validation scheme file (specifies which sample is part of each cv fold)
#   [7] fold:       cross-validation fold number
#   [8] number:     number of cross-validation repetitions
#   [9] save_name:  output file save name
#   [10] save_dir:   directory to save output file
# 
# Output:
#   [1] BGLR object as .RDS file
#   [2] Feature coefficients file
#   [3] R-squared values for each cross-validation repetition and of test set
#   [4] Predicted values for each cross-validation repetition and of test set
################################################################################

library(BGLR)
library(data.table)
library(parallel)

set.seed(20230426)

# Read in arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 10) {
    stop("Need 10 arguments: X_file Y_file feat_file test_file trait cvf_file fold number save_name save_dir", call.=FALSE)
} else {
    X_file <- args[1]
    Y_file <- args[2]
    feat_file <- args[3]
    test_file <- args[4]
    trait <- args[5]
    cvf_file <- args[6]
    fold <- as.numeric(args[7])
    number <- as.numeric(args[8])
    save_name <- args[9]
    save_dir <- args[10]
}

# Read in genotype data (-1, 0, 1 encoding) and standardize
print("Reading in data...")
if (feat_file != "all") {
    print("Pulling features to use...")
    FEAT <- scan(feat_file, what="character") # determine which features are included
    X <- fread(X_file, select=c("ID", FEAT), showProgress=TRUE) # subset genotype data 
    X <- as.matrix(X, rownames=1, colnames=1)
    Xs <- scale(X)/sqrt(ncol(X)) # scale and center
    Xs[1:5,1:5]
} else {
    X <- fread(X_file, showProgress=TRUE)
    X <- as.matrix(X, rownames=1, colnames=1)
    Xs <- scale(X)/sqrt(ncol(X)) # scale and center
    Xs[1:5,1:5]
}

# Read in remainind datasets
Y <- read.csv(Y_file, row.names=1) # phenotype data
Test <- scan(test_file, what="character") # test instances
cvs <- read.csv(cvf_file, row.names=1) # cross-validation fold and repetition assignments

# Process cross-validation file
cvs_all <- merge(Y, cvs, by="row.names", all.x=TRUE) # merge Y_file and cvs_file
rownames(cvs_all) <- cvs_all$Row.names # set row names to sample name 
cvs_all <- cvs_all[,(dim(Y)[2]+2):ncol(cvs_all)] # keep only cvs columns
cvs_all[is.na(cvs_all)] = 0 # samples in test file will be "NA", set to 0

# Make sure X and Y have the same order of rows as cvs_all
X <- X[rownames(cvs_all),]
Y <- Y[rownames(cvs_all),]

# Trait or traits to be modelled
if (trait == "all") {
    print("Modeling all traits...")
} else {
    Y <- Y[trait]
}

# Collect results
setwd(save_dir) # set output directory as working directory
file <- "RESULTS_BGLR_GBLUP.txt"
if (!file.exists(file)) {
    cat("Date", "Trait", "ID", "Alg", "NumInstances", "FeatureNum",
        "MSE_val", "MSE_val_sd", "MSE_val_se", "r2_val", "r2_val_sd",
        "r2_val_se", "PCC_val", "PCC_val_sd","PCC_val_se", "MSE_test",
        "MSE_test_sd", "MSE_test_se", "r2_test", "r2_test_sd", "r2_test_se",
        "PCC_test", "PCC_test_sd", "PCC_test_se\n", file=file, append=FALSE,
        sep="\t")
} else {message("RESULTS_BGLR_BLUPs.txt exists") }

# Evaluation metrics
mse <- function(preds, actual){ return(mean((actual-preds)^2)) } # mean squared error
se <- function(vector){ return(sd(vector)/length(vector)) } # standard error
r2_score <- function(preds, actual) {
	rss <- sum((preds - actual) ^ 2) # residual sum of squares
	tss <- sum((actual - mean(actual)) ^ 2) # total sum of squares
	return(1 - (rss/tss))
} # same as sklearn's r2_score function to compute coefficient of determination (R^2)

# Genomic BLUP
for (i in 1:length(Y)) { # loop through selected trait(s)
    print(sprintf("Modeling trait %s...", names(Y)[i]))
    pred_val <- c() # Predicted label for validation sets
    pred_test <- c() # Predicted label for test set
    for (k in 1:number) {
        print(sprintf("CV repetition number %i", k))
        tst <- cvs_all[,k] # column from cvs_all that specifies sample folds for this repetition
        Coeff <- c()
        y_test <- c() # predicted values of test set
        yhat <- data.frame(Y[i], yhat=0, row.names=rownames(Y)) # dataframe to hold predicted values
        for (j in 1:fold){
            print(sprintf("CV fold number %i", j))
            validation <- which(tst==j) # validation set for this fold
            training <- which(tst!=j & tst!=0) # training set is all other data excluding test set
            test <- which(tst==0) # testing set
            yNA <- Yy[,1] # label (dependent variable)
            yNA[validation] <- NA # mask validation sample values
            yNA[test] <- NA # mask test sample values
            
            # Build GBLUP model
            ETA=list(list(X=X[training,], model="BRR", saveEffects=TRUE))
            fit <- BGLR(y=Y[training,i], ETA=ETA, nIter=12000,
                burnIn=2000, saveAt = paste("gBLUP_", names(Y)[i], "_rep_", as.character(k), "_fold_", as.character(j), "_", sep=""),
                verbose=FALSE)
            print("Saving model...")
            saveRDS(fit, file=paste("gBLUP_", y, "_rep_", as.character(k), "_fold_", as.character(j), ".RDS", sep=""))

            # Extract results from model
            Coeff <- rbind(Coeff, fit$ETA[[1]]$b) # feature coefficients
            yhat$yhat[validation] <- fit$yHat[validation] # predicted labels for validation set
            yhat$yhat[test] <- fit$yHat[test] # predicted labels for test set
            yhat$yhat_sd[test] <- model$SD.yHat[test] # standard deviation of predicted labels
            y_test <- cbind(y_test, yhat$yhat) # collect predicted labels
        }
    }
    
    
    # Append results to file
    cat(Sys.time(), "\t" file=file, append=TRUE, sep="")
    cat(trait, "\t", file=file, append=TRUE, sep="")
    cat(varE, "\t", file=file, append=TRUE, sep="")
    cat(varU, "\t", file=file, append=TRUE, sep="")
    cat(h2_1, "\t", file=file, append=TRUE, sep="")
    cat(mean(r2), "\t", file=file, append=TRUE, sep="")
    cat(r2_sd, "\t", file=file, append=TRUE, sep="")
    
    list(noquote(Y), varE, varU, h2_1) # output

}
setwd("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_GBLUP/BGLR_GBLUP")
system.time(tmp <- mclapply(colnames(Y), gblup, mc.cores=35))

# write output to file
out <- do.call(rbind, tmp)
colnames(out) <- c("Y", "varE", "varB", "h2")
write.csv(out, "GBLUP_brr_training_results.csv", quote=FALSE, row.names=FALSE)