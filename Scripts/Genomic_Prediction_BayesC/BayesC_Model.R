##################################################################################################
# Description: BayesC Single Trait Model with Variable Selection 
# Not yet: within a cross-validation (cv) scheme
# Not yet: incorporate multi-trait capabilty
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
#
# Modified by: Kenia Segura Abá
##################################################################################################

# Load necessary packages
library(BGLR)
library(data.table)
library('parallel')

set.seed(42) # for reproducibility

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

# Arguments for debugging
# X_file <- "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv"
# Y_file <- "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv"
# test_file <- "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt"
# cvf_file <- "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/CVFs.csv"
# feat_file <- "all"
# trait <- "YPACETATE"
# fold <- 5
# number <- 10
# save_dir <- "/mnt/scratch/seguraab/yeast_project/yeast_BayesC_results"
# i <- 1; k <- 1; j <- 1

# Read in data
print("Reading in data...")
if (feat_file != "all") {
    print("Pulling features to use...")
    FEAT <- scan(feat_file, what="character") # determine which features are included
    X <- fread(X_file, select=c("ID", FEAT), showProgress=TRUE) # subset genotype data 
    X <- as.matrix(X, rownames=1, colnames=1)
    X[1:5,1:5]
} else {
    X <- fread(X_file, showProgress=TRUE)
    X <- as.matrix(X, rownames=1, colnames=1)
    X[1:5,1:5]
}

Y <- read.csv(Y_file, row.names=1)
Test <- scan(test_file, what="character")
cvs <- read.csv(cvf_file, row.names=1)

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

# 08/16/2022 Kenia: Added coefficient of determination (R^2) function
r2_score <- function(preds, actual) {
	# This function is comparable to sklearn's r2_score function
	# It computes the coefficient of determination (R^2)
	rss <- sum((preds - actual) ^ 2) # residual sum of squares
	tss <- sum((actual - mean(actual)) ^ 2) # total sum of squares
	return(1 - (rss/tss)) # return R^2 value
}

# BayesC model
setwd(save_dir)
start.time <- Sys.time() # start time
PCC_cv <- c() # CV Pearson Correlation Coefficient (PCC)
PCC_test <- c() # Test set PCC
R2_cv <- c() # CV performance R-squared
R2_test <- c() # Test set performance R-squared
Predict_validation <- c() # Predicted label for CV
Predict_test <- c() # Predicted label for test set
for (i in 1:length(Y)){
    print(sprintf("Modeling trait %s...", names(Y)[i]))
    corr_CV <- c() # PCC of cross-validation
    corr_test <- c() # PCC of test set
    Accuracy_CV <- c() # R-sq of cross-validation (performance)
    Accuracy_test <- c() # R-sq of test set (performance)
    Coef <- c() # feature coefficients
    pred_val <- c() # predicted value of validation
    pred_test <- c() # predicted value of test set
    for (k in 1:number) { # cross-validation repetitions
        print(sprintf("CV repetition number %i", k))
        tst <- cvs_all[,k] # column from cvs_all that specifies sample folds for this repetition
        Coeff <- c() # model coefficients for this repetition
        #~~~~~~~# Need to figure out what other model paramters I need to save
        y_test <- c() # predicted values of test set
        yhat <- data.frame(cbind(Y, yhat=0, yhat_sd = 0)) # dataframe to hold predicted values
        yhat$yhat <- as.numeric(yhat$yhat)
        yhat$yhat_sd <- as.numeric(yhat$yhat_sd)
        row.names(yhat) <- row.names(Y)
        for (j in 1:fold) { # cross-validion fold number
            print(sprintf("CV fold number %i", j))
            validation <- which(tst==j) # validation set for this fold
            training <- which(tst!=j & tst!=0) # training set is all other data excluding test set
            test <- which(tst==0) # testing set
            yNA <- Y[,i] # label (dependent variable)
            yNA[validation] <- NA # mask validation sample values
            yNA[test] <- NA # mask test sample values
            
            # Build BayesC model
            start <- Sys.time() # start time
            ETA <- list(list(X=X, model="BayesC", saveEffects=TRUE)) # regression function
            model <- BGLR(y=yNA, ETA=ETA, verbose=FALSE, nIter=12000, burnIn=2000, saveAt = paste("BayesC_", names(Y)[i], "_rep_", as.character(k), "_fold_", as.character(j), "_", sep="")) # about 12 minutes
            end <- Sys.time() # end time
            print("Saving model...")
            saveRDS(model, file=paste("BayesC_", names(Y)[i], "_rep_", as.character(k), "_fold_", as.character(j), ".RDS", sep=""))
            print(sprintf("Model elapsed time: %f", end-start))
            
            # Plot feature effects
            pdf(paste("BayesC_d_manhattan", names(Y)[i], "_rep_", as.character(k), "_fold_", as.character(j), ".pdf", sep=""))
            plot(model$ETA[[1]]$d, xlab="", ylab="")
            dev.off()
            pdf(paste("BayesC_b_manhattan", names(Y)[i], "_rep_", as.character(k), "_fold_", as.character(j), ".pdf", sep=""))
            plot(model$ETA[[1]]$b, xlab="", ylab="")
            dev.off()

            # Extract results from model
            Coeff <- rbind(Coeff, model$ETA[[1]]$d) # feature coefficients (what is b? the estimated coeffs, d is the posterior probs)
            yhat$yhat[validation] <- model$yHat[validation] # predicted labels for validation set
            yhat$yhat_sd[validation] <- model$SD.yHat[validation] # standard deviation of predicted labels
            yhat$yhat[test] <- model$yHat[test] # predicted labels for test set
            yhat$yhat_sd[test] <- model$SD.yHat[test] # standard deviation of predicted labels
            y_test <- cbind(y_test, yhat$yhat) # collect predicted labels
        }
        # Performance measures per repetition
        print("Calculating performance...")
        corr_cv <- cor(yhat[which(tst!=0), i], yhat$yhat[which(tst!=0)]) # PCC of validation
		corr_CV <- c(corr_CV, corr_cv)
        #Accuracy_CV <- c(Accuracy_CV, corr_cv^2) # R-sq of validation
        Accuracy_CV <- c(Accuracy_CV, r2_score(yhat[which(tst!=0), i], yhat$yhat[which(tst!=0)])) # 08/16/2022 Kenia: Added coefficient of determination (R^2) function
        
        y_test <- cbind(y_test, rowMeans(y_test)) # mean predicted values
        corr_Test <- cor(yhat[which(tst==0), i], y_test[which(tst==0), ncol(y_test)]) # PCC of test 
		corr_test <- c(corr_test, corr_Test)
        #Accuracy_test <- c(Accuracy_test, corr_Test^2) # R-sq of test
        Accuracy_test <- c(Accuracy_test, r2_score(yhat[which(tst==0), i], y_test[which(tst==0), ncol(y_test)])) # 08/16/2022 Kenia: Added coefficient of determination (R^2) function
        
        Coef <- rbind(Coef, colMeans(Coeff)) # mean feature coefficients
        pred_val <- cbind(pred_val, yhat$yhat[which(tst!=0)]) # predicted values of validation set
        pred_test <- cbind(pred_test, y_test[which(tst==0), ncol(y_test)]) # predicted values of test set
    }
    # Overall performance
    PCC_cv <- cbind(PCC_cv, corr_CV) # cross-validation PCC
	PCC_test <- cbind(PCC_test, corr_test) # test set PCC
    R2_cv <- cbind(R2_cv, Accuracy_CV) # cross-validation accuracy
    R2_test <- cbind(R2_test, Accuracy_test) # test set accuracy
    write.csv(Coef, paste("Posterior_", save_name, "_", names(Y)[i], ".csv", sep=""), row.names=FALSE, quote=FALSE) # coefficients
    colnames(pred_val) <- paste(names(Y)[i], "_", 1:number, sep="") # columns for predicted values for each repetition
    Predict_validation <- cbind(Predict_validation, pred_val) # validation predicted values
    colnames(pred_test) <- paste(names(Y)[i], "_", 1:number, sep="") # columns for predicted values for each repetition
    Predict_test <- cbind(Predict_test, pred_test) # test set predicted values
}
print("Complete.")
end.time <- Sys.time() # end time
print(sprintf("Total elapsed time: %f", end.time-start.time))

# Save results to files
print("Saving results files...")
colnames(PCC_cv) <- names(Y)
write.csv(PCC_cv,paste('PCC_cv_results_',save_name,'_',trait,'.csv',sep=''),row.names=FALSE,quote=FALSE)

colnames(PCC_test) <- names(Y)
write.csv(PCC_test,paste('PCC_test_results_',save_name,'_',trait,'.csv',sep=''),row.names=FALSE,quote=FALSE)

colnames(R2_cv) <- names(Y) #c(paste(names(Y), "_PCC", sep=""), paste(names(Y), "_R2", sep=""))
write.csv(R2_cv,paste("R2_cv_results_", save_name, "_", trait, ".csv", sep=""), row.names=FALSE, quote=FALSE)

colnames(R2_test) <- names(Y) #c(paste(names(Y), "_PCC", sep=""), paste(names(Y), "_R2", sep=""))
write.csv(R2_test,paste("R2_test_results_", save_name, "_", trait, ".csv", sep=""), row.names=FALSE, quote=FALSE)

rownames(Predict_validation) <- rownames(X)[which(tst!=0)]
write.csv(Predict_validation,paste("Predict_value_cv_", save_name, "_", trait, ".csv", sep=""), row.names=TRUE, quote=FALSE)

rownames(Predict_test) <- rownames(X)[test]
write.csv(Predict_test,paste("Predict_value_test_", save_name, "_", trait, ".csv", sep=""), row.names=TRUE, quote=FALSE)