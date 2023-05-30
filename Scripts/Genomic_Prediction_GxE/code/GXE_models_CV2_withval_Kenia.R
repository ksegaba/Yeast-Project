#!/usr/bin/env Rscript
#SBATCH --array=1-35
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --constraint=intel16|intel18
#SBATCH --output=../logs/%x_%A_%a
#-------------------------------------------------------------------------------
# Genotype-by-environment models
# Arguments:
#   [1]
#   [2]
#   [3]
#   [4] y/n to run models with genomic relationship matrix
# 

# Returns:
    model .RDS
    _ETA_int1_varB.dat
    _ETA_int2_varB.dat
    _mu.dat
    _varE.dat

# Written by Paulo Izquierdo
# Modified by Kenia Segura Abá
#-------------------------------------------------------------------------------

rm(list = ls())
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(BGLR))
suppressPackageStartupMessages(library(SFSI))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(Rfast))

setwd('/mnt/home/seguraab/Shiu_Lab/Project')


# Evaluation metrics
mse <- function(preds, actual){ return(mean((actual-preds)^2)) } # mean squared error
se <- function(vector){ return(sd(vector)/length(vector)) } # standard error
r2_score <- function(preds, actual) {
	rss <- sum((preds - actual) ^ 2) # residual sum of squares
	tss <- sum((actual - mean(actual)) ^ 2) # total sum of squares
	return(1 - (rss/tss))
} # same as sklearn's r2_score function to compute coefficient of determination (R^2)

`%notin%` <- Negate(`%in%`)

# Create file to save model results to
setwd(save_dir)  # set output directory as working directory
file <- "RESULTS_GxE.txt"
if (!file.exists(file)) {
    cat("Date", "RunTime", "Trait", "ID", "Alg", "NumInstances", "FeatureNum",
        "CVfold", "CV_rep", "MSE_val", "MSE_val_sd", "MSE_val_se", "r2_val",
        "r2_val_sd", "r2_val_se", "PCC_val", "PCC_val_sd","PCC_val_se",
        "MSE_test", "MSE_test_sd", "MSE_test_se", "r2_test", "r2_test_sd",
        "r2_test_se", "PCC_test", "PCC_test_sd", "PCC_test_se\n", file=file,
        append=FALSE, sep="\t")
} else {message("RESULTS_GxE.txt exists") }

#-------------------------------
# Read in data
#-------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 6) {
    stop("Need 4 arguments: X, feat, test, save", call.=FALSE)
} else {
    X_file <- args[1] # genetic marker feature table
    feat_file <- args[2] # file with list of features to include
    test_file <- args[3] # file with list of instances to include in the test set
    cvf_file <- args[4]
    fold <- args[5]
    reps <- args[6]
    calG <- args[7]
    prefix <- args[8] # prefix to add to output file names
    save_dir <- args[9]
}

print("Reading in data...")
if (feat_file != "all") {
    message("Pulling features to use...")
    FEAT <- scan(feat_file, what="character") # determine which features are included
    X <- fread(X_file, select=c("ID", FEAT), showProgress=TRUE) # subset feature table 
    X <- as.matrix(X, rownames=1, colnames=1)
    X[1:5,1:5]
} else {
    X <- fread(X_file, showProgress=TRUE)
    X <- as.matrix(X, rownames=1, colnames=1)
    X[1:5,1:5]
}

Y <- read.csv("Data/Peter_2018/pheno.csv", row.names=1) # phenotype matrix
Y_cor <- read.csv("Scripts/Data_Vis/Pheno_Figures/pheno_pairs_cor.csv") # phenotype matrix correlations
Test <- scan(test_file, what="character") # list of instances in test set
cvs <- read.csv(cvf_file, row.names=1) # cross-validation fold sample assignments

# Process cross-validation file
cvs_all <- merge(Y, cvs, by="row.names", all.x=TRUE) # merge Y_file and cvs_file
rownames(cvs_all) <- cvs_all$Row.names # set row names to sample name 
cvs_all <- cvs_all[,(dim(Y)[2]+2):ncol(cvs_all)] # keep only cvs columns
cvs_all[is.na(cvs_all)] = 0 # samples in test file will be "NA", set to 0
head(cvs_all)

# Make sure X and Y have the same order of rows as cvs_all
X <- X[rownames(cvs_all),]
Y <- Y[rownames(cvs_all),]

##### Genomic relationship
if calG == 'y'{
  G <- tcrossprod(scale(X))/ncol(X) # additive relationship
  D <- as.matrix(dist(X, method="euclidean"))^2 # euclidean distance
  D<-D/mean(D) # Gaussian kernel
  write.csv(G, paste(prefix, "_G.csv", sep=""))
  write.csv(D, paste(prefix, "_D.csv", sep=""))
}

#-------------------------------
# Get subsets of correlated traits
#-------------------------------
##### There are some strong correlations, but also negative correlation across enviroments. 
job <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
JOBS <- expand.grid(rep=1:20, trait=1:35)

trait <- as.vector(JOBS[job,"trait"])
rep <- as.vector(JOBS[job,"rep"])

trait_name <- colnames(Y[trait])

# for each trait, get indices of significantly correlated traits in Y_cor
index <- which(Y_cor$p.value < 0.05 & Y_cor$Env1==trait_name | Y_cor$Env2==trait_name)

# get the names of the top 4 correlated traits
trait_cor <- Y_cor[index,]
trait_cor <- trait_cor[order(trait_cor$PCC, decreasing=T),] # sort
correlated_traits <- c(unique(trait_cor$Env1[1:4]), unique(trait_cor$Env2[1:4]))
print(sprintf("%s is correlated with %s", trait_name,
      correlated_traits[which(correlated_traits!=trait_name)]))

#-------------------------------
# Training and testing sets 
#-------------------------------
YNA <- Y[correlated_traits] # subset the target trait and those it's correlated with
YNA[Test, trait_name] <- NA # mask only the test set values for the trait of interest

#---------------------------------------
# Gene-by-Environment Models
#---------------------------------------
set.seed(20230503)
nIter=32000
burnIn=3200

# scale and prepare the marker matrix
# X <- scale(X)/sqrt(ncol(X))
# X <- scale(X)/ncol(X)
X0 <- matrix(nrow=nrow(X),ncol=ncol(X),0) # a matrix full of zeros
X_main=rbind(X,X)
X_1=rbind(X,X0)
X_2=rbind(X0,X)

# Run models of all significant trait pairs with the target trait
for (idx in index) {
    y1 <- Y[Y_cor[idx,]$Env1]
    if (names(y1)==trait_name) y1[Test,] <- NA # mask test set if y1 is the target trait
    y2 <- Y[Y_cor[idx,]$Env2]
    if (names(y2)==trait_name) y2[Test,] <- NA
    y <- cbind(y1, y2) # env pair
    print(sprintf("Modeling trait pair %s, %s...", names(y)[1], names(y)[2]))

    # Cross-validation
    Coef_main <- c() # feature main coefficients
    Coef_env1 <- c() # feature environment 1 coefficients
    Coef_env2 <- c() # feature environment 2 coefficients
    pred_val <- c() # predicted value of validation
    pred_test <- c() # predicted value of test set
    Start <- Sys.time() # start time
    for (k in 1:reps){ # cross-validation repetitions
        print(sprintf("CV repetition number %i", k))
        rep <- cvs_all[,k] # column from cvs_all that specifies sample folds for this repetition
        Coeff_main <- c() # model main coefficients for this repetition
        Coeff_env1 <- c() # model env 1 coefficients for this repetition
        Coeff_env2 <- c() # model env 2 coefficients for this repetition
        yhat_val <- data.frame(Y[trait_name], yhat=0, row.names=row.names(Y)) # dataframe to hold predicted values
        yhat_test <- data.frame(Y[Test, trait_name], row.names=Test)
        colnames(yhat_test) <- colnames(Y[trait_name])
        for (j in 1:fold) { # cross-validion fold number
            print(sprintf("CV fold number %i", j))
            validation <- which(rep==j) # validation set for this fold
            training <- which(rep!=j & rep!=0) # training set is all other data excluding test set
            test <- which(rep==0) # testing set
            if (names(y)[1]==trait_name) y[validation, 1] <- NA # mask validation if y1 is the target trait
            if (names(y)[2]==trait_name) y[validation, 2] <- NA
            # model
            start <- Sys.time()
            ETA=list(main=list(X=X_main,model='BRR'),
                     int1=list(X=X_1, model='BRR'),
                     int2=list(X=X_2, model='BRR'), saveEffects=TRUE)
            fm <- BGLR(y=as.matrix(y), ETA=ETA, nIter=nIter, burnIn=burnIn, saveAt='GxE_')
            print("Saving model...")
            saveRDS(fm, file=paste("GxE_", trait_name, "_",
                    names(y)[which(names(y)!=trait_name)], "_rep_",
                    as.character(k), "_fold_", as.character(j), ".RDS", sep=""))
            end <- Sys.time() # end time
            print(sprintf("Model elapsed time: %f", end-start))
            
            # Extract results from model
            Coeff_main <- rbind(Coeff_main, fm$ETA$main$b) # feature main coefficients
            Coeff_env1 <- rbind(Coeff_env1, fm$ETA$int1$b) # feature env 1 coefficients
            Coeff_env2 <- rbind(Coeff_env2, fm$ETA$int2$b) # feature env 2 coefficients
            yhat_val$yhat[validation] <- fm$yHat[1:nrow(X)][validation] # predicted labels for validation set
            yhat_test[paste("rep",j,sep="_")] <- fm$yHat[1:nrow(X)][test] # collect predicted labels
        }
        Coef <- rbind(Coef, colMeans(Coeff)) # mean feature coefficients
        pred_val <- cbind(pred_val, yhat_val$yhat[which(yhat_val$yhat!=0)]) # predicted values of validation set
        pred_test <- cbind(pred_test, rowMeans(yhat_test[2:fold+1])) # predicted values of test set
    }
    # save average feature coefficients
    write.csv(Coef, paste("Coef_", save_name, "_", trait_name,
                          names(y)[which(names(y)!=trait_name)],
                          ".csv", sep=""), row.names=F, quote=F)
    
    # save predicted values
    write.csv(pred_val, paste("Predict_value_cv_", save_name, "_", trait_name,
                              names(y)[which(names(y)!=trait_name)],
                              ".csv", sep=""), row.names=F, quote=F)
    write.csv(pred_test, paste("Predict_value_test_", save_name, "_", trait_name,
                               names(y)[which(names(y)!=trait_name)],
                               ".csv", sep=""), row.names=F, quote=F)
    
    # save model performances
    print("Writing model results to file...")
    RunTime <- Sys.time()-Start
    ID <- paste(trait_name, names(y)[which(names(y)!=trait_name)], save_name, sep="_")
    NumInstances <- nrow(X)-length(Test)
    mse_val <- apply(pred_val, 2, mse, Y[which(rownames(Y) %notin% Test),i])
    r2_val <- apply(pred_val, 2, r2_score, Y[which(rownames(Y) %notin% Test),i])
    pcc_val <- apply(pred_val, 2, cor, Y[which(rownames(Y) %notin% Test),i])
    mse_test <- apply(pred_test, 2, mse, Y[which(rownames(Y) %in% Test),i])
    r2_test <- apply(pred_test, 2, r2_score, Y[which(rownames(Y) %in% Test),i])
    pcc_test <- apply(pred_test, 2, cor, Y[which(rownames(Y) %in% Test),i])
    cat(paste(Sys.time(), RunTime, names(Y)[i], ID, "Bayesian LASSO",
        NumInstances, ncol(X), fold, number, mean(mse_val), sd(mse_val), se(mse_val),
        mean(r2_val), sd(r2_val), se(r2_val), mean(pcc_val), sd(pcc_val),
        se(pcc_val), mean(mse_test), sd(mse_test), se(mse_test), mean(r2_test),
        sd(r2_test), se(r2_test), mean(pcc_test), sd(pcc_test), se(pcc_test),
        sep="\t"), file=file, append=T, sep="\n")
}