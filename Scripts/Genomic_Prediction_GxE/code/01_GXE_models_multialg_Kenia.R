#!/usr/bin/env/ Rscript
# SBATCH --reservation=glbrc
#SBATCH --array=1-350 # for 35 traits and 10 reps
#SBATCH --time=150:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#SBATCH --job-name=GxE_algs
#SBATCH --output=../logs/%x_%A_%a.out

system("module load GCC/9.3.0  OpenMPI/4.0.3  R/4.0.3")

#-------------------------------------------------------------------------------
# Multi-Environment and Multi-Trait models with either CV1 or CV2 train-test
# split scheme
# Arguments:
#   [1] X_file (str):    path to genetic marker feature file
#   [2] feat_file (str): "all" or path to file with list of markers to include
#   [3] test_file (str): file with list of instances to include in the test set
#   [4] prefix (str):    prefix to add to output file names
#   [5] save_dir (str):  path to output directory
#   [6] alg (str):       which algorithm to run (across, mxe, rnorm, or multi)
#   [7] reps (int):      number of times to repeat modeling
#   [8] which_cv (str):  cv1/cv2 scheme for training and testing split
#   [9] getG (str):      y/n to determine whether to calculate the genomic relationship matrix
# 
# Returns:
#   [1] Fitted model object as .RDS file for each rep
#   [2] Main marker coefficients file (all reps combined)
#   [3] Environment 1 marker coefficients file (all reps combined)
#   [4] Environment 2 marker coefficients file (all reps combined)
#   [5] Predicted test set values file for target environment (all reps combined)
#   [6] If which_cv = cv2, predicted test set values file for target environment (all reps combined)
#-------------------------------------------------------------------------------

rm(list = ls())
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(BGLR))

setwd('/mnt/home/seguraab/Shiu_Lab/Project')

# Debugging arguments
# X_file <- "Data/Peter_2018/ORFs_no_NA.csv"
# test_file <- "Data/Peter_2018/Test.txt"
# prefix <- "cno_cv1_TESTSCRIPT"
# save_dir <- "Scripts/Genomic_Prediction_GxE"
# reps <- 1
# nIter <- 100
# burnIn <- 10

#-----------------------------------
# Read in data
#-----------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 8) {
    stop("Need 9 arguments: X_file, feat_file, test_file, prefix, save_dir, alg, reps, which_cv, getG",
         call.=FALSE)
} else {
    X_file <- args[1] # genetic marker feature file
    feat_file <- args[2] # "all" or file with list of features to include
    test_file <- args[3] # file with list of instances to include in the test set
    prefix <- args[4] # prefix to add to output file names
    save_dir <- args[5] # path to output directory
    alg <- args[6] # which algorithm to run (across, mxe, rnorm, or multi)
    reps <- as.numeric(args[7]) # number of times to repeat modeling
    which_cv <- args[8] # cv1/cv2 scheme for training and testing split
    getG <- args[9] # y/n to determine whether to calculate the genomic relationship matrix
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
Test <- scan(test_file, what="character") # list of instances in test set
head(Y)

# Make sure Y and X have the same row order
X <- X[rownames(Y),]

#-----------------------------------
# Set output directory to save files
#-----------------------------------
setwd(save_dir)

#-----------------------------------
# Set-up job array variables
#-----------------------------------
job <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1")) # job ID number
JOBS <- expand.grid(rep=1:reps, trait=1:ncol(Y)) # total number of jobs is ncol(Y)*reps

trait <- as.vector(JOBS[job,"trait"]) # target trait index in Y
rep <- as.vector(JOBS[job,"rep"]) # repetition number

#-----------------------------------
# Get subsets of correlated traits
#-----------------------------------
# for each trait, get indices of significantly correlated traits in Y_cor
trait_name <- colnames(Y[trait])
print(sprintf("Get top four correlated traits for %s", trait_name))
Y_train <- Y[!(rownames(Y) %in% Test),]
Y_cor <- psych::corr.test(Y_train, method="pearson", adjust="holm") # pearson correlation
upper_tri <- reshape2::melt(upper.tri(Y_cor$r))
Y_cor <- merge(reshape2::melt(Y_cor$r), reshape2::melt(Y_cor$p), by=c("Var1", "Var2"))
Y_cor <- Y_cor[upper_tri$value,]
colnames(Y_cor) <- c("Env1", "Env2", "r", "p.value")
index <- which(Y_cor$p.value < 0.05 & Y_cor$Env1==trait_name | Y_cor$Env2==trait_name)
# write.csv(Y_cor, "Scripts/Data_Vis/Section_2/pheno_training_cor.csv", quote=F, index=F)

if (length(index) >= 4){
    y_cor <- Y_cor[index,] # subset correlated traits
    y_cor <- y_cor[order(y_cor$r, decreasing=T),] # sort
    y_cor <- unique(c(as.character(y_cor[1:4,"Env1"]), as.character(y_cor[1:4,"Env2"]))) # correlated traits
    y_cor <- y_cor[which(y_cor!=trait_name)] # names of correlated traits only
    y <- Y[,c(trait_name, y_cor)]
    n <- 5 # number of traits
} else {
    print(sprintf("%s doesn't have sufficient correlated traits", trait_name))
    y <- Y[,trait_name]
    n <- 1
}

#-----------------------------------
# Gene-by-Environment Models
#-----------------------------------
# scale X and calculate the genomic relationship matrix
if (getG=="y"){
    print("Standardizing marker matrix...")
    # X <- scale(X)/sqrt(ncol(X))
    X <- scale(X)/ncol(X)
    X[1:5,1:5]
}

if (which_cv=="cv1") y[Test,] <- NA # mask the test set values for all the traits
if (which_cv=="cv2") y[Test, trait_name] <- NA # mask the test set values for the target trait

# Run the models
set.seed(rep)
nIter=32000
burnIn=3200
start <- Sys.time() # start time
yNA <- as.vector(as.matrix(y)) # stack trait labels into one vector
yn <- do.call(rbind, replicate(n, as.matrix(y), simplify=F)) # to get rownames
names(yNA) <- rownames(yn)
if (alg=="across"){ # Across environments model
    Xn <- do.call(rbind, replicate(n, X, simplify=F)) # stack input matrix
    print(sum(rownames(Xn)==rownames(yNA))==nrow(Xn)) # sanity check
    envID <- factor(rep(colnames(y), each=nrow(y), levels=colnames(y)))
    ETA <- list(list(~envID-1,model="FIXED"), list(X=Xn, model='BRR'))
    fm <- BGLR(y=yNA, ETA=ETA, nIter=nIter, burnIn=burnIn)
    print("Saving model...")
    saveRDS(fm, file=paste(prefix, trait_name,
            paste(names(y)[which(names(y)!=trait_name)], collapse="_"),
            "rep", as.character(rep), "across_model.RDS", sep="_"))
    yHats <- matrix(fm$yHat, ncol=ncol(y))
    colnames(yHats) <- colnames(y)
    rownames(yHats) <- rownames(y)
    test_IDs <- matrix(fm$whichNa, ncol=ncol(y))
    yHats_test <- yHats[test_IDs[,1],] # predicted test set label
} else if (alg=="mxe"){ # Marker-by-environment model
    X0 <- matrix(nrow=nrow(X), ncol=ncol(X), 0) # a matrix full of zeros
    Xn <- do.call(rbind, replicate(n, X, simplify=F)) # main effects
    X_1 <- rbind(X, X0, X0, X0, X0) # environment-specific effects
    X_2 <- rbind(X0, X, X0, X0, X0)
    X_3 <- rbind(X0, X0, X, X0, X0)
    X_4 <- rbind(X0, X0, X0, X, X0)
    X_5 <- rbind(X0, X0, X0, X0, X)
    ETA <- list(main=list(X=Xn, model="BRR"), int1=list(X=X_1, model="BRR"),
                int2=list(X=X_2, model="BRR"), int3=list(X=X_3, model="BRR"),
                int4=list(X=X_4, model="BRR"), int5=list(X=X_5, model="BRR"))
    fm2 <- BGLR(y=yNA, ETA=ETA, nIter=nIter, burnIn=burnIn)
    print("Saving model...")
    saveRDS(fm2, file=paste(prefix, trait_name,
            paste(names(y)[which(names(y)!=trait_name)], collapse="_"),
            "rep", as.character(rep), "mxe_model.RDS", sep="_"))
    yHats <- matrix(fm2$yHat, ncol=ncol(y))
    colnames(yHats) <- colnames(y)
    rownames(yHats) <- rownames(y)
    test_IDs <- matrix(fm2$whichNa, ncol=ncol(y))
    yHats_test <- yHats[test_IDs[,1],] # predicted test set label
} else if (alg=="multi"){ # Multi-trait model
    ETA <- list(list(X=X,model="BRR"))
    fm3 <- Multitrait(y=as.matrix(y), ETA=ETA, nIter=nIter, burnIn=burnIn)
    print("Saving model...")
    saveRDS(fm3, file=paste(prefix, trait_name,
            paste(names(y)[which(names(y)!=trait_name)], collapse="_"),
            "rep", as.character(rep), "multitrait_model.RDS", sep="_"))
    yHats_test <- fm3$ETAHat[fm3$missing_records,] # predicted test set label
}
end <- Sys.time() # end time
print(sprintf("Model elapsed time: %f", end-start))

#-----------------------------------
# Write results to files
#-----------------------------------
# Extract coefficients from model and write to file
print(paste("Saving marker coefficients and predicted values files to", save_dir))
if (alg=="across"){
    file <- paste(prefix, trait_name, 
                  paste(names(y)[which(names(y)!=trait_name)], collapse="_"),
                  "coef_across.txt", sep="_") # marker effects same across all envs
    if (!file.exists(file)) {
        print(paste("Create marker coefficients file", file))
        cat(paste("Rep", paste(names(fm$ETA[[2]]$b), collapse="\t"), sep="\t"),
            file=file, append=F, sep="\n")
    } 
    print(paste("Writing marker coefficients to file for rep", rep))
    cat(paste(rep, paste(fm$ETA[[2]]$b, collapse="\t"), sep="\t"),
            file=file, append=T, sep="\n")
} else if (alg=="mxe"){
    file <- paste(prefix, trait_name, 
                  paste(names(y)[which(names(y)!=trait_name)], collapse="_"),
                  "coef_mxe.txt", sep="_") # main and specific env marker effects
    if (!file.exists(file)) {
        print(paste("Create marker coefficients file", file))
        cat(paste("Rep", "Type", paste(names(fm2$ETA$main$b), collapse="\t"),
            sep="\t"), file=file, append=F, sep="\n")
    }
    print(paste("Writing marker coefficients to file for rep", rep))
    message(paste(file,
        "exists. Will append main and interaction coefficients for rep", rep))
    cat(paste(rep, "main", paste(fm2$ETA$main$b, collapse="\t"), sep="\t"),
            file=file, append=T, sep="\n")
    cat(paste(rep, colnames(y)[1], paste(fm2$ETA$int1$b, collapse="\t"), sep="\t"),
            file=file, append=T, sep="\n")
    cat(paste(rep, colnames(y)[2], paste(fm2$ETA$int2$b, collapse="\t"), sep="\t"),
            file=file, append=T, sep="\n")
    cat(paste(rep, colnames(y)[3], paste(fm2$ETA$int3$b, collapse="\t"), sep="\t"),
            file=file, append=T, sep="\n")
    cat(paste(rep, colnames(y)[4], paste(fm2$ETA$int4$b, collapse="\t"), sep="\t"),
            file=file, append=T, sep="\n")
    cat(paste(rep, colnames(y)[5], paste(fm2$ETA$int5$b, collapse="\t"), sep="\t"),
            file=file, append=T, sep="\n")
} else if (alg=="multi"){
    file <- paste(prefix, trait_name, 
                  paste(names(y)[which(names(y)!=trait_name)], collapse="_"),
                  "coef_multitrait.txt", sep="_") # main marker coefs
    if (!file.exists(file)) {
        print(paste("Create marker coefficients file", file))
        cat(paste("Rep", paste(names(fm3$ETA[[1]]$beta), collapse="\t"), sep="\t"),
            file=file, append=F, sep="\n")
    }
    print(paste("Writing marker coefficients to file for rep", rep))
    cat(paste(rep, paste(fm3$ETA[[1]]$beta, collapse="\t"), sep="\t"),
            file=file, append=T, sep="\n")
}

# Write predicted values to file
file2 <- gsub("coef", "preds_test", file)
if (!file.exists(file2)) {
    print(paste("Create predicted values file", file))
    cat(paste("Rep", "Trait", paste(rownames(yHats_test), collapse="\t"),
        sep="\t"), file=file2, append=F, sep="\n")
} 
print(paste("Writing predicted values to file for rep", rep))
for (col in colnames(yHats_test)){
    cat(paste(rep, col, paste(yHats_test[,col], collapse="\t"),
        sep="\t"), file=file2, append=T, sep="\n")
}

