#!/usr/bin/env Rscript
#SBATCH --array=1-700 # for 35 traits and 20 reps
#SBATCH --time=15:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --constraint=[intel16|intel18]
#SBATCH --job-name=GxE
#SBATCH --output=../logs/%x_%A_%a.out

system("module load GCC/9.3.0  OpenMPI/4.0.3  R/4.0.3")

#-------------------------------------------------------------------------------
# Marker-by-environment models with either CV1 or CV2 train-test split scheme
# Arguments:
#   [1] X_file (str):    path to genetic marker feature file
#   [2] feat_file (str): "all" or path to file with list of markers to include
#   [3] test_file (str): file with list of instances to include in the test set
#   [4] prefix (str):    prefix to add to output file names
#   [5] save_dir (str):  path to output directory
#   [6] reps (int):      number of times to repeat modeling
#   [7] which_cv (str):  cv1/cv2 scheme for training and testing split
#   [8] scale (str):     y/n to determine whether to scale markers or not
# 
# Returns:
#   [1] Fitted model object as .RDS file for each rep
#   [2] Main marker coefficients file (all reps combined)
#   [3] Environment 1 marker coefficients file (all reps combined)
#   [4] Environment 2 marker coefficients file (all reps combined)
#   [5] Predicted test set values file for target environment (all reps combined)
#   [6] If which_cv = cv2, predicted test set values file for target environment (all reps combined)
# 
# Written by Kenia Segura Abá
#-------------------------------------------------------------------------------

rm(list = ls())
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(BGLR))

setwd('/mnt/home/seguraab/Shiu_Lab/Project')

#-----------------------------------
# Read in data
#-----------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 8) {
    stop("Need 8 arguments: X_file, feat_file, test_file, prefix, save_dir, reps, which_cv, scale", call.=FALSE)
} else {
    X_file <- args[1] # genetic marker feature file
    feat_file <- args[2] # "all" or file with list of features to include
    test_file <- args[3] # file with list of instances to include in the test set
    prefix <- args[4] # prefix to add to output file names
    save_dir <- args[5] # path to output directory
    reps <- as.numeric(args[6]) # number of times to repeat modeling
    which_cv <- args[7] # cv1/cv2 scheme for training and testing split
    scale <- args[8] # y/n to determine whether to scale markers or not
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
head(Y)
head(Y_cor)

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
print(sprintf("Modeling target trait %s, rep %s", trait_name, rep))
index <- which(Y_cor$p.value < 0.05 & Y_cor$Env1==trait_name | Y_cor$Env2==trait_name)

#-----------------------------------
# Gene-by-Environment Models
#-----------------------------------
set.seed(rep)

nIter=32000
burnIn=3200

# scale and prepare the marker matrix for marker-by-environment model
if (scale=="y"){
    print("Standardizing marker matrix...")
    # X <- scale(X)/sqrt(ncol(X))
    X <- scale(X)/ncol(X)
    X[1:5,1:5]
}

X0 <- matrix(nrow=nrow(X),ncol=ncol(X),0) # a matrix full of zeros
X_main=rbind(X,X)
X_1=rbind(X,X0)
X_2=rbind(X0,X)

# Run models of all significant trait pairs with the target trait
for (idx in index) {
    # Training and testing set split
    y1 <- Y[Y_cor[idx,]$Env1]
    y2 <- Y[Y_cor[idx,]$Env2]
    y <- cbind(y1, y2) # env pair
    print(sprintf("Modeling trait pair %s, %s...", names(y)[1], names(y)[2]))
    if (which_cv=="cv1") y[Test,] <- NA # mask the test set values for all the traits
    if (which_cv=="cv2") y[Test, trait_name] <- NA # mask the test set values for the target trait

    # Create files to save marker coefficients to
    file <- paste(prefix, trait_name, names(y)[which(names(y)!=trait_name)],
                "coef_main.txt", sep="_") # main marker coefs
    if (!file.exists(file)) {
        cat(paste("Rep", paste(colnames(X), collapse="\t"), sep="\t"),
            file=file, append=F, sep="\n")
    } else {message(paste(file, "exists. Will append main coefficients for rep",
                          rep))}
    
    file2 <- gsub("main", trait_name, file) # trait env 1 marker coefs
    if (!file.exists(file2)) {
        cat(paste("Rep", paste(colnames(X), collapse="\t"), sep="\t"),
            file=file2, append=F, sep="\n")
    } else {message(paste(file2, "exists. Will append", trait_name,
                          "coefficients for rep", rep))}
    
    file3 <- gsub("main", names(y)[which(names(y)!=trait_name)], file) # trait env 2 marker coefs
    if (!file.exists(file3)) {
        cat(paste("Rep", paste(colnames(X), collapse="\t"), sep="\t"),
            file=file3, append=F, sep="\n")
    } else {message(paste(file3, "exists. Will append",
                          names(y)[which(names(y)!=trait_name)],
                          "coefficients for rep", rep))}

    # Create files to save predicted test set values to
    file4 <- gsub("coef_main", paste("preds_test", trait_name, sep="_"), file)
    if (!file.exists(file4)) {
        cat(paste("Rep", paste(Test, collapse="\t"), sep="\t"),
            file=file4, append=F, sep="\n")
    } else {message(paste(file4, "exists. Will append predicted", trait_name,
                          "values for rep", rep))} # predictions of trait env 1

    if (which_cv=="cv1"){
        file5 <- gsub("coef_main", paste("preds_test",
                      names(y)[which(names(y)!=trait_name)], sep="_"), file)
        if (!file.exists(file5)) {
            cat(paste("Rep", paste(Test, collapse="\t"), sep="\t"),
                file=file5, append=F, sep="\n")
        } else {message(paste(file5, "exists. Will append predicted",
                              names(y)[which(names(y)!=trait_name)],
                              "values for rep", rep))} # predictions of trait env 2
    }

    # Run the models
    start <- Sys.time() # start time
    ETA=list(main=list(X=X_main,model='BRR'),
                int1=list(X=X_1, model='BRR'),
                int2=list(X=X_2, model='BRR'))
    fm <- BGLR(y=as.matrix(y), ETA=ETA, nIter=nIter, burnIn=burnIn)
    print("Saving model...")
    saveRDS(fm, file=paste(prefix, trait_name,
                        names(y)[which(names(y)!=trait_name)],
                        "rep", as.character(rep), "model.RDS", sep="_"))
    end <- Sys.time() # end time
    print(sprintf("Model elapsed time: %f", end-start))

    # Extract coefficients from model and write to file
    print("Saving marker coefficients to file...")
    cat(paste(rep, paste(fm$ETA$main$b, collapse="\t"), sep="\t"),
            file=file, append=T, sep="\n") # main marker coefficients
    if (names(y)[1]==trait_name){
        cat(paste(rep, paste(fm$ETA$int1$b, collapse="\t"), sep="\t"),
            file=file2, append=T, sep="\n") # trait env 1 marker coefficients
        cat(paste(rep, paste(fm$ETA$int2$b, collapse="\t"), sep="\t"),
            file=file3, append=T, sep="\n") # trait env 2 marker coefficients
    } else {
        cat(paste(rep, paste(fm$ETA$int1$b, collapse="\t"), sep="\t"),
            file=file3, append=T, sep="\n") # trait env 2 marker coefficients
        cat(paste(rep, paste(fm$ETA$int2$b, collapse="\t"), sep="\t"),
            file=file2, append=T, sep="\n") # trait env 1 marker coefficients
    }

    # Write predicted values to file
    if (names(y)[1]==trait_name){
        if (which_cv=="cv2"){
            cat(paste(rep, paste(fm$yHat[fm$whichNa], collapse="\t"), sep="\t"),
                file=file4, append=T, sep="\n")
        }
        if (which_cv=="cv1"){
            cat(paste(rep, paste(fm$yHat[fm$whichNa][1:length(Test)], collapse="\t"),
                      sep="\t"), file=file4, append=T, sep="\n")
            cat(paste(rep, paste(fm$yHat[fm$whichNa][length(Test)+1:length(Test)],
                      collapse="\t"), sep="\t"), file=file5, append=T, sep="\n")
        }
    } else {
        if (which_cv=="cv2"){
            cat(paste(rep, paste(fm$yHat[fm$whichNa], collapse="\t"), sep="\t"),
                file=file4, append=T, sep="\n")
        }
        if (which_cv=="cv1"){
            cat(paste(rep, paste(fm$yHat[fm$whichNa][length(Test)+1:length(Test)],
                      collapse="\t"), sep="\t"), file=file4, append=T, sep="\n")
            cat(paste(rep, paste(fm$yHat[fm$whichNa][1:length(Test)], collapse="\t"),
                      sep="\t"), file=file5, append=T, sep="\n")
            
        }
    }
}