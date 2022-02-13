##################################################################################################
# Description: Genomic BLUP for predicting phenotype within a cross-validation (cv) scheme
# 
# Note: To perform OBLUP and CBLUP as described in Li and Simianer 2018, use the open reading 
# frame (ORF) presence/absence matrix or the ORF copy number matrix as the input X_file.
#   
# Model: y = Xb + Zu + e
# y is a phenotype vector, b is a fixed effect estimate vector, u is random breeding value vector
# X, Z are incidence matrices for fixed and random effects, respectively, e is residuals vector
# 
# Arguments:
#   [1] X_file:     genotype matrix (default), orf presence/absence matrix, or orf copy number matrix
#   [2] Y_file:     phenotype matrix
#   [3] feat_file   selected features or "all" for all features in X_file
#   [4] test_file   file with samples in test set
#   [5] trait:      name of trait (column in phenotype matrix) or "all" for all traits in Y_file
#   [6] cvf_file:   cross-validation scheme file (specifies which sample is part of each cv fold)
#   [7] cv:         cross-validation fold number
#   [8] number:     number of cross-validation repetitions
#   [9] save_name:  output file save name
#   [10] orf:       y/n if X_file is orf presence/absence matrix
#   [11] cno:       y/n if X_file is orf copy number matrix
# 
# Output:
#   [1] Feature coefficients file
#   [2] R-squared values for each cross-validation repetition and of test set
#   [3] Predicted values for each cross-validation repetition and of test set
#
# Modified by: Kenia Segura Abá
##################################################################################################

# Load necessary pakages
library(rrBLUP)
library(data.table)
library(ggplot2)

set.seed(42) # for reproducibility

# Read in arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 11) {
    stop("Need 12 arguments: X_file Y_file feat_file test_file trait cvf_file fold number save_name", call.=FALSE)
} else {
    X_file <- args[1]
    Y_file <- args[2]
    feat_file <- args[3]
    test_file <- args[4]
    trait <- args[5]
    cvf_file <- args[6]
    cv <- as.numeric(args[7])
    number <- as.numeric(args[8])
    save_name <- args[9]
    orf <- args[10]
    cno <- args[11]
}

# Arguments for debugging
#setwd("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_GBLUP")
#X_file <- "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv"
#Y_file <- "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv"
#feat_file <- "/mnt/scratch/seguraab/yeast_project/yeast_rrBLUP_results/Markers_top2000.txt"
#trait <- "YPDCAFEIN40"
#test_file <- "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt"
#cv <- 5
#number <- 10
#cvf_file <- "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/CVFs.csv"
#save_name <- "rrBLUP_geno_Markers_top2000"

# Read in data
# if file is larger than 10Mb, using fread to read the file
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
Test <- scan(test_file, what='character')
cvs <- read.csv(cvf_file, row.names=1)

# Process cross-validation file
cvs_all <- merge(Y,cvs,by="row.names",all.x=TRUE)
rownames(cvs_all) <- cvs_all$Row.names
cvs_all <- cvs_all[,(dim(Y)[2]+2):ncol(cvs_all)]
cvs_all[is.na(cvs_all)] = 0

# make sure X and Y have the same order of rows as cvs_all
X <- X[rownames(cvs_all),]
Y <- Y[rownames(cvs_all),]

# Trait or traits to be modelled
if (trait == 'all') {
  print('Modeling all traits')
} else {
  Y <- Y[trait]
}

# Function to perform G/O/CBLUP within a cross-validation scheme
run_model <- function(G){
    start.time <- Sys.time() # start time
    PCC_cv <- c() # CV Pearson Correlation Coefficient (PCC)
    PCC_test <- c() # Test set PCC
    R2_cv <- c() # CV performance R-squared
    R2_test <- c() # Test set performance R-squared
    Predict_validation <- c() # Predicted label for CV
    Predict_test <- c() # Predicted label for test set
    for(i in 1:length(Y)){
        print(sprintf("Modeling trait %s...", names(Y)[i]))
        corr_CV <- c() # PCC of cross-validation
        corr_test <- c() # PCC of test set
        Accuracy_CV <- c() # R-sq of cross-validation (performance)
        Accuracy_test <- c() # R-sq of test set (performance)
        Coef <- c() # feature coefficients
        Error <- c() # residual error term (Ve) per trait
        Beta <- c() # fixed effects (β) per trait
        pred_val <- c() # predicted value of validation
        pred_test <- c() # predicted value of test set
        for(k in 1:number){
            print(sprintf("CV repetition number %i", k))
            tst <- cvs_all[,k] # column from cvs_all that specifies sample folds for this repetition
            Coeff <- c() # model coefficients for this repetition
            Errors <- c() # residual error term (Ve) per j-cv repetition
            Betas <- c() # fixed effects (β) per j-cv repetition
            y_test <- c() # predicted values of test set
            yhat <- data.frame(cbind(Y, yhat = 0)) # dataframe to hold predicted values
		    yhat$yhat <- as.numeric(yhat$yhat)
		    row.names(yhat) <- row.names(Y)
            for(j in 1:cv){
                print(sprintf("CV fold number %i", j))
                validation <- which(tst==j) # validation set for this fold
                training <- which(tst!=j & tst!=0) # training set is all other data excluding test set
                test <- which(tst==0) # testing set
                yNA <- Y[,i] # label (dependent variable)
                yNA[validation] <- NA # mask validation sample values
                yNA[test] <- NA # mask test sample values

                # Build G/O/CBLUP model
                start <- Sys.time()
                fit <- mixed.solve(y=Y[training,i], K=G[training,])
                end <- Sys.time() # end time
                print(sprintf("Model elapsed time: %f", end-start))
                
                # Extract results from model
                Coeff <- rbind(Coeff,fit$u) # feature coefficients
                effect_size <- as.matrix(fit$u)
                Errors <- rbind(Errors, fit$Ve) # residual error term (Ve) per cv fold
                Betas <- rbind(Betas, fit$beta) # fixed effects term (β) per cv fold
                yhat$yhat[validation] <- (as.matrix(G[validation,]) %*% effect_size)[,1] + c(fit$beta) # predicted labels for validation set
                yhat$yhat[test] <- (as.matrix(G[test,]) %*% effect_size)[,1] + c(fit$beta) # predicted labels for test set
                y_test <- cbind(y_test, yhat$yhat) # collect predicted labels
            }
            # Performance measures per repetition
            print("Calculating performance...")
            corr_cv <- cor(yhat[which(tst!=0),i], yhat$yhat[which(tst!=0)]) # PCC of validation column
            corr_CV <- c(corr_CV, corr_cv)
            Accuracy_CV <- c(Accuracy_CV,corr_cv^2) # R-sq of validation

            y_test <- cbind(y_test,rowMeans(y_test)) # mean predicted values
            corr_Test <- cor(yhat[which(tst==0),i], y_test[which(tst==0),ncol(y_test)]) # PCC of test
            corr_test <- c(corr_test, corr_Test)
            Accuracy_test <- c(Accuracy_test,corr_Test^2) # R-sq of test

            Coef <- rbind(Coef,colMeans(Coeff)) # mean feature coefficients
            Error <- rbind(Error, colMeans(Errors)) # mean residual errors across folds
            Beta <- rbind(Beta, colMeans(Betas)) # mean fixed effects across folds
            pred_val <- cbind(pred_val,yhat$yhat[which(tst!=0)]) # predicted values of validation set
            pred_test <- cbind(pred_test,y_test[which(tst==0),ncol(y_test)]) # predicted values of test set
        }
        # Overall performance
        PCC_cv <- cbind(PCC_cv, corr_CV) # cross-validation PCC
        PCC_test <- cbind(PCC_test, corr_test) # test set PCC
        R2_cv <- cbind(R2_cv,Accuracy_CV) # cross-validation R-sq
        R2_test <- cbind(R2_test,Accuracy_test) # test set R-sq
        write.csv(Coef,paste('Coef_',save_name,'_',names(Y)[i],'.csv',sep=''),row.names=FALSE,quote=FALSE) # Coefficients
        write.csv(Error,paste('Error_',save_name,'_',names(Y)[i],'.csv',sep=''),row.names=FALSE,quote=FALSE) # Save model average residual errors across folds
        write.csv(Beta,paste('Beta_',save_name,'_',names(Y)[i],'.csv',sep=''),row.names=FALSE,quote=FALSE) # Save model average fixed effects across folds
        colnames(pred_val) <- paste(names(Y)[i], "_", 1:number, sep="") # columns for predicted values for each repetition
        Predict_validation <- cbind(Predict_validation, pred_val) # validation predicted values
        colnames(pred_test) <- paste(names(Y)[i], "_", 1:number, sep="") # columns for predicted values for each repetition
        Predict_test <- cbind(Predict_test, pred_test) # test set predicted values
    }
    print("Complete.")
    end.time <- Sys.time() # end time
    print(sprintf("Total elapsed time: %f", end.time-start.time))

    # Save results to files
    colnames(PCC_cv) <- names(Y)
    colnames(PCC_test) <- names(Y)
    write.csv(PCC_cv,paste('PCC_cv_results_',save_name,'_',trait,'.csv',sep=''),row.names=FALSE,quote=FALSE)
    write.csv(PCC_test,paste('PCC_test_results_',save_name,'_',trait,'.csv',sep=''),row.names=FALSE,quote=FALSE)
    colnames(R2_cv) <- names(Y)
    colnames(R2_test) <- names(Y)
    write.csv(R2_cv,paste('R2_cv_results_',save_name,'_',trait,'.csv',sep=''),row.names=FALSE,quote=FALSE)
    write.csv(R2_test,paste('R2_test_results_',save_name,'_',trait,'.csv',sep=''),row.names=FALSE,quote=FALSE)
    rownames(Predict_validation) <- rownames(G)[which(tst!=0)]
    write.csv(Predict_validation,paste('Predict_value_cv_',save_name,'_',trait,'.csv',sep=''),row.names=TRUE,quote=FALSE)
    rownames(Predict_test) <- rownames(G)[test]
    write.csv(Predict_test,paste('Predict_value_test_',save_name,'_',trait,'.csv',sep=''),row.names=TRUE,quote=FALSE)
}

# Calculate covariance matrix and run the G/O/CBLUP
print("Calculating covariance matrix...")
if (orf == "y" & cno == "n") { # covariance of ORF presence/absence
    print("Performing OBLUP...")
    # OBLUP: O = WW'/sum(q(1-q)) 
    # W is ORF matrix centered as (0-q) and (1-q) for each entry
    q <- colSums(X)/nrow(X) # q is the frequency of each orf
    Q <- matrix(rep(q, nrow(X)), ncol=ncol(X), nrow=nrow(X), byrow=TRUE)
    rownames(Q) <- rownames(X) ; colnames(Q) <- colnames(X)
    W <- X - Q # W matrix
    sumq <- sum(q*(1-q)) # scaling factor
    O <- (W %*% t(W))/sumq
    # Perform OBLUP
    run_model(O)
} else if (cno == "y" & orf == "n") { # covariance of ORF copy number
    print("Performing CBLUP...")
    # CBLUP: C = SS'/f
    # S is copy number matrix with (b-u) entries where b is the copy number
    u <- colMeans(X) # u is mean copy number
    U <- matrix(rep(u, nrow(X)), ncol=ncol(X), nrow=nrow(X), byrow=TRUE)
    rownames(U) <- rownames(X) ; colnames(U) <- colnames(X)
    S <- X - U # S matrix
    f <- median(diag(S %*% t(S))) # f is the median of the diagonal of SS'
    C <- (S %*% t(S))/f # copy number
    # Perform CBLUP
    run_model(C)
} else { # covariance of genotypes
    print("Performing GBLUP...")
    # GBLUP: G = (ZZ')/(2*sum(p*(1-p))); VanRaden et al. 2008
    # Z is the MAF adjusted marker matrix with (0-2p), (1-2p), (2-2p) entries
    # genotypes are 0, 1, 2 for -1, 0, 1, respectively
    X[X==1] <- 2 # recode genotype matrix
    X[X==0] <- 1
    X[X==-1] <- 0
    write.csv(X, "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno_012.csv", quote=FALSE, row.names=FALSE)
    p <- colMeans(X)/2 # minor allele frequencies of each SNP
    # var(binomial) = n*theta(1-theta) ==> var(SNP) = 2p(1-p)
    sum2pq <- 2*sum(p*(1-p)) # scaling factor
    # Center X by mean allele frequencies
    P <- matrix(rep(2*p, nrow(X)), ncol=ncol(X), nrow=nrow(X), byrow=TRUE)
    rownames(P) <- rownames(X) ; colnames(P) <- colnames(X)
    Z <- X - P # Z matrix 
    # calculate covariance matrix G
    G <- (Z %*% t(Z)) / sum2pq
    # Perform GBLUP
    run_model(G)
    # plot of allele frequencies
    p <- as.data.frame(p)
    g <- ggplot(p2, aes(x=p)) +
            geom_histogram(bins=100, color="white") +
            xlab("Allele Frequencies") +
            theme_minimal()
    ggsave("SNP_allele_freq_histogram.pdf")
}
