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
# Model Output: Vu (additive genetic variance), Ve (residual variance), beta (intercept),
#               u (additive genetic values)
# Other stats: Vu^2/(Vu^2+Ve^2) as genomic heritability, Ve/Vu as lambda (ratio of variance components)
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
#   [12] snp:       y/n if X_file is snp genotype matrix
#   [13] cv:        y/n if cross-validation
# 
# Output:
#   [1] Breeding values (BV)
#   [2] Genomic heritability estimates (h2)
#   [3] Additive genetic variances (Va)
#   [4] Residual variances (Ve)
#   [5] Intercepts (Beta)
#   [3] R-squared values for each cross-validation repetition and of test set (R2)
#   [4] Predicted values for each cross-validation repetition and of test set 
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
if (length(args) < 12) {
    stop("Need 13 arguments: X_file Y_file feat_file test_file trait cvf_file fold number save_name orf cno snp cv", call.=FALSE)
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
    snp <- args[12]
    cv <- args[13]
}

# Arguments for debugging
#setwd("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_GBLUP")
#X_file <- "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv"
#Y_file <- "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv"
#feat_file <- "all"
#trait <- "YPDCAFEIN40"
#test_file <- "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt"
#cv <- 5
#number <- 10
#cvf_file <- "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/CVFs.csv"
#save_name <- "rrBLUP_geno_Markers_top2000"
#i=1;j=1;k=1
#y=Y[,trait]; W=X; maf=0.05

# Read in data
# if file is larger than 10Mb, using fread to read the file
print("Reading in data...")
if (feat_file != "all") {
    message("Pulling features to use...")
    FEAT <- scan(feat_file, what="character") # determine which features are included
    X <- fread(X_file, select=c("ID", FEAT), showProgress=TRUE) # subset genotype data 
    X <- as.matrix(X, rownames=1, colnames=1)
    X[1:5,1:5]
} else {
    X <- fread(X_file, showProgress=TRUE)
    X <- as.matrix(X, rownames=1, colnames=1)
    X[1:5,1:5]
}

Y <- read.csv(Y_file, row.names=1) # phenotype data

if (cv == "y"){ # Test file and Cross-validation file with fold assignments
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
}

# Trait or traits to be modelled
if (trait == 'all') {
  message('Modeling all traits')
} else {
  Y <- Y[trait]
}

# Function to perform G/O/CBLUP within a cross-validation scheme
cv_gblup <- function(G){ 
    for(i in 1:length(Y)){
        message(sprintf("Modeling trait %s...", names(Y)[i]))
        BV <- c() # breeding values
        Va <- c() # additive genetic variance term (Vu) per trait
        Ve <- c() # residual variance term (Ve) per trait
        Heritability <- c() # genomic heritability estimate
        Beta <- c() # intercepts (b) per trait
        for(k in 1:number){
            message(sprintf("CV repetition number %i", k))
            tst <- cvs_all[,k] # column from cvs_all that specifies sample folds for this repetition
            Bv <- c() # breeding values for this repetition
            va <- c() # additive genetic variance term (Vu) per j-cv repetition
            ve <- c() # residual variance term (Ve) per j-cv repetition
            Heritabilities <- c()  # genomic heritability estimate
            Betas <- c() # intercepts (b) per j-cv repetition
            for(j in 1:cv){
                message(sprintf("CV fold number %i", j))
                validation <- which(tst==j) # validation set for this fold
                training <- which(tst!=j & tst!=0) # training set is all other data excluding test set
                test <- which(tst==0) # testing set
                yNA <- Y[,i] # label (dependent variable)
                yNA[validation] <- NA # mask validation sample values
                yNA[test] <- NA # mask test sample values

                # Build G/O/CBLUP model
                system.time(fit <- mixed.solve(y=Y[training,i], K=G[training,training]))
                
                # Extract results from model
                Heritabilities <- rbind(Heritabilities, fit$Vu^2 / (fit$Vu^2 + fit$Ve^2)) # genomic heritability estimate
                Bv <- rbind(Bv,fit$u) # breeding values
                va <- rbind(va, fit$Vu) # additive genetic variance term (Vu) per cv fold
                ve <- rbind(ve, fit$Ve) # residual variance term (Ve) per cv fold
                Betas <- rbind(Betas, fit$beta) # intercept term (b) per cv fold
            }
            Heritability <- rbind(Heritability, colMeans(Heritabilities)) # genomic heritabilities
            BV <- rbind(BV,colMeans(Bv)) # mean breeding values
            Va <- rbind(Va,colMeans(va)) # additive genetic variances
            Ve <- rbind(Ve, colMeans(ve)) # mean residual variances across folds
            Beta <- rbind(Beta, colMeans(Betas)) # mean intercept across folds
        }
        # Write model parameters to file
        write.csv(Heritability,paste('h2_',save_name,'_',names(Y)[i],'.csv',sep=''),row.names=FALSE,quote=FALSE) # Genomic heritabilities
        write.csv(BV,paste('BV_',save_name,'_',names(Y)[i],'.csv',sep=''),row.names=FALSE,quote=FALSE) # Breeding values
        write.csv(Va,paste('Va_',save_name,'_',names(Y)[i],'.csv',sep=''),row.names=FALSE,quote=FALSE) # Additive genetic variances
        write.csv(Ve,paste('Ve_',save_name,'_',names(Y)[i],'.csv',sep=''),row.names=FALSE,quote=FALSE) # Save model average residual variances across folds
        write.csv(Beta,paste('Beta_',save_name,'_',names(Y)[i],'.csv',sep=''),row.names=FALSE,quote=FALSE) # Save model average intercept across folds
    }
    message("Complete.")
}

# Function to perform G/O/CBLUP (I may not need this function, and just add the fit mixed.solve line to the cv=='n' in the next section below.
g_blup <- function(W, y, maf=0.05, recode=FALSE){
    # Source: http://morotalab.org/apsc5984-2020/day29/day29a.html
    # Input: genotype matrix W, phenotype matrix y, minor allele frequency threshold, recode
    # Consideration: Does genotype encoding in W matter? W/out standardizing yes, w/standardizing the Vu, Ve, b, u barely change, but the heritability estimate does quite a bit.
    # # Try: -1, 0, 1 and 0, 1, 2
    if (recode == TRUE){
        W[W==1] <- 2 # recode genotype matrix
        W[W==0] <- 1
        W[W==-1] <- 0
        #write.csv(W, "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno_012.csv", quote=FALSE, row.names=FALSE)
    }
    # Genomic covariance matrix G
    p <- colMeans(W) # allele frequencies
    maf2 <- pmin(p, 1-p)
    index <- which(maf2 < maf) 
    W2 <- W[,-index] # remove genotypes with maf2 < maf
    p2 <- p[-index]
    W2 <- scale(W2, center=TRUE, scale=TRUE)
    G <- tcrossprod(W2)/(2 * sum(p2 * (1-p2))) # according to Li and Simianer 2018
    fit <- mixed.solve(y = y, K = G)
    return(fit)
}

# Collect results
file <- "RESULTS_BLUPs.csv"
if (!file.exists(file)) {
    cat("Date", "Model", "Trait", "PCC", "R-sq", "Intercept", "Vu", "Ve", "lambda", "heritability\n", file=file, append=FALSE, sep=",")
} else {message("RESULTS_BLUPs.csv exists") }

# Calculate covariance matrix and run the G/O/CBLUP
if (orf == "y") { # covariance of ORF presence/absence
    print("Calculating ORF covariance matrix...")
    # OBLUP: O = WW'/sum(q(1-q)) 
    # W is ORF matrix centered as (0-q) and (1-q) for each entry
    q <- colMeans(X) # q is the frequency of each orf
    Q <- matrix(rep(q, nrow(X)), ncol=ncol(X), nrow=nrow(X), byrow=TRUE)
    rownames(Q) <- rownames(X) ; colnames(Q) <- colnames(X)
    W <- X - Q # W matrix DO I NEED TO SCALE AND CENTER?
    W <- scale(W, center = TRUE, scale = FALSE) # standardize mean 0
    sumq <- sum(q*(1-q)) # scaling factor
    O <- (W %*% t(W))/sumq
    
    message("Performing OBLUP...")
    # Perform OBLUP
    if (cv == "y") {
        system.time(cv_gblup(O))
    } else {
        for (i in 1:length(Y)){
            system.time(fit <- mixed.solve(y = Y[,i], K = O))
            print("Saving results...")
            # Append results to file
            pcc <- cor(Y[,i], fit$u)
            cat(Sys.time(), file=file, append=TRUE, sep="")
            cat("OBLUP,", file=file, append=TRUE, sep="")
            cat(names(Y)[i], ",", file=file, append=TRUE, sep="")
            cat(pcc, ",", file=file, append=TRUE, sep="")
            cat(pcc^2, ",", file=file, append=TRUE, sep="")
            cat(fit$beta, ",", file=file, append=TRUE, sep="")
            cat(fit$Vu, ",", file=file, append=TRUE, sep="")
            cat(fit$Ve, ",", file=file, append=TRUE, sep="")
            cat(fit$Ve/fit$Vu, ",", file=file, append=TRUE, sep="")
            cat(fit$Vu^2 / (fit$Vu^2 + fit$Ve^2), "\n", file=file, append=TRUE, sep="")
            # Save predicted values
            write.csv(fit$u, paste(names(Y)[i], "OBLUP.csv", sep="_"))
        }
    }
} else if (cno == "y") { # covariance of ORF copy number
    print("Calculating ORF Copy Number covariance matrix...")
    # CBLUP: C = SS'/f
    # S is copy number matrix with (b-u) entries where b is the copy number
    u <- colMeans(X) # u is mean copy number
    U <- matrix(rep(u, nrow(X)), ncol=ncol(X), nrow=nrow(X), byrow=TRUE)
    rownames(U) <- rownames(X) ; colnames(U) <- colnames(X)
    S <- X - U # S matrix DO I NEED TO SCALE AND CENTER?
    S <- scale(S, center = TRUE, scale = FALSE) # standardize  mean 0
    f <- median(diag(S %*% t(S))) # f is the median of the diagonal of SS'
    C <- (S %*% t(S))/f # copy number
    print("Performing CBLUP...")
    # Perform CBLUP
    if (cv == "y") {
        print("Cross-validation...")
        system.time(cv_gblup(C))
    } else {
        for (i in 1:length(Y)){
            system.time(fit <- mixed.solve(y = Y[,i], K = C))
            print("Saving results...")
            # Append results to file
            pcc <- cor(Y[,i], fit$u)
            cat(Sys.time(), file=file, append=TRUE, sep="")
            cat("CBLUP,", file=file, append=TRUE, sep="")
            cat(names(Y)[i], ",", file=file, append=TRUE, sep="")
            cat(pcc, ",", file=file, append=TRUE, sep="")
            cat(pcc^2, ",", file=file, append=TRUE, sep="")
            cat(fit$beta, ",", file=file, append=TRUE, sep="")
            cat(fit$Vu, ",", file=file, append=TRUE, sep="")
            cat(fit$Ve, ",", file=file, append=TRUE, sep="")
            cat(fit$Ve/fit$Vu, ",", file=file, append=TRUE, sep="")
            cat(fit$Vu^2 / (fit$Vu^2 + fit$Ve^2), "\n", file=file, append=TRUE, sep="")
            # Save predicted values
            write.csv(fit$u, paste(names(Y)[i], "CBLUP.csv", sep="_"))
        }
    }
} else if (snp == "y") { # covariance of genotypes
    print("Calculating SNP covariance matrix...")
    # https://rpubs.com/amputz/GBLUP_and_ssGBLUP 
    # GBLUP: G = (ZZ')/(2*sum(p*(1-p))); VanRaden et al. 2008
    # Z is the MAF adjusted marker matrix with (0-2p), (1-2p), (2-2p) entries
    # genotypes are 0, 1, 2 for -1, 0, 1, respectively
    W=X; print("Recoding X_file from -1,0,1 to 0,1,2 encoding")
    W[W==1] <- 2 # recode genotype matrix
    W[W==0] <- 1
    W[W==-1] <- 0
    #sum(W[,1]==0)/750              genotype freq homozygous ref
    #[1] 0.076
    #sum(W[,1]==1)/750              genotype freq heterozygote
    #[1] 0.04533333
    #sum(W[,1]==2)/750              genotype freq heterozygous alt
    #[1] 0.8786667
    #0.8786667+0.04533333+0.076 = 1
    #0.076+(0.04533333/2)           allele freq p
    #[1] 0.09866667
    #0.8786667+(0.04533333/2)       allele freq q
    #[1] 0.9013334                  this matches colMeans/2
    p <- colMeans(W)/2 # allele frequencies of each SNP
    print(head(p^2)) #p2 genotype frequency
    print(head((1-p)^2)) # q2
    print(head(2*p*(1-p))) # 2pq
    print(head(sum((p^2),((1-p)^2),(2*p*(1-p))))) # should be 1
    ### Note: var(binomial) = n*theta(1-theta) ==> var(SNP) = 2p(1-p)
    sum2pq <- 2*sum(p*(1-p)) # scaling factor
    # Center X by mean allele frequencies
    P <- matrix(rep(2*p, nrow(W)), ncol=ncol(W), nrow=nrow(W), byrow=TRUE)
    rownames(P) <- rownames(W) ; colnames(P) <- colnames(W)
    Z <- W - P # Z matrix 
    Zs <- scale(Z, center = TRUE, scale = TRUE) # standardize mean 0
    # calculate covariance matrix G
    G <- (Zs %*% t(Zs)) / sum2pq
    # plot of allele frequencies
    p <- as.data.frame(p)
    g <- ggplot(p, aes(x=p)) +
            geom_histogram(bins=100, color="white") +
            xlab("Allele Frequencies") +
            theme_minimal()
    ggsave("SNP_allele_freq_histogram.pdf")
    print("Performing GBLUP...")
    # Perform GBLUP
    if (cv == "y") {
        print("Cross-validation...")
        system.time(cv_gblup(G))
    } else {
        for (i in 1:length(Y)){
            system.time(fit <- mixed.solve(y = Y[,i], K = G)) # GBLUP on 0, 1, 2 matrix
            #system.time(fit <- g_blup(X, y=Y[,i], recode=TRUE)) # 0,1,2 ENCODING
            print("Saving results...")
            # Append results to file
            pcc <- cor(Y[,i], fit$u)
            cat(Sys.time(), file=file, append=TRUE, sep="")
            cat("GBLUP,", file=file, append=TRUE, sep="")
            cat(names(Y)[i], ",", file=file, append=TRUE, sep="")
            cat(pcc, ",", file=file, append=TRUE, sep="")
            cat(pcc^2, ",", file=file, append=TRUE, sep="")
            cat(fit$beta, ",", file=file, append=TRUE, sep="")
            cat(fit$Vu, ",", file=file, append=TRUE, sep="")
            cat(fit$Ve, ",", file=file, append=TRUE, sep="")
            cat(fit$Ve/fit$Vu, ",", file=file, append=TRUE, sep="")
            cat(fit$Vu^2 / (fit$Vu^2 + fit$Ve^2), "\n", file=file, append=TRUE, sep="")
            # Save breeding values
            write.csv(fit$u, paste(names(Y)[i], "GBLUP.csv", sep="_"))
        }
    }
    # GBLUP with -1, 0, 1 encoding
    p <- colMeans(X) # minor allele frequencies of each SNP
    print(head(p^2)) #p2
    print(head((1-p)^2)) # q2
    print(head(2*p*(1-p))) # 2pq
    print(head(sum((p^2),((1-p)^2),(2*p*(1-p))))) # should be 1
    # var(binomial) = n*theta(1-theta) ==> var(SNP) = 2p(1-p)
    sum2pq <- 2*sum(p*(1-p)) # scaling factor
    # Center X by mean allele frequencies
    P <- matrix(rep(2*p, nrow(X)), ncol=ncol(X), nrow=nrow(X), byrow=TRUE)
    rownames(P) <- rownames(X) ; colnames(P) <- colnames(X)
    Z <- X - P # Z matrix 
    Zs <- scale(Z, center = TRUE, scale = TRUE) # standardize mean 0
    # calculate covariance matrix G
    G <- (Zs %*% t(Zs)) / sum2pq
    for (i in 1:length(Y)){
        system.time(fit <- mixed.solve(y = Y[,i], K = G)) # GBLUP on -1,0,1 matrix
        cat(Sys.time(), file=file, append=TRUE, sep="")
        cat("GBLUP-101,", file=file, append=TRUE, sep="")
        cat(names(Y)[i], ",", file=file, append=TRUE, sep="")
        cat(pcc, ",", file=file, append=TRUE, sep="")
        cat(pcc^2, ",", file=file, append=TRUE, sep="")
        cat(fit$beta, ",", file=file, append=TRUE, sep="")
        cat(fit$Vu, ",", file=file, append=TRUE, sep="")
        cat(fit$Ve, ",", file=file, append=TRUE, sep="")
        cat(fit$Ve/fit$Vu, ",", file=file, append=TRUE, sep="")
        cat(fit$Vu^2 / (fit$Vu^2 + fit$Ve^2), "\n", file=file, append=TRUE, sep="")
        write.csv(fit$u, paste(names(Y)[i], "GBLUP-101.csv", sep="_"))
    }
}