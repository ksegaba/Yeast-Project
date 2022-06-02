# BayesC Models for 35 Fitness Traits
library(BGLR)
library(data.table)
library(parallel)
library(ggplot2)
library(tidyr)
library(stringr)
library(matrixStats)

X <- fread("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv") # genotype
X <- as.matrix(X, rownames=1, colnames=1)
Xs <- scale(X, center=TRUE, scale=FALSE)
y <- read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv", row.names=1) # phenotype
Xs[1:5,1:5]; y[1:5,1:5]


### BayesC on each trait
ETA <- list(list(X=Xs, model="BayesC", counts=110, probIn=1/100, saveEffects=TRUE))
bayesc <- function(Y){ # Function to run BayesC on each trait
    print(Y)
    fit <- BGLR(y=as.matrix(y[Y]), ETA=ETA, saveAt=paste(Y, "bayesC_mapping_", sep="_"), nIter=7000, burnIn=4000, verbose=FALSE)

    # plot posterior probabilities
    pdf(paste(Y, "bayesC_mapping_posterior.pdf", sep="_"))
    plot(fit$ETA[[1]]$d, main=Y)
    dev.off()

    # collect output
    list(noquote(Y), var(y[Y]), fit$varE, 1-(fit$varE/var(y[Y])), max(fit$ETA[[1]]$d), 
         min(fit$ETA[[1]]$d), sum(fit$ETA[[1]]$d < 0.05), sum(fit$ETA[[1]]$d < 0.01))
}
system.time(tmp <- mclapply(colnames(y), bayesc, mc.cores=35))

# write output to file
out <- do.call(rbind, tmp) #### NEED TO ADD R2 AND PCC
colnames(out) <- c("Y", "varY", "varE", "h2", "max d", "min d", "d05", "d01"); head(out)
write.csv(out, "BayesC_mapping_results.csv", quote=FALSE, row.names=FALSE)


out <- read.csv("BayesC_mapping_results.csv", row.names=1)
#pdf("bayesC_heritability.pdf")
ggplot(out, aes(x=rownames(out), y=h2)) + geom_bar(stat="identity") + 
    theme_bw() + theme(axis.text.x = element_text(angle=90, vjust=0.5,
    hjust=1)) + xlab("Condition")
dev.off()

# Looking in which iteration the SNP was active in the model
setwd("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_BayesC/Results")
B=readBinMat("YPDBENOMYL500_bayesC_mapping_12000_ETA_1_b.bin")
dim(B)
colnames(B) <- colnames(Xs)[-1]
which(colMeans(B!=0)>0.5)
# integer(0)
which(colMeans(B!=0)>0.3)
# [1] 11154

grep("_33476", colnames(X))
# [1] 14625 49345
colnames(X)[14625]
# [1] "chromosome4_334761"

pdf("test.pdf")
plot(B[,11154])
dev.off()

pdf("test.pdf") # manhattan plot
plot(colMeans(B!=0))
abline(v=11154) # SNP 11154
dev.off()




### BayesC 12000 iterations
ETA <- list(list(X=Xs, model="BayesC", count=110, probIn=1/100, saveEffects=TRUE))
bayesc <- function(Y){ # Function to run BayesC on each trait
    print(Y)
    fit <- BGLR(y=as.matrix(y[Y]), ETA=ETA, nIter=12000, burnIn=2000, saveAt=paste(Y, "bayesC_mapping_12000_", sep="_"), verbose=FALSE)

    # plot posterior probabilities
    pdf(paste(Y, "bayesC_mapping_posterior.pdf", sep="_"))
    plot(fit$ETA[[1]]$d, main=Y)
    dev.off()

    # collect output
    list(noquote(Y), var(y[Y]), fit$varE, 1-(fit$varE/var(y[Y])), max(fit$ETA[[1]]$d), 
         min(fit$ETA[[1]]$d), sum(fit$ETA[[1]]$d < 0.05), sum(fit$ETA[[1]]$d < 0.01))
}
system.time(tmp <- mclapply(colnames(y), bayesc, mc.cores=35))
# write output to file
out <- do.call(rbind, tmp)
colnames(out) <- c("Y", "varY", "varE", "h2", "max d", "min d", "d05", "d01"); out
write.csv(out, "BayesC_mapping_results_12000iter.csv", quote=FALSE, row.names=FALSE)




#### BayesC 12000 iterations training and testing set (for lab meeting)
test <- read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/Test.txt", header=F)
X_train <- X[!(rownames(X) %in% test$V1),] ; Xs_train <- scale(X_train, center=TRUE, scale=FALSE)
X_test <- X[test$V1,] ; Xs_test <- scale(X_test, center=TRUE, scale=FALSE) # idk if test is scaled when computing genotypic values
y_train <- y[!(rownames(y) %in% test$V1),]
y_test <- y[test$V1,]

ETA <- list(list(X=Xs_train, model="BayesC", count=110, probIn=1/100, saveEffects=TRUE))
bayesc <- function(Y){ # Function to run BayesC on each trait
    print(Y)
    fit <- BGLR(y=as.matrix(y_train[Y]), ETA=ETA, nIter=12000, burnIn=2000, saveAt=paste(Y, "training_bayesC_mapping_12000_", sep="_"), verbose=FALSE)

    # plot posterior probabilities
    pdf(paste(Y, "training_bayesC_mapping_posterior.pdf", sep="_"))
    plot(fit$ETA[[1]]$d, main=Y)
    dev.off()

    # collect output
    list(noquote(Y), var(y_train[Y]), fit$varE, 1-(fit$varE/var(y_train[Y])), max(fit$ETA[[1]]$d), 
         min(fit$ETA[[1]]$d), sum(fit$ETA[[1]]$d < 0.05), sum(fit$ETA[[1]]$d < 0.01))
}
system.time(tmp <- mclapply(colnames(y_train), bayesc, mc.cores=35))

# write output to file
out <- do.call(rbind, tmp)
colnames(out) <- c("Y", "varY", "varE", "h2", "max d", "min d", "d05", "d01"); out
write.csv(out, "BayesC_mapping_training_results_12000iter.csv", quote=FALSE, row.names=FALSE)

#### Calculate genotypic value and R2, PCC
setwd("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_BayesC/Results")
files <- list.files(pattern="[A-Z0-9]+_training_bayesC_mapping_12000_ETA_1_b.bin", 
                    full.names=TRUE, recursive=FALSE)

# Collect results
file <- "RESULTS_BayesC.csv"
if (!file.exists(file)) {
    cat("Date", "Trait", "yHat_test_SD", "PCC_test_mean", "PCC_test_SD", "R2_test_mean", "R2_test_SD\n", file=file, append=FALSE, sep=",")
} else {message("RESULTS_BLUPs.csv exists") }

cols <- c("Trait", "yHat_test_SD", "PCC_test_mean", "PCC_test_SD", "R2_test_mean", "R2_test_SD")
results <- data.frame(matrix(0, nrow=35, ncol=6, dimnames=list(colnames(y), cols)))
for (f in files){
    B=readBinMat(f) # marker effects for each model iteration
    colnames(B) <- colnames(X)
    G <- data.frame(matrix(nrow=125, ncol=2000)) # predicted genotypic values for each model iteration of the test set
    rownames(G) <- rownames(Xs_test)
    name <- str_extract(f, "[A-Z0-9]+_training_bayesC_mapping_12000_ETA_1_b.bin")
    trait <- gsub("_training_bayesC_mapping_12000_ETA_1_b.bin", "", name) # phenotype name (same as in y dataframe)
    print(trait)
    for (j in 1:nrow(B)){ # calculate predicted genotypic values
        g <- Xs_test %*% B[j,colnames(X_test)] ### IS X_test SCALED WHEN DOING THIS CALCULATION? YES, YOU SCALE TRAINING AND TESTING AFTER SPLITTING
        G[j] <- g
    }
    write.csv(G, paste("Predicted_values_test", trait, "bayesc_12000iter.csv", sep="_"), quote=F, row.names=F)
    sd <- mean(sqrt(rowVars(as.matrix(G)))) # mean standard deviation of yHats (predicted genotypic values) across model iterations for each instance
    pcc <- cor(G, y_test[trait], method="pearson") # Pearson's correlation coefficient
    pcc_sd <- sd(pcc)
    r2 <- pcc^2 # R-squared
    r2_sd <- sd(r2)

    # Append results to file
    cat(Sys.time(), file=file, append=TRUE, sep="")
    cat(trait, ",", file=file, append=TRUE, sep="")
    cat(sd, ",", file=file, append=TRUE, sep="")
    cat(mean(pcc), ",", file=file, append=TRUE, sep="")
    cat(pcc_sd, ",", file=file, append=TRUE, sep="")
    cat(mean(r2), ",", file=file, append=TRUE, sep="")
    cat(r2_sd, "\n", file=file, append=TRUE, sep="")

    results[trait,] <- c(trait, sd, mean(pcc), pcc_sd, mean(r2), r2_sd) # save metrics
}
write.csv(results, "results_bayesc.csv", quote=F, row.names=F)