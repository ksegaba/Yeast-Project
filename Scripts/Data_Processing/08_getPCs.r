#!/usr/bin/Rscript
#####################################
# Uses a genotype file and outputs PCs and PC plots
#
# Arguments: [1] id (i.e. wheat_599_CIMMYT)
#
# Plots include: - Percent of variance explained by the PCs
#                - Plots of PCs: PC vs PC (PC 1-5) 
#                - Plot the G matrix diagonals (related to population structure)
#
# Written by: Christina Azodi
# Original: 4.24.17
# Modified by：Peipei Wang, 02/25/2021, to save the top 5 PCs as a matrix 
# Modified by: Kenia Segura Aba, 05/10/2021, to perform rrBLUP in one condition
#####################################

# Removes all existing variables from the workspace
rm(list=ls())

# Compute PCs for descriptive stats and for later use in modeling 
args = commandArgs(TRUE)

#### Testing parameters

####

# Input data
# Note from Kenia: I used the filtered phenotype and genotype matrices (750 isolates out of 1011) 
# from 04_filter_diploid.py and 04_filter_pheno.py
X_file <- args[1] # your genetic matrix, e.g., geno.csv
Y_file <- args[2] # your phenotypic matrix, e.g., pheno.csv

X <- read.csv(X_file, row.names=1)
Y <- read.csv(Y_file, row.names=1)

# Create output folder if not already present

#source('../parameters/parameters.r')
# get rid of markers with no variation among individuals
X_tem <- c()
col_names <- c()
for(i in 1:ncol(X)){
	if(length(unique(X[,i])) > 1){
		X_tem <- cbind(X_tem,X[,i])
		col_names <- c(col_names,colnames(X)[i])
		}
	}
rownames(X_tem) <- rownames(X)
colnames(X_tem) <- col_names
X = X_tem	

#### Calculate the G matrix
X=scale(X) # Centers (subtract the column means) and Scales (dividing the centered columns by their stdev)
G=tcrossprod(X) # Take the cross product X transpose
G=G/mean(diag(G))


#### Calculate eigenvalues
EVD=eigen(G)
rownames(EVD$vectors)=rownames(G)
save(EVD,file='EVD.RData')
write.csv(EVD$vectors,"PCA_matrix.csv",row.names=T,quote=F)
write.csv(EVD$vectors[,1:5],"PCA5_geno.csv",row.names=T,quote=F)
#### Plot proportion of variance explained ####
pdf('VarExplained.pdf')
var_exp <- (EVD$values)/nrow(X)
cum_exp <- list(var_exp[1])

for(i in 2:length(var_exp)){
  cum_exp[[i]] <- as.numeric(cum_exp[i-1]) + as.numeric(var_exp[i])
}

plot(1:nrow(X), cum_exp, xlab='# PCs', ylab='% Variance Explained')
abline(h= c(0.5, 0.9), col='blue')

dev.off()
 
#### Plot the G matrix diagonals ####
pdf('diag.pdf')
hist(diag(G), xlab = 'Diagonal Element (G-Matrix)',breaks=20)
dev.off()


#### Plots the PCs ####
pdf('PCs.pdf')

par(mfrow=c(3,4))
for(i in 1:4){
  for(j in (i+1):5){
    plot(x=EVD$vectors[,i],y=EVD$vectors[,j],
         main=paste0('PC-',i,' Vs. PC-',j),
         xlab=paste0('PC-',i),
         ylab=paste0('PC-',j),
         cex=.1 #,col=as.integer(factor(SUBJECTS$line))
    )
    #print(c(i,j))
  }
}
dev.off()

#### Export % X variance and Yi variance explained by the top 5 PCs ####

XY_var_exp <- cbind('X', sum(var_exp[0:5]))
cols <- names(Y)

# First let's try one condition
y=Y[,1]
fm = lm(y~EVD$vectors[,1:5])
XY_var_exp <- rbind(XY_var_exp, cbind(names(Y)[1], summary(fm)$adj.r.squared))
write.table(XY_var_exp, "VarExplainedOneCondition.csv", sep=",", row.names=FALSE, col.names=FALSE)

# Then, let's try all 35 conditions
for(i in 1:length(Y)){
  y=Y[, cols[i]]
  fm = lm(y~EVD$vectors[,1:5])
  XY_var_exp <- rbind(XY_var_exp, cbind(names(Y)[i], summary(fm)$adj.r.squared))
}
  
write.table(XY_var_exp, 'VarExplained.csv', sep=',', row.names=FALSE, col.names=FALSE)

quit(save='no')
