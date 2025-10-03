# Description: Calculate kinship and the marker-based narrow-sense heritability using R package 'heritability'
# Author: Kenia Segura Abá

# Install and load necessary packages
#install.packages("heritability")
library(heritability)
library(rrBLUP)
library(data.table)

# Add arguments
args = commandArgs(trailingOnly=TRUE)
X_file = args[1] # genotype matrix, e.g., geno.csv
Y_file = args[2] # phenotype matrix, e.g., pheno.csv
feat_file = args[3] # selected features file or "all" for all the markers in the genetic matrix
trait = args[4]

# Load phenotype data
Y <- read.csv(Y_file, row.names=1)

# Load genotype data
if (feat_file != 'all'){ # Subset genotype data if feat_file is not all
  print('Pulling features to use...')
  FEAT <- scan(feat_file, what='character') 
  X <- fread(X_file,select=c('ID',FEAT))
  X <- as.matrix(X,rownames=1) # set ID to row names
} else{
	X <- as.matrix(fread(X_file),rownames=1)
	}

# Make sure X and Y have the same row order
X <- X[rownames(Y),]

# Calculate kinship matrix
K <- A.mat(X)
write.csv(K, paste('Kinship_',trait,'.csv',sep=''))

# Calculate narrow-sense heritability
h2 <- marker_h2(data.vector=Y[[trait]],geno.vector=rownames(X),K=K)
write.csv(h2, paste('h2_heritability_',trait,'.csv',sep=''))