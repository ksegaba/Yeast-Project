# Description: Run sERRBLUP on diploid S. cerevisiae yeast isolate biallelic SNPs data from Peter 2018

args = commandArgs(trailingOnly = T)
trait <- args[1]
K <- args[2]

### 0. Install Necessary packages
#unzip("EpiGP-v0.1.0.zip") ; list.files()
#devtools::install_local("evojgani-EpiGP-c5761a1/miraculix_0.9.6.tar.gz")
#devtools::install_local("evojgani-EpiGP-c5761a1/RandomFieldsUtils_0.5.9.tar.gz")
#devtools::install_local("evojgani-EpiGP-c5761a1/pkg/")
#install.packages("BGLR")

### 1. Load packages and upload data
library(EpiGP)
#library(BGLR) # contains example data

#data(wheat) # example genetic marker and phenotype data
#gi <- read.csv("Shiu_Lab/Project/Data File S1. Raw genetic interaction datasets: Pair-wise interaction format/SGA_combined.txt",
#                sep = "\t", header=T) # Costanzo 2016 GI scores and fitness data
geno <- read.csv("~/Shiu_Lab/Project/geno.csv", header=T, row.names=1) # Peter 2018 biallelic SNPs data
pheno <- read.csv("~/Shiu_Lab/Project/pheno.csv", header=T, row.names=1) # Peter 2018 diploid S. cerevisiae isolate fitness data

#### 2. Recode marker matrix 
##m <- Recodemarker(wheat.X) # recode to 2, 1, 0 if m == 1, 0, or -1, respectively
#rownames(m) <- rownames(wheat.Y)

geno_sub <- geno[,1:30000]
m_costanzo <- Recodemarker(geno_sub)
print("m_costanzo") ; m_costanzo[1:5,1:5]

### 3. Calculate ERRBLUP relationship matrix
#G_ERRBLUP <- Gall(m) # generate relationship matrix based on all pairwise SNP interactions
#rownames(G_ERRBLUP) <- rownames(m)

G_ERRBLUP_costanzo <- Gall(m_costanzo)
rownames(G_ERRBLUP_costanzo) <- rownames(m_costanzo)
print("G_ERRBLUP_costanzo") ; G_ERRBLUP_costanzo[1:5,1:5]

### 3. Afterwards, we consider a subset of phenotypic values as a training set and do the
#phenotype prediction based on ERRBLUP; we have done this by 5-fold cross-validation
#with 5 replicates. As an example, here we randomly select 60 lines as the test set and
#the remaining lines are considered as the training set as follows:
#Pheno <- wheat.Y[,1] # subset 1 phenotype treatment
#Pheno_train <- Pheno[1:round(4*length(Pheno)/5)] # create training set
#ERRBLUP_pred <- ERRBLUP(Pheno_train, G_ERRBLUP) # predicted test set phenotypes

pheno2 <- as.matrix(pheno[,trait])
rownames(pheno2) <- rownames(pheno) ; colnames(pheno2) <- trait
pheno_train <- pheno2[1:round(4*length(pheno2)/5)]
ERRBLUP_pred_costanzo <- ERRBLUP(pheno_train, G_ERRBLUP_costanzo)
head(ERRBLUP_pred_costanzo)

### 4. Then, interaction effects and interaction effect variances are calculated as follows:
#t_hat <- SNP_effect(m, Pheno_train, G_ERRBLUP)
#sigma_hat <- SNP_var(m, Pheno_train, t_hat)

t_hat_costanzo <- SNP_effect(m, pheno_train, G_ERRBLUP_costanzo)
sigma_hat_costanzo <- SNP_var(m, pheno_train, t_hat_costanzo)
head(t_hat_costanzo) ; head(sigma_hat_costanzo)

### 5. Based on the desired proportion of interactions, the sERRBLUP relationship matrix is
#calculated for both selections based on estimated effects and estimated effect variances.
#k <- 10
#Gtop_effect <- Gtop(m, Pheno, t_hat, k)
#Gtop_var <- Gtop(m, Pheno, sigma_hat, k)

k <- K
Gtop_effect_costanzo <- Gtop(m, pheno2, t_hat_costanzo, k)
Gtop_var_costanzo <- Gtop(m, pheno2, sigma_hat_costanzo, k)
head(Gtop_effect_costanzo) ; head(Gtop_var_costanzo)

### Finally, sERRBLUP is performed based on both approaches as follows:
#sERRBLUP_effect <- sERRBLUP(Pheno_train, G_ERRBLUP, Gtop_effect)
#sERRBLUP_var <- sERRBLUP(Pheno_train, G_ERRBLUP, Gtop_var)
#sERRBLUP_effect <- as.matrix(sERRBLUP_effect) ; sERRBLUP_var <- as.matrix(sERRBLUP_var)
#rownames(sERRBLUP_effect) <- rownames(m) ; rownames(sERRBLUP_var) <- rownames(m)

sERRBLUP_effect_costanzo <- sERRBLUP(pheno_train, G_ERRBLUP_costanzo, Gtop_effect_costanzo)
sERRBLUP_var_costanzo <- sERRBLUP(pheno_train, G_ERRBLUP_costanzo, Gtop_var_costanzo)
sERRBLUP_effect_costanzo <- as.matrix(sERRBLUP_effect_costanzo) ; sERRBLUP_var_costanzo <- as.matrix(sERRBLUP_var_costanzo)
rownames(sERRBLUP_effect_costanzo) <- rownames(m_costanzo) ; rownames(sERRBLUP_var_costanzo) <- rownames(m_costanzo)
head(sERRBLUP_effect_costanzo) ; head(sERRBLUP_var_costanzo)

to_save <- cbind(sERRBLUP_effect_costanzo, sERRBLUP_var_costanzo) ; colnames(to_save) <- c("Effect", "Var")
write.csv(to_save, paste("sERRBLUP_eff_var_",trait,".csv", sep=""), row.names=T, quote=F)


