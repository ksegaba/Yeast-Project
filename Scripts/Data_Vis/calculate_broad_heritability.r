# Description: Estimate the narrow- and broad-sense heritability of diploid S. cerevisiae isolates using biallelic SNPs data (Peter 2018)

### 0. Install necessary packages
#install.packages('sommer')


### 1. Load packages and data
library(sommer)

geno <- read.csv("~/Shiu_Lab/Project/geno.csv", header=T, row.names=1) # Peter 2018 diploid S. cerevisiae isolate biallelic SNPs data
pheno <- read.csv("~/Shiu_Lab/Project/pheno.csv", header=T, row.names=1) # Peter 2018 diploid S. cerevisiae isolate fitness data

### 2. Compute relationship matrices
A <- A.mat(as.matrix(geno)) # additive relationship matrix
D <- D.mat(as.matrix(geno)) # dominance relationship matrix
E <- E.mat(as.matrix(geno)) # epistatic relationship matrix

### 3. Fit the mixed model 
traits <- colnames(pheno)
results <- data.frame(Conditions=traits,h2=0,h2_SE=0,H2_ADE=0,H2_ADE_SE=0, H2_AD=0, H2_AD_SE=0)
for (i in 1:length(traits)){
  p <- pheno[i] # subset pheno matrix
  p$ID <- rownames(geno)
  p$IDD <- rownames(geno)
  p$IDE <- rownames(geno)
  colnames(p) <- c("trt", "ID", "IDD", "IDE")
  head(p)
  
  ### 4. Estimate the genomic heritability (narrow and broad)
  ans.ADE <- mmer(fixed=trt~1, random=~vs(ID,Gu=A) + vs(IDD,Gu=D), rcov=~units, data=p, verbose = FALSE)
  summary(ans.ADE)$varcomp
  h2 <- vpredict(ans.ADE, h2 ~ (V1) / (V1+V3))
  results$h2[i] <- h2$Estimate[1] # narrow sense heritability
  results$h2_SE[i] <- h2$SE[1]
  H2 <- vpredict(ans.ADE, h2 ~ (V1+V2)/(V1+V2+V3)) # broad sense heritability
  results$H2_AD[i] <- H2$Estimate[1]
  results$H2_AD_SE[i] <- H2$SE[1]

  #ans.ADE2 <- mmer(trt~1, random=~vs(ID,Gu=A) + vs(IDD,Gu=D) + vs(IDE,Gu=E), rcov=~units, data=p, verbose = FALSE)
  #summary(ans.ADE2)
  #H2_ADE <- vpredict(ans.ADE2, h2 ~ (V1+V2+V3) / (V1+V2+V3+V4))
  #results$H2_ADE[i] <- H2_ADE$Estimate[1] # broad sense heritability (w/epistasis)
  #results$H2_ADE_SE[i] <- H2_ADE$SE[1]
}

### 5. Write to file
write.csv(results, "/mnt/scratch/seguraab/yeast_project/Heritability_h2_H2_sommer.csv")

### References
#G C (2016). Genome assisted prediction of quantitative traits using the R package sommer. PLoS ONE, 11, 1-15.
#Peter, J., De Chiara, M., Friedrich, A., Yue, J. X., Pflieger, D., Bergström, A., ... & Schacherer, J. (2018). Genome evolution across 1,011 Saccharomyces cerevisiae isolates. Nature, 556(7701), 339-344.


