################################################################################
# Gene set enrichment analysis of top random forest features in each condition #
################################################################################

library(dplyr)
library(stringr)
library(parallel)
library(qvalue) #for q value transformation

setwd("/mnt/home/seguraab/Shiu_Lab/Project")

# SNP features mapped to genes file
genes <- read.csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes.txt",
                  sep=",", header=FALSE)

# Loop through SHAP files for each condition and add gene information
get_genes <- function(f) {
    df <- read.delim(f, sep="\t") # read in shap values
    df2 <- t(df[2:ncol(df)]) # transpose
    colnames(df2) <- df$ID
    out <- merge(genes, df2, by.x="V1", by.y="row.names") # add gene information
    colnames(out)[1:4] <- c("snp", "chr", "pos", "gene")
    name <- str_extract(f, "SHAP_values_sorted_[A-Z0-9]+_training[^txt]") # extract file name
    write.csv(out, paste("Scripts/Genomic_Prediction_RF/SHAP/Genes_", 
              substr(name,1,nchar(name)-1), ".csv", sep=""), quote=FALSE, 
              row.names=FALSE)
}
dir <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP" # path to SHAP files
files <- list.files(path=dir, pattern="SHAP_values_sorted_Y", full.names=TRUE, 
                    recursive=FALSE)
mclapply(X=files, FUN=get_genes, mc.cores=4) # get genes

#################################################
#        Gene set enrichment analysis           #
# contingency table for each GENE:              #
#             |SNP in gene | SNP not in GENE    #
#-------------------------------------------    #
# Top SNPs    |            |                    #
# Not top SNPs|            |                    #
#################################################

# top feature files for each condition
dir <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP" # path to SHAP files
files <- list.files(path=dir, pattern="SHAP_values_sorted_[A-Z0-9]+_training.csv", 
                    full.names=TRUE, recursive=FALSE)

# gene set enrichment on each file 
for (f in files){
    df <- read.csv(f, header=TRUE) # read in top feature set
    bg <- genes[!(genes$V1 %in% df$snp),] # remove top from background set
    print(nrow(df))
    
    # contingency table
    cols <- c("Genes", "Top_SNP_in_gene", "Top_SNP_not_in_gene", 
                    "SNP_not_top_in_gene", "SNP_not_top_not_in_gene", 
                    "p.val", "odds ratio")
    contingency <- data.frame(matrix(0, nrow=length(unique(genes$V4)), ncol=7, 
                        dimnames=list(unique(genes$V4), cols)))
    contingency[,"Genes"] <- rownames(contingency)

    # fill in contingency table for each gene
    for (g in unique(genes$V4)){
        a <- nrow(df[which(df$gene==g),]) # SNPs in top features and in the gene
        b <- nrow(df[which(df$gene!=g),]) # SNPs in top features and not in the gene
        c <- nrow(bg[which(bg$V4==g),]) # SNPs not top features and in the gene
        d <- nrow(bg[which(bg$V4!=g),]) # SNPs not top features and not in the gene
        tbl <- matrix(as.numeric(contingency[g,2:5]), ncol=2, byrow=TRUE) # gene contingency table
        res <- fisher.test(tbl, alternative="two.sided") # fisher's exact test
        pval <- as.numeric(res$p.value) # p-value
        odds <- as.numeric(res$estimate) # odds ratio
        contingency[g,] <- c(g, a, b, c, d, pval, odds)
    }
    # Calculate q-values
    sub <- contingency[which(contingency[,"Top_SNP_in_gene"]!=0),] # don't include rows with p_val=1
    qobj <- qvalue(sub[,"p_val"])# q-value and FDR
    contingency <- cbind(contingency, qobj$qvalues)
    contingency <- cbind(contingency, qobj$lfdr)

    # check with p.adjust BH method
    qvals <- p.adjust(, method="BH")
    name <- str_extract(f, "SHAP_values_sorted_[A-Z0-9]+_training.csv") # extract file name
    write.csv(contingency, paste("Scripts/Genomic_Prediction_RF/SHAP/Enrichment_Genes_", 
              name, sep=""), quote=FALSE, row.names=FALSE)
}

# ADD GO column
Enrichment <- function(k,n,C,G){ return((k / C) / (n/ G)) }

# Yeast GO term file
#go <- read.csv("Data/yeast_GO/genes_go_sgd.tsv", sep="\t")
#newdf <- go %>% group_by(Gene) %>% summarise_at("GO", toString)
#write.csv(newdf, "Data/yeast_GO/genes_go_sgd_collapsed.csv", quote=FALSE, 
#          row.names=FALSE)


### Interpretation: pvalue = 1 means snp is only in negative set, so remove before qvalue transformation
# q value transformation looks at pvalue dist which is saying multiple pvalues (multiple tests) contributing to a hypothesis,
# the transformation corrects for multiple testing. If many pvalues are 1, then it skews the dist and they have to be removed
