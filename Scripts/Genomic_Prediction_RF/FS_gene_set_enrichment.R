################################################################################
# Gene set enrichment analysis of top random forest features in each condition #
#             * Separate pathway and GO term enrichment analyses *             #
################################################################################

library(tidyr)
library(dplyr)
library(stringr)
library(parallel)
library(qvalue)

setwd("/mnt/home/seguraab/Shiu_Lab/Project")

# 20210421 S288C GO term file
go <- read.csv("Data/yeast_GO/genes_go_sgd.tsv", sep="\t")
# > length(unique(go$GO))
# [1] 5649

# 20210421 S288C Pathway file
pwys <- read.delim("Data/S288C_reference_genome_R64-3-1_20210421/All-instances-of-Genes-in-Saccharomyces-cerevisiae-S288c.txt", skip=1, header=TRUE)

# SNP features mapped to genes file
genes <- read.csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes.txt",
                  sep=",", header=FALSE)
genes <- left_join(genes, pwys, by=c("V4"="Genes")) # add pathway annotations then unnest them
genes$Pathways.of.gene <- genes$Pathways.of.gene %>% replace_na("none") # replace NAs
genes$Pathways.of.gene[genes$Pathways.of.gene==""] <- "none" # replace ""
genes <- genes[!(genes$V4=="intergenic"),] # drop intergenic genes
genes2 <- genes %>% dplyr::mutate(Pathways.of.gene=strsplit(as.character(Pathways.of.gene), " // ")) %>% unnest(Pathways.of.gene)
#write.csv(genes2, "Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_pwys.csv", quote=F, row.names=F)
#write.csv(unique(genes2$Pathways.of.gene), "Data/Peter_2018/biallelic_snps_diploid_S288C_genes_pwys_unique.csv", quote=F, row.names=F)

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
# contingency table for each PATHWAY:           #
#               | Gene in top | Gene not top    #
#-------------------------------------------    #
# In pathway    |             |                 #
# Not in pathway|             |                 #
#################################################

# top feature files for each condition
dir <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP" # path to SHAP files
files <- list.files(path=dir, pattern="Genes_SHAP_values_sorted_[A-Z0-9]+_training.csv", 
                    full.names=TRUE, recursive=FALSE)

# pathway enrichment on each file 
#for (f in files){
enrichment <- function(f, scores, name, save_name){ 
    #name <- str_extract(f, "SHAP_values_sorted_[A-Z0-9]+_training.csv") # extract file name
    sprintf("File %s", name)
    top <- read.csv(f, header=TRUE) # read in top feature set
    
    # Collapse rows by max SNP importance score
    scores <- read.csv(scores, sep='\t') # RF importance score
    top <- top[1:4] # keep only snp, chr, pos, gene columns
    top_scores <- merge(top, scores, by.x="snp", by.y="X") # merge top and scores
    top_scores_max <- top_scores %>% group_by(gene) %>% summarise_all(max)   # collapse

    #bg <- genes[!(genes$V1 %in% top_$snp),] # remove top from background set based on snps
    #### should I be dropping based on genes instead of snps? that way I don't have overlapping genes in top and bg
    bg <- genes[!(genes$V4 %in% top_scores_max$gene),]  # based on genes
    bg$Pathways.of.gene <- bg$Pathways.of.gene %>% replace_na("none") # replace NAs
    bg$Pathways.of.gene[bg$Pathways.of.gene==""] <- "none" # replace ""
    bg2 <- bg %>% dplyr::mutate(Pathways.of.gene=strsplit(as.character(Pathways.of.gene), " // ")) %>% unnest(Pathways.of.gene, keep_empty=TRUE)
    sprintf("   Top genes: %s", nrow(top_scores_max))
    sprintf("   Genes not in top: %s", length(unique(bg$V4)))

    # drop intergenic snps
    top_scores_max <- top_scores_max[!(top_scores_max$gene=="intergenic"),]
    bg2 <- bg2[!(bg2$V4=="intergenic"),]

    # concatenate pathway information
    top2 <- left_join(top_scores_max, genes, by=c("snp"="V1"))
    top2$Pathways.of.gene <- top2$Pathways.of.gene %>% replace_na("none") # replace NAs
    top2$Pathways.of.gene[top2$Pathways.of.gene==""] <- "none" # replace ""
    top2 <- top2 %>% dplyr::mutate(Pathways.of.gene=strsplit(as.character(Pathways.of.gene), " // ")) %>% unnest(Pathways.of.gene, keep_empty=TRUE)
    #save_name <- str_extract(scores, "[A-Z0-9]+_exp_rf_[0-9]+_imp")
    write.csv(top2, paste("Scripts/Genomic_Prediction_RF/SHAP/Pathways_mapped_", save_name, "_max_collapsed.csv", sep=""), quote=FALSE, row.names=FALSE)

    # contingency table
    cols <- c("Pathways", "Genes_top_in_pathway", "Genes_not_top_in_pathway", 
                    "Genes_top_not_in_pathway", "Genes_not_top_not_in_pathway", 
                    "p.val", "odds ratio", "qvalues", "lfdr")
    contingency <- data.frame(matrix(0, nrow=length(unique(genes2$Pathways.of.gene)), ncol=9, 
                        dimnames=list(unique(genes2$Pathways.of.gene), cols)))
    contingency[,"Pathways"] <- rownames(contingency)

    # fill in contingency table for each gene
    print("   Doing gene set enrichment analysis...")
    for (pathway in unique(genes2$Pathways.of.gene)){
        a <- length(unique(top2[which(top2$Pathways.of.gene==pathway),]$V4)) # Genes in top features and in `pathway`
        b <- length(unique(bg2[which(bg2$Pathways.of.gene==pathway),]$V4)) # Genes not in top features and in `pathway`
        c <- length(unique(top2$V4)) - a # Genes in top feature and not in `pathway`
        d <- length(unique(bg2$V4)) - b # Genes not in top features and not in `pathway`
        tbl <- matrix(c(a, b, c, d), ncol=2, byrow=TRUE) # gene contingency table
        res <- fisher.test(tbl, alternative="two.sided") # fisher's exact test
        contingency[pathway,] <- c(pathway, a, b, c, d, res$p.value, res$estimate, NA, NA)
    }
    # Calculate q-values
    print("   Calculating q values...")
    contingency$p.val <- as.numeric(contingency$p.val)
    sub <- contingency[!(contingency$p.val==1),] # don't include rows with p_val=1
    #name <- gsub(".csv", ".tsv", name) # rename file
    if (dim(sub)[1] > 1){
        qvals <- qvalue(as.numeric(sub$p.val), lambda=0)
        sub$qvalues <- qvals$qvalues
        sub$lfdr <- qvals$lfdr
        out <- rbind(sub, contingency[which(contingency$p.val==1),])
        write.csv(out, paste("Scripts/Genomic_Prediction_RF/SHAP/Pathways_Enriched_", 
                  name, sep=""), quote=FALSE, row.names=FALSE)
    } else {
        write.csv(sub, paste("Scripts/Genomic_Prediction_RF/SHAP/Pathways_Enriched_", 
                   name, sep=""), quote=FALSE, row.names=FALSE)
    }
}

enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPACETATE_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPACETATE_exp_rf_16_imp", "SHAP_values_sorted_YPACETATE_training.csv", "YPACETATE_exp_rf_16_imp")            
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPD14_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPD14_exp_rf_256_imp", "SHAP_values_sorted_YPD14_training.csv", "YPD14_exp_rf_256_imp")                
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPD40_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPD40_exp_rf_1024_imp", "SHAP_values_sorted_YPD40_training.csv", "YPD40_exp_rf_1024_imp")                
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPD42_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPD42_exp_rf_512_imp", "SHAP_values_sorted_YPD42_training.csv", "YPD42_exp_rf_512_imp")                
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPD6AU_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPD6AU_exp_rf_1024_imp", "SHAP_values_sorted_YPD6AU_training.csv", "YPD6AU_exp_rf_1024_imp")               
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDANISO10_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDANISO10_exp_rf_128_imp", "SHAP_values_sorted_YPDANISO10_training.csv", "YPDANISO10_exp_rf_128_imp")           
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDANISO20_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDANISO20_exp_rf_1024_imp", "SHAP_values_sorted_YPDANISO20_training.csv", "YPDANISO20_exp_rf_1024_imp")           
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDANISO50_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDANISO50_exp_rf_256_imp", "SHAP_values_sorted_YPDANISO50_training.csv", "YPDANISO50_exp_rf_256_imp")           
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDBENOMYL200_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDBENOMYL200_exp_rf_512_imp", "SHAP_values_sorted_YPDBENOMYL200_training.csv", "YPDBENOMYL200_exp_rf_512_imp")        
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDBENOMYL500_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDBENOMYL500_exp_rf_1024_imp", "SHAP_values_sorted_YPDBENOMYL500_training.csv", "YPDBENOMYL500_exp_rf_1024_imp")        
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDCAFEIN40_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDCAFEIN40_exp_rf_512_imp", "SHAP_values_sorted_YPDCAFEIN40_training.csv", "YPDCAFEIN40_exp_rf_512_imp")          
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDCAFEIN50_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDCAFEIN50_exp_rf_1024_imp", "SHAP_values_sorted_YPDCAFEIN50_training.csv", "YPDCAFEIN50_exp_rf_1024_imp")          
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDCHX05_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDCHX05_exp_rf_1024_imp", "SHAP_values_sorted_YPDCHX05_training.csv", "YPDCHX05_exp_rf_1024_imp")             
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDCHX1_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDCHX1_exp_rf_512_imp", "SHAP_values_sorted_YPDCHX1_training.csv", "YPDCHX1_exp_rf_512_imp")              
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDCUSO410MM_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDCUSO410MM_exp_rf_1024_imp", "SHAP_values_sorted_YPDCUSO410MM_training.csv", "YPDCUSO410MM_exp_rf_1024_imp")         
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDDMSO_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDDMSO_exp_rf_256_imp", "SHAP_values_sorted_YPDDMSO_training.csv", "YPDDMSO_exp_rf_256_imp")              
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDETOH_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDETOH_exp_rf_32_imp", "SHAP_values_sorted_YPDETOH_training.csv", "YPDETOH_exp_rf_32_imp")              
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDFLUCONAZOLE_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDFLUCONAZOLE_exp_rf_256_imp", "SHAP_values_sorted_YPDFLUCONAZOLE_training.csv", "YPDFLUCONAZOLE_exp_rf_256_imp")       
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDFORMAMIDE4_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDFORMAMIDE4_exp_rf_256_imp", "SHAP_values_sorted_YPDFORMAMIDE4_training.csv", "YPDFORMAMIDE4_exp_rf_256_imp")        
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDFORMAMIDE5_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDFORMAMIDE5_exp_rf_128_imp", "SHAP_values_sorted_YPDFORMAMIDE5_training.csv", "YPDFORMAMIDE5_exp_rf_128_imp")        
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDHU_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDHU_exp_rf_128_imp", "SHAP_values_sorted_YPDHU_training.csv", "YPDHU_exp_rf_128_imp")                
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDKCL2M_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDKCL2M_exp_rf_256_imp", "SHAP_values_sorted_YPDKCL2M_training.csv", "YPDKCL2M_exp_rf_256_imp")             
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDLICL250MM_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDLICL250MM_exp_rf_128_imp", "SHAP_values_sorted_YPDLICL250MM_training.csv", "YPDLICL250MM_exp_rf_128_imp")         
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDMV_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDMV_exp_rf_256_imp", "SHAP_values_sorted_YPDMV_training.csv", "YPDMV_exp_rf_256_imp")                
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDNACL15M_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDNACL15M_exp_rf_512_imp", "SHAP_values_sorted_YPDNACL15M_training.csv", "YPDNACL15M_exp_rf_512_imp")           
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDNACL1M_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDNACL1M_exp_rf_1024_imp", "SHAP_values_sorted_YPDNACL1M_training.csv", "YPDNACL1M_exp_rf_1024_imp")            
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDNYSTATIN_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDNYSTATIN_exp_rf_16_imp", "SHAP_values_sorted_YPDNYSTATIN_training.csv", "YPDNYSTATIN_exp_rf_16_imp")          
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDSDS_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDSDS_exp_rf_1024_imp", "SHAP_values_sorted_YPDSDS_training.csv", "YPDSDS_exp_rf_1024_imp")               
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPDSODIUMMETAARSENITE_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPDSODIUMMETAARSENITE_exp_rf_256_imp", "SHAP_values_sorted_YPDSODIUMMETAARSENITE_training.csv", "YPDSODIUMMETAARSENITE_exp_rf_256_imp")
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPETHANOL_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPETHANOL_exp_rf_1024_imp", "SHAP_values_sorted_YPETHANOL_training.csv", "YPETHANOL_exp_rf_1024_imp")            
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPGALACTOSE_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPGALACTOSE_exp_rf_128_imp", "SHAP_values_sorted_YPGALACTOSE_training.csv", "YPGALACTOSE_exp_rf_128_imp")          
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPGLYCEROL_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPGLYCEROL_exp_rf_32_imp", "SHAP_values_sorted_YPGLYCEROL_training.csv", "YPGLYCEROL_exp_rf_32_imp")           
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPRIBOSE_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPRIBOSE_exp_rf_128_imp", "SHAP_values_sorted_YPRIBOSE_training.csv", "YPRIBOSE_exp_rf_128_imp")             
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPSORBITOL_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPSORBITOL_exp_rf_256_imp", "SHAP_values_sorted_YPSORBITOL_training.csv", "YPSORBITOL_exp_rf_256_imp")           
enrichment("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/Genes_SHAP_values_sorted_YPXYLOSE_training.csv", "/mnt/scratch/seguraab/yeast_project/yeast_rf_results/YPXYLOSE_exp_rf_256_imp", "SHAP_values_sorted_YPXYLOSE_training.csv", "YPXYLOSE_exp_rf_256_imp") 


# GO term enrichment
#Enrichment <- function(k,n,C,G){ return((k / C) / (n/ G)) }




### Interpretation: pvalue = 1 means snp is only in negative set, so remove before qvalue transformation
# q value transformation looks at pvalue dist which is saying multiple pvalues (multiple tests) contributing to a hypothesis,
# the transformation corrects for multiple testing. If many pvalues are 1, then it skews the dist and they have to be removed
