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

# Filter ORF features mapped to genes file
blastx <- read.delim("Data/S288C_reference_genome_R64-3-1_20210421/S288C_orf_peter_blastx.txt", header=TRUE, skip=1)
strict <- blastx %>% group_by(qacc) %>% 
    dplyr::filter(pident == max(pident) & pident >= 95 & evalue == min(evalue)) %>% 
    as.data.frame()
loose <- blastx %>% group_by(qacc) %>% 
    dplyr::filter(pident == max(pident) & pident >= 30 & evalue == min(evalue) & evalue <= 1e-5) %>% 
    as.data.frame()

# find duplicates
qdupes_strict <- as.data.frame(strict[duplicated(strict[,c('qacc')]) | duplicated(strict[,c('qacc')], fromLast=T),])
qdupes_loose <- as.data.frame(loose[duplicated(loose[,c('qacc')]) | duplicated(loose[,c('qacc')], fromLast=T),])

# Loop through SHAP files for each condition and add gene information
get_genes <- function(f) {
    df <- read.delim(f, sep="\t") # read in shap valueshead
    df2 <-  filter(df, if_all(everything(), ~ .x != 0)) # only keep features with SHAP values
    #df2$orfs <- rownames(df2)
    print(intersect(loose$qacc, df2$X))
    #out <- inner_join(strict, df2, by=c('qacc'='orfs')) # add gene information
    #colnames(out)[1:4] <- c("snp", "chr", "pos", "gene")
    #name <- str_extract(f, "SHAP_values_sorted_[A-Z0-9]+_[a-z]_training[^txt]") # extract file name
    #write.csv(out, paste("Scripts/Genomic_Prediction_RF/SHAP/ORFs/Genes_", 
    #          substr(name,1,nchar(name)-1), ".csv", sep=""), quote=FALSE, 
    #          row.names=FALSE)
}
dir <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/ORFs" # path to SHAP files
files <- list.files(path=dir, pattern="SHAP_values_sorted_average_Y", full.names=TRUE, 
                    recursive=FALSE)
#mclapply(X=files, FUN=get_genes, mc.cores=4) # get genes
lapply(files, get_genes)

##########################################################
#              Gene set enrichment analysis              #
# contingency table for each PATHWAY:                    #
#               | ORF in top | ORF not top               #
#-----------------------------------------               #
# In GO term    |            |                           #
# Not in GO term|            |                           #
##########################################################

# top feature files for each condition
dir <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/ORFs" # path to SHAP files
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
