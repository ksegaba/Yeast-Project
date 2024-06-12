################################################################################
# Pathway enrichment (overrepresentation) analysis of top random forest 
# features in each condition       
################################################################################

rm(list=ls())

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(GSEABase))

############ Read pathway data (doesn't have pathway descriptions) #############
setwd("/mnt/home/seguraab/Shiu_Lab/Co-function/Data/MetaCyc/")
pwy <- read.csv("All-genes-pathways-S288c_pivoted.txt")

########################### Map pathways to ORFs ###########################
setwd("/mnt/home/seguraab/Shiu_Lab/Project/")
genes <- read.delim("Data/Peter_2018/final_map_orf_to_gene.txt", sep="\t", header=1)
pwy2 <- left_join(pwy, genes, by=c("Accession.1"="gene")) # merge
with_pwy <- pwy2[pwy2$Pathways.of.gene!="",] # genes with pathway annotations
with_pwy <- with_pwy[which(!is.na(with_pwy$orf)),]
no_pwy <- pwy2[pwy2$Pathways.of.gene=="",] # genes with no pathway annotations
write.csv(pwy2, 
    "Data/Peter_2018/ORFs_and_S288C_genes_pwy_all.csv", quote=F, row.names=F)
write.csv(with_pwy, 
    "Data/Peter_2018/ORFs_and_S288C_genes_pwy.csv", quote=F, row.names=F)
write.csv(no_pwy, 
    "Data/Peter_2018/ORFs_and_S288C_genes_no_pwy.csv", quote=F, row.names=F)

################################################################################
#                         Pathway Enrichment Analysis                          #
#                    contingency table for each GO term:                       #
#                                  | Gene in top | Gene not top                #
#                   -----------------------------------------                  #
#                    In pathway    |             |                             #
#                    Not in pathway|             |                             #
################################################################################
### Map all ORFs to pathways - prep for background set
all_orfs <- t(read.csv("Data/Peter_2018/ORFs_pres_abs.csv", nrows=1, header=F))
all_orfs <- all_orfs[-1,] # drop "ID"
all_orfs <- gsub(".", "-", all_orfs, fixed=T) # replace periods with dashes
all_orfs <- gsub("^X", "", all_orfs) # remove leading X
all_orfs <- left_join(data.frame(all_orfs), genes, by=c("all_orfs"="orf")) # add gene information
all_orfs <- left_join(all_orfs, pwy[,c("Accession.1", "Pathways.of.gene")],
                      by=c("gene"="Accession.1"), relationship="many-to-many") # add PWY information
colnames(all_orfs) <- c("orf", "gene", "organism", "pathway")
remove(pwy)

enrichment <- function(k, n, C, G){ 
    # determine direction of enrichment
    # if >= 1: + (overrepresented)
    # if < 1: - (underrepresented)
    # k: number of genes in cluster with GO
    # n: total number of genes in cluster
    # C: total number of genes (in cluster + background) with GO
    # G: total number of genes (in cluster + background)
    return((k/C)/(n/G))
}

ora <- function(all_orfs, top, bg, path){
    # Overrepresentation Analysis
    # all_orfs: dataframe of all genes and pathway annotations
    # top: dataframe of genes of interest
    # bg: dataframe of genes in background set
    # path: file path and name to save as

    # create contingency table
    cols <- c("PWY", "Gene_top_has_PWY", "Gene_not_top_has_PWY", "Gene_top_no_PWY",
              "Gene_not_top_no_PWY", "direction", "p.val", "odds ratio", "qvalues")
    contingency <- data.frame(matrix(nrow=1, ncol=9))
    colnames(contingency) <- cols

    # fill in contingency table for each gene
    print("   Running ORA...")
    for (pwy in unique(all_orfs$pathway)){
        a <- length(unique(top[which(top$pathway==pwy),]$gene)) # Genes in top features and have `pwy`
        b <- length(unique(bg[which(bg$pathway==pwy),]$gene)) # Genes not in top features and have `pwy`
        c <- length(unique(top$gene)) - a # Genes in top features and do not have `pwy`
        d <- length(unique(bg$gene)) - b # Genes not in top features and do not have `pwy`
        tbl <- matrix(c(a, b, c, d), ncol=2, byrow=TRUE) # gene contingency table
        res <- fisher.test(tbl, alternative="two.sided") # fisher's exact test
        if (a+b!=0){
            if(enrichment(k=a, n=a+c, C=a+b, G=a+b+c+d) >= 1) direction = '+' else direction = '-' # direction of enrichment
            contingency <- rbind(contingency, list(pwy, a, b, c, d, direction, res$p.value, res$estimate, 'NA'))
        }
    }
    contingency <- contingency[!is.na(contingency$PWY),] # drop first row with NAs
    contingency <- contingency[contingency$PWY!="",]

    # Calculate q-values
    sub <- contingency[!(contingency$p.val==1),] # don't include rows with p_val=1
    
    # Calculate q-values & get heatmap of significant pathways
    if (nrow(sub)!=0){
        print("   Calculating q values...")
        qvals <- p.adjust(sub$p.val, method="BH")
        sub$qvalues <- qvals

        # save contingency table
        sub <- sub[order(sub$qvalues),]
        write.table(sub, path, sep="\t", quote=F, row.names=F)
    }
}

pwy_enrichment <- function(f){
    ### ORA of top gene features
    # read in top gene feature file
    top <- read.delim(f, sep="\t") # read in shap values
    
    # add pathway information
    top <- left_join(top, all_orfs[,c("orf", "pathway")], by=c("orf"="orf"))

    # name of importance score file
    name <- str_extract(f, "[A-Z0-9]+_[a-z]+_[0-9]+_imp") # extract file name
    print(paste("File", name, sep=" "))
    
    ### Step 2: calculate the enrichment score
    bg <- all_orfs[!(all_orfs$orf %in% top$orf),] # remove top from background set
    bg <- bg[!duplicated(bg),] # remove duplicates
    print(paste("   Top ORFs: ", length(unique(top$orf)), sep=""))
    print(paste("   ORFs not in top: ", length(unique(bg$orf)), sep=""))
    print(paste("   Top genes: ", length(unique(top$gene)), sep=""))
    print(paste("   Genes not in top: ", length(unique(bg$gene)), sep=""))
    print(paste("   Total number of ORFs is correct: ", length(unique(top$orf))+length(unique(bg$orf))==length(unique(all_orfs$orf))))
    
    ## Overrepresentation Analysis
    path <- gsub("Genes_", "ORA_PWY_Genes_", f)
    path <- gsub("GO_Enrichment", "PWY_Enrichment", path)
    ora(all_orfs, top, bg, path)
}


# Read in top features' (FS) average SHAP values files
dir <- "Scripts/Genomic_Prediction_RF/GO_Enrichment/ORFs_fs"
orf_files <- list.files(path=dir, pattern="^Genes_[A-Z0-9]+_orf", full.names=T)
cnv_files <- list.files(path=dir, pattern="^Genes_[A-Z0-9]+_cno", full.names=T)

mclapply(X=orf_files, FUN=pwy_enrichment, mc.cores=35)
mclapply(X=cnv_files, FUN=pwy_enrichment, mc.cores=35)
