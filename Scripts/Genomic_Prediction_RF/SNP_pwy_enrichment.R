################################################################################
# Pathway enrichment (overrepresentation) analysis of top random forest 
# features in each condition       
################################################################################

rm(list=ls())

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(GSEABase))
suppressPackageStartupMessages(library(topGO)) 
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(pvclust))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(RColorBrewer))

########### Reshape pathway data (doesn't have pathway descriptions) ###########
setwd("/mnt/home/seguraab/Shiu_Lab/Co-function/Data/MetaCyc/")
pwy <- read.csv("All-genes-pathways-S288c.txt", sep="\t")
pwy2 <- pwy %>% separate_rows(Pathways.of.gene, sep="[/]+") %>% as.data.frame()
write.csv(pwy2, "All-genes-pathways-S288c_pivoted.txt", quote=F, row.names=F)

########################### Map pathways to SNPs ###########################
setwd("/mnt/home/seguraab/Shiu_Lab/Project/")
genes <- read.csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes.txt",
                  header=F)
colnames(genes) <- c("snp", "chr", "pos", "gene")
genes <- genes[!(genes$Accession.1=="intergenic"),] # drop intergenic snps
pwy2 <- left_join(pwy2, genes, by=c("Accession.1"="gene")) # merge
with_pwy <- pwy2[pwy2$Pathways.of.gene!="",] # genes with pathway annotations
sum(is.na(with_pwy$snp)) # check for no missing data
with_pwy <- with_pwy[which(!is.na(with_pwy$snp)),]
no_pwy <- pwy2[!pwy2$Pathways.of.gene!="",] # genes with no pathway annotations
write.csv(pwy2, 
    "Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_pwy_all.csv", 
    quote=F, row.names=F)
write.csv(with_pwy, 
    "Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_pwy.csv", 
    quote=F, row.names=F)
write.csv(no_pwy, 
    "Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_no_pwy.csv", 
    quote=F, row.names=F)

################################################################################
#                         Pathway Enrichment Analysis                          #
#                    contingency table for each GO term:                       #
#                                  | Gene in top | Gene not top                #
#                   -----------------------------------------                  #
#                    In pathway    |             |                             #
#                    Not in pathway|             |                             #
################################################################################

get_genes <- function(f, baseline=F) {
    # add the genes and pathway annotations to each SHAP file
    df <- read.delim(f, sep="\t") # read in shap values
    colnames(df) <- c("SNP", "shap_value")
    out <- right_join(with_pwy, df, by=c("snp"="SNP")) # add gene & pwy information
    name <- str_extract(f, "SHAP_values_sorted_average_[A-Z0-9]+_[0-9]+_training[^txt]") # extract file name
    name <- gsub("SHAP", "PWY_SHAP", name)
    path <- as.character("Scripts/Genomic_Prediction_RF/SHAP/SNPs/fs")
    write.csv(out, paste(file.path(path, name), "csv", sep=""), quote=F, row.names=F)
    return(out)
}

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

heatmap2 <- function(toplot, rownames, path){
    toplot$logP <- 0
    toplot <- toplot[order(toplot$qvalues, decreasing=T),]
    for(i in 1:nrow(toplot)){ # take the log of the q-values
        if(toplot$direction[i] == '-') toplot$logP[i] <- log10(toplot$qvalues[i])
        if(toplot$direction[i] == '+') toplot$logP[i] <- -log10(toplot$qvalues[i])
    }
    if (nrow(toplot[toplot$logP < -10,])!=0) toplot[toplot$logP < -10,]$logP <- -10
    if (nrow(toplot[toplot$logP > 10,])!=0) toplot[toplot$logP > 10,]$logP <- 10
    
    print("   Generating heatmap...")
    path <- gsub('ORA_qval_vs_lfdr', 'ORA_hm', path)
    plotp <- as.matrix(toplot$logP)
    rownames(plotp) <- rownames
    rdbu_r <- rev(brewer.pal(n=3, "RdBu")) # reversed color palette
    col_fun = colorRamp2(c(-round(max(plotp), 2), 0, round(max(plotp), 2)), rdbu_r)
    pdf(path) # file to save
    p <- Heatmap(plotp, 
        col=col_fun,
        show_row_dend=F,
        row_names_side="left",
        row_names_gp = gpar(fontsize = 9),
        name="log10(qvalues)",
        width=ncol(plotp)*unit(2, "cm"),
        heatmap_legend_param=list(title="log10(qval)", at=c(-round(max(plotp), 2), 0, round(max(plotp), 2)))
    )
    draw(p)
    dev.off()
    path <- gsub("pdf", "tsv", path)
    write.table(toplot, path, sep="\t", quote=F, row.names=T)
}

ora <- function(with_pwy, top, bg, path){
    # Overrepresentation Analysis
    # with_pwy: dataframe of all genes and GO annotations
    # top: dataframe of genes of interest
    # bg: dataframe of genes in background set
    # path: file path and name to save as

    # create contingency table
    cols <- c("PWY", "Gene_top_has_PWY", "Gene_not_top_has_PWY", "Gene_top_no_PWY",
              "Gene_not_top_no_PWY", "direction", "p.val", "odds ratio", "qvalues", "lfdr")
    contingency <- data.frame(matrix(nrow=1, ncol=10))
    colnames(contingency) <- cols

    # fill in contingency table for each gene
    print("   Running ORA...")
    for (pathway in unique(with_pwy$Pathways.of.gene)){
        a <- length(unique(top[which(top$Pathways.of.gene==pathway),]$Accession.1)) # Genes in top features and have `pwy`
        b <- length(unique(bg[which(bg$Pathways.of.gene==pathway),]$Accession.1)) # Genes not in top features and have `pwy`
        c <- length(unique(top$Accession.1)) - a # Genes in top features and do not have `pwy`
        d <- length(unique(bg$Accession.1)) - b # Genes not in top features and do not have `pwy`
        tbl <- matrix(c(a, b, c, d), ncol=2, byrow=TRUE) # gene contingency table
        res <- fisher.test(tbl, alternative="two.sided") # fisher's exact test
        if (a+b!=0){
            if(enrichment(k=a, n=a+c, C=a+b, G=a+b+c+d) >= 1) direction = '+' else direction = '-' # direction of enrichment
            contingency <- rbind(contingency, c(pathway, a, b, c, d, direction, res$p.value, res$estimate, 'NA', 'NA'))
        }
    }
    contingency <- contingency[!is.na(contingency$PWY),] # drop first row with NAs
    
    # Calculate q-values
    contingency$p.val <- as.numeric(contingency$p.val)
    sub <- contingency[!(contingency$p.val==1),] # don't include rows with p_val=1
    
    # add biological process, cellular component, and molecular function info
    if (nrow(sub)!=0){
        print("   Calculating q values...")
        qvals <- qvalue(sub$p.val, lambda=0)
        sub$qvalues <- qvals$qvalues
        sub$lfdr <- qvals$lfdr
        # save contingency table
        sub <- sub[order(sub$qvalues),]
        write.table(sub, paste(path, "tsv", sep=""), sep="\t", quote=F, row.names=F)

        # plot qvalues and lfdr 
        path <- gsub("ORA", "ORA_qval_vs_lfdr", path)
        ggplot(sub, aes(x=qvalues, y=lfdr)) + geom_point() + 
            geom_smooth(method="lm", color="black") +
            labs(x="q-values", y="local FDR") + theme_bw()
        ggsave(paste(path, "pdf", sep=""))

        # Generate heatmaps for significant pathways
        toplot <- sub[which(sub$lfdr < 0.05),] # sig GO terms # sub$qvalues < 0.05 & 
        print(nrow(toplot))
        if (nrow(toplot) != 0) heatmap2(toplot, toplot$PWY, paste(path, "pdf", sep=""))
    }
}

pwy_enrichment <- function(f){
    ### ORA of top gene features
    # add gene information to top snp file
    top <- get_genes(f)

    # drop snps without pathway annotations (some of these may be intergenic)
    top <- top[which(!is.na(top$Pathways.of.gene)),]
    
    # name of importance score file
    name <- str_extract(f, "[A-Z0-9]+_[0-9]+[^_training.csv]") # extract file name
    print(paste("File", name, sep=" "))
    
    ### Step 2: calculate the enrichment score
    bg <- with_pwy[!(with_pwy$Accession.1 %in% top$Accession.1),] # remove top from background set
    bg <- bg[c(2,5:6)] # keep only Accession.1, Product, Pathways.of.gene
    bg <- bg[!duplicated(bg),] # remove duplicates
    print(paste("   Top Genes: ", length(unique(top$Accession.1)), sep=""))
    print(paste("   Genes not in top: ", length(unique(bg$Accession.1)), sep=""))
    print(paste("   Total number of genes is correct: ", length(unique(top$Accession.1))+length(unique(bg$Accession.1))==length(unique(with_pwy$Accession.1))))
    
    ## Overrepresentation Analysis
    path <- paste("Scripts/Genomic_Prediction_RF/SHAP/SNPs/fs/PWY_ORA_", str_extract(f,
        "SHAP_values_sorted_average_[A-Z0-9]+_[0-9]+_training[^txt]"), sep="")
    ora(with_pwy, top, bg, path)
}


# Read in top features' (FS) average SHAP values files
dir <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/SNPs/fs" # path to FS average SHAP files
files <- list.files(path=dir, pattern="^SHAP_values_sorted_average_Y", 
                    full.names=TRUE, recursive=FALSE)

mclapply(X=files, FUN=pwy_enrichment, mc.cores=35) # match go to orfs
# lapply(files, go_enrichment)
