################################################################################
# Go enrichment analysis of top random forest features in each condition       #
################################################################################

rm(list=ls())

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(GSEABase))

setwd("/mnt/home/seguraab/Shiu_Lab/Project")

#### Evidence code-filtered GO term annotation file (see SNP_gene_set_enrichment.R)
go <- read.csv("Data/yeast_GO/sgd_GO_BP.csv")
length(unique(go$GO.ID)) # [1] 4909

#### ORF to gene map file
map <- read.delim("Data/Peter_2018/final_map_orf_to_gene.txt", sep="\t", header=1)

#### Get RF (after FS) feature importance score files
rf_res_orf <- read.csv("Results/RESULTS_RF_ORFs_FS.txt", sep="\t", header=1) # RF FS model results
rf_res_cnv <- read.csv("Results/RESULTS_RF_CNVs_FS.txt", sep="\t", header=1)
dir <- "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results/fs" # directory  to search
orf_files <- unlist(lapply(paste(rf_res_orf$ID, "imp", sep="_"), 
                    function(x){list.files(path=dir, pattern=x, full.names=T)}))
cnv_files <- unlist(lapply(paste(rf_res_cnv$ID, "imp", sep="_"),
                    function(x){list.files(path=dir, pattern=x, full.names=T)}))


### Map ORFs to genes & GO terms
get_genes <- function(f){
    # add the genes and GO annotations to each file
    df <- read.delim(f, sep="\t") # read file
    df$X <- gsub(".", "-", df$X, fixed=T) # replace periods with dashes
    df$X <- gsub("^X", "", df$X) # remove leading X
    out <- right_join(map, df, by=c("orf"="X")) # add gene information
    out <- left_join(out, go[,c('Gene', 'GO.ID', 'BP', 'CC', 'MF')], by=c("gene"="Gene")) # add GO information
    out <- out[order(out$mean_imp, decreasing=T),] # order by importance score
    name <- str_extract(f, "[A-Z0-9]+_[a-z]+_[0-9]+_imp") # extract file name
    name <- paste("Genes_", name, sep="")
    path <- as.character("Scripts/Genomic_Prediction_RF/GO_Enrichment/ORFs_fs")
    # write.table(out, file.path(path, name), "tsv", sep="."), sep="\t", quote=F, row.names=F)
}
mclapply(orf_files, get_genes, mc.cores=35)
mclapply(cnv_files, get_genes, mc.cores=35)

################################################################################
#                           GO Enrichment Analysis                             #
#                    contingency table for each GO term:                       #
#                                  | ORF in top | ORF not top                  #
#                   -----------------------------------------                  #
#                    In GO term    |            |                              #
#                    Not in GO term|            |                              #
################################################################################
### Get files with genes
dir <- "Scripts/Genomic_Prediction_RF/GO_Enrichment/ORFs_fs"
orf_files <- list.files(path=dir, pattern="Genes_[A-Z0-9]+_orf", full.names=T)
cnv_files <- list.files(path=dir, pattern="Genes_[A-Z0-9]+_cno", full.names=T)

### Map all ORFs to GO terms - prep for background set
all_orfs <- t(read.csv("Data/Peter_2018/ORFs_pres_abs.csv", nrows=1, header=F))
all_orfs <- all_orfs[-1,] # drop "ID"
all_orfs <- gsub(".", "-", all_orfs, fixed=T) # replace periods with dashes
all_orfs <- gsub("^X", "", all_orfs) # remove leading X
all_orfs <- left_join(data.frame(all_orfs), map, by=c("all_orfs"="orf")) # add gene information
all_orfs <- left_join(all_orfs, go[,c('Gene', 'GO.ID', 'BP', 'CC', 'MF')],
                      by=c("gene"="Gene"), relationship="many-to-many") # add GO information
colnames(all_orfs) <- c("orf", "gene", "organism", "GO", "BP", "CC", "MF")
remove(go)

### Enrichment analysis functions
enrichment <- function(k, n, C, G){ 
    # Determine the direction of enrichment
    # if >= 1: + (overrepresented)
    # if < 1: - (underrepresented)
    # k: number of ORFs in cluster with GO
    # n: total number of ORFs in cluster
    # C: total number of ORFs (in cluster + background) with GO
    # G: total number of ORFs (in cluster + background)
    return((k/C)/(n/G))
}

ora <- function(f, all_orfs, path){
    # Overrepresentation Analysis
    # all_orfs: dataframe of all ORFs and GO annotations
    # top2: dataframe of ORFs of interest
    # path: file path and name to save as

    # read in top ORF feature file
    top <- read.delim(f, sep="\t")
    colnames(top) <- c("orf", "gene", "organism", "mean_imp", "GO", "BP", "CC", "MF")
    
    # get background set
    bg <- all_orfs[!(all_orfs$orf %in% top$orf),] # remove top from background set
    print(paste("   Top ORFs: ", length(unique(top$orf)), sep=""))
    print(paste("   ORFs not in top: ", length(unique(bg$orf)), sep=""))
    print(paste("   Top genes: ", length(unique(top$gene)), sep=""))
    print(paste("   Genes not in top: ", length(unique(bg$gene)), sep=""))
    print(paste("   Total number of ORFs is correct: ", length(unique(top$orf))+length(unique(bg$orf))==length(unique(all_orfs$orf))))
    
    # make contingency table for overrepresntatino analysis
    cols <- c("GO", "ORF_top_has_GO", "ORF_not_top_has_GO", "ORF_top_no_GO",
              "ORF_not_top_no_GO", "direction", "p.val", "odds ratio", "qvalues")
    contingency <- data.frame(matrix(nrow=1, ncol=9))
    colnames(contingency) <- cols

    # fill in contingency table for each gene
    print("   Running ORA...")
    for (go in unique(all_orfs$GO)){
        if (!is.na(go)){
            a <- length(unique(top[which(top$GO==go),]$orf)) # ORFs in top features and have `go`
            b <- length(unique(bg[which(bg$GO==go),]$orf)) # ORFs not in top features and have `go`
            c <- length(unique(top$orf)) - a # ORFs in top feature and do not have `go`
            d <- length(unique(bg$orf)) - b # ORFs not in top features and do not have `go`
            tbl <- matrix(c(a, b, c, d), ncol=2, byrow=TRUE) # gene contingency table
            res <- fisher.test(tbl, alternative="two.sided") # fisher's exact test
            if (a+b!=0){
                if(enrichment(k=a, n=a+c, C=a+b, G=a+b+c+d) >= 1) direction = '+' else direction = '-' # overrepresentation analysis
                contingency <- rbind(contingency, list(go, a, b, c, d, direction, res$p.value, res$estimate, NA))
            }
        }
    }
    contingency <- contingency[!is.na(contingency$GO),] # drop first row with NAs

    # Calculate q-values
    sub <- contingency[!(contingency$p.val==1),] # don't include rows with p_val=1
    
    # add biological process, cellular component, and molecular function info
    if (nrow(sub)!=0){
        print("   Grabbing BP, CC, and MF info...")
        sub$BP <- "" # biological process
        sub$CC <- "" # cellular component
        sub$MF <- "" # molecular function
        for(i in 1:nrow(sub)){
            tryCatch({
                if(!is.null(getGOTerm(sub[i,1])$BP[1])) sub[i,10] <- getGOTerm(sub[i,1])$BP[1]
                if(!is.null(getGOTerm(sub[i,1])$CC[1])) sub[i,11] <- getGOTerm(sub[i,1])$CC[1]
                if(!is.null(getGOTerm(sub[i,1])$MF[1])) sub[i,12] <- getGOTerm(sub[i,1])$MF[1]
            }, error = function(e){print(paste("no GO for ", sub[i,1])); NaN},
                finally = {})
        }
        print("   Calculating q values...")
        qvals <- p.adjust(sub$p.val, method="BH")
        sub$qvalues <- qvals
        
        # save contingency table
        sub <- sub[order(sub$qvalues),]
        write.table(sub, path, sep="\t", quote=F, row.names=F)
    }
}


go_enrichment <- function(f){
    # Run GO Enrichment
    save <-  gsub("Genes_", "ORA_Genes_", f)
    ora(f, all_orfs, save)
}

mclapply(orf_files, go_enrichment, mc.cores=35)
mclapply(cnv_files, go_enrichment, mc.cores=35)