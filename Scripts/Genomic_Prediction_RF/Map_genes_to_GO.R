#
# Map genes represented in SNP features to GO terms
#
# Kenia E. Segura Abá

library('GOstats')
library('GSEABase')
library(tidyverse)
library('parallel')

# Go term data
go <- read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/yeast_GO/genes_go_sgd.tsv", sep="\t")
newdf <- go %>% group_by(Gene) %>% summarise_at("GO", toString)
write.csv(newdf, "/mnt/home/seguraab/Shiu_Lab/Project/Data/yeast_GO/genes_go_sgd_collapsed.csv", quote=FALSE, row.names=FALSE)

# SNP features mapped to genes
genes <- read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes.txt", sep=",", header=FALSE)

# Feature importance score datasets for each condition
dir <- "/mnt/home/seguraab/Shiu_Lab/Project/yeast_RF_results/SNPs_as_Features"
files <- list.files(path=dir, pattern="_rf_baseline_imp", full.names=TRUE, recursive=FALSE)
x=files[1]
# Function to btain biological process (BP), cell component(CC) and, molecular function (MF) info
get_GO <- function(dat){
    for(i in 1:nrow(dat)){
        tryCatch( {
            if(!is.null(getGOTerm(dat[i,"GO"])$BP[1])) 
            {dat[i,"BP"] <- getGOTerm(dat[i,"GO"])$BP[1]}
            if(!is.null(getGOTerm(dat[i,"GO"])$CC[1])) 
            {dat[i,"CC"] <- getGOTerm(dat[i,"GO"])$CC[1]}
            if(!is.null(getGOTerm(dat[i,"GO"])$MF[1])) 
            {dat[i,"MF"] <- getGOTerm(dat[i,"GO"])$MF[1]}
        },
        error = function(e) {print(paste("no GO for ",dat[i,"GO"]));NaN}, finally = {})
    }
}

# Function to combine elements in a list into a string
list_to_str <- function(multi_lvl_list){
    res <- c()
    for (i in 1:length(multi_lvl_list)){ 
        res[i] <- paste(multi_lvl_list[[i]], collapse=", ")
    } 
    return(res) 
}

# Function to merge data and call get_GO
merge_data <- function(x){
    # Read in files
    f <- read.table(x, header=TRUE, sep="\t")
    name <- str_extract(x, pattern="[:alnum:]+_rf_baseline_imp") # file name
    trait <- strsplit(name, "_")[[1]][1] # trait name
    # Merge datasets
    out <- right_join(x=genes, y=f, by=c("V1"="X")) # by feature name
    out2 <- left_join(x=out, y=newdf, by=c("V4"="Gene")) # by gene name
    test <- strsplit(out2$GO, ", ") # separate GO terms in each cell
    test <- lapply(test, unique) # remove duplicate GO terms
    out2$GO <- list_to_str(test) # replace duplicate GO terms with unique GO terms
    out3 <- left_join(x=out, y=go, by=c("V4"="Gene"))
    out3 <- distinct(out3) # remove duplicate rows
    colnames(out2) <- c("Feature", "Chr", "Pos", "Gene", "mean_imp", "GO")
    colnames(out3) <- c("Feature", "Chr", "Pos", "Gene", "mean_imp", "GO", "SGD ID", "GO Description", "UniProt ID")
    # Obtain biological process (BP), cell component(CC) and, molecular function (MF) info
    out3 <- cbind(out3,"BP"="")
    out3 <- cbind(out3,"CC"="")
    out3 <- cbind(out3,"MF"="")
    out3 <- get_GO(out3)
    write.table(out2, paste(dir, "/", trait, "_rf_baseline_go_imp", sep=""), quote=FALSE, row.names=FALSE)
    write.table(out3, paste(dir, "/", trait, "_rf_baseline_go_expanded_imp", sep=""), quote=FALSE, row.names=FALSE)
}

# Add GO term/gene information to feature importance score data
system.time(tmp <- mclapply(X=files, FUN=merge_data, mc.cores=35))