################################################################################
# Go enrichment analysis of top random forest features in each condition       #
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
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(RColorBrewer))

setwd("/mnt/home/seguraab/Shiu_Lab/Project")

#### 2022-09-20T19:04 S288C GO term annotation file from GO database
go <- fread("Data/yeast_GO/sgd.gaf", sep="\t", skip=35, header=F, fill=T)
go <- as.data.frame(go)
colnames(go) <- c("DB", "DB.ID", "Protein", "Qualifier", "GO.ID", "DF.Ref",
    "Evidence", "With.or.From", "Aspect", "Gene.Description", "Gene.Synonym",
    "Type", "Taxon", "Date", "Assigned.By", "Annotation.Extension",
    "Gene.Product") # GAF Format Column Names

# Biological process go term gene counts
tmp <- aggregate(go$Gene.Synonym, by=list(go$GO.ID), FUN=length)
tmp[tmp$Group.1=='GO:0008150',]
#         Group.1    x
# 1899 GO:0008150 1075

# Filter GO terms by evidence code
exp <- c("IDA","IPI", "IMP", "IGI", "IEP", "HDA", "HMP", "HGI", "HEP") # EXP & HTP evidence codes
go <- subset(go, Evidence %in% exp)

# Extract gene name from Gene.Synonym
unique(go$Type)
Gene <- str_extract(go$Gene.Synonym, "[A-Z0-9-]+|")
go <- cbind(go, Gene)

# Extract protein complexes
# tmp <- go[grep("protein_complex", go$Type),]
# cpx <- str_extract(tmp$Gene.Synonym, "[A-Za-z0-9[:punct:]\\+\\) ]+|")
# tmp <- cbind(tmp, cpx)
# colnames(tmp)[18] <- "Gene" # rename to Gene
# Gene <- as.data.frame(Gene) # convert to dataframe
# row.names(Gene) <- row.names(go) # set correct indices
# Gene[as.integer(row.names(tmp)),] # this should match with go, but it doesn't
# # combine with go, how?? indicies between Gene and go are not matching even after setting them

# Add extra GO term information
print("   Grabbing BP, CC, and MF info...")
go$BP <- "" # biological process
go$CC <- "" # cellular component
go$MF <- "" # molecular function
for(i in 1:nrow(go)){
    tryCatch({
        if(!is.null(getGOTerm(go[i,5])$BP[1])) go[i,19] <- getGOTerm(go[i,5])$BP[1]
        if(!is.null(getGOTerm(go[i,5])$CC[1])) go[i,20] <- getGOTerm(go[i,5])$CC[1]
        if(!is.null(getGOTerm(go[i,5])$MF[1])) go[i,21] <- getGOTerm(go[i,5])$MF[1]
    }, error = function(e){print(paste("no GO for ", go[i,5])); NaN},
        finally = {})
}
# write.csv(go, "Data/yeast_GO/sgd_GO_BP.csv", quote=F, row.names=F)

########################### Map GO terms to SNPs ###########################
genes <- read.csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes.txt",
                  sep=",", header=T)
genes <- genes[!(genes$gene=="intergenic"),] # drop intergenic snps
out <- left_join(genes, go, by=c("gene"="Gene"))
with_go <- out[complete.cases(out),] # genes with go 
no_go <- out[is.na(out$GO.ID),] # genes with no go
write.table(with_go, 
    "Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_go.csv", sep="\t", 
    quote=F, row.names=F)
write.table(no_go, 
    "Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_no_go.csv", sep="\t", 
    quote=F, row.names=F)

################################################################################
#                           GO Enrichment Analysis                             #
#                    contingency table for each GO term:                       #
#                                  | ORF in top | ORF not top                  #
#                   -----------------------------------------                  #
#                    In GO term    |            |                              #
#                    Not in GO term|            |                              #
################################################################################

get_genes <- function(f, baseline=F) {
    # add the genes and GO annotations to each SHAP file
    df <- read.delim(f, sep="\t") # read in shap values
    colnames(df) <- c("SNP", "shap_value")
    out <- right_join(with_go, df, by=c("snp"="SNP")) # add gene & go information
    name <- str_extract(f, "SHAP_values_sorted_average_[A-Z0-9]+_[0-9]+_training[^txt]") # extract file name
    name <- gsub("SHAP", "GO_SHAP", name)
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

heatmap <- function(toplot, path){
    if (nrow(toplot) > 1){
        print("   Generating heatmap...")
        plotp <- as.matrix(cbind(c(10,-10), toplot$logP)) # for heatmap
        rownames(plotp) <- toplot$BP
        path <- gsub('ORA', 'ORA_hm', path); path <- gsub('training.', 'training.pdf', path)
        pdf(path) # file to save
        heatmap.2(plotp,
            Rowv=F, # turn off row clustering
            Colv=F, # turn off column clustering
            dendrogram="none", # remove dendrogram
            col = colorRampPalette(c("blue","white","red"))(21), # color palette
            trace="none", # remove trace lines
            margins=c(5,15), # adjust margins (col, row)
            cellnote=cbind(round(as.numeric(toplot$lfdr), 2), round(as.numeric(toplot$odds.ratio), 2)), # cell labels
            notecol="black",
            cexRow=0.8,# row label font size
            labCol=F, # remove column labels
            keysize=0.9, # legend size
            density.info="none", # remove distribution from legend
            key.title="log(q-value)", # legend title
            main=str_extract(path, "[A-Z0-9]+_[a-z]+_[0-9]+[^_training.pdf]"), # figure title
            lmat=rbind(4:3,2:1), # relative position of each element
            lhei=c(0.8,2), # dimensions of display cell heights
            lwid=c(0.8,2))
        dev.off()
        path <- gsub(".pdf", ".tsv", path)
        write.table(toplot, path, sep="\t", quote=F, row.names=T)
    }
}

heatmap2 <- function(toplot, path){
    plotp <- as.matrix(toplot$logP)
    rownames(plotp) <- toplot$BP
    path <- gsub('ORA', 'ORA_hm2', path); path <- gsub('training.', 'training.pdf', path)
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
    path <- gsub(".pdf", ".tsv", path)
    write.table(toplot, path, sep="\t", quote=F, row.names=T)
}

make_contingency <- function(all_orfs){
    cols <- c("GO", "ORF_top_has_GO", "ORF_not_top_has_GO", "ORF_top_no_GO",
              "ORF_not_top_no_GO", "direction", "p.val", "odds ratio", "qvalues", "lfdr")
    contingency <- data.frame(matrix(0, nrow=length(unique(all_orfs$GO.ID)),
        ncol=10, dimnames=list(unique(all_orfs$GO.ID), cols)))
    contingency$GO <- rownames(contingency)

    print("   Grabbing BP, CC, and MF info...")
    contingency$BP <- "" # biological process
    contingency$CC <- "" # cellular component
    contingency$MF <- "" # molecular function
    for(i in 1:nrow(contingency)){
        tryCatch({
            if(!is.null(getGOTerm(contingency[i,1])$BP[1])) contingency[i,11] <- getGOTerm(contingency[i,1])$BP[1]
            if(!is.null(getGOTerm(contingency[i,1])$CC[1])) contingency[i,12] <- getGOTerm(contingency[i,1])$CC[1]
            if(!is.null(getGOTerm(contingency[i,1])$MF[1])) contingency[i,13] <- getGOTerm(contingency[i,1])$MF[1]
        }, error = function(e){print(paste("no GO for ", contingency[i,1])); NaN},
            finally = {})
    }
    return(contingency)
}

ora <- function(with_go, top, bg, contingency, path){
    # Overrepresentation Analysis
    # with_go: dataframe of all genes and GO annotations
    # top2: dataframe of genes of interest
    # bg: dataframe of genes in background set
    # path: file path and name to save as

    # fill in contingency table for each gene
    print("   Running ORA...")
    for (go in unique(with_go$GO)){
        a <- length(unique(top2[which(top2$GO.ID==go),]$gene)) # Genes in top features and have `go`
        b <- length(unique(bg[which(bg$GO.ID==go),]$gene)) # Genes not in top features and have `go`
        c <- length(unique(top2$gene)) - a # Genes in top features and do not have `go`
        d <- length(unique(bg$gene)) - b # Genes not in top features and do not have `go`
        tbl <- matrix(c(a, b, c, d), ncol=2, byrow=TRUE) # gene contingency table
        res <- fisher.test(tbl, alternative="two.sided") # fisher's exact test
        if(enrichment(k=a, n=a+c, C=a+b, G=a+b+c+d) >= 1) direction = '+' else direction = '-' # direction of enrichment
        contingency[go,] <- c(go, a, b, c, d, direction, res$p.value, res$estimate, NA, NA)
    }
    
    # Calculate q-values
    print("   Calculating q values...")
    contingency$p.val <- as.numeric(contingency$p.val)
    sub <- contingency[!(contingency$p.val==1),] # don't include rows with p_val=1
    if (dim(sub)[1] > 1){
        qvals <- qvalue(sub$p.val, lambda=0)
        sub$qvalues <- qvals$qvalues
        sub$lfdr <- qvals$lfdr
        # save contingency table
        out <- rbind(sub, contingency[which(contingency$p.val==1),])
        out <- out[order(out$qvalues),]
        write.table(out, paste(path, "tsv", sep=""), sep="\t", quote=F, row.names=F)

        # Generate heatmaps for significant go terms        
        toplot <- sub[(sub$BP != '' & sub$lfdr < 0.05),] # sig GO terms # sub$qvalues < 0.05 & 
        print(nrow(toplot))
        if (dim(toplot)[1] != 0){
            toplot$logP <- 0
            toplot <- toplot[order(toplot$qvalues, decreasing=T),]
            for(i in 1:nrow(toplot)){ # take the log of the q-values
                if(toplot$direction[i] == '-') toplot$logP[i] <- log10(toplot$qvalues[i])
                if(toplot$direction[i] == '+') toplot$logP[i] <- -log10(toplot$qvalues[i])
            }
            if (nrow(toplot[toplot$logP < -10,])!=0) toplot[toplot$logP < -10,]$logP <- -10
            if (nrow(toplot[toplot$logP > 10,])!=0) toplot[toplot$logP > 10,]$logP <- 10
            print("   Generating heatmap...")
            heatmap2(toplot, path)
        }
        #plotfdr <- as.matrix(cbind(toplot$lfdr, c(1))
    } else {
        write.table(contingency, path, sep="\t", quote=F, row.names=F)
    }
}

go_enrichment <- function(f){
    ### ORA and GSEA of top gene features
    # add gene information to top snp file
    top <- get_genes(f)
    
    # name of importance score file
    name <- str_extract(f, "[A-Z0-9]+_[0-9]+[^_training.csv]") # extract file name
    name <- gsub("_", "_exp_rf_", name)
    print(paste("File", name, sep=" "))
    
    ### Step 1. add mean importance scores to rank genes for GSEA
    # scratch <- "/mnt/gs21/scratch/seguraab/yeast_project/yeast_rf_results"
    # imp <- read.delim(paste(scratch, "/", name, "_imp", sep=""))
    # imp <- imp[order(imp$mean_imp, decreasing=T),] # highest to lowest
    # top2 <- left_join(top, imp, by=c("snp"="X")) # add importance scores
    # top2 <- top2[complete.cases(top2),] # drop rows with no match in with_go
    # # aggregate all snps within the same genes and take the max shap and imp score
    # top2 <- top2 %>% group_by(gene) %>% 
    #     mutate(max_shap=max(shap_value), sd_mean_shap=sd(shap_value), 
    #         max_mean_imp=max(mean_imp), sd_mean_mean_imp=sd(mean_imp)) %>%
    #     as.data.frame()
    # # drop snps and remove duplicates rows
    # top2 <- top2[c(4,9,22:30)] # keep only genes
    # top2 <- top2[!duplicated(top2),] # remove duplicated rows
    # top2 <- top2[order(top2$max_mean_imp, decreasing=T),] # reorder
    # top2$rank <- seq(1:nrow(top2)) # add gene rank THIS IS NOT CORRECT ******** 
    # name <- str_extract(f, "SHAP_values_sorted_average_[A-Z0-9]+_[0-9]+_training[^txt]") # extract file name
    # write.table(top2, 
    #     paste("Scripts/Genomic_Prediction_RF/SHAP/SNPs/fs/GO_genes_", name,"tsv",
    #     sep=""), sep="\t", quote=FALSE, row.names=FALSE)
    
    ### Step 2: calculate the enrichment score
    bg <- with_go[!(with_go$gene %in% top$gene),] # remove top from background set
    bg <- bg[c(4,9,22:24)] # keep only gene, GO.ID, BP, CC, MF
    bg <- bg[!duplicated(bg),] # remove duplicates
    print(paste("   Top Genes: ", length(unique(top$gene)), sep=""))
    print(paste("   Genes not in top: ", length(unique(bg$gene)), sep=""))
    
    ## Overrepresentation Analysis
    save <- paste("Scripts/Genomic_Prediction_RF/SHAP/SNPs/fs/ORA_", str_extract(f,
        "SHAP_values_sorted_average_[A-Z0-9]+_[0-9]+_training[^txt]"), sep="")
    ora(with_go, top, bg, save)
}

# GO mapped to genes
with_go <- read.csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_go.csv", sep="\t")

# Read in top features' (FS) average SHAP values files
dir <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/SNPs/fs" # path to FS average SHAP files
files <- list.files(path=dir, pattern="SHAP_values_sorted_average_Y", 
                    full.names=TRUE, recursive=FALSE)

mclapply(X=files, FUN=go_enrichment, mc.cores=35) # match go to orfs
#lapply(files, go_enrichment)


### Generate a large heatmap combining all SHAP values GO enrichment
