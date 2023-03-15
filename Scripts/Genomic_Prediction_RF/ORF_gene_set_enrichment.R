################################################################################
# Go enrichment analysis of top random forest features in each condition       #
################################################################################

rm(list=ls())

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

setwd("/mnt/home/seguraab/Shiu_Lab/Project")

#### Evidence code-filtered GO term annotation file (SNP_gene_set_enrichment.R)
go <- read.csv("Data/yeast_GO/sgd_GO_BP.csv")
# length(unique(go$GO.ID))
# [1] 4909

########################### Read Blast Output Files ###########################
#### Filter ORF features mapped to S288C genes file
s288c_blastx <- read.delim("Data/S288C_reference_genome_R64-3-1_20210421/S288C_orf_peter_blastx.txt", header=TRUE, skip=1)
strict <- s288c_blastx %>% group_by(qacc) %>% 
    dplyr::filter(pident == max(pident) & pident >= 95 & evalue == min(evalue)) %>% 
    as.data.frame()
strict$gene <- strict$sacc
loose <- s288c_blastx %>% group_by(qacc) %>% 
    dplyr::filter(pident == max(pident) & pident >= 30 & evalue == min(evalue) & evalue <= 1e-5) %>% 
    as.data.frame()

# find duplicates
qdupes_strict <- as.data.frame(strict[duplicated(strict[,c('qacc')]) | duplicated(strict[,c('qacc')], fromLast=T),])
qdupes_loose <- as.data.frame(loose[duplicated(loose[,c('qacc')]) | duplicated(loose[,c('qacc')], fromLast=T),])

#### Filter ORF features mapped to Arabidopsis thaliana genes
tair10_blastx <- read.delim("Data/Arabidopsis_Genome_TAIR10.1/peter_orf_TAIR10_blastx.txt", header=TRUE)
tair10_strict <- tair10_blastx %>% group_by(qacc) %>% 
    dplyr::filter(pident == max(pident) & pident >= 95 & evalue == min(evalue)) %>% 
    as.data.frame()
tair10_strict$gene <- tair10_strict$sacc

#### Filter ORF features mapped to Drosophila melanogaster genes
iso1mt_blastx <- read.delim("Data/Drosophila_Genome_R6_ISO1MT/peter_orf_ISO1MT_blastx.txt", header=TRUE)
iso1mt_strict <- iso1mt_blastx %>% group_by(qacc) %>% 
    dplyr::filter(pident == max(pident) & pident >= 95 & evalue == min(evalue)) %>% 
    as.data.frame()
iso1mt_strict$gene <- iso1mt_strict$sacc

#### Filter ORF features mapped to Human genes
human_blastx <- read.delim("Data/Human_Genome_GRCh38.p14/peter_orf_human_blastx.txt", header=TRUE)
human_strict <- human_blastx %>% group_by(qacc) %>% 
    dplyr::filter(pident == max(pident) & pident >= 95 & evalue == min(evalue)) %>% 
    as.data.frame()
human_strict$gene <- human_strict$sacc

#### Filter ORF features mapped to Neurospora crossa genes
nc12_blastx <- read.delim("Data/Neurospora_OR74A_Genome_NC12/peter_orf_NC12_blastx.txt", header=TRUE)
nc12_strict <- nc12_blastx %>% group_by(qacc) %>% 
    dplyr::filter(pident == max(pident) & pident >= 95 & evalue == min(evalue)) %>% 
    as.data.frame()
nc12_strict$gene <- nc12_strict$sacc

#### Filter ORF features mapped to non-redundant database
nr_blastx <- read.delim("Data/BLAST_nr_db/peter_orf_nr_blastx.txt", header=TRUE)
nr_strict <- nr_blastx %>% group_by(qacc) %>% 
    dplyr::filter(pident == max(pident) & pident >= 95 & evalue == min(evalue)) %>% 
    as.data.frame()
# write.csv(unique(nr_strict$sacc), "Data/BLAST_nr_db/nr_genes_for_sgd.txt",
#     quote=FALSE, row.names=FALSE)

### SGD and NCBI nr_strict gene matches
sgd <- read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/BLAST_nr_db/nr_genes_sgd_matches.csv", header=1)
ncbi <- read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/BLAST_nr_db/nr_genes_ncbi_matches.csv", header=1)
colnames(sgd)[6] <- "sacc" # rename "secondaryIdentifier"
colnames(ncbi)[2] <- "sacc" # rename "Symbol"
nr_strict_matches <- rbind.fill(sgd, ncbi)

# get original ORF qacc values
tmp <- nr_strict[(nr_strict$sacc %in% nr_strict_matches$input) | (nr_strict$sacc %in% nr_strict_matches$Protein),]
nr_strict_matches <- right_join(tmp, nr_strict_matches, by=c("sacc"="input"))
colnames(nr_strict_matches)[19] <- "gene"# rename sacc.y
#unique(nr_strict_matches$Organism) # may need to get GO terms if available

####################### Match Annotations to SHAP Files #######################
### Ensure ORF IDs match in the blast and shap files
match_id <- function(df) {
    df$qacc <- gsub("-", ".", df$qacc)
    df$qacc <- sub("^", "X", df$qacc)
    return(df)
}
strict <- match_id(strict)
tair10_strict <- match_id(tair10_strict)
iso1mt_strict <- match_id(iso1mt_strict)
human_strict <- match_id(human_strict)
nc12_strict <- match_id(nc12_strict)
nr_strict_matches <- match_id(nr_strict_matches)


### Loop through SHAP files for each condition and add gene information
get_genes <- function(f) {
    df <- read.delim(f, sep="\t") # read in shap files
    colnames(df) <- c("X", "shap_value")
    #df2 <-  filter(df, if_all(everything(), ~ .x != 0)) # only keep features with nonzero SHAP values
    
    # Match to S288C genes
    df2 <- right_join(strict, df, by=c("qacc"="X"))
    missing <- df2[is.na(df2$sacc),c(1,17)] # ORFs that did not match to S288C genes
    df2 <- df2[!is.na(df2$sacc),] # matches only

    # Match missing ORFs to genes in other species
    if (dim(missing)[1]!=0){
        # Match to Arabidopsis thaliana genes
        match_tair10 <- right_join(tair10_strict, missing, by="qacc")
        missing2 <- match_tair10[is.na(match_tair10$gene),c(1,14)]
        match_tair10 <- match_tair10[!is.na(match_tair10$gene),]

        # Match to Drosophila melanogaster genes
        match_iso1mt <- right_join(iso1mt_strict, missing2, by="qacc")
        missing2 <- match_iso1mt[is.na(match_iso1mt$gene),c(1,14)]
        match_iso1mt <- match_iso1mt[!is.na(match_iso1mt$gene),]

        # Match to Human genes
        match_human <- right_join(human_strict, missing2, by="qacc")
        missing2 <- match_human[is.na(match_human$gene),c(1,14)]
        match_human <- match_human[!is.na(match_human$gene),]

        # Match to Neurospora crossa genes
        match_nc12 <- right_join(nc12_strict, missing2, by="qacc")
        missing2 <- match_nc12[is.na(match_nc12$gene),c(1,14)]
        match_nc12 <- match_nc12[!is.na(match_nc12$gene),]

        # Match to non-redundant database genes
        match_nr <- right_join(nr_strict_matches, missing2, by="qacc")
        missing2 <- match_nr[is.na(match_nr$gene),c(1,30)]
        match_nr <- match_nr[!is.na(match_nr$gene),]
        
        # Combine matches to df2
        all_match <- rbind.fill(match_tair10, 
            rbind.fill(match_iso1mt, 
            rbind.fill(match_human, 
            rbind.fill(match_nc12, 
            rbind.fill(match_nr, df2[!is.na(df2$gene),])))))
        return(all_match)
    } else {
        return(df2)
    }
}

### Add the GO terms to each SHAP file
get_go <- function(f, baseline=F){
    # Map ORFs to genes and extract file name
    genes <- get_genes(f)
    name <- str_extract(f, "SHAP_values_sorted_average_[A-Z0-9]+_[a-z]+_[0-9]+_training[^txt]")
    
    # Match to S288C GO annotations
    go_match <- left_join(genes, go, by=c('gene'='Gene'))
    go_match_sub <- go_match[c("qacc", "sacc", "gene", "shap_value", "GO.ID", "BP")]

    # collapse rows for each qacc so that all the go terms and genes are in one cell.
    if (baseline==F){
        write.table(go_match_sub, 
            paste("Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs/GO_", name,"tsv",
            sep=""), sep="\t", quote=FALSE, row.names=FALSE)
    } else  
        return(go_match_sub)
}

# Read in top features' (FS) average SHAP values files
dir <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs" 
files <- list.files(path=dir, pattern="^SHAP_values_sorted_average_Y", full.names=TRUE, 
                    recursive=FALSE)
mclapply(X=files, FUN=get_go, mc.cores=4) # match go to orfs

# Read in list of ORFs
dir <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/ORFs/"
# path <- paste(dir,
#     "baseline/SHAP_values_sorted_average_YPACETATE_orf_training.txt", sep="")
# somehow during debugging the YPACETE orf and cno files were rewritten
path <- paste(dir, 
    "baseline/SHAP_values_sorted_average_YPD6AU_orf_training.txt", sep="")
all_orfs <- get_go(path, baseline=T)
map <- all_orfs[1:3]
map <- map[!duplicated(map),]
# write.table(map, "Data/Peter_2018/ORFs_gene_map.tsv", sep="\t", row.names=F, quote=F)

# remove ORF-gene pairs with missing GO annotations
all_orfs <- all_orfs[!duplicated(all_orfs),] # remove duplicate rows
unique(all_orfs[which(is.na(all_orfs$GO.ID)),]$gene) # 778 unique genes with no GO terms
dim(all_orfs[which(is.na(all_orfs$GO.ID)),]) # [1] 830   6
# write.table(all_orfs, "Data/Peter_2018/all_orfs_with_go.tsv", sep="\t", row.names=F, quote=F)

# save the genes with no GO terms and manually search on SGD database (see 2023_Lab_Notebook)
write.table(unique(all_orfs[which(is.na(all_orfs$GO.ID)),]$gene), "Data/Peter_2018/ORFs_no_GO_gene_list.txt", row.names=F, quote=F)

# Read in SGD results and combine with all_orfs
sgd <- read.csv("Data/Peter_2018/ORFs_no_GO_gene_list_sgd.tsv", sep="\t")
# Filter GO terms by evidence code
exp <- c("IDA","IPI", "IMP", "IGI", "IEP", "HDA", "HMP", "HGI", "HEP") # EXP & HTP evidence codes
dim(sgd) # [1] 15867     8
sgd <- subset(sgd, Evidence.Code %in% exp)
dim(sgd) # [1] 2544    8
length(unique(sgd$gene)) # [1] 332
dim(all_orfs[which(all_orfs$gene %in% sgd$gene),]) # [1] 364   6
tmp <- left_join(all_orfs[which(all_orfs$gene %in% sgd$gene),], sgd[c("gene", "GO.ID")], by="gene")
sum(is.na(tmp$GO.ID.x)) # [1] 2833
sum(is.na(tmp$GO.ID.y)) # [1] 0
tmp <- tmp[-5] # drop "GO.ID.x"
tmp <- tmp[c("qacc", "sacc", "gene", "shap_value","GO.ID.y", "BP")] # reorder
colnames(tmp) <- c("qacc", "sacc", "gene", "shap_value","GO.ID", "BP")
all_orfs2 <- rbind.fill(all_orfs, tmp)
dim(all_orfs2) # [1] 34578     6
all_orfs2 <- all_orfs2[!duplicated(all_orfs2),]
dim(all_orfs2) # [1] 34066     6
length(unique(tmp$qacc)) # [1] 362
sum(is.na(all_orfs2$GO.ID)) # [1] 830 (this means that these missings were kept, and tmp just concatenated those we found without replacing NAs
length(unique(all_orfs$gene)) # [1] 5758
length(unique(all_orfs2$gene)) # [1] 5758
length(unique(all_orfs2$GO.ID)) # [1] 4710
length(unique(all_orfs$GO.ID)) # [1] 4490
`%notin%` = Negate(`%in%`)
dim(all_orfs2[which(all_orfs2$gene %in% sgd$gene & is.na(all_orfs2$GO.ID)),]) # [1] 364   6  delete these!
dim(all_orfs2[which(all_orfs2$gene %in% sgd$gene & !is.na(all_orfs2$GO.ID)),]) # [1] 2321    6 using head you can see the 364 might be there
lookfor <- all_orfs2[which(all_orfs2$gene %in% sgd$gene & is.na(all_orfs2$GO.ID)),]$gene
which(all_orfs2[which(all_orfs2$gene %in% sgd$gene & is.na(all_orfs2$GO.ID)),]$gene %in% all_orfs2[which(all_orfs2$gene %in% sgd$gene & !is.na(all_orfs2$GO.ID)),]$gene) # all 364 have GO.IDs
all_orfs2 <- all_orfs2[-which(all_orfs2$gene %in% sgd$gene & is.na(all_orfs2$GO.ID)),] # drop the rows
dim(all_orfs2) # [1] 33702     6 (34066-33702=364)

# save ORFs mapped to genes and GO terms
write.table(all_orfs2[-4], "Data/Peter_2018/ORFs_GO_genes.tsv", sep="\t", quote=F, row.names=F) 
all_orfs <- read.csv("Data/Peter_2018/ORFs_GO_genes.tsv", sep="\t")
with_go <- all_orfs[!is.na(all_orfs$GO.ID),] # should I not do this? Should genes with no GO should be included in the backgroundset for GO enrichment? Thilanka says no since we don't know what the GO term could be.
sum(is.na(with_go$GO.ID)) # 0
sum(is.na(with_go$gene)) # 0, BP has missing values
write.table(with_go, "Data/Peter_2018/ORFs_with_GO.tsv", sep="\t", row.names=F, quote=F) # save only ORFs with GO annotations

################################################################################
#                           GO Enrichment Analysis                             #
#                    contingency table for each GO term:                       #
#                                  | ORF in top | ORF not top                  #
#                   -----------------------------------------                  #
#                    In GO term    |            |                              #
#                    Not in GO term|            |                              #
################################################################################

# all ORFs and GO annotations
with_go <- read.csv('Data/Peter_2018/ORFs_with_GO.tsv', sep="\t", header=T) 
colnames(with_go) <- c("qacc", "sacc", "gene", "GO", "BP")

get_genes <- function(f) {
    # add the genes and pathway annotations to each SHAP file
    df <- read.delim(f, sep="\t") # read in shap values
    colnames(df) <- c("ORF", "shap_value")
    out <- right_join(with_go, df, by=c("qacc"="ORF")) # add gene & pwy information
    name <- str_extract(f, "SHAP_values_sorted_average_Y[A-Z0-9_a-z]+_training[^txt]") # extract file name
    name <- gsub("SHAP", "GO_SHAP", name)
    path <- as.character("Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs")
    write.csv(out, paste(file.path(path, name), "csv", sep=""), quote=F, row.names=F)
    return(out)
}

enrichment <- function(k, n, C, G){ 
    # determine direction of enrichment
    # if >= 1: + (overrepresented)
    # if < 1: - (underrepresented)
    # k: number of ORFs in cluster with GO
    # n: total number of ORFs in cluster
    # C: total number of ORFs (in cluster + background) with GO
    # G: total number of ORFs (in cluster + background)
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

heatmap2 <- function(toplot, rownames, path){
    toplot$logP <- 0
    toplot <- toplot[order(toplot$qvalues),]

    for(i in 1:nrow(toplot)){ # take the log of the q-values
        if(toplot$direction[i] == '-') toplot$logP[i] <- log10(toplot$qvalues[i])
        if(toplot$direction[i] == '+') toplot$logP[i] <- -log10(toplot$qvalues[i])
    }
    if (nrow(toplot[toplot$logP < -10,])!=0) toplot[toplot$logP < -10,]$logP <- -10
    if (nrow(toplot[toplot$logP > 10,])!=0) toplot[toplot$logP > 10,]$logP <- 10
    
    print("   Generating heatmap...")
    path <- gsub('ORA', 'ORA_hm2', path)
    plotp <- as.matrix(toplot$logP)
    rownames(plotp) <- rownames
    palette_func <- colorRampPalette(c('red', 'white', 'blue'))
    palette <- rev(palette_func(3)) # reversed color palette
    col_fun = colorRamp2(c(-round(max(plotp), 2), 0, round(max(plotp), 2)), palette)
    pdf(path) # file to save
    p <- Heatmap(plotp, 
        col=col_fun,
        show_row_dend=F,
        row_names_side="left",
        row_names_gp = gpar(fontsize = 9),
        name="log10(qvalues)",
        width=ncol(plotp)*unit(2, "cm"),
        heatmap_legend_param=list(title="log10(qval)", at=c(-round(max(plotp), 2), 0, round(max(plotp), 2))))
    draw(p)
    dev.off()
    path <- gsub("pdf", "tsv", path)
    write.table(toplot, path, sep="\t", quote=F, row.names=T)
}

ora <- function(with_go, top, bg, path){
    # Overrepresentation Analysis
    # with_go: dataframe of all ORFs and GO annotations
    # top2: dataframe of ORFs of interest
    # bg: dataframe of ORFs in background set
    # path: file path and name to save as

    # make contingency table for overrepresntatino analysis
    cols <- c("GO", "ORF_top_has_GO", "ORF_not_top_has_GO", "ORF_top_no_GO",
              "ORF_not_top_no_GO", "direction", "p.val", "odds ratio", "qvalues", "lfdr")
    contingency <- data.frame(matrix(nrow=1, ncol=10))
    colnames(contingency) <- cols

    # fill in contingency table for each gene
    print("   Running ORA...")
    for (go in unique(with_go$GO)){
        if (!is.na(go)){
            a <- length(unique(top[which(top$GO==go),]$qacc)) # ORFs in top features and have `go`
            b <- length(unique(bg[which(bg$GO==go),]$qacc)) # ORFs not in top features and have `go`
            c <- length(unique(top$qacc)) - a # ORFs in top feature and do not have `go`
            d <- length(unique(bg$qacc)) - b # ORFs not in top features and do not have `go`
            tbl <- matrix(c(a, b, c, d), ncol=2, byrow=TRUE) # gene contingency table
            res <- fisher.test(tbl, alternative="two.sided") # fisher's exact test
            if (a+b!=0){
                if(enrichment(k=a, n=a+c, C=a+b, G=a+b+c+d) >= 1) direction = '+' else direction = '-' # overrepresentation analysis
                contingency <- rbind(contingency, c(go, a, b, c, d, direction, res$p.value, res$estimate, NA, NA))
            }
        }
    }
    contingency <- contingency[!is.na(contingency$GO),] # drop first row with NAs

    # Calculate q-values
    contingency$p.val <- as.numeric(contingency$p.val) 
    sub <- contingency[!(contingency$p.val==1),] # don't include rows with p_val=1
    
    # add biological process, cellular component, and molecular function info
    if (nrow(sub)!=0){
        print("   Grabbing BP, CC, and MF info...")
        sub$BP <- "" # biological process
        sub$CC <- "" # cellular component
        sub$MF <- "" # molecular function
        for(i in 1:nrow(sub)){
            tryCatch({
                if(!is.null(getGOTerm(sub[i,1])$BP[1])) sub[i,11] <- getGOTerm(sub[i,1])$BP[1]
                if(!is.null(getGOTerm(sub[i,1])$CC[1])) sub[i,12] <- getGOTerm(sub[i,1])$CC[1]
                if(!is.null(getGOTerm(sub[i,1])$MF[1])) sub[i,13] <- getGOTerm(sub[i,1])$MF[1]
            }, error = function(e){print(paste("no GO for ", sub[i,1])); NaN},
                finally = {})
        }
        print("   Calculating q values...")
        qvals <- qvalue(sub$p.val, lambda=0)
        sub$qvalues <- qvals$qvalues
        sub$lfdr <- qvals$lfdr
        # save contingency table
        sub <- sub[order(sub$qvalues),]
        write.table(sub, paste(path, "tsv", sep=""), sep="\t", quote=F, row.names=F)

        # plot qvalues and lfdr 
        path2 <- gsub("ORA", "ORA_qval_vs_lfdr", path)
        ggplot(sub, aes(x=qvalues, y=lfdr)) + geom_point() + 
            geom_smooth(method="lm", color="black") +
            labs(x="q-values", y="local FDR") + theme_bw()
        ggsave(paste(path2, "pdf", sep=""))

        # Generate heatmaps for significant go terms        
        # Biological process
        toplot <- sub[(sub$BP != '' & sub$lfdr < 0.05),] # sig GO terms # sub$qvalues < 0.05 & 
        print(nrow(toplot))
        if (dim(toplot)[1] != 0) heatmap2(toplot, toplot$BP, paste(path, "BP.pdf", sep=""))
        # Cellular Component
        toplot <- sub[(sub$CC != '' & sub$lfdr < 0.05),] # sig GO terms # sub$qvalues < 0.05 & 
        print(nrow(toplot))
        if (dim(toplot)[1] != 0) heatmap2(toplot, toplot$CC, paste(path, "CC.pdf", sep=""))
        # Molecular Function
        toplot <- sub[(sub$MF != '' & sub$lfdr < 0.05),] # sig GO terms # sub$qvalues < 0.05 & 
        print(nrow(toplot))
        if (dim(toplot)[1] != 0) heatmap2(toplot, toplot$MF, paste(path, "MF.pdf", sep=""))
    }
}

topORFs <- function(group){return(group!=0)} # function to filter ORFs for topGO

gsea <- function(top2){
    ### Gene Set Enrichment Analysis with TopGO
    # make topGOdata object
    allORFs <- top2[c(1,7)] # ORF identifiers + importance scores
    allORFs <- rbind.fill(allORFs, bg)[1:2] # add background set
    allORFs <- allORFs[!duplicated(allORFs),] # drop duplicates
    allORFs[is.na(allORFs$mean_imp),]$mean_imp <- 0 # set background ORFs to 0
    allORFs$group <- ifelse(allORFs$mean_imp > 0, 1, 0) # assign groups (0 is bg, 1 is top)
    orfList <- setNames(allORFs$group, as.character(allORFs$qacc)) # named vector
    # Note: some top2 ORFs have an imp score of 0, these are now part of the background set
    
    # Note: How to incorporate orf ranks when imp score is 0 for most?
    GOdata <- new('topGOdata', 
        ontology="BP", # ontology of interest
        allGenes=orfList, # named vector of ORF identifiers and group assignments
        geneSel=topORFs, # function to get interesting ORFs
        annot=annFUN.file , # function to map ORFs to GO terms
        file="Data/Peter_2018/ORFs_GO_genes.map" # ORF GO mappings
    )
    print(GOdata)

    # run enrichment analysis
    classic.fisher <- runTest(GOdata, algorithm="classic", statistic="fisher") # Classic fisher's exact test
    weight01.fisher <- runTest(GOdata, algorithm="weight01", statistic="fisher") # Weighted fisher's exact test
    elim.ks <- runTest(GOdata, algorithm="elim", statistic="ks") # Kolmogorov-Smirnov statistic

    # results
    all_res <- GenTable(GOdata, classicFisher=classic.fisher,
        weightFisher=weight01.fisher, ks = elim.ks,
        orderBy='classicFisher', ranksOf='classicFisher',
        topNodes=length(usedGO(object=GOdata)), # number of nodes
        numChar=1000) # avoid truncation of string
    all_res <- left_join(all_res, go2orfs, by=c("GO.ID"="GO")) # add ORFs
    save <- paste("Scripts/Genomic_Prediction_RF/SHAP/topGO_", str_extract(f,
            "SHAP_values_sorted_average_[A-Z0-9]+_[a-z]+_[0-9]_training.[^csv]"),
            "tsv", sep="")
    write.csv(all_res, save, sep="\t", quote=F, row.names=F)
}

go_enrichment <- function(f){
    # ORA and GSEA of top ORF features
    
    # read in top ORF feature file
    top <- get_genes(f)
    
    # drop orfs without GO annotations
    top <- top[which(!is.na(top$GO)),]
    
    # name of importance score file
    name <- str_extract(f, "[A-Z0-9]+_[a-z]+_[0-9]+[^_training.csv]") # extract file name
    print(paste("File", name, sep=" "))
    
    ### Step 1. add mean importance scores to rank orfs (for GSEA)
    # scratch <- "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results"
    # imp <- read.delim(paste(scratch, "/", name, "_imp", sep=""))
    # imp <- imp[order(imp$mean_imp, decreasing=T),] # highest to lowest
    # imp$rank <- seq(1:nrow(imp)) # add rank
    # top2 <- left_join(top, imp, by=c("qacc"="X")) # add importance scores
    # top2 <- top2[!duplicated(top2),] # remove duplicated rows
    # top2 <- top2[order(top2$mean_imp, decreasing=T),] # reorder

    ### Step 2: calculate the enrichment score
    bg <- with_go[!(with_go$qacc %in% top$qacc),] # remove top from background set
    #bg$GO <- bg$GO %>% replace_na("none") # replace NAs
    #bg$GO.Description[bg$GO.Description==""] <- "none" # replace ""
    print(paste("   Top ORFs: ", length(unique(top$qacc)), sep=""))
    print(paste("   ORFs not in top: ", length(unique(bg$qacc)), sep=""))
    print(paste("   Top genes: ", length(unique(top$gene)), sep=""))
    print(paste("   Genes not in top: ", length(unique(bg$gene)), sep=""))
    print(paste("   Total number of ORFs is correct: ", length(unique(top$qacc))+length(unique(bg$qacc))==length(unique(with_go$qacc))))
    
    ## Overrepresentation Analysis
    path <- paste("Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs/ORA_", str_extract(f,
        "SHAP_values_sorted_average_[A-Z0-9]+_[a-z]+_[0-9]+_training[^txt]"), sep="")
    ora(with_go, top, bg, path)

    ## Gene Set Enrichment Analysis
    #gsea(top2)
}

## create topGO mapping file
# map <- all_orfs %>% group_by(qacc) %>% dplyr::summarize(GO = toString(GO))
# write.table(map, "Data/Peter_2018/ORFs_GO_genes.map", sep="\t", quote=F, row.names=F, col.names=F)

# group ORFs by GO
go2orfs <- with_go %>% group_by(GO) %>% dplyr::summarize(ORFs=toString(qacc))
write.table(go2orfs, "Data/Peter_2018/ORFs_with_GO_reduced_map.tsv", sep="\t", quote=F, row.names=F)

## top features' (FS) average SHAP values files w/GO annotations for each env
dir <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs" # path to FS average SHAP files
files <- list.files(path=dir, pattern="^SHAP_values_sorted_average_Y", 
                    full.names=TRUE, recursive=FALSE)
mclapply(X=files, FUN=go_enrichment, mc.cores=35) # match go to orfs

##################### Heatmaps combining all environments #####################
envs <- c("YPACETATE", "YPD6AU", "YPD14", "YPD40", "YPD42", "YPDANISO10", "YPDANISO20", 
    "YPDANISO50", "YPDBENOMYL200", "YPDBENOMYL500", "YPDCAFEIN40", "YPDCAFEIN50", 
    "YPDCHX1", "YPDCHX05", "YPDCUSO410MM", "YPDDMSO", "YPDETOH", "YPDFLUCONAZOLE", 
    "YPDFORMAMIDE4", "YPDFORMAMIDE5", "YPDHU", "YPDKCL2M", "YPDLICL250MM", "YPDMV", 
    "YPDNACL1M", "YPDNACL15M", "YPDNYSTATIN", "YPDSDS", "YPDSODIUMMETAARSENITE", 
    "YPETHANOL", "YPGALACTOSE", "YPGLYCEROL", "YPRIBOSE", "YPSORBITOL", "YPXYLOSE")
ORF_feats <- c(250, 1250, 250, 500, 250, 500, 500, 250, 750, 250, 500, 500, 250, 
    250, 500, 500, 250, 250, 250, 250, 250, 250, 750, 250, 250, 250, 500, 250, 
    250, 250, 500, 250, 250, 250, 250)
CNO_feats <- c(250, 500, 250, 250, 250, 500, 250, 250, 500, 250, 250, 250, 500, 
    250, 500, 500, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 
    250, 250, 250, 250, 500, 250, 250)
genes <- read.csv(
    "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs_GO_genes.tsv", sep="\t")
genes <- genes[c(1,3)] # keep only orfs and genes
genes <- genes[!duplicated(genes),] # orfs mapped to genes
colnames(genes) <- c("orf", "gene")
go <- read.csv(
    "Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs/ORA_SHAP_values_sorted_average_YPACETATE_cno_250_training.tsv", sep="\t")
go <- go[c(1,11:13)]

##### Feature importance scores of most predictive models for each environment
orf_files <- c()
cno_files <- c()
for (i in 1:length(envs)){
    orf_files[i] <- paste(envs[i], "orf", ORF_feats[i], "imp", sep="_")
    cno_files[i] <- paste(envs[i], "cno", CNO_feats[i], "imp", sep="_")
}

dir <- "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results"
min_max <- function(x){ # min-max normalization
    return((x-min(x))/(max(x)-min(x)))
}
combine_data <- function(type="orf", dir=dir, files=orf_files){
    master_imp_orf <- read.csv(file.path(dir, files[1]), sep="\t")
    colnames(master_imp_orf) <- c(type, envs[1])
    for (i in 2:length(files)){
        df2 <- read.csv(file.path(dir, files[i]), sep="\t")
        colnames(df2) <- c(type, envs[i])
        master_imp_orf <- full_join(master_imp_orf, df2, by=type)
    }
    return(master_imp_orf)
}
combine_top5 <- function(type="orf", dir=dir, files=orf_files){
    df <- read.csv(file.path(dir, files[1]), sep="\t")
    colnames(df) <- c(type, envs[1])
    df <- right_join(genes, df, by="orf") # add genes
    df <- df[order(df[,3], decreasing=T),] # sort by importance score
    df <- df[complete.cases(df),] # remove orfs with no gene match
    df <- aggregate(df[,2:3], by=list(df$gene), median)[,c(1,3)] # take median of orf imp scores for those that mapped to the same gene
    colnames(df) <- c("gene", envs[1])
    df[,2] <- min_max(df[,2]) # min-max normalization
    df <- df[order(df[,2], decreasing=T),] # sort
    df <- df[1:5,] # subset top 5 orfs
    for (i in 2:length(files)){
        df2 <- read.csv(file.path(dir, files[i]), sep="\t")
        colnames(df2) <- c(type, envs[i])
        df2 <- right_join(genes, df2, by="orf") # add genes
        df2 <- df2[order(df2[,3], decreasing=T),] # sort by importance score
        df2 <- df2[complete.cases(df2),] # remove orfs with no gene match
        df2 <- aggregate(df2[,2:3], by=list(df2$gene), median)[,c(1,3)] # take median of orf imp scores for those that mapped to the same gene
        colnames(df2) <- c("gene", envs[i])
        df2[,2] <- min_max(df2[,2]) # min-max normalization
        df2 <- df2[order(df2[,2], decreasing=T),] # sort
        df2 <- df2[1:5,] # subset top 5 orfs
        df <- full_join(df, df2, by="gene") # combine to master df
    }
    return(df)
}
plot_hm <- function(toplot, col_fun, n=2, rl, s=10, h, w=35*unit(0.5, "cm"), name, at){
    p <- Heatmap(toplot, 
    col=col_fun(n),
    na_col="grey",
    border_gp = gpar(col = "black"),
    cluster_rows=F,
    show_row_dend=F,
    cluster_columns=F,
    show_column_dend=F,
    row_names_side="left",
    row_labels=rl,
    row_names_gp = gpar(fontsize = s),
    name=name,
    height=h, 
    width=w,
    heatmap_legend_param=list(title=name, at=at))
    return(p)
}

master_imp_orf <- combine_data(type="orf", dir=dir, files=orf_files) # ORF presence/absence feature importance scores
master_imp_orf <- right_join(genes, master_imp_orf, by="orf") # add genes
master_imp_orf_top20 <- combine_top20(type="orf", dir=dir, files=orf_files)
master_imp_cno <- combine_data(type="orf", dir=dir, files=cno_files) # ORF copy number feature importance scores
master_imp_cno <- right_join(genes, master_imp_cno, by="orf") # add genes
master_imp_cno_top20 <- combine_top20(type="orf", dir=dir, files=cno_files)
#summary(master_imp_orf)
#summary(master_imp_cno)

col_fun <- colorRampPalette(c("white", "red"))
toplot <- as.matrix(master_imp_orf[3:37]) # ORF presence/absence
path <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs/orf_imp_all.pdf"
pdf(path, height=unit(163, "in"), width=unit(11, "in")) # file to save
p <- plot_hm(toplot, col_fun, n=2, rl=master_imp_orf$gene, s=5,
    h=2004*unit(0.2, "cm"), w=35*unit(0.5, "cm"),
    name="Importance Score", at=c(0,0.1))
draw(p)
dev.off()

toplot <- as.matrix(master_imp_orf_top20[2:36])# ORF presence/absence top 20 genes
path <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs/orf_imp_top20.pdf"
pdf(path, height=unit(35, "in"), width=unit(11, "in")) # file to save
p <- plot_hm(toplot, col_fun, n=2, rl=master_imp_orf_top20$gene, s=10,
    h=160*unit(0.5, "cm"), w=35*unit(0.5, "cm"),
    name="Importance Score", at=c(0,0.1))
draw(p)
dev.off()
# max(unique(master_imp_orf$YPETHANOL)[3:216])
# [1] 0.1178526

toplot <- as.matrix(master_imp_cno[3:37]) # ORF copy number
path <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs/cno_imp_all.pdf"
pdf(path, height=unit(163, "in"), width=unit(11, "in")) # file to save
p <- plot_hm(toplot, col_fun, n=2, rl=master_imp_orf$gene, s=5,
    h=2004*unit(0.2, "cm"), w=35*unit(0.5, "cm"),
    name="Importance Score", at=c(0,0.1))
draw(p)
dev.off()

toplot <- as.matrix(master_imp_cno_top20[2:36])# ORF copy number top 20 genes
path <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs/cno_imp_top20.pdf"
pdf(path, height=unit(35, "in"), width=unit(11, "in")) # file to save
p <- plot_hm(toplot, col_fun, n=2, rl=master_imp_orf_top20$gene, s=10,
    h=160*unit(0.5, "cm"), w=35*unit(0.5, "cm"),
    name="Importance Score", at=c(0,0.1))
draw(p)
dev.off()
# max(unique(master_imp_cno[17])[3:499,])
# [1] 0.5575024

##### SHAP values of most predictive models for each environment
dir <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs" # path to FS average SHAP files
orf_files <- list.files(path=dir, pattern="^SHAP_values_sorted_average_Y[A-Z0-9]+_orf", 
    full.names=FALSE, recursive=FALSE)
cno_files <- list.files(path=dir, pattern="^SHAP_values_sorted_average_Y[A-Z0-9]+_cno", 
    full.names=FALSE, recursive=FALSE)
master_shap_orf <- combine_data(type="orf", dir=dir, files=orf_files) # ORF presence/absence
master_shap_orf <- right_join(genes, master_shap_orf, by="orf") # add genes
master_shap_orf_top5 <- combine_top5(type="orf", dir=dir, files=orf_files) # top 5 genes
master_shap_orf_top5[is.na(master_shap_orf_top5)] <- 0 # replace missing values with 0
write.csv(master_shap_orf, file.path(dir, "orf_shap_all.csv"), quote=F, row.names=F)
write.csv(master_shap_orf_top5, file.path(dir, "orf_shap_top5_normalized.csv"), quote=F, row.names=F)
master_shap_cno <- combine_data(type="orf", dir=dir, files=cno_files) # ORF copy number
master_shap_cno <- right_join(genes, master_shap_cno, by="orf") # add genes
master_shap_cno_top5 <- combine_top5(type="orf", dir=dir, files=cno_files) # top 5 genes
master_shap_cno_top5[is.na(master_shap_cno_top5)] <- 0 # replace missing values with 0
write.csv(master_shap_cno, file.path(dir, "cnv_shap_all.csv"), quote=F, row.names=F)
write.csv(master_shap_cno_top5, file.path(dir, "cnv_shap_top5_normalized.csv"), quote=F, row.names=F)
#summary(master_shap_orf)
#summary(master_shap_cno)

col_fun <- colorRampPalette(c("white", "red"))
#col_fun <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))
toplot <- as.matrix(master_shap_orf[3:37]) # ORF presence/absence
path <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs/orf_shap_all.pdf"
pdf(path, height=unit(163, "in"), width=unit(11, "in")) # file to save
p <- plot_hm(toplot, col_fun, n=2, rl=master_shap_orf$gene, s=5,
    h=2004*unit(0.2, "cm"), w=35*unit(0.5, "cm"),
    name="SHAP Value", at=c(0,1))
draw(p)
dev.off()

toplot <- as.matrix(master_shap_orf_top5[2:36])# ORF presence/absence top 5 genes
path <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs/orf_shap_top5_normalized.pdf"
pdf(path, height=unit(17, "in"), width=unit(11, "in")) # file to save
p <- plot_hm(toplot, col_fun, n=2, rl=master_shap_orf_top5$gene, s=10,
    h=65*unit(0.5, "cm"), w=35*unit(0.5, "cm"),
    name="Normalized SHAP Value", at=c(0,1))
draw(p)
dev.off()
# max(unique(master_shap_orf$YPDCUSO410MM)[4:421])
# [1] 0.191616

toplot <- as.matrix(master_shap_cno[3:37]) # ORF copy number
path <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs/cno_shap_all.pdf"
pdf(path, height=unit(125, "in"), width=unit(11, "in")) # file to save
p <- plot_hm(toplot, col_fun, n=2, rl=master_shap_cno$gene, s=5,
    h=1437*unit(0.2, "cm"), w=35*unit(0.5, "cm"),
    name="SHAP Value", at=c(0,1))
draw(p)
dev.off()

toplot <- as.matrix(master_shap_cno_top5[2:36])# ORF copy number top 5 genes
path <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs/cno_shap_top5_normalized.pdf"
pdf(path, height=unit(18, "in"), width=unit(11, "in")) # file to save
p <- plot_hm(toplot, col_fun, n=2, rl=master_shap_cno_top5$gene, s=10,
    h=71*unit(0.5, "cm"), w=35*unit(0.5, "cm"),
    name="Normalized SHAP Value", at=c(0,1))
draw(p)
dev.off()
# max(unique(master_shap_cno$YPDCUSO410MM)[3:501])
# [1] 0.2948732

##### GO Enrichment
dir <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs" # path to GO Enrichment files
orf_files <- list.files(path=dir, pattern="ORA_SHAP_values_sorted_average_[A-Z0-9]+_orf",
    full.names=FALSE, recursive=FALSE)
cno_files <- list.files(path=dir, pattern="ORA_SHAP_values_sorted_average_[A-Z0-9]+_cno",
    full.names=FALSE, recursive=FALSE)

combine_data <- function(dir=dir, files=orf_files){
    df <- read.csv(file.path(dir, files[1]), sep="\t")
    df <- df[which(df$lfdr < 0.05),] # keep significant GO terms
    df <- df[,c(1,8)] # keep GO and odds ratio
    colnames(df) <- c("GO", envs[1])
    for (i in 2:length(files)) {
        df2 <- read.csv(file.path(dir, files[i]), sep="\t")
        df2 <- df2[which(df2$lfdr < 0.05),]
        df2 <- df2[,c(1,8)]
        colnames(df2) <- c("GO", envs[i])
        df <- full_join(df, df2, by="GO")
    }
    return(df)
}

master_go_orf <- combine_data(dir=dir, files=orf_files) # ORF presence/absence
# add BP, CC, MF
master_go_orf$BP <- ""
master_go_orf$CC <- ""
master_go_orf$MF <- ""
for(i in 1:nrow(master_go_orf)){
    tryCatch({
        if(!is.null(getGOTerm(master_go_orf$GO[i])$BP[1])) master_go_orf$BP[i] <- getGOTerm(master_go_orf$GO[i])$BP[1]
        if(!is.null(getGOTerm(master_go_orf$GO[i])$CC[1])) master_go_orf$CC[i] <- getGOTerm(master_go_orf$GO[i])$CC[1]
        if(!is.null(getGOTerm(master_go_orf$GO[i])$MF[1])) master_go_orf$MF[i] <- getGOTerm(master_go_orf$GO[i])$MF[1]
    }, error = function(e){print(paste("no GO for ", master_go_orf$GO[i])); NaN},
        finally = {})
}

master_go_cno <- combine_data(dir=dir, files=cno_files) # ORF copy number
master_go_cno$BP <- ""
master_go_cno$CC <- ""
master_go_cno$MF <- ""
for(i in 1:nrow(master_go_cno)){
    tryCatch({
        if(!is.null(getGOTerm(master_go_cno$GO[i])$BP[1])) master_go_cno$BP[i] <- getGOTerm(master_go_cno$GO[i])$BP[1]
        if(!is.null(getGOTerm(master_go_cno$GO[i])$CC[1])) master_go_cno$CC[i] <- getGOTerm(master_go_cno$GO[i])$CC[1]
        if(!is.null(getGOTerm(master_go_cno$GO[i])$MF[1])) master_go_cno$MF[i] <- getGOTerm(master_go_cno$GO[i])$MF[1]
    }, error = function(e){print(paste("no GO for ", master_go_cno$GO[i])); NaN},
        finally = {})
}

toplot <- as.matrix(master_go_orf[2:36]) # ORF presence/absence Log Odds Ratio
toplot[toplot>10] <- 10 # to avoid Inf issue
toplot[is.na(toplot)] <- 0 # no mising values
col_fun <- colorRampPalette(c("white", "red"))  
path <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs/orf_go_all.pdf"
pdf(path, height=unit(8, "in"), width=unit(11, "in"))# file to save
p <- plot_hm(toplot, col_fun, n=2, rl=master_go_orf$BP, s=10,
    h=7*unit(0.5, "cm"), w=35*unit(0.5, "cm"),
    name="Odds Ratio", at=c(0,10))
draw(p)
dev.off()

toplot <- as.matrix(master_go_cno[2:36]) # ORF copy number Log Odds Ratio
toplot[toplot>10] <- 10 # to avoid Inf issue
toplot[is.na(toplot)] <- 0 # no mising values
col_fun <- colorRampPalette(c("white", "red"))  
path <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs/cno_go_all.pdf"
pdf(path, height=unit(8, "in"), width=unit(11, "in"))# file to save
p <- plot_hm(toplot, col_fun, n=2, rl=master_go_cno$BP, s=10,
    h=7*unit(0.5, "cm"), w=35*unit(0.5, "cm"),
    name="Odds Ratio", at=c(0,10))
draw(p)
dev.off()

################################################################################
# Plotting q-values (cannot plot infinite odds ratios)
files <- list.files(path=dir, pattern="ORA_SHAP_values_sorted_average_[A-Z0-9]+_[a-z]+_[0-9]+_training.tsv",
    full.names=FALSE, recursive=FALSE)
df <- read.csv(files[1], sep="\t")
df_qval <- df[,c(1,6,8,9,10)]
df_qval <- subset(df_qval, df_qval$lfdr < 0.05) # & df_qval$BP != "")
env <- str_extract(files[1], "YP[A-Z0-9]+_[a-z]+[^_]")
colnames(df_qval) <- c("GO", paste(env, "direction", sep="_"), paste(env, "OR", sep="_"), paste(env, "qval", sep="_"), paste(env, "fdr", sep="_"))

for (f in files[2:70]) {
    df2 <- read.csv(f, sep="\t")
    df2_qval <- df2[,c(1,6,8,9,10)]
    df2_qval <- subset(df2_qval, df2_qval$lfdr < 0.05)
    env <- str_extract(f, "YP[A-Z0-9]+_[a-z]+[^_]")
    colnames(df2_qval) <- c("GO", paste(env, "direction", sep="_"), paste(env, "OR", sep="_"), paste(env, "qval", sep="_"), paste(env, "fdr", sep="_"))
    df_qval <- full_join(df_qval, df2_qval, by="GO")
}

df_qval$BP <- ""
df_qval$CC <- ""
df_qval$MF <- ""
for(i in 1:nrow(df_qval)){
        tryCatch({
            if(!is.null(getGOTerm(df_qval[i,1])$BP[1])) df_qval$BP[i] <- getGOTerm(df_qval$GO[i])$BP[1]
            if(!is.null(getGOTerm(df_qval[i,1])$CC[1])) df_qval$CC[i] <- getGOTerm(df_qval$GO[i])$CC[1]
            if(!is.null(getGOTerm(df_qval[i,1])$MF[1])) df_qval$MF[i] <- getGOTerm(df_qval$GO[i])$MF[1]
        }, error = function(e){print(paste("no GO for ", df_qval[i,1])); NaN},
            finally = {})
}

orf_qval <- df_qval[grep("_orf", colnames(df_qval))]
orf_qval <- cbind(df_qval$BP, orf_qval)
orf_qval_sub <- subset(orf_qval, orf_qval[,1] != "")
rownames(orf_qval_sub) <- orf_qval_sub[,1]
cno_qval <- df_qval[grep("_cno", colnames(df_qval))]
cno_qval <- cbind(df_qval$BP, cno_qval)
cno_qval_sub <- subset(orf_qval, cno_qval[,1] != "")
rownames(cno_qval_sub) <- cno_qval_sub[,1]

# plot qvalues
orf_qvals <- orf_qval_sub[grep("_OR", colnames(orf_qval_sub))]
orf_qvals <- orf_qvals[noquote(order(colnames(orf_qvals)))]
rdbu <- c('red', 'blue')
path <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs/orf_go_all.pdf"
heatmap2(as.matrix(orf_qvals), path, "SHAP value", n=3, rdbu, as.numeric(1), row=F, col=F)

palette_func <- colorRampPalette(rdbu)
palette <- rev(palette_func(3)) # reversed color palette
col_fun = colorRamp2(c(-5,0,5), palette)

p <- Heatmap(log(as.matrix(orf_qvals[-1])), 
    col=col_fun,
    na_col="grey",
    cluster_rows=F,
    show_row_dend=F,
    cluster_columns=F,
    show_column_dend=F,
    row_names_side="left",
    row_names_gp = gpar(fontsize = 10),
    name="q-value",
    width=35*unit(2, "cm"),
    #heatmap_width=unit(15, "cm"),
    heatmap_legend_param=list(title="Log Odds Ratio", at=c(-5,0,5)))

pdf(path, height=unit(7, "in"), width=unit(30, "in"))# file to save
draw(p)
dev.off()



cno_qvals <- cno_qval_sub[grep("_OR", colnames(cno_qval_sub))]
cno_qvals <- cno_qvals[noquote(order(colnames(cno_qvals)))]
rdbu <- c('red', 'blue')
path <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs/cno_go_all.pdf"
heatmap2(as.matrix(orf_qvals), path, "SHAP value", n=3, rdbu, as.numeric(1), row=F, col=F)

p <- Heatmap(log(as.matrix(cno_qvals[-1])), 
    col=col_fun,
    na_col="grey",
    cluster_rows=F,
    show_row_dend=F,
    cluster_columns=F,
    show_column_dend=F,
    row_names_side="left",
    row_names_gp = gpar(fontsize = 10),
    name="q-value",
    width=35*unit(2, "cm"),
    #heatmap_width=unit(30, "in"),
    heatmap_legend_param=list(title="q-value", at=c(-5, 0, 5)))

pdf(path, height=unit(7, "in"), width=unit(30, "in")) # file to save
draw(p)
dev.off()