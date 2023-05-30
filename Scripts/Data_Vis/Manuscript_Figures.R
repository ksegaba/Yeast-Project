################################################################################
# Chapter 1 Manuscript Figures
# Table of Contents:
# line : Fitness
# line : Kinship
# line : Fitness vs Kinship
# line : Linkage Disequilibrium Decay
# line : ORF presence/absence and copy number variations
# line : Isolate fitness correlatinos vs ORF dataset correlations
# line : Kinship and ORF dataset correlations
# line : TASSEL5 Principal Component Analysis (PCA) plots
# line : Random Forest (RF) prediction performances (baseline using all features)
# line : Prediction performances of all algorithms (using RF feature selection features)
# line : RF feature selection (FS) curves
# line : RF after FS: Test set performance comparisons
# line : Heritability & RF FS performance comparison
# line : Baseline RF gene importance rank density scatter plots comparing data types
# line : GO Enrichment (RF after FS models)
# line : Pathway Enrichment (RF after FS models)
# line : Gene RF importance scores (after FS) comparisons across environments
# line : Fitness Factors impacting RF performance
# line : RF performance comparisons using different test sets (baseline using all features and randomized label)
# line : Multi-Output RF prediction performances (baseline using all features)

rm(list=ls())
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(gridtext))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpmisc))

setwd("/mnt/home/seguraab/Shiu_Lab/Project")

# Isolate growth condition labels
cond <- c("YPACETATE", "YPD14", "YPD40", "YPD42", "YPD6AU", "YPDANISO10", 
        "YPDANISO20", "YPDANISO50", "YPDBENOMYL200", "YPDBENOMYL500", "YPDCAFEIN40", 
        "YPDCAFEIN50", "YPDCHX05", "YPDCHX1", "YPDCUSO410MM", "YPDDMSO", "YPDETOH", 
        "YPDFLUCONAZOLE", "YPDFORMAMIDE4", "YPDFORMAMIDE5", "YPDHU", "YPDKCL2M", 
        "YPDLICL250MM", "YPDMV", "YPDNACL15M", "YPDNACL1M", "YPDNYSTATIN", "YPDSDS", 
        "YPDSODIUMMETAARSENITE", "YPETHANOL", "YPGALACTOSE", "YPRIBOSE", "YPGLYCEROL", 
        "YPXYLOSE", "YPSORBITOL")
new_cond <- c("YP Acetate 2%", "YPD 14ºC", "YPD 40ºC", "YPD 42ºC", "YPD 6-Azauracile 600 µg/ml",
        "YPD Anisomycin 10 µg/ml", "YPD Anisomycin 20 µg/ml", "YPD Anisomycin 50 µg/ml",
        "YPD Benomyl 200 µg/ml", "YPD Benomyl 500 µg/ml", "YPD Caffeine 40 mM", "YPD Caffeine 50 mM",
        "YPD Cycloheximide 0.5 µg/ml", "YPD Cycloheximide 1 µg/ml", "YPD CuSO4 10 mM", "YPD DMSO 6%",
        "YPD Ethanol 15%", "YPD Fluconazole 20 µg/ml", "YPD Formamide 4%", "YPD Formamide 5%",
        "YPD Hydroxyurea 30 mg/ml", "YPD KCL 2 M", "YPD LiCl 250 mM", "YPD Methylviologen 20 mM",
        "YPD NaCl 1.5 M", "YPD NaCl 1 M", "YPD Nystatin 10 µg/ml", "YPD SDS 0.2%", 
        "YPD Sodium metaarsenite 2.5 mM", "YP Ethanol 2%", "YP Galactose 2%", "YP Ribose 2%",
        "YP Glycerol 2%", "YP Xylose 2%", "YP Sorbitol 2%")
conds <- as.data.frame(cbind(cond, new_cond))

################################################################################
# FITNESS
################################################################################
pheno <- read.csv("Data/Peter_2018/pheno.csv", row.names=1) # fitness data
pheno <- pheno[,cond[order(cond)]] # re-order columns
pheno <- reshape2::melt(pheno) # pivot longer

# Boxplot of fitness distributions in each environment
ggplot(pheno, aes(x=variable, y=value)) + theme_bw(8) +
        geom_boxplot(outlier.shape=NA, color="#5B9BD5") + 
        theme(axis.text.x=element_text(color="black", size=9, angle=45, hjust=1)) +
        theme(axis.text.y=element_text(color="black", size=9))
ggsave("Scripts/Data_Vis/fitness_boxplot.pdf", width=8, height=4, device="pdf", useDingbats=FALSE)

# Fitness correlations among replicates (ensure data quality)
reps <- read.csv("Data/Peter_2018/1002_pheno_all_conditions_4Rep_40h.csv")
reps <- reps[!is.na(reps$Growth.ratio),]
length(unique(reps$Strain)) # only have replicate info for 382 strains :(
stats <- reps %>% group_by(Condition, Strain) %>% 
        dplyr::summarize(min=min(Growth.ratio), Q1=quantile(Growth.ratio, 0.25), 
                median=median(Growth.ratio), mean=mean(Growth.ratio), 
                sd=sd(Growth.ratio), Q3=quantile(Growth.ratio, 0.75), 
                max=max(Growth.ratio))
IDs <- read.csv("Data/Peter_2018/IDcorrespondance.txt", sep="\t")
stats <- right_join(IDs, stats, by=c("YJS_name"="Strain"))
pheno <- read.csv("Data/Peter_2018/pheno.csv")
stats <- stats[stats$Standardized_name %in% pheno$ID, ] # keep only diploids
length(unique(stats$Standardized_name)) # only have replicate info for 124 diploids
write.csv(stats, "Data/Peter_2018/1002_pheno_all_conditions_4Rep_40h_only_diploids.csv", quote=F, row.names=F)
summary(stats[4:10]) # interested in the sd, these should be small
ggplot() + geom_point(aes(x=sd, y=mean, group=Condition, alpha=0.1), stats)
ggsave("Data/Peter_2018/1002_pheno_all_conditions_4Rep_40h_only_diploids_stdev.pdf")

# FITNESS CORRELATIONS

pCorEnvs <- read.csv("Data/Peter_2018/pheno_corr_envs.csv", row.names=1) # Correlation between conditions
colnames(pCorEnvs) <- rownames(pCorEnvs)
pCorIso <- read.csv("Data/Peter_2018/pheno_corr_isolates.csv", row.names=1) # Correlation between isolates across all conditions
pCorIso[pCorIso < 0] <- 0

pdf("Scripts/Data_Vis/pheno_corr_isolates.pdf") # v1, v2 [set <0 to 0] (kinship order), v3 (not kinship order, changed color, set neg to 0)
hm2 <- heatmap.2(as.matrix(pCorIso),
                col = colorRampPalette(c("blue","yellow"))(21),
                trace="none",
                dendrogram="none",
                #Rowv=hm$rowDendrogram, # same order as kinship heatmap
                #Colv=hm$rowDendrogram,
                labRow = FALSE,
                labCol = FALSE,
                cex.axis=1,
                cexRow=1,
                cexCol=1,
                density.info = "none",
                keysize=1,
                key.title="PCC",
                key.par = list(cex=1, mar=c(3,1,3,0)), # bottom, left, top, right
                margins=c(1,1)) # column, row
dev.off()

# pdf("Scripts/Data_Vis/pheno_corr_envs_v2.pdf")
# hm3 <- heatmap.2(as.matrix(pCorEnvs),
#                 col = colorRampPalette(c("blue","white","red"))(21),
#                 #hclustfun = function(x) hclust(x, method="complete"),
#                 trace = "none",
#                 dendrogram = "none",
#                 notecex = 1,
#                 cex.axis = 1,
#                 cexRow = 0.6,
#                 cexCol = 0.6,
#                 density.info = "none",
#                 keysize = 1,
#                 key.title = "PCC",
#                 key.par = list(cex=1, mar=c(3,1,3,0)),
#                 margins = c(11,11))
# dev.off()

# heatmap of pCorEnvs
rdbu_r <- rev(brewer.pal(n=11, "RdBu")) # reversed color palette
col_fun = colorRamp2(seq(-1,1,.2), rdbu_r) #same as cm
#cm = ColorMapping(name="PCC", colors=rdbu_r, levels=seq(-1,1,.2)) # color mapping
pdf("Scripts/Data_Vis/pheno_corr_envs_v3.pdf", height=9, width=9)
hm3 <- ComplexHeatmap::Heatmap(as.matrix(pCorEnvs),
        col=col_fun, show_row_dend=F, show_column_dend=F,
        heatmap_legend_param=list(title="PCC", at=seq(-1,1,.2), color_bar="continuous"))
draw(hm3)
dev.off()

# Binned pCorEnv heatmap
pCorEnvs_binned <- apply(pCorEnvs, 2, cut, seq(-1,1,.2)) # bins on each column
rownames(pCorEnvs_binned) <- rownames(pCorEnvs) # set rownames
pCorEnvs_binned_melt <- reshape2::melt(pCorEnvs_binned) # pivot longer
pCorEnvs_melt <- reshape2::melt(as.matrix(pCorEnvs)) # pivot cor data longer
pCorEnvs_binned <- merge(pCorEnvs_melt, pCorEnvs_binned_melt, by=c("Var1", "Var2")) # add bin info to cor data
pCorEnvs_binned_av <- pCorEnvs_binned %>% group_by(Var1, Var2, value.y) %>% summarize(bin_mean = mean(value.x))  # bin average
pCorEnvs_binned_av2 <- reshape2::dcast(pCorEnvs_binned_av, formula=Var1~Var2, value.var="bin_mean") # pivot binned cor data wider
rownames(pCorEnvs_binned_av2) <- pCorEnvs_binned_av2$Var1 # set rownames
pCorEnvs_binned_av2 <- pCorEnvs_binned_av2[-1] # drop Var1

# draw heatmap
#col_fun = colorRamp2(seq(-1,1,.2), rdbu_r) #same as cm
pdf("Scripts/Data_Vis/pheno_corr_envs_binned.pdf", height=9, width=9)
hm3 <- ComplexHeatmap::Heatmap(as.matrix(pCorEnvs_binned_av2), 
        col=color_mapping_legend(cm@colors),
        heatmap_legend_param=list(title="PCC", at=seq(-1,1,.2), 
                color_bar="discrete"))
draw(hm3)
dev.off()
# the binned and original heatmaps look the same

# log of PCC values
pCorEnvs_log <- log(pCorEnvs)
pdf("Scripts/Data_Vis/pheno_corr_envs_log.pdf", height=9, width=9)
hm3 <- ComplexHeatmap::Heatmap(as.matrix(pCorEnvs_log),
        heatmap_legend_param=list(title="log(PCC)", color_bar="discrete"))
draw(hm3)
dev.off()

################################################################################
# KINSHIP
################################################################################
kin <- fread("Data/Peter_2018/geno_transposed.csv_ordered_kinship.txt", skip=3)
kin <- as.matrix(kin, rownames=1, colnames=1)
colnames(kin) <- rownames(kin)
write.csv(kin, "Data/Peter_2018/kinship.csv", quote=F, row.names=T)

pdf("Scripts/Data_Vis/kinship.pdf")
hm <- heatmap.2(kin, 
                col = colorRampPalette(c("blue", "yellow"))(21),
                trace = "none",
                dendrogram = "none",
                Rowv = TRUE,
                Colv = TRUE,
                notecex = 1,
                cex.axis = 1,
                cexRow = 1,
                cexCol = 1,
                density.info = "none",
                keysize=1,
                key.title = "Kinship",
                key.par = list(cex=1, mar=c(3,1,3,0)), # bottom, left, top, right
                margins=c(1,1))
dev.off()

################################################################################
# FITNESS VS KINSHIP
################################################################################
add_bins <- function(mat, step) {
        to_bin <- reshape2::melt(mat) # reshape kinship matrix
        #quants <- quantile(i_bin$value, seq(0,1,0.05)) # pCor quantiles
        #min <- min(to_bin$value) - step
        max <- max(to_bin$value) + step
        breaks <- seq(min(to_bin$value), max, step) # regular bins
        to_bin$breaks <- cut(to_bin$value, breaks=breaks, include.lowest=T) # regular bins
        to_bin$breaks <- as.character(lapply(strsplit(as.character(to_bin$breaks), split=","), head, n=1))
        to_bin$breaks <- gsub("\\(", "", to_bin$breaks)
        to_bin$breaks <- gsub("\\[", "", to_bin$breaks)
        to_bin$breaks <- as.numeric(to_bin$breaks)
        return (to_bin)
}

kin <- read.csv("Scripts/Data_Vis/kinship.csv", header=T)
k_bin <- add_bins(kin, 0.5) # reshape kinship matrix
colnames(k_bin) <- c("Var1", "Var2", "value", "breaks")

pCorIso <- read.csv("Data/Peter_2018/pheno_corr_isolates.csv", header=T, row.names=1) # Correlation between isolates across all conditions
pCorIso <- as.matrix(pCorIso) # original, no values have been reset to 0
i_bin <- add_bins(pCorIso, 0.12) # reshape fitness across isolates correlation matrix

k_i_bin <- merge(k_bin, i_bin, by=c("Var1","Var2")) # merge dataframes
colnames(k_i_bin) <- c("Isolate 1", "Isolate 2", "Kinship", "Kinship bin", "pCor", "pCor bin")
# > min(k_i_bin$pCor)
# [1] -0.4334015
# > max(k_i_bin$pCor)
# [1] 1
# > mean(k_i_bin$pCor)
# [1] 0.5916347
k.i.cor <- cor.test(k_i_bin$Kinship, k_i_bin$pCor, method="spearman") # Spearman correlation
k_i_bin_count <- aggregate(cbind(count = `Isolate 2`) ~ `Kinship bin`, # Kinship bin counts
        data=k_i_bin, FUN=function(x){NROW(x)})
k_i_bin_median <- k_i_bin %>% group_by(`Kinship bin`) %>% summarise(median=median(pCor)) # median pCor per kinship bin 

get_density <- function(x, y, ...) { # density
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

k_i_bin$density <- get_density(k_i_bin$Kinship, k_i_bin$pCor, n=100)

# 5% and 95% quantiles for each kinship bin
#quants <- as.data.frame(do.call("rbind", tapply(k_i_bin$pCor, k_i_bin$`Kinship bin`, quantile, c(0.05, 0.95))))
#names(quants) <- c("x5","x95")
#quants$`Kinship bin` <- as.numeric(row.names(quants))

# Density scatter plot
#ggplot(k_i_bin, aes(x=Kinship, y=pCor, color=density)) + geom_point(size=0.2, alpha=0.2) +
#        geom_smooth(method="lm", formula=y~x, se=FALSE) +
#        scale_color_viridis(direction=-1) + theme_bw(8)
#ggsave("Scripts/Data_Vis/kinpCorIso_Density.pdf", width=5, height=4, device="pdf", useDingbats=FALSE)

ggplot(k_i_bin, aes(x=Kinship, y=pCor, color=density)) + geom_point(size=0.2, alpha=0.2) +
        scale_color_viridis(direction=-1) + theme_bw(base_size=12) + 
        theme(axis.text.x=element_text(color="black"), axis.text.y = element_text(color="black")) +
        geom_line(data=k_i_bin_median, aes(x=`Kinship bin`, y=median, group=1), color="black")
        #annotate(geom="text", x=2, y=-0.3, label=paste("Spearman's rank correlation = ", round(k.i.cor, digits=2),sep=""))
ggsave("Scripts/Data_Vis/kinpCorIso_Density_median_v3.pdf", width=5, height=4, device="pdf", useDingbats=FALSE)

#ggplot(k_i_bin, aes(x=`Kinship bin`, y=pCor, color=density)) + geom_point(size=0.2, alpha=0.2) +
#        scale_color_viridis(direction=-1) + theme_bw(8) + 
#        geom_line(data=k_i_bin_median, aes(x=`Kinship bin`, y=median, group=1), color="red")
#ggsave("Scripts/Data_Vis/kinpCorIso_Density_median_bin.pdf", width=5, height=4, device="pdf", useDingbats=FALSE)

# Binned plot with quantiles
#ggplot(data = k_i_bin, aes(x=`Kinship bin`, y=pCor, group=`Kinship bin`))  + geom_violin() +
#        geom_boxplot(outlier.shape=NA) + theme_bw(10) +
#        geom_line(data=quants, aes(x=`Kinship bin`, y=x95, group = 1), color="red") +
#        geom_line(data=quants, aes(x=`Kinship bin`, y=x5, group = 1), color="blue")
#ggsave("Scripts/Data_Vis/kinpCorIso_Distribution.pdf", width=4, height = 4, device="pdf", useDingbats=FALSE)

# Bin counts plot
ggplot(k_i_bin_count, aes(x=`Kinship bin`, y=log(count))) + theme_bw(10) +geom_bar(stat="identity", position= "dodge")
ggsave("Scripts/Data_Vis/kinpCorIso_BinCount_v3.pdf", width=3.42, height = 2.42, device="pdf", useDingbats=FALSE)

################################################################################
# LINKAGE DISEQUILIBRIUM DECAY
################################################################################
ld <- fread("Data/Peter_2018/geno_LD_window_50.txt") # Window size 50, TASSEL5 output
ld_sub <- ld[which(ld$Dist_bp != "N/A" & ld$`R^2` != "NaN"), c("Dist_bp", "R^2", "DPrime", "pDiseq")] # remove rows with missing values
ld_sub$Dist_bp <- as.numeric(ld_sub$Dist_bp)
ld_sub$bin <- cut(ld_sub$Dist_bp, breaks=c(seq(1, max(ld_sub$Dist_bp), by=50))) # bin by distance
ld_sub$bin <- as.character(lapply(strsplit(as.character(ld_sub$bin), split=","),head, n=1)) # rename bins
ld_sub$bin <- gsub("\\(", "", ld_sub$bin)

# plot median R2 in each 50 bp bin
out <- ld_sub %>% group_by(bin) %>% summarise(med_R2=median(`R^2`), sd_R2=sd(`R^2`))
pdf("Scripts/Data_Vis/ld_w50_binned.pdf")
out %>% ggplot(aes(x=as.numeric(bin), y=med_R2)) + geom_point() + theme_bw() +
        xlab("Distance (bp)") + ylab(parse(text="R^2")) + 
        geom_smooth()
dev.off()

# LD decay plot of first 5000 bp
out <- right_join(out, ld_sub, by="bin") # add Dist_bp info
out_sub <- out[out$Dist_bp <= 5000,]
out_sub <- out_sub[!duplicated(out_sub[1:3]),]
pdf("Scripts/Data_Vis/ld_w50_binned_first_5000bp.pdf")
out_sub %>% ggplot(aes(x=as.numeric(bin), y=med_R2)) + geom_point() + theme_bw() +
        xlab("Distance (bp)") + ylab(parse(text="R^2")) + geom_smooth()
dev.off()

 # save to file
write.table(out, "Data/Peter_2018/geno_LD_window_50_binned.tsv", sep="\t", quote=F, row.names=F)

################################################################################
# Allele frequencies
################################################################################
X <- fread("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv")
X <- as.matrix(X, rownames=1, colnames=1)
X = X+1

################################################################################
# ORF PRESENCE/ABSENCE & COPY NUMBER VARIATIONS
################################################################################
orf <- fread("Data/Peter_2018/ORFs_pres_abs.csv") # ORF presence/absence
orf <- as.matrix(orf, rownames=1, colnames=1)
oCor <- cor(t(orf)) # between isolates
# > min(oCor)
# [1] 0.6016447
# > max(oCor)
# [1] 1
# > mean(oCor)
# [1] 0.9396035
oCor[oCor < 0.9] <- 0.9

cno <- fread("Data/Peter_2018/ORFs_no_NA.csv") # ORF copy number variants
cno <- as.matrix(cno, rownames=1, colnames=1)
cCor <- cor(t(cno))
# > min(cCor)
# [1] 0.2334361
# > max(cCor)
# [1] 1
# > mean(cCor)
# [1] 0.8875162
cCor[cCor < 0.8] <- 0.8

pdf("Scripts/Data_Vis/orf_corr_isolates.pdf") # kinship order (v1, v2 [reset < 0.9 or 0.8]), v3 does not follow kinship order, v4 (orf and cno same order, not kinship order)
hm5 <- heatmap.2(as.matrix(oCor),
        col = colorRampPalette(c("blue","yellow"))(21),
        trace="none",
        dendrogram="none",
        #Rowv=hm$rowDendrogram, # same order as kinship heatmap
        #Colv=hm$rowDendrogram,
        labRow = "Isolates",
        labCol = "Isolates",
        cexRow = 1,
        cexCol = 1,
        density.info = "none",
        keysize = 1,
        key.title = "PCC",
        key.par = list(cex=1, mar=c(3,1,3,0)), # bottom, left, top, right
        margins = c(1,1)) # column, row
dev.off()

pdf("Scripts/Data_Vis/cno_corr_isolates.pdf")
heatmap.2(as.matrix(cCor),
        col = colorRampPalette(c("blue","yellow"))(21),
        trace="none",
        dendrogram="none",
        Rowv=hm5$rowDendrogram, # same order as v3 of orf_corr_isolates
        Colv=hm5$rowDendrogram,
        labRow = "Isolates",
        labCol = "Isolates",
        cexRow = 1,
        cexCol = 1,
        density.info = "none",
        keysize = 1,
        key.title = "PCC",
        key.par = list(cex=1, mar=c(3,1,3,0)), # bottom, left, top, right
        margins = c(1,1)) # column, row
dev.off()

################################################################################
# ISOLATE FITNESS CORRELATION VS ORF DATASET CORRELATIONS
################################################################################
orf <- fread("Data/Peter_2018/ORFs_pres_abs.csv") # ORF presence/absence
orf <- as.matrix(orf, rownames=1, colnames=1)
oCor <- cor(t(orf)) # between isolates
o_bin <- add_bins(oCor, 0.05)
cno <- fread("Data/Peter_2018/ORFs_no_NA.csv") # ORF copy number variants
cno <- as.matrix(cno, rownames=1, colnames=1)
cCor <- cor(t(cno))
c_bin <- add_bins(cCor, 0.05)

# oCor vs pCor
o_i_bin <- merge(o_bin, i_bin, by=c("Var1", "Var2")) # merge dataframes
colnames(o_i_bin) <- c("Isolate 1", "Isolate 2", "oCor", "oCor bin", "pCor", "pCor bin")
# > min(o_i_bin$oCor)
# [1] 0.6016447
# > max(o_i_bin$oCor)
# [1] 1
# > mean(o_i_bin$oCor)
# [1] 0.9396035
o.i.cor <- cor.test(o_i_bin$oCor, o_i_bin$pCor, method="spearman") # Spearman correlation
o_i_bin_count <- aggregate(cbind(count = `Isolate 2`) ~ `oCor bin`, # oCor bin counts
        data=o_i_bin, FUN=function(x){NROW(x)})
o_i_bin_median <- o_i_bin %>% group_by(`oCor bin`) %>% summarise(median=median(pCor)) # median pCor per oCor bin
o_i_bin$density <- get_density(o_i_bin$oCor, o_i_bin$pCor, n=100)
ggplot(o_i_bin, aes(x=oCor, y=pCor, color=density)) + geom_point(size=0.2, alpha=0.2) +
        scale_color_viridis(direction=-1) + theme_bw(base_size=12) + 
        theme(axis.text.x=element_text(color="black"), axis.text.y = element_text(color="black")) +
        geom_line(data=o_i_bin_median, aes(x=`oCor bin`, y=median, group=1), color="black")
        #annotate(geom="text", x=2, y=-0.3, label=paste("Spearman's rank correlation = ", round(o.i.cor[[4]][[1]], digits=2),sep=""))
ggsave("Scripts/Data_Vis/oCorpCorIso_Density_median.pdf", width=5, height=4, device="pdf", useDingbats=FALSE)

# cCor vs pCor
c_i_bin <- merge(c_bin, i_bin, by=c("Var1", "Var2")) # merge dataframes
colnames(c_i_bin) <- c("Isolate 1", "Isolate 2", "cCor", "cCor bin", "pCor", "pCor bin")
# > min(c_i_bin$cCor)
# [1] 0.2334361
# > max(c_i_bin$cCor)
# [1] 1
# > mean(c_i_bin$cCor)
# [1] 0.8875162
c.i.cor <- cor.test(c_i_bin$cCor, c_i_bin$pCor, method="spearman") # Spearman correlation
c_i_bin_count <- aggregate(cbind(count = `Isolate 2`) ~ `cCor bin`, # cCor bin counts
        data=c_i_bin, FUN=function(x){NROW(x)})
c_i_bin_median <- c_i_bin %>% group_by(`cCor bin`) %>% summarise(median=median(pCor)) # median pCor per cCor bin
c_i_bin$density <- get_density(c_i_bin$cCor, c_i_bin$pCor, n=100)
ggplot(c_i_bin, aes(x=cCor, y=pCor, color=density)) + geom_point(size=0.2, alpha=0.2) +
        scale_color_viridis(direction=-1) + theme_bw(base_size=12) + 
        theme(axis.text.x=element_text(color="black"), axis.text.y = element_text(color="black")) +
        geom_line(data=c_i_bin_median, aes(x=`cCor bin`, y=median, group=1), color="black")
        #annotate(geom="text", x=2, y=-0.3, label=paste("Spearman's rank correlation = ", round(c.i.cor[[4]][[1]], digits=2),sep=""))
ggsave("Scripts/Data_Vis/cCorpCorIso_Density_median.pdf", width=5, height=4, device="pdf", useDingbats=FALSE)

################################################################################
# KINSHIP VS ORF DATASET CORRELATIONS
################################################################################
# Kinship vs oCor
k_o_bin <- merge(k_bin, o_bin, by=c("Var1", "Var2")) # merge dataframes
colnames(k_o_bin) <- c("Isolate 1", "Isolate 2", "Kinship", "Kinship bin", "oCor", "oCor bin")
k.o.cor <- cor.test(k_o_bin$Kinship, k_o_bin$oCor, method="spearman") # Spearman correlation
k_o_bin_median <- k_o_bin %>% group_by(`Kinship bin`) %>% summarise(median=median(oCor)) # median oCor per Kinship bin
k_o_bin$density <- get_density(k_o_bin$Kinship, k_o_bin$oCor, n=100)
ggplot(k_o_bin, aes(x=Kinship, y=oCor, color=density)) + geom_point(size=0.2, alpha=0.2) +
        scale_color_viridis(direction=-1) + theme_bw(base_size=12) + 
        theme(axis.text.x=element_text(color="black"), axis.text.y = element_text(color="black")) +
        geom_line(data=k_o_bin_median, aes(x=`Kinship bin`, y=median, group=1), color="black")
        #annotate(geom="text", x=2, y=-0.3, label=paste("Spearman's rank correlation = ", round(k.o.cor[[4]][[1]], digits=2),sep=""))
ggsave("Scripts/Data_Vis/kinoCor_Density_median.pdf", width=5, height=4, device="pdf", useDingbats=FALSE)

# Kinship vs cCor
k_c_bin <- merge(k_bin, c_bin, by=c("Var1", "Var2")) # merge dataframes
colnames(k_c_bin) <- c("Isolate 1", "Isolate 2", "Kinship", "Kinship bin", "cCor", "cCor bin")
k.c.cor <- cor.test(k_c_bin$Kinship, k_c_bin$cCor, method="spearman") # Spearman correlation
k_c_bin_median <- k_c_bin %>% group_by(`Kinship bin`) %>% summarise(median=median(cCor)) # median cCor per Kinship bin
k_c_bin$density <- get_density(k_c_bin$Kinship, k_c_bin$cCor, n=100)
ggplot(k_c_bin, aes(x=Kinship, y=cCor, color=density)) + geom_point(size=0.2, alpha=0.2) +
        scale_color_viridis(direction=-1) + theme_bw(base_size=12) + 
        theme(axis.text.x=element_text(color="black"), axis.text.y = element_text(color="black")) +
        geom_line(data=k_c_bin_median, aes(x=`Kinship bin`, y=median, group=1), color="black")
        #annotate(geom="text", x=2, y=-0.3, label=paste("Spearman's rank correlation = ", round(k.c.cor[[4]][[1]], digits=2),sep=""))
ggsave("Scripts/Data_Vis/kincCor_Density_median.pdf", width=5, height=4, device="pdf", useDingbats=FALSE)

# oCor vs cCor
o_c_bin <- merge(o_bin, c_bin, by=c("Var1", "Var2")) # merge dataframes
colnames(o_c_bin) <- c("Isolate 1", "Isolate 2", "oCor", "oCor bin", "cCor", "cCor bin")
o.c.cor <- cor.test(o_c_bin$oCor, o_c_bin$cCor, method="spearman") # Spearman correlation
o_c_bin_median <- o_c_bin %>% group_by(`oCor bin`) %>% summarise(median=median(cCor)) # median cCor per oCor bin
o_c_bin$density <- get_density(o_c_bin$oCor, o_c_bin$cCor, n=100)
ggplot(o_c_bin, aes(x=oCor, y=cCor, color=density)) + geom_point(size=0.2, alpha=0.2) +
        scale_color_viridis(direction=-1) + theme_bw(base_size=12) + 
        theme(axis.text.x=element_text(color="black"), axis.text.y = element_text(color="black")) +
        geom_line(data=o_c_bin_median, aes(x=`oCor bin`, y=median, group=1), color="black")
        #annotate(geom="text", x=2, y=-0.3, label=paste("Spearman's rank correlation = ", round(o.c.cor[[4]][[1]], digits=2),sep=""))
ggsave("Scripts/Data_Vis/oCorcCor_Density_median.pdf", width=5, height=4, device="pdf", useDingbats=FALSE)

################################################################################
# TASSEL5 PCA PLOTS
################################################################################
pcs <- read.table("Data/Peter_2018/geno_PCs1.csv", skip=2, sep=",", header=T) # PCs
meta <- read_excel("Data/Peter_2018/Peter_2018_Supplementary_Tables.xls", 
        sheet=1, range=cell_rows(4:1015)) # Isolate metadata
pcs <- right_join(meta[c(2,4,5)], pcs, by=c("Standardized name"="Taxa")) # add strain origin

pdf("Scripts/Data_Vis/PCA_tassel_eco_1v2.pdf", height=6, width=9)
pcs %>% ggplot(aes(x=PC1, y=PC2, col=`Ecological origins`)) + geom_point()
dev.off()
pdf("Scripts/Data_Vis/PCA_tassel_eco_1v3.pdf", height=6, width=9)
pcs %>% ggplot(aes(x=PC1, y=PC3, col=`Ecological origins`)) + geom_point()
dev.off()
pdf("Scripts/Data_Vis/PCA_tassel_eco_1v4.pdf", height=6, width=9)
pcs %>% ggplot(aes(x=PC1, y=PC4, col=`Ecological origins`)) + geom_point()
dev.off()
pdf("Scripts/Data_Vis/PCA_tassel_eco_1v5.pdf", height=6, width=9)
pcs %>% ggplot(aes(x=PC1, y=PC5, col=`Ecological origins`)) + geom_point()
dev.off()
pdf("Scripts/Data_Vis/PCA_tassel_geo_1v2.pdf", height=6, width=6)
pcs %>% ggplot(aes(x=PC1, y=PC2, col=`Geographical origins`)) + geom_point() + theme(legend.position="none")
dev.off()
pdf("Scripts/Data_Vis/PCA_tassel_geo_1v3.pdf", height=6, width=6)
pcs %>% ggplot(aes(x=PC1, y=PC3, col=`Geographical origins`)) + geom_point() + theme(legend.position="none")
dev.off()
pdf("Scripts/Data_Vis/PCA_tassel_geo_1v4.pdf", height=6, width=6)
pcs %>% ggplot(aes(x=PC1, y=PC4, col=`Geographical origins`)) + geom_point() + theme(legend.position="none")
dev.off()
pdf("Scripts/Data_Vis/PCA_tassel_geo_1v5.pdf", height=6, width=6)
pcs %>% ggplot(aes(x=PC1, y=PC5, col=`Geographical origins`)) + geom_point() + theme(legend.position="none")
dev.off()

################################################################################
# BASELINE RANDOM FOREST PREDICTION PERFORMANCE (figures in excel)
################################################################################
## PC Models
# These results were originally with the SNP results but I took them out and moved them to a new file
# path <- "/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results"
# rf_fs <- read.csv(file.path(path, "RESULTS_reg.txt"), sep="\t")
# rf_pc <- rf_fs[grep("PCs", rf_fs$ID),]# population structure results
# rf_pc <- right_join(conds, rf_pc, by=c("cond"="Y")) # add full condition names
# write.table(rf_pc, "Results/RESULTS_RF_PCs.txt", sep="\t", row.names=F, quote=F)
path <- "/mnt/home/seguraab/Shiu_Lab/Project/Results/RESULTS_RF_PCs.txt"
rf_pc <- read.csv(path, sep="\t")
rf_pc <- rf_pc[order(rf_pc$r2_test, decreasing=T),]
cond_order <- rf_pc$cond
write.table(rf_pc, "Results/RESULTS_RF_PCs_sorted.txt", sep="\t", row.names=F, quote=F)

## SNP Models
path <- "/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results"
rf_base <- read.csv(file.path(path, "baseline/RESULTS_reg_baseline.txt"), sep="\t")
rf_base$cond <- gsub("_rf_baseline", "", rf_base$ID)
rf_base <- right_join(conds, rf_base, by="cond") # add full condition names
rownames(rf_base) <- rf_base$cond
rf_base <- rf_base[cond_order,]
#  rownames(rf_base)==cond_order
write.table(rf_base, "Results/RESULTS_RF_SNPs_baseline.txt", sep="\t", row.names=F, quote=F)

## ORF Copy Number Models
path <- "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results"
rf_orf <- read.csv(file.path(path, "RESULTS_reg.txt"), sep="\t")
rf_orf <- rf_orf[grep("2022-0[4-6]-[0-2][5-9]",rf_orf$DateTime),] # extract correct runs by date
rf_cnv <- rf_orf[grep("_cno_baseline", rf_orf$ID),] # ORF copy number models
which(duplicated(rf_cnv$ID)) # YPDDANISO20 is duplicated (almost exact numbers due to the nature of the algorithm)
rf_cnv <- rf_cnv[!duplicated(rf_cnv$ID),] # remove the 2022-06-09 run; keep april date for consistency purposes
rf_cnv <- right_join(conds, rf_cnv, by=c("cond"="Y")) # add full condition names
rownames(rf_cnv) <- rf_cnv$cond
rf_cnv <- rf_cnv[cond_order,]
# rownames(rf_cnv)==cond_order
write.table(rf_cnv, "Results/RESULTS_RF_CNVs_baseline.txt", sep="\t", row.names=F, quote=F)

## ORF Presence/Absence Models
rf_orf <- rf_orf[grep("_orf_baseline", rf_orf$ID),] # ORF presence/absence
# rf_orf[which(duplicated(rf_orf$ID)),] # 8 traits duplicated
rf_orf <- rf_orf[!duplicated(rf_orf$ID),] # remove duplicates that are later by date
rf_orf <- right_join(conds, rf_orf, by=c("cond"="Y")) # add full condition names & remove duplicate rows
rownames(rf_orf) <- rf_orf$cond
rf_orf <- rf_orf[cond_order,]
write.table(rf_orf, "Results/RESULTS_RF_ORFs_baseline.txt", sep="\t", row.names=F, quote=F)

################################################################################
# PREDICTION PERFORMANCES OF ALL ALGORITHMS (USING RF FEATURE SELECTION FEATURES)
################################################################################
## SNP
# read in results files for each algorithm
dir <- "/mnt/gs21/scratch/seguraab/yeast_project/"
rf_snp <- read.csv("Results/RESULTS_RF_SNPs_FS.txt", sep="\t") # Random Forest results
xgb_snp <- read.csv(paste(
        dir,"SNP_yeast_XGBoost_results/fs/RESULTS_xgboost.txt", sep=""), sep="\t") # XGBoost results
bl_snp <- read.csv(paste(
        dir, "SNP_yeast_BL_results/fs/RESULTS_BL.txt", sep=""), sep="\t") # Bayesian LASSO results
bayesc_snp_files <- list.files() # BayesC results
gblup_snp_files <- list.files() # gBLUP results
rrblup_snp_files_val <- paste(dir,
        "SNP_yeast_rrBLUP_results/R2_cv_results_rrBLUP_top", rf_snp$FeatureNum,
        "_", rf_snp$cond, ".csv", sep="") # rrBLUP validation results files
rrblup_snp_files_test <- paste(dir,
        "SNP_yeast_rrBLUP_results/R2_test_results_rrBLUP_top", rf_snp$FeatureNum,
        "_", rf_snp$cond, ".csv", sep="") # rrBLUP test results files
rrblup_snp_preds_val <- paste(dir,
        "SNP_yeast_rrBLUP_results/Predict_value_cv_rrBLUP_top", rf_snp$FeatureNum,
        "_", rf_snp$cond, ".csv", sep="") # rrBLUP predicted labels of validation set
rrblup_snp_preds_test <- paste(dir,
        "SNP_yeast_rrBLUP_results/Predict_value_test_rrBLUP_top", rf_snp$FeatureNum,
        "_", rf_snp$cond, ".csv", sep="") # rrBLUP predicted labels of test set

# add missing columns to xgb_snp
xgb_snp <- left_join(
        rf_snp[,c("cond", "new_cond", "Tag", "NumInstances", "FeatureNum", "CVfold", "CV_rep")],
        xgb_snp, by=c("cond"="Trait"))
colnames(xgb_snp) <- c("cond", "new_cond", "Tag", "NumInstances", "FeatureNum",
        "CVfold", "CV_rep", "DateTime", "RunTime", "Data", "MSE_val", "MSE_val_sd",
        "RMSE_val", "RMSE_val_sd", "EVS_val", "EVS_val_sd","r2_val", "r2_val_sd",
        "PCC_val", "PCC_val_sd", "MSE_test", "MSE_test_sd", "RMSE_test",
        "RMSE_test_sd", "EVS_test", "EVS_test_sd", "r2_test", "r2_test_sd",
        "PCC_test", "PCC_test_sd") # set same as rf_snp
write.table(xgb_snp, "Results/RESULTS_XGB_SNPs_FS.txt", sep="\t", quote=F, row.names=F)

# combine rf, xgb, and bl
all <- as.data.frame(rbind.fill(rf_snp, xgb_snp))
all <- dplyr::select(all,-c("DateTime", "RunTime", "Tag", "Data", "EVS_val",
        "EVS_val_sd", "EVS_val_se", "EVS_test", "EVS_test_sd", "EVS_test_se",
        "RMSE_val", "RMSE_val_sd", "RMSE_test", "RMSE_test_sd")) # drop unnecessary columns and fill in missing info
all[is.na(all$ID), "ID"] <- gsub("rf", "xgboost", all[!is.na(all$ID), "ID"])
all[is.na(all$Alg), "Alg"] <- "XGBoost"
all[is.na(all$MSE_val_se), "MSE_val_se"] <- all[all$Alg=="XGBoost", "MSE_val_sd"]/sqrt(10)
all[is.na(all$r2_val_se), "r2_val_se"] <- all[all$Alg=="XGBoost", "r2_val_sd"]/sqrt(10)
all[is.na(all$PCC_val_se), "PCC_val_se"] <- all[all$Alg=="XGBoost", "PCC_val_sd"]/sqrt(10)
all[is.na(all$MSE_test_se), "MSE_test_se"] <- all[all$Alg=="XGBoost", "MSE_test_sd"]/sqrt(10)
all[is.na(all$r2_test_se), "r2_test_se"] <- all[all$Alg=="XGBoost", "r2_test_sd"]/sqrt(10)
all[is.na(all$PCC_test_se), "PCC_test_se"] <- all[all$Alg=="XGBoost", "PCC_test_sd"]/sqrt(10)
bl_snp <- left_join(conds, bl_snp, by=c("cond"="Trait"))
bl_snp <- tibble::add_column(bl_snp, CVfold = 5, .before = "MSE_val")
bl_snp <- tibble::add_column(bl_snp, CV_rep = 10, .before = "MSE_val")
write.table(bl_snp, "Results/RESULTS_BL_SNPs_FS.txt", sep="\t", quote=F, row.names=F)
all <- as.data.frame(rbind.fill(all, dplyr::select(bl_snp, -c("Date", "RunTime"))))

# combine rrblup_snp_files_val and _test results
rrblup_snp <- data.frame(matrix(nrow=1, ncol=length(colnames(all)))) # rrBLUP results
colnames(rrblup_snp) <- colnames(all)
mse <- function(preds, actual){ return(mean((actual-preds)^2)) }
se <- function(vector){ return(sd(vector)/length(vector)) }
pheno <- read.csv("Data/Peter_2018/pheno.csv", row.names=1)
for (i in 1:length(rrblup_snp_files_val)){
        r2_val <- read.csv(rrblup_snp_files_val[i])
        r2_test <- read.csv(rrblup_snp_files_test[i])
        cond <- colnames(r2_val)
        preds_val <- read.csv(rrblup_snp_preds_val[i], row.names=1)
        preds_test <- read.csv(rrblup_snp_preds_test[i], row.names=1)
        actual_val <- pheno[rownames(preds_val), cond]
        actual_test <- pheno[rownames(preds_test), cond]
        mse_val <- apply(preds_val, 2, mse, actual_val)
        mse_test <- apply(preds_test, 2, mse, actual_test)
        pcc_val <- apply(preds_val, 2, cor, actual_val)
        pcc_test <- apply(preds_test, 2, cor, actual_test)
        rrblup_snp <- rbind(rrblup_snp, list(cond, "", "", "rrBLUP", 625,
                0, 5, 10, mean(mse_val), sd(mse_val), se(mse_val),
                mean(as.matrix(r2_val)), sd(as.matrix(r2_val)),
                se(as.matrix(r2_val)), mean(pcc_val), sd(pcc_val), se(pcc_val),
                mean(mse_test), sd(mse_test), se(mse_test),
                mean(as.matrix(r2_test)), sd(as.matrix(r2_test)),
                se(as.matrix(r2_val)), mean(pcc_test), sd(pcc_test),
                se(pcc_test)))
}
rrblup_snp <- rrblup_snp[-1,] # drop first row of NAs
rrblup_snp$new_cond <- ifelse(rrblup_snp$cond==rf_snp$cond, rf_snp$new_cond, NA)
rrblup_snp$FeatureNum <- ifelse(rrblup_snp$cond==rf_snp$cond, rf_snp$FeatureNum, NA)
rrblup_snp$ID <- paste(rrblup_snp$cond, "rrblup", rrblup_snp$FeatureNum, sep="_")
write.table(rrblup_snp, "Results/RESULTS_rrBLUP_SNPs_FS.txt", sep="\t", quote=F, row.names=F)

# combine bayesc_snp_files results

# combine gblup_snp_files results

# combine them into one data matrix
all <- as.data.frame(rbind.fill(all, rrblup_snp)) # combine all and rrblup
all <- as.data.frame(rbind.fill(all, bayesc_snp) # combine all and bayesc
all <- as.data.frame(rbind.fill(all, gblup_snp) # combine all and gblup
write.table(all, "Results/RESULTS_ALL_ALG_FS.txt", sep="\t", quote=F, row.names=F)

all <- as.data.frame(rbind.fill(all, bl_snp[-1,])) # combine all and bayesian lasso
# fill in missing info
all[which(all$Alg=="Bayesian LASSO"), "new_cond"] <- all[which(all$Alg=="RF"), "new_cond"]
all[which(all$Alg=="Bayesian LASSO"), "ID"] <- gsub("rf", "bl", all[which(all$Alg=="RF"), "ID"])
all[which(all$Alg=="Bayesian LASSO"), "FeatureNum"] <- all[which(all$Alg=="RF"), "FeatureNum"]
all <- as.data.frame(rbind.fill(all, bl_snp) # combine all and bayesian lasso
all <- as.data.frame(rbind.fill(all, gblup_snp) # combine all and gblup
all <- as.data.frame(rbind.fill(all, bayesc_snp) # combine all and bayesc
write.table(all, "Results/RESULTS_ALL_ALG_FS.txt", sep="\t", quote=F, row.names=F)

## ORF Copy Number



## ORF Presence/Absence


################################################################################
# RANDOM FOREST FEATURE SELECTION CURVES & PERFORMANCES (figures in excel)
################################################################################
plot_fs <- function(fs, save){
        for (i in 1:length(cond)){
                sprintf("%s",cond[i]) ; sprintf("%s",new_cond[i])
                df <- fs[grep(cond[i],fs$ID),] # subset data corresponding to condition
                labels <- data.frame(FeatureNum=df[which(df$r2_val==max(df$r2_val)), "FeatureNum"], 
                        r2_val = df[which(df$r2_val==max(df$r2_val)), "r2_val"],
                        text = df[which(df$r2_val==max(df$r2_val)), "FeatureNum"]) #paste0("Max value at FeatureNum = ", 
                                #df[which(df$r2_val==max(df$r2_val)), "FeatureNum"], " and r2_val = ", max(df$r2_val)))
                ggplot(df) + 
                        geom_line(aes(x=FeatureNum, y=r2_test, color="r2_test")) + theme_bw(8) +
                        labs(title=new_cond[i], x="Features", y="Performance") + 
                        geom_line(aes(x=FeatureNum, y=r2_val, color="r2_val")) +
                        scale_color_manual(name="set", values=c("r2_test"="black", "r2_val"="red")) +
                        geom_text(data=labels, aes(x=FeatureNum, y=r2_val, label=text)) + 
                        annotate("segment",
                                x=df[which(df$r2_val==max(df$r2_val)), "FeatureNum"],
                                xend=df[which(df$r2_val==max(df$r2_val)), "FeatureNum"],
                                y=df[which(df$r2_val==max(df$r2_val)), "r2_val"]-0.01,
                                yend=df[which(df$r2_val==max(df$r2_val)), "r2_val"],
                                arrow=arrow(type = "closed", length = unit(0.02, "npc")), color="blue")
                ggsave(paste("Scripts/Data_Vis/",cond[i],save, sep=""), 
                        width=7, height=4, device="pdf", useDingbats=FALSE)
        }
}
get_opt_feat <- function(fs){
        out <- c()
        for (i in 1:length(cond)){
                df <- fs[grep(cond[i],fs$ID),] # subset data corresponding to condition
                #df <- right_join(conds, df, by=c("cond"="Y")) 
                feat <- df[df$r2_val==max(df$r2_val),]
                out <- rbind(out, feat)
        }
        return(out)
}

## PC models
path <- "/mnt/home/seguraab/Shiu_Lab/Project/Results/RESULTS_RF_PCs_sorted.txt"
rf_pc <- read.csv(path, sep="\t")
cond_order <- rf_pc$cond

## SNP Models
path <- "/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results/fs"
rf <- read.csv(file.path(path, "RESULTS_reg.txt"), sep="\t")
rf_fs <- rf[grep("rf_[0-9]+", rf$ID),] # SNP feature selection results
rf_fs <- right_join(conds, rf_fs, by=c("cond"="Y")) 
plot_fs(rf_fs, "_FS.pdf") # 2 to 40k features
rf_exp_fs <- rf_fs[grep("exp", rf_fs$ID),] 
plot_fs(rf_exp_fs, "_exp_FS.pdf") # 2 to 2048 features
rf_fs <- rbind(rf_fs, rf_exp_fs)
out <- get_opt_feat(rf_fs) # optimal models
sum(duplicated(out$cond)) # check for duplicate rows (0)
out <- out[!duplicated(out),] # remove duplicates
rownames(out) <- out$cond
out <- out[cond_order,] # set to PC env order
write.table(out, "Results/RESULTS_RF_SNPs_FS.txt", sep="\t", row.names=F, quote=F)

## ORF Copy Number Models
path <- "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results"
rf_orf <- read.csv(file.path(path, "RESULTS_reg.txt"), sep="\t")
rf_cnv_fs <- rf_orf[grep("_cno_[^b]", rf_orf$ID),]
plot_fs(rf_cnv_fs, "_cnv_FS.pdf")
out <- get_opt_feat(rf_cnv_fs) # optimal models
out <- right_join(conds, out, by=c("cond"="Y"))
rownames(out) <- out$cond
out <- out[cond_order,] # same order as PCs
write.table(out, "Results/RESULTS_RF_CNVs_FS.txt", sep="\t", row.names=F, quote=F)

## ORF Presence/Absence Models
rf_orf_fs <- rf_orf[grep("_orf_[0-9]+", rf_orf$ID),]
plot_fs(rf_orf_fs, "_orf_FS.pdf")
out <- get_opt_feat(rf_orf_fs) # optimal models
out <- right_join(conds, out, by=c("cond"="Y"))
rownames(out) <- out$cond
out <- out[cond_order,]
write.table(out, "Results/RESULTS_RF_ORFs_FS.txt", sep="\t", row.names=F, quote=F)

################################################################################
# RANDOM FOREST AFTER FEATURE SELECTION: TEST SET PERFORMANCE COMPARISONS
################################################################################
pcs <- read.table("Results/RESULTS_RF_PCs_sorted.txt", sep="\t", header=T, row.names=1)
snp <- read.table("Results/RESULTS_RF_SNPs_FS.txt", sep="\t", header=T, row.names=1)
orf <- read.table("Results/RESULTS_RF_ORFs_FS.txt", sep="\t", header=T, row.names=1)
cnv <- read.table("Results/RESULTS_RF_CNVs_FS.txt", sep="\t", header=T, row.names=1)

# Bar plots with test, validation and PC performances for SNPs, ORFs, and CNVs
rownames(pcs)==rownames(snp) # ensure all environments are in the same order
rownames(pcs)==rownames(orf)
rownames(pcs)==rownames(cnv)

##### added 05/08/2023 (I vaguely remember the rest of the code outside of this)
all <- cbind(pcs$r2_test, snp$r2_test, orf$r2_test, cnv$r2_test)
rownames(all) <- rownames(pcs)
colnames(all) <- c("PCs", "SNPs", "ORFs", "CNVs")
write.table(all, "Results/RESULTS_RF_ALL_FS.txt", sep="\t", quote=F)
all <- reshape2::melt(all)
all <- left_join(all, conds, by=c("Var1"="cond"))
all <- all[order(all$value, decreasing=T),]
h2 <- read.csv("Results/Heritability_h2_H2_sommer.csv") # narrow-sense heritability
h2 <- left_join(h2, conds, by=c("Conditions"="cond"))
ggplot(all, aes(x=reorder(new_cond, -value), y=Var2, fill=value)) +
        geom_tile(aes(fill=value)) +
        geom_text(aes(label=round(value, 2)), size=3, color="white") +
        scale_fill_viridis(discrete=F) + theme_minimal(base_size=12) +
        theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, color="black", face="bold"),
              axis.text.y=element_text(color="black", face="bold.italic"))
ggsave("Scripts/Data_Vis/RF_FS_All_Data_Types_R2_Test.pdf", width=13, height=4.5)
hm <- ggplot(all, aes(x=reorder(new_cond, -value), y=Var2, fill=value)) +
        geom_tile(aes(fill=value)) +
        geom_text(aes(label=round(value, 2)), size=8, color="white") +
        scale_fill_viridis(discrete=F) + theme_minimal(base_size=22) +
        theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, color="black", face="bold"),
              axis.text.y=element_text(color="black", face="bold.italic"),
              text=element_text(size=30))
ggMarginal(hm, data=h2, x=reorder(new_cond, -h2), y=h2, type="histogram")
grid.arrange(bar, hm, nrow = 2, ncol = 1, widths=c(35), heights=c(11, 4))
ggsave("Scripts/Data_Vis/RF_FS_All_Data_Types_R2_Test_for_poster.pdf", width=35, height=11)
#####

ggplot() + geom_bar(aes(x=new_cond,y=r2_val), snp, fill="#90CE4F", stat="identity") +
        geom_errorbar(aes(x=new_cond, ymin=r2_val-r2_val_sd, 
                ymax=r2_val+r2_val_sd), snp, width=0.2) +
        geom_point(aes(x=new_cond, y=r2_test), pcs, size=2.5, col="#702EA0") +
        geom_tile(aes(x=new_cond, y=r2_test, width=.9, height=.005), snp, fill="black") +
        theme_bw(10) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("Scripts/Data_Vis/snp_RF_performance_after_fs.pdf", width=8, height=4, units="in")
dev.off()

ggplot() + geom_bar(aes(x=new_cond,y=r2_val), orf, fill="#FFC000", stat="identity") +
        geom_errorbar(aes(x=new_cond, ymin=r2_val-r2_val_sd, 
                ymax=r2_val+r2_val_sd), orf, width=0.2) +
        geom_point(aes(x=pcs$new_cond, y=r2_test), pcs, size=2.5, col="#702EA0") +
        geom_tile(aes(x=new_cond, y=r2_test, width=.9, height=.005), orf, fill="black") +
        theme_bw(10) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("Scripts/Data_Vis/orf_RF_performance_after_fs.pdf", width=8, height=4, units="in")
dev.off()

ggplot() + geom_bar(aes(x=new_cond,y=r2_val), cnv, fill="#5B9BD5", stat="identity") +
        geom_errorbar(aes(x=new_cond, ymin=r2_val-r2_val_sd, 
                ymax=r2_val+r2_val_sd), cnv, width=0.2) +
        geom_point(aes(x=pcs$new_cond, y=r2_test), pcs, size=2.5, col="#702EA0") +
        geom_tile(aes(x=new_cond, y=r2_test, width=.9, height=.005), cnv, fill="black") +
        theme_bw(10) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("Scripts/Data_Vis/cnv_RF_performance_after_fs.pdf", width=8, height=4, units="in")
dev.off()

# Scatter plots comparing test performances between data types
# Combine performances
test <- list(pcs[c(1,31)], snp[c(1,31)], orf[c(1,31)], cnv[c(1,31)]) %>% 
        purrr::reduce(left_join, by="cond")
colnames(test) <- c("cond", "PCs", "SNPs", "ORFs", "CNVs")
test <- test[order(test$cond),] # ribose and glycerol out of order

# Correlations of test set performances between data types
rows <- test$cond
test.cor <- cor(as.matrix(t(test[-1]))) # pearson's correlation
colnames(test.cor) <- rownames(test.cor) <-  rows

# Heatmap of test set performance correlations
rdbu_r <- rev(brewer.pal(n=11, "RdBu")) # reversed color palette
col_fun <- colorRamp2(seq(-1,1,.2), rdbu_r) #same as cm
pdf("Scripts/Data_Vis/RF_test_correlations.pdf", height=9, width=9)
hm <- Heatmap(test.cor, col=col_fun, show_row_dend=F, show_column_dend=F,
        heatmap_legend_param=list(title="PCC", at=seq(-1,1,.2), color_bar="continuous"))
draw(hm)
dev.off()

# Differences in performance between data types
test$`SNPs-PCs` <- test$SNPs-test$PCs
test$`ORFs-PCs` <- test$ORFs-test$PCs
test$`CNVs-PCs` <- test$CNVs-test$PCs
test$`ORFs-SNPs` <- test$ORFs-test$SNPs
test$`CNVs-SNPs` <- test$CNVs-test$SNPs
test$`CNVs-ORFs` <- test$CNVs-test$ORFs

# Percent performance increase or decrease between data types
test$`SNPs/PCs` <- (test$SNPs/test$PCs-1)*100 # should I not convert to percent?
test$`ORFs/PCs` <- (test$ORFs/test$PCs-1)*100
test$`CNVs/PCs` <- (test$CNVs/test$PCs-1)*100
test$`ORFs/SNPs` <- (test$ORFs/test$SNPs-1)*100
test$`CNVs/SNPs` <- (test$CNVs/test$SNPs-1)*100
test$`CNVs/ORFs` <- (test$CNVs/test$ORFs-1)*100
test$`log_SNPs/PCs` <- ifelse(test$`SNPs/PCs`< 0, -log10(abs(test$`SNPs/PCs`)), log10(test$`SNPs/PCs`))
test$`log_ORFs/PCs` <- ifelse(test$`ORFs/PCs` < 0, -log10(abs(test$`ORFs/PCs`)), log10(test$`ORFs/PCs`))
test$`log_CNVs/PCs` <- ifelse(test$`CNVs/PCs` < 0, -log10(abs(test$`CNVs/PCs`)), log10(test$`CNVs/PCs`))
test$`log_ORFs/SNPs` <- ifelse(test$`ORFs/SNPs` < 0, -log10(abs(test$`ORFs/SNPs`)), log10(test$`ORFs/SNPs`))
test$`log_CNVs/SNPs` <- ifelse(test$`CNVs/SNPs` < 0, -log10(abs(test$`CNVs/SNPs`)), log10(test$`CNVs/SNPs`))
test$`log_CNVs/ORFs` <- ifelse(test$`CNVs/ORFs` < 0, -log10(abs(test$`CNVs/ORFs`)), log10(test$`CNVs/ORFs`))
# write.table(test, "Results/RESULTS_RF_FS_ALL.txt", sep="\t", quote=F, row.names=F)

# Heatmap of performance comparisons
rownames(test) <- test$cond
rdbu_r <- rev(brewer.pal(n=11, "RdBu")) # reversed color palette
col_fun <- colorRamp2(seq(-4,4,.8), rdbu_r) #same as cm
# Log of percent difference in performance
pdf("Scripts/Data_Vis/RF_FS_test_log_percent_difference.pdf", height=9, width=9)
hm <- Heatmap(as.matrix(test[18:23]), col=col_fun, show_row_dend=F, show_column_dend=F,
        heatmap_legend_param=list(title="Log of percent difference", at=seq(-4,4,.8), color_bar="continuous"))
draw(hm)
dev.off()
# Subtracted difference in performance
col_fun <- colorRamp2(seq(-0.4,0.4,0.08), rdbu_r) #same as cm
pdf("Scripts/Data_Vis/RF_FS_test_subtracted_difference.pdf", height=9, width=9)
hm <- Heatmap(as.matrix(test[6:11]), col=col_fun, show_row_dend=F, show_column_dend=F,
        heatmap_legend_param=list(title="R^2 Difference", at=seq(-0.4,0.4,0.08), color_bar="continuous"))
draw(hm)
dev.off()

# Plot feature types against each other (The function is not working)
plot_lm <- function(df, save){
        formula <- y ~ x
        pdf(save)
        p <- ggplot(df, aes(x, y)) +
                geom_point() +
                stat_quant_line(formula = formula, quantiles = 0.5) +
                stat_quant_eq(formula = formula, quantiles = 0.5)
        dev.off()
        return(p)
}
df <- cbind(pcs$r2_test, snp$r2_test); colnames(df) <- c("x", "y")
plot_lm(as.data.frame(df), "Results/PCs_v_SNPs_RF_FS.pdf")
df <- cbind(pcs$r2_test, orf$r2_test); colnames(df) <- c("x", "y")
plot_lm(as.data.frame(df), "Results/PCs_v_ORFs_RF_FS.pdf")
df <- cbind(pcs$r2_test, cnv$r2_test); colnames(df) <- c("x", "y")
plot_lm(as.data.frame(df), "Results/PCs_v_CNVs_RF_FS.pdf")
df <- cbind(snp$r2_test, orf$r2_test); colnames(df) <- c("x", "y")
plot_lm(as.data.frame(df), "Results/SNPs_v_ORFs_RF_FS.pdf")
df <- cbind(snp$r2_test, cnv$r2_test); colnames(df) <- c("x", "y")
plot_lm(as.data.frame(df), "Results/SNPs_v_CNVs_RF_FS.pdf")


################################################################################
# HERITABILITY & RF FS PERFORMANCE COMPARISON
################################################################################
h2 <- read.csv("Results/Heritability_h2_H2_sommer.csv") # heritability data
h2 <- merge(conds, h2, by.x="cond", by.y="Conditions")
h2 <- h2[order(h2$h2, descending=True),] # sort in ascending order
h2$new_cond <- factor(h2$new_cond, levels = h2$new_cond)
test <- read.csv("Results/RESULTS_RF_FS_ALL.txt", sep="\t") # performance data
h2 <- merge(h2, test, by.x="cond", by.y="cond")

# Plot narrow-sense heritability
ggplot(h2, aes(x=new_cond, y=h2, fill="#F8766D")) + theme_bw(8) +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_errorbar(aes(ymin=h2-h2_SE, ymax=h2+h2_SE), width=0.2,
                position=position_dodge(0.9)) +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        theme(axis.text.y = element_text(size=9, color="black"))
ggsave("Scripts/Data_Vis/h2_environments.pdf", width=10, height=4, device="pdf", useDingbats=FALSE)

# Plot broad-sense heritability
ggplot(h2, aes(x=new_cond, y=H2_AD, fill="#F8766D")) + theme_bw(8) +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_errorbar(aes(ymin=H2_AD-H2_AD_SE, ymax=H2_AD+H2_AD_SE), width=0.2,
                position=position_dodge(0.9)) +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        theme(axis.text.y = element_text(size=9, color="black"))
ggsave("Scripts/Data_Vis/H2_AD_environments.pdf", width=10, height=4, device="pdf", useDingbats=FALSE)

# Narrow-sense heritability vs RF FS performance
ggplot(h2, aes(x=SNPs, y=h2, col=new_cond, shape=new_cond)) + geom_point() + 
        scale_shape_manual(values=1:nlevels(h2$new_cond)/2) +
        theme_bw(8) + geom_point(size=3) + 
        #geom_smooth(mapping=aes(SNPs, h2), method="lm", formula=h2~SNPs, stat="smooth") +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        theme(axis.text.y = element_text(size=9, color="black"))
ggsave("Scripts/Data_Vis/h2_vs_test_R2_RF_FS.pdf", width=8, height=5, device="pdf", useDingbats=FALSE)

ggplot(h2, aes(x=SNPs, y=h2)) + geom_point() + 
        theme_bw(8) + geom_point(size=3) + 
        geom_smooth(method="lm") +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        theme(axis.text.y = element_text(size=9, color="black"))
ggsave("Scripts/Data_Vis/h2_vs_test_R2_RF_FS_lm.pdf", width=5, height=5, device="pdf", useDingbats=FALSE)

################################################################################
# BASELINE RF GENE IMPORTANCE RANK DENSITY SCATTER PLOTS COMPARING DATA TYPES
################################################################################
# paths to feature importance score files
dir <- "/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results/baseline"
snp_files <- unlist(lapply(cond, 
        function(x) list.files(dir, paste(x, "rf_baseline_imp", sep="_"),
        full.names=T)))
dir <- "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results/baseline"
orf_files <- unlist(lapply(cond, 
        function(x) list.files(dir, paste(x, "orf_baseline_imp", sep="_"),
        full.names=T)))
cnv_files <- unlist(lapply(cond, 
        function(x) list.files(dir, paste(x, "cno_baseline_imp", sep="_"),
        full.names=T)))

# read marker to gene map files
map_snps <- read.csv(
        "Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes.txt", header=F)
map_orfs <- read.csv("Data/Peter_2018/ORFs_gene_map.tsv", sep="\t")

# function to calculate point density
get_density <- function(x, y, ...){
        dens <- MASS::kde2d(x, y, ...) # 2D kernel density estimation
        ix <- findInterval(x, dens$x)
        iy <- findInterval(y, dens$y)
        ii <- cbind(ix, iy)
        return(dens$z[ii])
}

plot_density <- function(data, x, y, density, corr, px, py, xlab, ylab, title){
        # Plot dashed lines for two percentile cutoffs per omic type (4 total lines)
        p <- ggplot(data) + geom_point(aes(x=x, y=y, color=density),
                                       shape=16, alpha=0.75, size=1) + 
                scale_color_viridis() + theme_bw() +
                geom_hline(yintercept=py, color="blue", linetype="dashed") + 
                geom_vline(xintercept=px, color="blue", linetype="dashed") +
                ggtitle(paste(title, "; rho =", round(corr$estimate, 3), "; P-value =",
                              round(corr$p.value, 3))) +
                xlab(xlab) + ylab(ylab) + 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"))
    return (p)
}

plot_contour <- function(data, x, y, density, corr, px, py, xlab, ylab, title){
        # Plot dashed lines for two percentile cutoffs per omic type (4 total lines)
        p <- ggplot(data) + 
                geom_point(aes(x=x, y=y), shape=16, alpha=0.75, size=1) + 
                geom_density_2d_filled(aes(x=x, y=y)) + theme_bw() +
                ggtitle(paste(title, "; rho =", round(corr$estimate, 3),
                        "; P-value =", round(corr$p.value, 3))) +
                scale_x_continuous(limits = c(0, 1)) + # for normalized ranks
                scale_y_continuous(limits = c(0, 1)) +
                # scale_x_continuous(limits = c(0, 4)) + # for log10 ranks
                # scale_y_continuous(limits = c(0, 4)) +
                xlab(xlab) + ylab(ylab) + 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"))
        return(p)
}

# create a snp vs orf, snp vs cnv, orf vs cnv plot for each environment
correlations <- c("Comparison", "Environment", "rho", "p-value")
snp_ranks <- unique(map_snps["V4"])
snp_imps <- unique(map_snps["V4"])
orf_ranks <- unique(map_orfs["gene"])
orf_imps <- unique(map_orfs["gene"])
cnv_ranks <- unique(map_orfs["gene"])
cnv_imps <- unique(map_orfs["gene"])
source("http://peterhaschke.com/Code/multiplot.R") #load multiplot function
plots_rnorm_sxo <- vector('list', length(cond)) # list of plots, SNPs vs ORFs
plots_rnorm_sxc <- vector('list', length(cond)) # SNPs vs CNVs
plots_rnorm_oxc <- vector('list', length(cond)) # ORFs vs CNVs
i <- 1
for (env in cond){
        print(i)
        # read feature importance score files
        snp <- read.csv(snp_files[grep(env, snp_files, fixed=T)], sep="\t")
        orf <- read.csv(orf_files[grep(env, orf_files, fixed=T)], sep="\t")
        cnv <- read.csv(cnv_files[grep(env, cnv_files, fixed=T)], sep="\t")
        
        # map markers to genes
        snp <- left_join(snp, map_snps[c("V1", "V4")], by=c("X"="V1"))
        orf <- left_join(orf, map_orfs[c("qacc", "gene")], by=c("X"="qacc"))
        cnv <- left_join(cnv, map_orfs[c("qacc", "gene")], by=c("X"="qacc"))

        # get max importance score per gene
        snp <- snp %>% group_by(V4) %>% dplyr::summarise(max=max(mean_imp)) %>% as.data.frame()
        orf <- orf %>% group_by(gene) %>% dplyr::summarise(max=max(mean_imp)) %>% as.data.frame()
        cnv <- cnv %>% group_by(gene) %>% dplyr::summarise(max=max(mean_imp)) %>% as.data.frame()

        # drop rows with missing annotations
        snp <- stats::na.omit(snp)
        snp <- snp[which(snp$V4!="intergenic"),]
        orf <- stats::na.omit(orf)
        cnv <- stats::na.omit(cnv)

        # rank genes
        snp <- snp[order(snp$max, decreasing=T),]
        snp$rank <- 1:nrow(snp)
        orf <- orf[order(orf$max, decreasing=T),]
        orf$rank <- 1:nrow(orf)
        cnv <- cnv[order(cnv$max, decreasing=T),]
        cnv$rank <- 1:nrow(cnv)

        # normalize ranks
        snp$rank_norm <- 1 - ((snp$rank - 1)/(nrow(snp) - 1))
        orf$rank_norm <- 1 - ((orf$rank - 1)/(nrow(orf) - 1))
        cnv$rank_norm <- 1 - ((cnv$rank - 1)/(nrow(cnv) - 1))

        # log of the rank
        snp$log_rank <- log10(snp$rank)
        orf$log_rank <- log10(orf$rank)
        cnv$log_rank <- log10(cnv$rank)

        # take the intersection of the 3 data types
        gene_ranks <- full_join(orf[c("gene", "rank")],
                snp[c("V4", "rank")], by=c("gene"="V4"))
        gene_ranks <- full_join(cnv[c("gene", "rank")], gene_ranks, by="gene")
        colnames(gene_ranks) <- c("gene", "CNVs", "ORFs", "SNPs")
        gene_ranks <- stats::na.omit(gene_ranks)

        gene_rnorm <- full_join(orf[c("gene", "rank_norm")],
                snp[c("V4", "rank_norm")], by=c("gene"="V4"))
        gene_rnorm <- full_join(cnv[c("gene", "rank_norm")], gene_rnorm, by="gene")
        colnames(gene_rnorm) <- c("gene", "CNVs", "ORFs", "SNPs")
        gene_rnorm <- stats::na.omit(gene_rnorm)

        gene_rlog <- full_join(orf[c("gene", "log_rank")],
                snp[c("V4", "log_rank")], by=c("gene"="V4"))
        gene_rlog <- full_join(cnv[c("gene", "log_rank")], gene_rlog, by="gene")
        colnames(gene_rlog) <- c("gene", "CNVs", "ORFs", "SNPs")
        gene_rlog <- stats::na.omit(gene_rlog)

        # compare gene ranks and importance scores between conditions
        snp_ranks <- left_join(snp_ranks, snp[c("V4", "rank")], by="V4")
        snp_imps <- left_join(snp_imps, snp[c("V4", "max")], by="V4")
        orf_ranks <- left_join(orf_ranks, orf[c("gene", "rank")], by="gene")
        orf_imps <- left_join(orf_imps, orf[c("gene", "max")], by="gene")
        cnv_ranks <- left_join(cnv_ranks, cnv[c("gene", "rank")], by="gene")
        cnv_imps <- left_join(cnv_imps, cnv[c("gene", "max")], by="gene")

        # create table of rank correlations between data types
        sxo_corr = cor.test(gene_ranks$ORFs, gene_ranks$SNPs, method='spearman')
        sxc_corr = cor.test(gene_ranks$CNVs, gene_ranks$SNPs, method='spearman')
        oxc_corr = cor.test(gene_ranks$CNVs, gene_ranks$ORFs, method='spearman')
        correlations <- rbind(correlations, list("SNPs vs ORFs", env, sxo_corr$estimate, sxo_corr$p.value))
        correlations <- rbind(correlations, list("SNPs vs CNVs", env, sxc_corr$estimate, sxc_corr$p.value))
        correlations <- rbind(correlations, list("ORFs vs CNVs", env, oxc_corr$estimate, oxc_corr$p.value))

        # calculate density of normalized ranks & max importance scores per pair of data types
        # gene_ranks$`SNPs vs ORFs` <- get_density(gene_ranks$ORFs, gene_ranks$SNPs, n=100)
        # gene_ranks$`SNPs vs CNVs` <- get_density(gene_ranks$CNVs, gene_ranks$SNPs, n=100)
        # gene_ranks$`ORFs vs CNVs` <- get_density(gene_ranks$CNVs, gene_ranks$ORFs, n=100)
        
        gene_rnorm$`SNPs vs ORFs` <- get_density(gene_rnorm$ORFs, gene_rnorm$SNPs, n=100)
        gene_rnorm$`SNPs vs CNVs` <- get_density(gene_rnorm$CNVs, gene_rnorm$SNPs, n=100)
        gene_rnorm$`ORFs vs CNVs` <- get_density(gene_rnorm$CNVs, gene_rnorm$ORFs, n=100)

        gene_rlog$`SNPs vs ORFs` <- get_density(gene_rlog$ORFs, gene_rlog$SNPs, n=100)
        gene_rlog$`SNPs vs CNVs` <- get_density(gene_rlog$CNVs, gene_rlog$SNPs, n=100)
        gene_rlog$`ORFs vs CNVs` <- get_density(gene_rlog$CNVs, gene_rlog$ORFs, n=100)
        
        # contour density plot of gene ranks
        # p1 <- plot_contour(data=gene_ranks, x=gene_ranks$ORFs,
        #         y=gene_ranks$SNPs, density=gene_ranks$`SNPs vs ORFs`,
        #         corr=sxo_corr,
        #         px=quantile(gene_ranks$ORFs, 0.95), xlab="ORFs Gene Rank",
        #         py=quantile(gene_ranks$SNPs, 0.95), ylab="SNPs Gene Rank", title=env)
        # plots_rank_sxo[[i]] <- p1
        # ggsave(plot=p1, paste("Scripts/Data_Vis/SNPs_vs_ORFs_gene_rank_contour_", env, ".pdf", sep=""))
        # # graphics.off()
        # p2 <- plot_contour(data=gene_ranks, x=gene_ranks$CNVs,
        #         y=gene_ranks$SNPs, density=gene_ranks$`SNPs vs CNVs`,
        #         corr=sxc_corr,
        #         px=quantile(gene_ranks$CNVs, 0.95), xlab="CNVs Gene Rank",
        #         py=quantile(gene_ranks$SNPs, 0.95), ylab="SNPs Gene Rank", title=env)
        # plots_rank_sxc[[i]] <- p2
        # ggsave(plot=p2, paste("Scripts/Data_Vis/SNPs_vs_CNVs_gene_rank_contour_", env, ".pdf", sep=""))
        # # graphics.off()
        # p3 <- plot_contour(data=gene_ranks, x=gene_ranks$CNVs,
        #         y=gene_ranks$ORFs, density=gene_ranks$`ORFs vs CNVs`,
        #         corr=oxc_corr,
        #         px=quantile(gene_ranks$CNVs, 0.95), xlab="CNVs Gene Rank",
        #         py=quantile(gene_ranks$CNVs, 0.95), ylab="ORFs Gene Rank", title=env)
        # plots_rank_oxc[[i]] <- p3
        # ggsave(plot=p3, paste("Scripts/Data_Vis/ORFs_vs_CNVs_gene_rank_contour_", env, ".pdf", sep=""))
        # graphics.off()
        
        # density plot of normalized gene ranks
        p4 <- plot_contour(data=gene_rnorm, x=gene_rnorm$ORFs,
                y=gene_rnorm$SNPs, density=gene_rnorm$`SNPs vs ORFs`,
                corr=cor.test(gene_rnorm$ORFs, gene_rnorm$SNPs, method='spearman'),
                px=quantile(gene_rnorm$ORFs, 0.95), xlab="ORFs Normalized Gene Rank",
                py=quantile(gene_rnorm$SNPs, 0.95), ylab="SNPs Normalized Gene Rank", title=env)
        plots_rnorm_sxo[[i]] <- p4
        ggsave(plot=p4, paste("Scripts/Data_Vis/SNPs_vs_ORFs_gene_rnorm_contour_", env, ".pdf", sep=""), width=8)
        p5 <- plot_contour(data=gene_rnorm, x=gene_rnorm$CNVs,
                y=gene_rnorm$SNPs, density=gene_rnorm$`SNPs vs CNVs`,
                corr=cor.test(gene_rnorm$CNVs, gene_rnorm$SNPs, method='spearman'),
                px=quantile(gene_rnorm$CNVs, 0.95), xlab="CNVs Normalized Gene Rank",
                py=quantile(gene_rnorm$SNPs, 0.95), ylab="SNPs Normalized Gene Rank", title=env)
        plots_rnorm_sxc[[i]] <- p5
        ggsave(plot=p5, paste("Scripts/Data_Vis/SNPs_vs_CNVs_gene_rnorm_contour_", env, ".pdf", sep=""), width=8)
        p6 <- plot_contour(data=gene_rnorm, x=gene_rnorm$CNVs,
                y=gene_rnorm$ORFs, density=gene_rnorm$`ORFs vs CNVs`,
                corr=cor.test(gene_rnorm$CNVs, gene_rnorm$ORFs, method='spearman'),
                px=quantile(gene_rnorm$CNVs, 0.95), xlab="CNVs Normalized Gene Rank",
                py=quantile(gene_rnorm$CNVs, 0.95), ylab="ORFs Normalized Gene Rank", title=env)
        plots_rnorm_oxc[[i]] <- p6
        ggsave(plot=p6, paste("Scripts/Data_Vis/ORFs_vs_CNVs_gene_rnorm_contour_", env, ".pdf", sep=""), width=8)

        # density plot of log gene ranks
        # pdf(paste("Scripts/Data_Vis/SNPs_vs_ORFs_gene_rlog_contour_", env, ".pdf", sep=""), width=8)
        # p8 <- plot_contour(data=gene_rlog, x=gene_rlog$ORFs,
        #         y=gene_rlog$SNPs, density=gene_rlog$`SNPs vs ORFs`,
        #         corr=cor.test(gene_rlog$ORFs, gene_rlog$SNPs, method='spearman'),
        #         px=quantile(gene_rlog$ORFs, 0.95), xlab="ORFs Log10 Gene Rank",
        #         py=quantile(gene_rlog$SNPs, 0.95), ylab="SNPs Log10 Gene Rank", title=env)
        # print(p8)
        # dev.off()
        # pdf(paste("Scripts/Data_Vis/SNPs_vs_CNVs_gene_rlog_contour_", env, ".pdf", sep=""), width=8)
        # p9 <- plot_contour(data=gene_rlog, x=gene_rlog$CNVs,
        #         y=gene_rlog$SNPs, density=gene_rlog$`SNPs vs CNVs`,
        #         corr=cor.test(gene_rlog$CNVs, gene_rlog$SNPs, method='spearman'),
        #         px=quantile(gene_rlog$CNVs, 0.95), xlab="CNVs Log10 Gene Rank",
        #         py=quantile(gene_rlog$SNPs, 0.95), ylab="SNPs Log10 Gene Rank", title=env)
        # print(p9)
        # dev.off()
        # pdf(paste("Scripts/Data_Vis/ORFs_vs_CNVs_gene_rlog_contour_", env, ".pdf", sep=""), width=8)
        # p10 <- plot_contour(data=gene_rlog, x=gene_rlog$CNVs,
        #         y=gene_rlog$ORFs, density=gene_rlog$`ORFs vs CNVs`,
        #         corr=cor.test(gene_rlog$CNVs, gene_rlog$ORFs, method='spearman'),
        #         px=quantile(gene_rlog$CNVs, 0.95), xlab="CNVs Log10 Gene Rank",
        #         py=quantile(gene_rlog$CNVs, 0.95), ylab="ORFs Log10 Gene Rank", title=env)
        # print(p10)
        # dev.off()

        i <- i + 1
        # rm(snp, orf, cnv, gene_ranks, gene_rnorm, gene_rlog, gene_imps)
        if (i!=35) rm(snp, orf, cnv, gene_ranks, gene_rnorm, gene_rlog)
}
correlations <- data.frame(correlations)
write.table(as.matrix(correlations), "Scripts/Data_Vis/Gene_rank_rho_between_data_types.tsv", sep="\t", quote=F, row.names=F)
# Plot contour plots between data types
pdf("Scripts/Data_Vis/SNPs_vs_ORFs_gene_rank_norm_contour_ALL.pdf", height=20, width=20) # normalized gene ranks
multiplot(plotlist=plots_rnorm_sxo, cols = 5)
dev.off()
pdf("Scripts/Data_Vis/SNPs_vs_CNVs_gene_rank_norm_contour_ALL.pdf", height=20, width=20)
multiplot(plotlist=plots_rnorm_sxc, cols = 5)
dev.off()
pdf("Scripts/Data_Vis/ORFs_vs_CNVs_gene_rank_norm_contour_ALL.pdf", height=20, width=20)
multiplot(plotlist=plots_rnorm_oxc, cols = 5)
dev.off()

# Calculate correlations between environments per data type
colnames(snp_ranks) <- c("gene", cond) # reset column names
colnames(snp_imps) <- c("gene", cond)
colnames(orf_ranks) <- c("gene", cond)
colnames(orf_imps) <- c("gene", cond)
colnames(cnv_ranks) <- c("gene", cond)
colnames(cnv_imps) <- c("gene", cond)
snp_ranks <- snp_ranks[which(snp_ranks$gene!="intergenic"),] # remove intergenic snps
snp_imps <- snp_imps[which(snp_imps$gene!="intergenic"),]
snp_ranks_cor <- reshape2::melt(cor(snp_ranks[-1], method="spearman")) # calculate spearman's rho
snp_ranks_cor <- snp_ranks_cor[!duplicated(t(apply(snp_ranks_cor[1:2], 1, sort))),] # keep only unique pairs
snp_imps_cor <- reshape2::melt(cor(snp_imps[-1], method="spearman"))
snp_imps_cor <- snp_imps_cor[!duplicated(t(apply(snp_imps_cor[1:2], 1, sort))),]
orf_ranks_cor <- reshape2::melt(cor(orf_ranks[-1], method="spearman"))
orf_ranks_cor <- orf_ranks_cor[!duplicated(t(apply(orf_ranks_cor[1:2], 1, sort))),]
orf_imps_cor <- reshape2::melt(cor(orf_imps[-1], method="spearman"))
orf_imps_cor <- orf_imps_cor[!duplicated(t(apply(orf_imps_cor[1:2], 1, sort))),]
cnv_ranks_cor <- reshape2::melt(cor(cnv_ranks[-1], method="spearman"))
cnv_ranks_cor <- cnv_ranks_cor[!duplicated(t(apply(cnv_ranks_cor[1:2], 1, sort))),]
cnv_imps_cor <- reshape2::melt(cor(cnv_imps[-1], method="spearman"))
cnv_imps_cor <- cnv_imps_cor[!duplicated(t(apply(cnv_imps_cor[1:2], 1, sort))),]
write.table(snp_ranks_cor, "Scripts/Data_Vis/SNP_gene_ranks_rho.tsv", sep="\t", quote=F, row.names=F)
write.table(snp_imps_cor, "Scripts/Data_Vis/SNP_gene_imps_rho.tsv", sep="\t", quote=F, row.names=F)
write.table(orf_ranks_cor, "Scripts/Data_Vis/ORF_gene_ranks_rho.tsv", sep="\t", quote=F, row.names=F)
write.table(orf_imps_cor, "Scripts/Data_Vis/ORF_gene_imps_rho.tsv", sep="\t", quote=F, row.names=F)
write.table(cnv_ranks_cor, "Scripts/Data_Vis/CNV_gene_ranks_rho.tsv", sep="\t", quote=F, row.names=F)
write.table(cnv_imps_cor, "Scripts/Data_Vis/CNV_gene_imps_rho.tsv", sep="\t", quote=F, row.names=F)


plot_cor_hm <- function(df, save){
        ggplot(df, aes(Var1, Var2, fill=round(value, 2))) +
        geom_tile(color="white") +
        scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0,
                limit=c(-1,1), space="Lab", name="Spearman's\nrho") +
        theme_bw() + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
                axis.title.x=element_blank(), axis.title.y=element_blank(),
                panel.grid.major=element_blank(), panel.border=element_blank(),
                panel.background=element_blank(), axis.ticks=element_blank()) +
        coord_fixed()
        ggsave(save)
}
plot_cor_hm(snp_ranks_cor, "Scripts/Data_Vis/SNP_gene_ranks_rho.pdf")
plot_cor_hm(snp_imps_cor, "Scripts/Data_Vis/SNP_gene_imps_rho.pdf")
plot_cor_hm(orf_ranks_cor, "Scripts/Data_Vis/ORF_gene_ranks_rho.pdf")
plot_cor_hm(orf_imps_cor, "Scripts/Data_Vis/ORF_gene_imps_rho.pdf")
plot_cor_hm(cnv_ranks_cor, "Scripts/Data_Vis/CNV_gene_ranks_rho.pdf")
plot_cor_hm(cnv_imps_cor, "Scripts/Data_Vis/CNV_gene_imps_rho.pdf")

################################################################################
# GO ENRICHMENT (RF MODELS AFTER FS)
################################################################################
## SNPs
dir <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/SNPs/fs" # path to ORA results files
files <- list.files(path=dir, pattern="^ORA_SHAP_", full.names=T, recursive=F)

# read in the first file
trait <- str_extract(files[1], "YP[A-Z0-9]+[^_]")
df <- read.csv(files[1], sep="\t")
df$Env <- trait # add environment to keep track where the GO term is enriched
df <- df[df$lfdr < 0.05,] # filter

# combine with remaining environments
for (i in 2:length(files)){
    trait <- str_extract(files[i], "YP[A-Z0-9]+[^_]") # environment name
    df2 <- read.csv(files[i], sep="\t")
    df2$Env <- trait
    df2 <- df2[df2$lfdr < 0.05,]
    df <- rbind.fill(df, df2)
}

# take the log of the q-values
df$logQ <- 0
df <- df[order(df$qvalues, decreasing=F),]
for(i in 1:nrow(df)){ 
        if(df$direction[i] == '-') df$logQ[i] <- log10(df$qvalues[i])
        if(df$direction[i] == '+') df$logQ[i] <- -log10(df$qvalues[i])
}
if (nrow(df[df$logQ < -10,])!=0) df[df$loggQ < -10,]$logQ <- -10 # set limits for extreme values
if (nrow(df[df$logQ > 10,])!=0) df[df$logQ > 10,]$logQ <- 10

# reshape
df_logQ <- df %>% group_by(Env) %>% 
        pivot_wider(id_cols=GO, names_from=Env, values_from=logQ) %>% 
        as.data.frame()
df_logQ <- left_join(df_logQ, df[c('GO','BP','CC','MF')], by='GO')
df_logQ <- df_logQ[!duplicated(df_logQ),]
# save and manually add GO:0017056 MF and GO:0005762 CC using amiGO 2, delete extra duplicated rows
write.table(df_logQ, 
        "Scripts/Data_Vis/SNPs_GO_logQ.tsv", sep="\t", quote=F, row.names=F)

df_logQ <- read.csv("Scripts/Data_Vis/SNP_Figures/SNPs_GO_logQ.tsv", sep="\t") # manually cleaned
df_logQ[is.na(df_logQ)] <- 0 # set insignificant q-values to zero
df_logQ <- df_logQ[which(df_logQ$YPD6AU >= 0),] # drop rows with negative values
df_logQ <- df_logQ %>% dplyr::select(which(apply(df_logQ, 2, n_distinct) > 1)) # drop columns with only zeroes
new_order <- c("GO", sort(colnames(df_logQ[2:8])), "BP", "CC", "MF") # sort columns
df_logQ <- df_logQ[,new_order]

# get genes annotated with these go terms
path <- "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_go.tsv"
genes <- read.csv(path, sep="\t")
df_logQ <- left_join(df_logQ, genes[c(4,9)], by=c("GO"="GO.ID"))
df_logQ <- df_logQ[!duplicated(df_logQ),] # drop duplicate rows due to many snps per gene
write.csv(df_logQ, "Scripts/Data_Vis/SNPs_GO_BP_CC_MF_logQ_genes.csv", quote=F, row.names=F)

# separate into GO categories
BP <- df_logQ[which(df_logQ$BP!=""),-c(1,10,11)]
rownames(BP) <- BP$BP
BP <- BP[-8]

CC <- df_logQ[which(df_logQ$CC!=""),-c(1,9,11)]
rownames(CC) <- CC$CC
CC <- CC[-8]

MF <- df_logQ[which(df_logQ$MF!=""),-c(1,9,10)]
rownames(MF) <- MF$MF
MF <- MF[-8]

# heatmaps (red=enriched, blue=underrepresented)
plot_hm <- function(df, range){
        # Function to plot heatmap
        p <- Heatmap(as.matrix(df), 
                col=col_fun,
                na_col="grey",
                cluster_rows=F, 
                show_row_dend=F,
                cluster_columns=F,
                show_column_dend=F,
                row_names_side="left",
                row_names_gp = gpar(fontsize=8),
                column_names_gp = gpar(fontize=8),
                border_gp = gpar(col="black", lty=1),
                #rect_gp = gpar(col="black", lty=1),
                name="log10(q-values)",
                width=ncol(df)*unit(0.5, "cm"),
                height=nrow(df)*unit(0.5, "cm"),
                heatmap_legend_param=list(title="log10(qval)", 
                at=range))
        return (p)
}
palette_func <- colorRampPalette(c("red", "white"))
palette <- rev(palette_func(2)) 
col_fun <- colorRamp2(c(0, max(df_logQ[2:8])), palette)
pdf("Scripts/Data_Vis/SNPs_GO_BP_logQ.pdf", height=unit(9, "in"), width=unit(10, "in"))
p1 <- plot_hm(BP, range=c(0, max(df_logQ[2:8])))
p2 <- plot_hm(CC, range=c(0, max(df_logQ[2:8])))
p3 <- plot_hm(MF, range=c(0, max(df_logQ[2:8])))
draw(p1 %v% p2 %v% p3)
dev.off()

## ORFs
# see /mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/ORF_gene_set_enrichment.R
# line 442
dir <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs" # path to ORA results files
files <- list.files(path=dir, pattern="^ORA_SHAP_", full.names=T, recursive=F)

# separate into orf and cnv data types
files_orf <- c()
files_cnv <- c()
for (i in 1:length(files)){
        if (grepl("_orf_", files[i])) files_orf <- c(files_orf, files[i])
        if (grepl("_cno_", files[i])) files_cnv <- c(files_cnv, files[i])
}

combine <- function(files, path){
        # function to combine files in files list into a dataframe
        # read in the first file
        trait <- str_extract(files[1], "YP[A-Z0-9]+[^_]")
        df <- read.csv(files[1], sep="\t")
        df$Env <- trait # add environment to keep track where the GO term is enriched
        df <- df[df$lfdr < 0.05,] # filter

        # combine with remaining environments
        for (i in 2:length(files)){
                trait <- str_extract(files[i], "YP[A-Z0-9]+[^_]") # environment name
                df2 <- read.csv(files[i], sep="\t")
                df2$Env <- trait
                df2 <- df2[df2$lfdr < 0.05,]
                df <- rbind.fill(df, df2)
        }

        # take the log of the q-values
        df$logQ <- 0
        df <- df[order(df$qvalues, decreasing=F),]
        for(i in 1:nrow(df)){ 
                if(df$direction[i] == '-') df$logQ[i] <- log10(df$qvalues[i])
                if(df$direction[i] == '+') df$logQ[i] <- -log10(df$qvalues[i])
        }
        if (nrow(df[df$logQ < -10,])!=0) df[df$loggQ < -10,]$logQ <- -10 # set limits for extreme values
        if (nrow(df[df$logQ > 10,])!=0) df[df$logQ > 10,]$logQ <- 10

        # reshape
        df_logQ <- df %>% group_by(Env) %>% 
                pivot_wider(id_cols=GO, names_from=Env, values_from=logQ) %>% 
                as.data.frame()
        df_logQ <- left_join(df_logQ, df[c('GO','BP','CC','MF')], by='GO')
        df_logQ <- df_logQ[!duplicated(df_logQ),]
        # save and manually add missing descriptions of GO terms using amiGO 2, delete extra duplicated rows
        write.table(df_logQ, path, sep="\t", quote=F, row.names=F)
}

df_logQ_orf <- combine(files_orf, "Scripts/Data_Vis/ORFs_pres_abs_GO_logQ.tsv")
df_logQ_cnv <- combine(files_cnv, "Scripts/Data_Vis/ORFs_copy_num_GO_logQ.tsv")

# heatmaps for ORFs presence/absence
df_logQ <- read.csv("Scripts/Data_Vis/ORF_CNO_Figures/ORFs_pres_abs_GO_logQ.tsv", sep="\t") # manually cleaned
df_logQ[is.na(df_logQ)] <- 0 # set insignificant q-values to zero
df_logQ <- df_logQ[which(df_logQ$YPDANISO10 >= 0),] # drop rows with negative values (underreprsented GO terms)
df_logQ <- df_logQ %>% dplyr::select(which(apply(df_logQ, 2, n_distinct) > 1)) # drop columns with only zeroes
new_order <- c("GO", sort(colnames(df_logQ[2:32])), "BP", "CC", "MF") # sort columns
df_logQ <- df_logQ[,new_order]

# separate into GO categories
BP <- df_logQ[which(df_logQ$BP!=""),-c(1,34,35)]
rownames(BP) <- BP$BP
BP <- BP[-32]

CC <- df_logQ[which(df_logQ$CC!=""),-c(1,33,35)]
rownames(CC) <- CC$CC
CC <- CC[-32]

MF <- df_logQ[which(df_logQ$MF!=""),-c(1,33,34)]
rownames(MF) <- MF$MF
MF <- MF[-32]

# heatmaps (red=enriched, blue=underrepresented)
palette_func <- colorRampPalette(c("red", "white"))
palette <- rev(palette_func(2)) 
col_fun <- colorRamp2(c(0, max(df_logQ[2:32])), palette)
pdf("Scripts/Data_Vis/ORFs_pres_abs_GO_BP_CC_MF_logQ.pdf", height=unit(10, "in"), width=unit(15, "in"))
p1 <- plot_hm(BP, range=c(0, max(df_logQ[2:32])))
p2 <- plot_hm(CC, range=c(0, max(df_logQ[2:32])))
p3 <- plot_hm(MF, range=c(0, max(df_logQ[2:32])))
draw(p1 %v% p2 %v% p3)
dev.off()

## heatmaps for CNVs
df_logQ <- read.csv("Scripts/Data_Vis/ORF_CNO_Figures/ORFs_copy_num_GO_logQ.tsv", sep="\t") # manually cleaned
df_logQ[is.na(df_logQ)] <- 0 # set insignificant q-values to zero
df_logQ <- df_logQ[which(df_logQ$YPDANISO10 >= 0),] # drop rows with negative values (underreprsented GO terms)
df_logQ <- df_logQ %>% dplyr::select(which(apply(df_logQ, 2, n_distinct) > 1)) # drop columns with only zeroes
new_order <- c("GO", sort(colnames(df_logQ[2:30])), "BP", "CC", "MF") # sort columns
df_logQ <- df_logQ[,new_order]

# separate into GO categories
BP <- df_logQ[which(df_logQ$BP!=""),-c(1,32,33)]
rownames(BP) <- BP$BP
BP <- BP[-30]

CC <- df_logQ[which(df_logQ$CC!=""),-c(1,31,33)]
rownames(CC) <- CC$CC
CC <- CC[-30]

MF <- df_logQ[which(df_logQ$MF!=""),-c(1,31,32)]
rownames(MF) <- MF$MF
MF <- MF[-30]

# heatmaps (red=enriched, blue=underrepresented)
palette_func <- colorRampPalette(c("red", "white"))
palette <- rev(palette_func(2)) 
col_fun <- colorRamp2(c(0, max(df_logQ[2:35])), palette)
pdf("Scripts/Data_Vis/ORFs_copy_num_GO_BP_CC_MF_logQ.pdf", height=unit(10, "in"), width=unit(15, "in"))
p1 <- plot_hm(BP, range=c(0, max(df_logQ[2:30])))
p2 <- plot_hm(CC, range=c(0, max(df_logQ[2:30])))
p3 <- plot_hm(MF, range=c(0, max(df_logQ[2:30])))
draw(p1 %v% p2 %v% p3)
dev.off()

################################################################################
# PATHWAY ENRICHMENT (RF MODELS AFTER FS)
################################################################################
## SNPs
dir <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/SNPs/fs" # path to ORA results files
files <- list.files(path=dir, pattern="PWY_ORA_SHAP_", full.names=T, recursive=F)
path <- "/mnt/home/seguraab/Shiu_Lab/Co-function/Data/MetaCyc/All-pathways-S288c_descriptions.txt"
info <- read.csv(path, sep="\t")# pathway descriptions

combine <- function(files, path){
        # read in the first file
        trait <- str_extract(files[1], "YP[A-Z0-9]+[^_]")
        df <- read.csv(files[1], sep="\t")
        df$Env <- trait # add environment to keep track where the GO term is enriched
        df <- df[df$lfdr < 0.05,] # filter by local FDR

        # combine with remaining environments
        for (i in 2:length(files)){
        trait <- str_extract(files[i], "YP[A-Z0-9]+[^_]") # environment name
        df2 <- read.csv(files[i], sep="\t")
        df2$Env <- trait
        df2 <- df2[df2$lfdr < 0.05,]
        df <- rbind.fill(df, df2)
        }

        # take the log of the q-values
        df$logQ <- 0
        df <- df[order(df$qvalues, decreasing=F),]
        for(i in 1:nrow(df)){ 
                if(df$direction[i] == '-') df$logQ[i] <- log10(df$qvalues[i])
                if(df$direction[i] == '+') df$logQ[i] <- -log10(df$qvalues[i])
        }
        if (nrow(df[df$logQ < -10,])!=0) df[df$loggQ < -10,]$logQ <- -10 # set limits for extreme values
        if (nrow(df[df$logQ > 10,])!=0) df[df$logQ > 10,]$logQ <- 10

        # reshape
        df_logQ <- df %>% group_by(Env) %>% 
                pivot_wider(id_cols=PWY, names_from=Env, values_from=logQ) %>% 
                as.data.frame()
        df_logQ$PWY <- gsub(" ", "", df_logQ$PWY) # remove extra white spaces in df_logQ
        df_logQ <- left_join(df_logQ, info[c('Pathways','Object.ID')], by=c('PWY'='Object.ID'))
        # df_logQ <- df_logQ[!duplicated(df_logQ),]
        write.table(df_logQ, path, sep="\t", quote=F, row.names=F)
        return(df_logQ)
}

df_logQ <- combine(files, "Scripts/Data_Vis/SNPs_PWY_logQ.tsv")
# df_logQ[is.na(df_logQ)] <- 0 # set insignificant q-values to zero
# df_logQ <- df_logQ[which(df_logQ$YPD6AU >= 0),] # drop rows with negative values
# df_logQ <- df_logQ %>% dplyr::select(which(apply(df_logQ, 2, n_distinct) > 1)) # drop columns with only zeroes

# heatmap of significant pathways
palette_func <- colorRampPalette(c("red", "white"))
palette <- rev(palette_func(2)) 
col_fun <- colorRamp2(c(0, max(df_logQ[2])), palette)
rownames(df_logQ) <- df_logQ$Pathways
pdf("Scripts/Data_Vis/SNPs_PWY_logQ.pdf", height=unit(9, "in"), width=unit(10, "in"))
p <- plot_hm(df_logQ[2], range=c(0, max(df_logQ[2]))) # function in snp go enrichment
draw(p)
dev.off()

# get genes annotated with these pathways
path <- "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_pwy.csv"
genes <- read.csv(path)
genes$Pathways.of.gene <- gsub(" ", "", genes$Pathways.of.gene)
df_logQ <- left_join(df_logQ, genes[c(2,6)], by=c("PWY"="Pathways.of.gene"))
df_logQ <- df_logQ[!duplicated(df_logQ),] # drop duplicate rows due to many snps per gene
write.csv(df_logQ, "Scripts/Data_Vis/SNPs_PWY_logQ_genes.csv", quote=F, row.names=F)

## ORF presence/absence & ORF copy number
dir <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/ORFs/fs" # path to ORA results files
files <- list.files(path=dir, pattern="PWY_ORA_SHAP_", full.names=T, recursive=F)
path <- "/mnt/home/seguraab/Shiu_Lab/Co-function/Data/MetaCyc/All-pathways-S288c_descriptions.txt"
info <- read.csv(path, sep="\t")# pathway descriptions

files_orf <- c()
files_cnv <- c()
for (i in 1:length(files)){
        if (grepl("_orf_", files[i])) files_orf <- c(files_orf, files[i])
        if (grepl("_cno_", files[i])) files_cnv <- c(files_cnv, files[i])
}

df_logQ_orf <- combine(files_orf, "Scripts/Data_Vis/ORFs_pres_abs_PWY_logQ.tsv") # function on line 996
df_logQ_cnv <- combine(files_cnv, "Scripts/Data_Vis/ORFs_copy_num_PWY_logQ.tsv")

df_logQ_orf[is.na(df_logQ_orf)] <- 0 # set insignificant q-values to zero
df_logQ_cnv[is.na(df_logQ_cnv)] <- 0
df_logQ_orf <- df_logQ_orf %>% dplyr::select(which(apply(df_logQ_orf, 2, n_distinct) > 1)) # drop columns with only zeroes
df_logQ_cnv <- df_logQ_cnv %>% dplyr::select(which(apply(df_logQ_cnv, 2, n_distinct) > 1))
new_order <- c("PWY", sort(colnames(df_logQ_orf[2:23])), "Pathways") # sort columns
df_logQ_orf <- df_logQ_orf[,new_order]
new_order <- c("PWY", sort(colnames(df_logQ_cnv[2:26])), "Pathways")
df_logQ_cnv <- df_logQ_cnv[,new_order]

# heatmaps of significant pathways
palette_func <- colorRampPalette(c("red", "white"))
palette <- rev(palette_func(2))
col_fun <- colorRamp2(c(0, max(df_logQ_orf[2:23])), palette)
rownames(df_logQ_orf) <- df_logQ_orf$Pathways
pdf("Scripts/Data_Vis/ORFs_pres_abs_PWY_logQ.pdf", height=unit(9, "in"), width=unit(10, "in"))
p <- plot_hm(df_logQ_orf[2:23], range=c(0, max(df_logQ_orf[2:23]))) # function on line 818
draw(p)
dev.off()

col_fun <- colorRamp2(c(0, max(df_logQ_cnv[2:26])), palette)
rownames(df_logQ_cnv) <- df_logQ_cnv$Pathways
pdf("Scripts/Data_Vis/ORFs_copy_num_PWY_logQ.pdf", height=unit(9, "in"), width=unit(10, "in"))
p <- plot_hm(df_logQ_cnv[2:26], range=c(0, max(df_logQ_cnv[2:26]))) # function on line 818
draw(p)
dev.off()

# get genes annotated with these pathways
path <- "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs_and_S288C_genes_pwy.csv"
genes <- read.csv(path)
genes$Pathways.of.gene <- gsub(" ", "", genes$Pathways.of.gene)
df_logQ_orf <- left_join(df_logQ_orf, genes[c(2,6)], by=c("PWY"="Pathways.of.gene"))
df_logQ_orf <- df_logQ_orf[!duplicated(df_logQ_orf),] # drop duplicate rows due to several orf gene matches
write.csv(df_logQ_orf, "Scripts/Data_Vis/ORFs_pres_abs_PWY_logQ_genes.csv", quote=F, row.names=F)
df_logQ_cnv <- left_join(df_logQ_cnv, genes[c(2,6)], by=c("PWY"="Pathways.of.gene"))
df_logQ_cnv <- df_logQ_cnv[!duplicated(df_logQ_cnv),]
write.csv(df_logQ_cnv, "Scripts/Data_Vis/ORFs_copy_num_PWY_logQ_genes.csv", quote=F, row.names=F)

################################################################################
# GENE RF IMPORTANCE SCORES (AFTER FS) COMPARISONS ACROSS ENVIRONMENTS
################################################################################
## SNPs
snp <- read.csv("Results/RESULTS_RF_SNPs_FS.txt", sep="\t", header=T)
map <- read.csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes.txt", header=F) # snp to gene map
colnames(map) <- c("snp", "chr", "pos", "gene")

# Figure: how many snps per gene
counts <- map %>% group_by(gene) %>% tally()
ggplot(counts, aes(n)) + geom_histogram(bins=100) + xlim(c(0,100))
ggsave("Scripts/Data_Vis/SNPs_per_gene_counts.pdf")
dev.off()
table(counts$n) # most genes have about 40 or less snps

# Combine feature importance files to make a heatmap
dir <- "/mnt/scratch/seguraab/yeast_project/SNP_yeast_RF_results/fs/" # data directory
files <- paste(dir, snp$ID, "_imp", sep="") # path to feature selection importance score files
save_names <- paste("Scripts/Data_Vis/", snp$ID, "_imp.tsv", sep="") # new file names
# read in the first file
imp_df <- read.csv(files[1], sep="\t")
imp_df <- imp_df[order(imp_df$mean_imp, decreasing=T),] # sort by mean importance score
# add gene information and calculate max score per gene
imp_df <- left_join(imp_df, map[,c(1,4)], by=c("X"="snp"))
imp_df <- imp_df %>% group_by(gene) %>% dplyr::summarise(count=n(), max=max(mean_imp)) %>% as.data.frame()
write.table(imp_df, save_names[1], sep="\t", quote=F, row.names=F) # save to file before subsetting columns
# bin the gene level max importance scores by percentiles
imp_df_sub <- imp_df_sub[,c('gene', 'max')] # subset
imp_df_sub <- imp_df_sub[which(imp_df_sub$gene!="intergenic"),] # drop intergenic snps
imp_df_sub <- imp_df_sub[!duplicated(imp_df_sub),] # remove duplicate genes
imp_df_sub$gene_bin <- cut(imp_df_sub$max,
        breaks=quantile(imp_df_sub$max, c(0, 0.95, 0.99, 1)),
        labels=c("(0,95]", "(95,99]", "(99,100]"), include.lowest=T)
imp_df_sub <- imp_df_sub[order(imp_df_sub$max, decreasing=T),] # re-order by gene level mean imp
# subset imp_df_sub for combining with the rest of the environments
imp_df_sub2 <- imp_df_sub[,c("gene", "gene_bin")] # keep only gene and bin columns
colnames(imp_df_sub2) <- c("gene", str_extract(files[1],"(?<=/)[A-Z0-9]+(?=_rf|_exp)")) # add env name
imp_df_sub3 <- imp_df_sub[,c("gene", "max")] # keep only gene and max mean_imp columns
colnames(imp_df_sub3) <- c("gene", str_extract(files[1],"(?<=/)[A-Z0-9]+(?=_rf|_exp)"))
imp_df_sub2_top10 <- imp_df_sub2[1:10,] # take the top 10 genes
imp_df_sub3_top10 <- imp_df_sub3[1:10,]
# combine dataframes
for (i in 2:length(files)){
        df <- read.csv(files[i], sep="\t")
        df <- df[order(df$mean_imp, decreasing=T),] # sort by mean importance score
        # add gene information and calculate mean and standard error of mean scores per gene
        df <- left_join(df, map[,c(1,4)], by=c("X"="snp"))
        df <- df %>% group_by(gene) %>% dplyr::summarise(count=n(), max=max(mean_imp)) %>% as.data.frame()
        write.table(df, save_names[i], sep="\t", quote=F, row.names=F) # save to file before subsetting columns
        # bin the gene level max importance scores by percentiles
        df_sub <- df[,c('gene', 'max')] # subset
        df_sub <- df_sub[which(df_sub$gene!="intergenic"),] # drop intergenic snps
        df_sub <- df_sub[!duplicated(df_sub),] # remove duplicate genes
        df_sub$gene_bin <- cut(df_sub$max,
                breaks=quantile(df_sub$max, c(0, 0.95, 0.99, 1)),
                labels=c("(0,95]", "(95,99]", "(99,100]"), include.lowest=T)
        df_sub <- df_sub[order(df_sub$max, decreasing=T),] # re-order by gene level mean imp
        # subset imp_df_sub for combining with the rest of the environments
        df_sub2 <- df_sub[,c("gene", "gene_bin")] # keep only gene and bin columns
        colnames(df_sub2) <- c("gene", str_extract(files[i],"(?<=/)[A-Z0-9]+(?=_rf|_exp)"))
        df_sub3 <- df_sub[,c("gene", "max")] # keep only gene and max mean_imp columns
        colnames(df_sub3) <- c("gene", str_extract(files[i],"(?<=/)[A-Z0-9]+(?=_rf|_exp)"))
        # combine with imp_df_sub
        imp_df_sub2 <- full_join(imp_df_sub2, df_sub2, by="gene") # percentiles
        imp_df_sub3 <- full_join(imp_df_sub3, df_sub3, by="gene") # importance scores
        # take the top 10 genes
        df_sub2_top10 <- df_sub2[1:10,]
        df_sub3_top10 <- df_sub3[1:10,]
        imp_df_sub2_top10 <- full_join(imp_df_sub2_top10, df_sub2_top10, by="gene") # combine with imp_df_sub_top10
        imp_df_sub3_top10 <- full_join(imp_df_sub3_top10, df_sub3_top10, by="gene")
        
}
# save combined data
write.table(imp_df_sub2, "Scripts/Data_Vis/SNPs_imp_percentile_all_sorted.tsv", sep="\t", quote=F, row.names=F)
write.table(imp_df_sub3, "Scripts/Data_Vis/SNPs_imp_max_all_sorted.tsv", sep="\t", quote=F, row.names=F)
write.table(imp_df_sub2_top10, "Scripts/Data_Vis/SNPs_imp_percentile_top10_sorted.tsv", sep="\t", quote=F, row.names=F)
write.table(imp_df_sub3_top10, "Scripts/Data_Vis/SNPs_imp_max_top10_sorted.tsv", sep="\t", quote=F, row.names=F)

# heatmap of importance scores for all the top genes in each environment
fac_cols <- sapply(imp_df_sub2, is.factor) # identify all factor columns
imp_df_sub2[fac_cols] <- lapply(imp_df_sub2[fac_cols], as.character) # convert all factors to characters
imp_df_sub2[is.na(imp_df_sub2)] <- "0" # replace string values to plot heatmap
imp_df_sub2[imp_df_sub2=="(0,95]"] <- "1"
imp_df_sub2[imp_df_sub2=="(95,99]"] <- "2"
imp_df_sub2[imp_df_sub2=="(99,100]"] <- "3"
rownames(imp_df_sub2) <- imp_df_sub2$gene # set rownames
imp_df_sub2 <- imp_df_sub2[,-1] # remove gene column
imp_df_sub2_num <- matrix(as.numeric(as.matrix(imp_df_sub2)), ncol = ncol(imp_df_sub2)) # convert to numeric matrix
rownames(imp_df_sub2_num) <- rownames(imp_df_sub2)
colnames(imp_df_sub2_num) <- colnames(imp_df_sub2)
col_fun = colorRamp2(c(0, 1, 2, 3), c("#fff8fe", "#b898cf", "#ce0b87", "#6b2700"))
# clustered heatmap
p <- Heatmap(imp_df_sub2_num,
        col=col_fun,
        # na_col="grey",
        cluster_rows=T, 
        show_row_dend=F,
        cluster_columns=T,
        show_column_dend=F,
        show_row_names=F,
        show_column_names=T,
        column_names_side="top",
        column_names_gp = gpar(fontize=7),
        border_gp = gpar(col="black", lty=1),
        use_raster=T,
        name="Importance percentiles",
        heatmap_legend_param=list(title="Importance percentiles",
                color_bar = "discrete", 
                labels = gt_render(c("missing", "(0,95]", "(95,99]", "(99,100]"))))
pdf("Scripts/Data_Vis/SNPs_imp_percentile_all_clustered.pdf", height=unit(8.75, "in"), width=unit(7.5, "in"))
draw(p)
dev.off()

# sorted heatmap
p <- Heatmap(imp_df_sub2_num,
        col=col_fun,
        # na_col="grey",
        cluster_rows=F, 
        show_row_dend=F,
        cluster_columns=F,
        show_column_dend=F,
        show_row_names=F,
        show_column_names=T,
        column_names_side="top",
        column_names_gp = gpar(fontize=7),
        border_gp = gpar(col="black", lty=1),
        use_raster=T,
        name="Importance percentiles",
        heatmap_legend_param=list(title="Importance percentiles",
                color_bar = "discrete", 
                labels = gt_render(c("missing", "(0,95]", "(95,99]", "(99,100]"))))
pdf("Scripts/Data_Vis/SNPs_imp_percentile_all_sorted.pdf", height=unit(8.75, "in"), width=unit(7.5, "in"))
draw(p)
dev.off()

# heatmap of importance scores for top 10 genes in each environment
fac_cols <- sapply(imp_df_sub2_top10, is.factor) # identify all factor columns
imp_df_sub2_top10[fac_cols] <- lapply(imp_df_sub2_top10[fac_cols], as.character) # convert all factors to characters
imp_df_sub2_top10[is.na(imp_df_sub2_top10)] <- "0" # replace string values to plot heatmap
imp_df_sub2_top10[imp_df_sub2_top10=="(0,95]"] <- "1"
imp_df_sub2_top10[imp_df_sub2_top10=="(95,99]"] <- "2"
imp_df_sub2_top10[imp_df_sub2_top10=="(99,100]"] <- "3"
rownames(imp_df_sub2_top10) <- imp_df_sub2_top10$gene # set rownames
imp_df_sub2_top10 <- imp_df_sub2_top10[,-1] # remove gene column
imp_df_sub2_top10_num <- matrix(as.numeric(as.matrix(imp_df_sub2_top10)), ncol = ncol(imp_df_sub2_top10)) # convert to numeric matrix
rownames(imp_df_sub2_top10_num) <- rownames(imp_df_sub2_top10)
colnames(imp_df_sub2_top10_num) <- colnames(imp_df_sub2_top10)
col_fun = colorRamp2(c(0, 1, 2, 3), c("#fff8fe", "#b898cf", "#ce0b87", "#6b2700"))
# clustered heatmap
p <- Heatmap(imp_df_sub2_top10_num,
        col=col_fun,
        # na_col="grey",
        cluster_rows=T, 
        show_row_dend=F,
        cluster_columns=T,
        show_column_dend=F,
        show_row_names=F,
        show_column_names=T,
        column_names_side="top",
        column_names_gp = gpar(fontize=7),
        border_gp = gpar(col="black", lty=1),
        use_raster=T,
        name="Importance percentiles",
        heatmap_legend_param=list(title="Importance percentiles",
                color_bar = "discrete", 
                labels = gt_render(c("missing", "(0,95]", "(95,99]", "(99,100]"))))
pdf("Scripts/Data_Vis/SNPs_imp_percentile_top10_clustered.pdf", height=unit(8.75, "in"), width=unit(7.5, "in"))
draw(p)
dev.off()

# sorted heatmap
p <- Heatmap(imp_df_sub2_top10_num,
        col=col_fun,
        # na_col="grey",
        cluster_rows=F, 
        show_row_dend=F,
        cluster_columns=F,
        show_column_dend=F,
        show_row_names=F,
        show_column_names=T,
        column_names_side="top",
        column_names_gp = gpar(fontize=7),
        border_gp = gpar(col="black", lty=1),
        use_raster=T,
        name="Importance percentiles",
        heatmap_legend_param=list(title="Importance percentiles",
                color_bar = "discrete", 
                labels = gt_render(c("missing", "(0,95]", "(95,99]", "(99,100]"))))
pdf("Scripts/Data_Vis/SNPs_imp_percentile_top10_sorted.pdf", height=unit(8.75, "in"), width=unit(7.5, "in"))
draw(p)
dev.off()

# Histogram of the number of environments each gene is a top predictor for
imp_df_sub2_num[imp_df_sub2_num > 1] <- 1 # set all non-zero values to 1
num_envs <- rowSums(imp_df_sub2_num)# get the number of envs for each gene
write.table(num_envs, "Scripts/Data_Vis/SNPs_top_genes_num_env.tsv", sep="\t", quote=F)
num_envs <- num_envs[which(names(num_envs)!="intergenic")] # remove intergenic category
ggplot(as.data.frame(num_envs), aes(x=num_envs)) + geom_histogram(bins=35) + 
        xlab("Number of Environments") + ylab("Number of Genes") + 
        scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
        theme_bw()
ggsave("Scripts/Data_Vis/SNPs_top_genes_num_env.pdf", height=unit(3.5, "in"), width=unit(3.5, "in"))
dev.off()

############################# ORF presence/absence #############################
orf <- read.csv("Results/RESULTS_RF_ORFs_FS.txt", sep="\t", header=T)
map <- read.csv("Data/Peter_2018/ORFs_GO_genes.tsv", sep="\t")
map <- map[c("qacc", "gene")]
map <- map[!duplicated(map),]
colnames(map) <- c("orf", "gene")

# Combine feature importance files to make a heatmap
dir <- "/mnt/scratch/seguraab/yeast_project/ORF_yeast_RF_results/fs/" # data directory
files <- paste(dir, orf$ID, "_imp", sep="") # path to feature selection importance score files
save_names <- paste("Scripts/Data_Vis/", orf$ID, "_imp.tsv", sep="") # new file names

# read in the first file
imp_df <- read.csv(files[1], sep="\t")
imp_df <- imp_df[order(imp_df$mean_imp, decreasing=T),] # sort by mean importance score
# add gene information and calculate max score per gene
imp_df <- left_join(imp_df, map, by=c("X"="orf"))
imp_df <- imp_df %>% group_by(gene) %>% dplyr::summarise(max=max(mean_imp)) %>% as.data.frame()
write.table(imp_df, save_names[1], sep="\t", quote=F, row.names=F) # save to file before subsetting columns
# bin the gene level max importance scores by percentiles
imp_df <- imp_df[!duplicated(imp_df),] # remove duplicate genes
imp_df$gene_bin <- cut(imp_df$max,
        breaks=quantile(imp_df$max, c(0, 0.95, 0.99, 1)),
        labels=c("(0,95]", "(95,99]", "(99,100]"), include.lowest=T)
imp_df <- imp_df[order(imp_df$max, decreasing=T),] # re-order by gene level mean imp
# subset imp_df for combining with the rest of the environments
imp_df2 <- imp_df[,c("gene", "gene_bin")] # keep only gene and bin columns
colnames(imp_df2) <- c("gene", str_extract(files[1],"(?<=/)[A-Z0-9]+(?=_orf|_exp)")) # add env name
imp_df3 <- imp_df[,c("gene", "max")] # keep only gene and max mean_imp columns
colnames(imp_df3) <- c("gene", str_extract(files[1],"(?<=/)[A-Z0-9]+(?=_orf|_exp)"))
imp_df2_top10 <- imp_df2[1:10,] # take the top 10 genes
imp_df3_top10 <- imp_df3[1:10,]
# combine dataframes
for (i in 2:length(files)){
        df <- read.csv(files[i], sep="\t")
        df <- df[order(df$mean_imp, decreasing=T),] # sort by mean importance score
        # add gene information and calculate mean and standard error of mean scores per gene
        df <- left_join(df, map, by=c("X"="orf"))
        df <- df %>% group_by(gene) %>% dplyr::summarise(max=max(mean_imp)) %>% as.data.frame()
        # bin the gene level max importance scores by percentiles
        df <- df[!duplicated(df),] # remove duplicate genes
        df$gene_bin <- cut(df$max,
                breaks=quantile(df$max, c(0, 0.95, 0.99, 1)),
                labels=c("(0,95]", "(95,99]", "(99,100]"), include.lowest=T)
        df <- df[order(df$max, decreasing=T),] # re-order by gene level mean imp
        # subset imp_df for combining with the rest of the environments
        df2 <- df[,c("gene", "gene_bin")] # keep only gene and bin columns
        colnames(df2) <- c("gene", str_extract(files[i],"(?<=/)[A-Z0-9]+(?=_orf|_exp)"))
        df3 <- df[,c("gene", "max")] # keep only gene and max mean_imp columns
        colnames(df3) <- c("gene", str_extract(files[i],"(?<=/)[A-Z0-9]+(?=_orf|_exp)"))
        # combine with imp_df
        imp_df2 <- full_join(imp_df2, df2, by="gene") # percentiles
        imp_df3 <- full_join(imp_df3, df3, by="gene") # importance scores
        # take the top 10 genes
        df2_top10 <- df2[1:10,]
        df3_top10 <- df3[1:10,]
        imp_df2_top10 <- full_join(imp_df2_top10, df2_top10, by="gene") # combine with imp_df_top10
        imp_df3_top10 <- full_join(imp_df3_top10, df3_top10, by="gene")
        
}
# save combined data
write.table(imp_df2, "Scripts/Data_Vis/ORFs_imp_percentile_all_sorted.tsv", sep="\t", quote=F, row.names=F)
write.table(imp_df3, "Scripts/Data_Vis/ORFs_imp_max_all_sorted.tsv", sep="\t", quote=F, row.names=F)
write.table(imp_df2_top10, "Scripts/Data_Vis/ORFs_imp_percentile_top10_sorted.tsv", sep="\t", quote=F, row.names=F)
write.table(imp_df3_top10, "Scripts/Data_Vis/ORFs_imp_max_top10_sorted.tsv", sep="\t", quote=F, row.names=F)

# heatmap of importance scores for all the top genes in each environment
fac_cols <- sapply(imp_df2, is.factor) # identify all factor columns
imp_df2[fac_cols] <- lapply(imp_df2[fac_cols], as.character) # convert all factors to characters
imp_df2[is.na(imp_df2)] <- "0" # replace string values to plot heatmap
imp_df2[imp_df2=="(0,95]"] <- "1"
imp_df2[imp_df2=="(95,99]"] <- "2"
imp_df2[imp_df2=="(99,100]"] <- "3"
rownames(imp_df2) <- imp_df2$gene # set rownames
imp_df2 <- imp_df2[,-1] # remove gene column
imp_df2_num <- matrix(as.numeric(as.matrix(imp_df2)), ncol = ncol(imp_df2)) # convert to numeric matrix
rownames(imp_df2_num) <- rownames(imp_df2)
colnames(imp_df2_num) <- colnames(imp_df2)
col_fun = colorRamp2(c(0, 1, 2, 3), c("#fff8fe", "#b898cf", "#ce0b87", "#6b2700"))
# clustered heatmap
p <- Heatmap(imp_df2_num,
        col=col_fun,
        # na_col="grey",
        cluster_rows=T, 
        show_row_dend=F,
        cluster_columns=T,
        show_column_dend=F,
        show_row_names=F,
        show_column_names=T,
        column_names_side="top",
        column_names_gp = gpar(fontize=7),
        border_gp = gpar(col="black", lty=1),
        use_raster=T,
        name="Importance percentiles",
        heatmap_legend_param=list(title="Importance percentiles",
                color_bar = "discrete", 
                labels = gt_render(c("missing", "(0,95]", "(95,99]", "(99,100]"))))
pdf("Scripts/Data_Vis/ORFs_imp_percentile_all_clustered.pdf", height=unit(8.75, "in"), width=unit(7.5, "in"))
draw(p)
dev.off()

# sorted heatmap
p <- Heatmap(imp_df2_num,
        col=col_fun,
        # na_col="grey",
        cluster_rows=F, 
        show_row_dend=F,
        cluster_columns=F,
        show_column_dend=F,
        show_row_names=F,
        show_column_names=T,
        column_names_side="top",
        column_names_gp = gpar(fontize=7),
        border_gp = gpar(col="black", lty=1),
        use_raster=T,
        name="Importance percentiles",
        heatmap_legend_param=list(title="Importance percentiles",
                color_bar = "discrete", 
                labels = gt_render(c("missing", "(0,95]", "(95,99]", "(99,100]"))))
pdf("Scripts/Data_Vis/ORFs_imp_percentile_all_sorted.pdf", height=unit(8.75, "in"), width=unit(7.5, "in"))
draw(p)
dev.off()

# heatmap of importance scores for top 10 genes in each environment
fac_cols <- sapply(imp_df2_top10, is.factor) # identify all factor columns
imp_df2_top10[fac_cols] <- lapply(imp_df2_top10[fac_cols], as.character) # convert all factors to characters
imp_df2_top10[is.na(imp_df2_top10)] <- "0" # replace string values to plot heatmap
imp_df2_top10[imp_df2_top10=="(0,95]"] <- "1"
imp_df2_top10[imp_df2_top10=="(95,99]"] <- "2"
imp_df2_top10[imp_df2_top10=="(99,100]"] <- "3"
rownames(imp_df2_top10) <- imp_df2_top10$gene # set rownames
imp_df2_top10 <- imp_df2_top10[,-1] # remove gene column
imp_df2_top10_num <- matrix(as.numeric(as.matrix(imp_df2_top10)), ncol = ncol(imp_df2_top10)) # convert to numeric matrix
rownames(imp_df2_top10_num) <- rownames(imp_df2_top10)
colnames(imp_df2_top10_num) <- colnames(imp_df2_top10)
col_fun = colorRamp2(c(0, 1, 2, 3), c("#fff8fe", "#b898cf", "#ce0b87", "#6b2700"))
# clustered heatmap
p <- Heatmap(imp_df2_top10_num,
        col=col_fun,
        # na_col="grey",
        cluster_rows=T, 
        show_row_dend=F,
        cluster_columns=T,
        show_column_dend=F,
        show_row_names=F,
        show_column_names=T,
        column_names_side="top",
        column_names_gp = gpar(fontize=7),
        border_gp = gpar(col="black", lty=1),
        use_raster=T,
        name="Importance percentiles",
        heatmap_legend_param=list(title="Importance percentiles",
                color_bar = "discrete", 
                labels = gt_render(c("missing", "(0,95]", "(95,99]", "(99,100]"))))
pdf("Scripts/Data_Vis/ORFs_imp_percentile_top10_clustered.pdf", height=unit(8.75, "in"), width=unit(7.5, "in"))
draw(p)
dev.off()

# sorted heatmap
p <- Heatmap(imp_df2_top10_num,
        col=col_fun,
        # na_col="grey",
        cluster_rows=F, 
        show_row_dend=F,
        cluster_columns=F,
        show_column_dend=F,
        show_row_names=F,
        show_column_names=T,
        column_names_side="top",
        column_names_gp = gpar(fontize=7),
        border_gp = gpar(col="black", lty=1),
        use_raster=T,
        name="Importance percentiles",
        heatmap_legend_param=list(title="Importance percentiles",
                color_bar = "discrete", 
                labels = gt_render(c("missing", "(0,95]", "(95,99]", "(99,100]"))))
pdf("Scripts/Data_Vis/ORFs_imp_percentile_top10_sorted.pdf", height=unit(8.75, "in"), width=unit(7.5, "in"))
draw(p)
dev.off()

# Histogram of the number of environments each gene is a top predictor for
imp_df2_num[imp_df2_num > 1] <- 1 # set all non-zero values to 1
num_envs <- rowSums(imp_df2_num)# get the number of envs for each gene
write.table(num_envs, "Scripts/Data_Vis/ORFs_top_genes_num_env.tsv", sep="\t", quote=F)
ggplot(as.data.frame(num_envs), aes(x=num_envs)) + geom_histogram(bins=35) + 
        xlab("Number of Environments") + ylab("Number of Genes") + 
        scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
        theme_bw()
ggsave("Scripts/Data_Vis/ORFs_top_genes_num_env.pdf", height=unit(3.5, "in"), width=unit(3.5, "in"))
dev.off()

############################### ORF copy number ################################
orf <- read.csv("Results/RESULTS_RF_CNVs_FS.txt", sep="\t", header=T)
map <- read.csv("Data/Peter_2018/ORFs_GO_genes.tsv", sep="\t")
map <- map[c("qacc", "gene")]
map <- map[!duplicated(map),]
colnames(map) <- c("orf", "gene")

# Combine feature importance files to make a heatmap
dir <- "/mnt/scratch/seguraab/yeast_project/ORF_yeast_RF_results/fs/" # data directory
files <- paste(dir, orf$ID, "_imp", sep="") # path to feature selection importance score files
save_names <- paste("Scripts/Data_Vis/", orf$ID, "_imp.tsv", sep="") # new file names

# read in the first file
imp_df <- read.csv(files[1], sep="\t")
imp_df <- imp_df[order(imp_df$mean_imp, decreasing=T),] # sort by mean importance score
# add gene information and calculate max score per gene
imp_df <- left_join(imp_df, map, by=c("X"="orf"))
imp_df <- imp_df %>% group_by(gene) %>% dplyr::summarise(max=max(mean_imp)) %>% as.data.frame()
write.table(imp_df, save_names[1], sep="\t", quote=F, row.names=F) # save to file before subsetting columns
# bin the gene level max importance scores by percentiles
imp_df <- imp_df[!duplicated(imp_df),] # remove duplicate genes
imp_df$gene_bin <- cut(imp_df$max,
        breaks=quantile(imp_df$max, c(0, 0.95, 0.99, 1)),
        labels=c("(0,95]", "(95,99]", "(99,100]"), include.lowest=T)
imp_df <- imp_df[order(imp_df$max, decreasing=T),] # re-order by gene level mean imp
# subset imp_df for combining with the rest of the environments
imp_df2 <- imp_df[,c("gene", "gene_bin")] # keep only gene and bin columns
colnames(imp_df2) <- c("gene", str_extract(files[1],"(?<=/)[A-Z0-9]+(?=_cno|_exp)")) # add env name
imp_df3 <- imp_df[,c("gene", "max")] # keep only gene and max mean_imp columns
colnames(imp_df3) <- c("gene", str_extract(files[1],"(?<=/)[A-Z0-9]+(?=_cno|_exp)"))
imp_df2_top10 <- imp_df2[1:10,] # take the top 10 genes
imp_df3_top10 <- imp_df3[1:10,]
# combine dataframes
for (i in 2:length(files)){
        df <- read.csv(files[i], sep="\t")
        df <- df[order(df$mean_imp, decreasing=T),] # sort by mean importance score
        # add gene information and calculate mean and standard error of mean scores per gene
        df <- left_join(df, map, by=c("X"="orf"))
        df <- df %>% group_by(gene) %>% dplyr::summarise(max=max(mean_imp)) %>% as.data.frame()
        # bin the gene level max importance scores by percentiles
        df <- df[!duplicated(df),] # remove duplicate genes
        df$gene_bin <- cut(df$max,
                breaks=quantile(df$max, c(0, 0.95, 0.99, 1)),
                labels=c("(0,95]", "(95,99]", "(99,100]"), include.lowest=T)
        df <- df[order(df$max, decreasing=T),] # re-order by gene level mean imp
        # subset imp_df for combining with the rest of the environments
        df2 <- df[,c("gene", "gene_bin")] # keep only gene and bin columns
        colnames(df2) <- c("gene", str_extract(files[i],"(?<=/)[A-Z0-9]+(?=_cno|_exp)"))
        df3 <- df[,c("gene", "max")] # keep only gene and max mean_imp columns
        colnames(df3) <- c("gene", str_extract(files[i],"(?<=/)[A-Z0-9]+(?=_cno|_exp)"))
        # combine with imp_df
        imp_df2 <- full_join(imp_df2, df2, by="gene") # percentiles
        imp_df3 <- full_join(imp_df3, df3, by="gene") # importance scores
        # take the top 10 genes
        df2_top10 <- df2[1:10,]
        df3_top10 <- df3[1:10,]
        imp_df2_top10 <- full_join(imp_df2_top10, df2_top10, by="gene") # combine with imp_df_top10
        imp_df3_top10 <- full_join(imp_df3_top10, df3_top10, by="gene")
        
}
# save combined data
write.table(imp_df2, "Scripts/Data_Vis/CNVs_imp_percentile_all_sorted.tsv", sep="\t", quote=F, row.names=F)
write.table(imp_df3, "Scripts/Data_Vis/CNVs_imp_max_all_sorted.tsv", sep="\t", quote=F, row.names=F)
write.table(imp_df2_top10, "Scripts/Data_Vis/CNVs_imp_percentile_top10_sorted.tsv", sep="\t", quote=F, row.names=F)
write.table(imp_df3_top10, "Scripts/Data_Vis/CNVs_imp_max_top10_sorted.tsv", sep="\t", quote=F, row.names=F)

# heatmap of importance scores for all the top genes in each environment
fac_cols <- sapply(imp_df2, is.factor) # identify all factor columns
imp_df2[fac_cols] <- lapply(imp_df2[fac_cols], as.character) # convert all factors to characters
imp_df2[is.na(imp_df2)] <- "0" # replace string values to plot heatmap
imp_df2[imp_df2=="(0,95]"] <- "1"
imp_df2[imp_df2=="(95,99]"] <- "2"
imp_df2[imp_df2=="(99,100]"] <- "3"
rownames(imp_df2) <- imp_df2$gene # set rownames
imp_df2 <- imp_df2[,-1] # remove gene column
imp_df2_num <- matrix(as.numeric(as.matrix(imp_df2)), ncol = ncol(imp_df2)) # convert to numeric matrix
rownames(imp_df2_num) <- rownames(imp_df2)
colnames(imp_df2_num) <- colnames(imp_df2)
col_fun = colorRamp2(c(0, 1, 2, 3), c("#fff8fe", "#b898cf", "#ce0b87", "#6b2700"))
# clustered heatmap
p <- Heatmap(imp_df2_num,
        col=col_fun,
        # na_col="grey",
        cluster_rows=T, 
        show_row_dend=F,
        cluster_columns=T,
        show_column_dend=F,
        show_row_names=F,
        show_column_names=T,
        column_names_side="top",
        column_names_gp = gpar(fontize=7),
        border_gp = gpar(col="black", lty=1),
        use_raster=T,
        name="Importance percentiles",
        heatmap_legend_param=list(title="Importance percentiles",
                color_bar = "discrete", 
                labels = gt_render(c("missing", "(0,95]", "(95,99]", "(99,100]"))))
pdf("Scripts/Data_Vis/CNVs_imp_percentile_all_clustered.pdf", height=unit(8.75, "in"), width=unit(7.5, "in"))
draw(p)
dev.off()

# sorted heatmap
p <- Heatmap(imp_df2_num,
        col=col_fun,
        # na_col="grey",
        cluster_rows=F, 
        show_row_dend=F,
        cluster_columns=F,
        show_column_dend=F,
        show_row_names=F,
        show_column_names=T,
        column_names_side="top",
        column_names_gp = gpar(fontize=7),
        border_gp = gpar(col="black", lty=1),
        use_raster=T,
        name="Importance percentiles",
        heatmap_legend_param=list(title="Importance percentiles",
                color_bar = "discrete", 
                labels = gt_render(c("missing", "(0,95]", "(95,99]", "(99,100]"))))
pdf("Scripts/Data_Vis/CNVs_imp_percentile_all_sorted.pdf", height=unit(8.75, "in"), width=unit(7.5, "in"))
draw(p)
dev.off()

# heatmap of importance scores for top 10 genes in each environment
fac_cols <- sapply(imp_df2_top10, is.factor) # identify all factor columns
imp_df2_top10[fac_cols] <- lapply(imp_df2_top10[fac_cols], as.character) # convert all factors to characters
imp_df2_top10[is.na(imp_df2_top10)] <- "0" # replace string values to plot heatmap
imp_df2_top10[imp_df2_top10=="(0,95]"] <- "1"
imp_df2_top10[imp_df2_top10=="(95,99]"] <- "2"
imp_df2_top10[imp_df2_top10=="(99,100]"] <- "3"
rownames(imp_df2_top10) <- imp_df2_top10$gene # set rownames
imp_df2_top10 <- imp_df2_top10[,-1] # remove gene column
imp_df2_top10_num <- matrix(as.numeric(as.matrix(imp_df2_top10)), ncol = ncol(imp_df2_top10)) # convert to numeric matrix
rownames(imp_df2_top10_num) <- rownames(imp_df2_top10)
colnames(imp_df2_top10_num) <- colnames(imp_df2_top10)
col_fun = colorRamp2(c(0, 1, 2, 3), c("#fff8fe", "#b898cf", "#ce0b87", "#6b2700"))
# clustered heatmap
p <- Heatmap(imp_df2_top10_num,
        col=col_fun,
        # na_col="grey",
        cluster_rows=T, 
        show_row_dend=F,
        cluster_columns=T,
        show_column_dend=F,
        show_row_names=F,
        show_column_names=T,
        column_names_side="top",
        column_names_gp = gpar(fontize=7),
        border_gp = gpar(col="black", lty=1),
        use_raster=T,
        name="Importance percentiles",
        heatmap_legend_param=list(title="Importance percentiles",
                color_bar = "discrete", 
                labels = gt_render(c("missing", "(0,95]", "(95,99]", "(99,100]"))))
pdf("Scripts/Data_Vis/CNVs_imp_percentile_top10_clustered.pdf", height=unit(8.75, "in"), width=unit(7.5, "in"))
draw(p)
dev.off()

# sorted heatmap
p <- Heatmap(imp_df2_top10_num,
        col=col_fun,
        # na_col="grey",
        cluster_rows=F, 
        show_row_dend=F,
        cluster_columns=F,
        show_column_dend=F,
        show_row_names=F,
        show_column_names=T,
        column_names_side="top",
        column_names_gp = gpar(fontize=7),
        border_gp = gpar(col="black", lty=1),
        use_raster=T,
        name="Importance percentiles",
        heatmap_legend_param=list(title="Importance percentiles",
                color_bar = "discrete", 
                labels = gt_render(c("missing", "(0,95]", "(95,99]", "(99,100]"))))
pdf("Scripts/Data_Vis/CNVs_imp_percentile_top10_sorted.pdf", height=unit(8.75, "in"), width=unit(7.5, "in"))
draw(p)
dev.off()

# Histogram of the number of environments each gene is a top predictor for
imp_df2_num[imp_df2_num > 1] <- 1 # set all non-zero values to 1
num_envs <- rowSums(imp_df2_num)# get the number of envs for each gene
write.table(num_envs, "Scripts/Data_Vis/CNVs_top_genes_num_env.tsv", sep="\t", quote=F)
ggplot(as.data.frame(num_envs), aes(x=num_envs)) + geom_histogram(bins=35) + 
        xlab("Number of Environments") + ylab("Number of Genes") + 
        scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
        theme_bw()
ggsave("Scripts/Data_Vis/CNVs_top_genes_num_env.pdf", height=unit(3.5, "in"), width=unit(3.5, "in"))
dev.off()

################################################################################
# FITNESS FACTORS IMPACTING RF PERFORMANCE
################################################################################
# Scripts: Scripts/Genomic_Prediction_RF/performance_factors_regression_plots.py
#          Scripts/Genomic_Prediction_RF/performance_factors_regression.py
## SNPs
# from linear model regression, mean and Q3 had the highest performances
# see ../Genomic_Prediction_RF/snps_RF_FS_ols_env_results.txt
pheno <- read.csv("Data/Peter_2018/pheno.csv", row.names=1) # fitness data
rf <- read.table("Results/RESULTS_RF_SNPs_FS.txt", sep="\t", header=T)
performance <- rf[c(1,30)] # select only cond and r2_test
rownames(performance) <- performance$cond # set cond as row names
performance <- performance[order(performance$r2_test), ] # sort by increasing performance
pheno <- pheno[performance$cond] # re-order columns in pheno
pheno <- reshape2::melt(pheno) # pivot longer
pheno$variable <- as.factor(pheno$variable) # create factor for plotting
# make violin plot
mean <- mean(pheno[pheno$variable=="YPDCAFEIN40",]$value) # mean fitness of YPDCAFEIN40
q3 <- quantile(pheno[pheno$variable=="YPDCAFEIN40",]$value, 0.75) # 75th percentile of YPDCAFEIN40
ggplot(pheno, aes(x=variable, y=value)) + 
        geom_hline(yintercept=mean, linetype="dashed", color="red") + 
        geom_hline(yintercept=q3, linetype="dashed", color="blue") + 
        geom_violin() + geom_boxplot(width=0.1) + theme_bw(8) +
        theme(axis.text.x=element_text(color="black", size=9, angle=45, hjust=1)) +
        theme(axis.text.y=element_text(color="black", size=9)) + ylab("Fitness") +
        xlab("Environment")
ggsave("Scripts/Data_Vis/pheno_violin.pdf", width=10, height=4, device="pdf", useDingbats=FALSE)
dev.off()


## ORFs presence/absence



## ORFs copy number


################################################################################
# BASELINE RF PERFORMANCE COMPARISONS USING DIFFERENT TEST SETS
################################################################################
## PCs
# evaluated on original test set
pcs_base <- read.csv("Results/RESULTS_RF_PCs_sorted.txt", sep="\t")
pcs_base <- pcs_base[c("ID", "cond", "r2_test")]
colnames(pcs_base) <- c("ID", "Y", "r2_test")
pcs_base$group <- "Original Test"
# evaluated on the five additional test sets and/or randomized label
pcs <- read.csv("/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results/RESULTS_reg_PCs_splits.txt", sep="\t")
pcs <- pcs[c("ID", "Y", "r2_test")]
pcs$group <- ""
pcs[which(endsWith(pcs$ID, "base_rand")),]$group <- "Original Test (random)"
pcs[which(endsWith(pcs$ID, "base_rand_split1")),]$group <- "Test 1 (random)"
pcs[which(endsWith(pcs$ID, "base_rand_split2")),]$group <- "Test 2 (random)"
pcs[which(endsWith(pcs$ID, "base_rand_split3")),]$group <- "Test 3 (random)"
pcs[which(endsWith(pcs$ID, "base_rand_split4")),]$group <- "Test 4 (random)"
pcs[which(endsWith(pcs$ID, "base_rand_split5")),]$group <- "Test 5 (random)"
pcs[which(endsWith(pcs$ID, "baseline_split1")),]$group <- "Test 1"
pcs[which(endsWith(pcs$ID, "baseline_split2")),]$group <- "Test 2"
pcs[which(endsWith(pcs$ID, "baseline_split3")),]$group <- "Test 3"
pcs[which(endsWith(pcs$ID, "baseline_split4")),]$group <- "Test 4"
pcs[which(endsWith(pcs$ID, "baseline_split5")),]$group <- "Test 5"
pcs <- rbind.fill(pcs, pcs_base)
pcs <- pcs[order(pcs$group),]
ggplot(pcs, aes(x=group, y=Y, fill=r2_test)) + geom_tile() +
        geom_text(aes(label=round(r2_test, 2)), size=2.5, color="white") +
        theme_minimal(base_size=8) + ylab("Conditions") + xlab("Test Set")
ggsave("Scripts/Data_Vis/Test_comparisons_PCs.pdf", width=10, height=9)

## SNPs
# evaluated on original test set
snps_base <- read.csv("Results/RESULTS_RF_SNPs_baseline.txt", sep="\t")
snps_base <- snps_base[c("ID", "Y", "r2_test")]
snps_base$group <- "Original Test"
# evaluated on the five additional test sets and/or randomized label
snps <- read.csv("/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results/RESULTS_reg_SNPs_splits.txt", sep="\t")
snps <- snps[c("ID", "Y", "r2_test")]
snps$group <- ""
snps[which(endsWith(snps$ID, "baseline")),]$group <- "Original Test"
snps[which(endsWith(snps$ID, "base_rand")),]$group <- "Original Test (random)"
snps[which(endsWith(snps$ID, "base_rand_split1")),]$group <- "Test 1 (random)"
snps[which(endsWith(snps$ID, "base_rand_split2")),]$group <- "Test 2 (random)"
snps[which(endsWith(snps$ID, "base_rand_split3")),]$group <- "Test 3 (random)"
snps[which(endsWith(snps$ID, "base_rand_split4")),]$group <- "Test 4 (random)"
snps[which(endsWith(snps$ID, "base_rand_split5")),]$group <- "Test 5 (random)"
snps[which(endsWith(snps$ID, "baseline_split1")),]$group <- "Test 1"
snps[which(endsWith(snps$ID, "baseline_split2")),]$group <- "Test 2"
snps[which(endsWith(snps$ID, "baseline_split3")),]$group <- "Test 3"
snps[which(endsWith(snps$ID, "baseline_split4")),]$group <- "Test 4"
snps[which(endsWith(snps$ID, "baseline_split5")),]$group <- "Test 5"
snps <- rbind.fill(snps, snps_base)
snps <- snps[order(snps$group),]
ggplot(snps, aes(x=group, y=Y, fill=r2_test)) + geom_tile() +
        geom_text(aes(label=round(r2_test, 2)), size=2.5, color="white") +
        theme_minimal(base_size=8) + ylab("Conditions") + xlab("Test Set")
ggsave("Scripts/Data_Vis/Test_comparisons_SNPs.pdf", width=10, height=9)

## ORFs presence/absence
orfs_base <- read.csv("Results/RESULTS_RF_ORFs_baseline.txt", sep="\t")
orfs_base <- orfs_base[c("ID", "cond", "r2_test")]
colnames(orfs_base) <- c("ID", "Y", "r2_test")
orfs <- read.csv("/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results/RESULTS_reg.txt", sep="\t")
orfs <- orfs[c("ID", "Y", "r2_test")]
orfs$group = ""
orfs[which(endsWith(orfs$ID, "base_rand")),]$group <- "Original Test (random)" # these are all missing, need to run
orfs[which(endsWith(orfs$ID, "base_rand_split1")),]$group <- "Test 1 (random)"
orfs[which(endsWith(orfs$ID, "base_rand_split2")),]$group <- "Test 2 (random)"
orfs[which(endsWith(orfs$ID, "base_rand_split3")),]$group <- "Test 3 (random)"
orfs[which(endsWith(orfs$ID, "base_rand_split4")),]$group <- "Test 4 (random)"
orfs[which(endsWith(orfs$ID, "base_rand_split5")),]$group <- "Test 5 (random)"
orfs[which(endsWith(orfs$ID, "baseline_split1")),]$group <- "Test 1"
orfs[which(endsWith(orfs$ID, "baseline_split2")),]$group <- "Test 2"
orfs[which(endsWith(orfs$ID, "baseline_split3")),]$group <- "Test 3"
orfs[which(endsWith(orfs$ID, "baseline_split4")),]$group <- "Test 4"
orfs[which(endsWith(orfs$ID, "baseline_split5")),]$group <- "Test 5"
orfs <- orfs[which(orfs$group!=""),] # don't include FS results
orfs <- rbind.fill(orfs, orfs_base)
orfs[which(endsWith(orfs$ID, "baseline")),]$group <- "Original Test"
cnvs <- orfs[grep("cno", orfs$ID),] # ORFs copy number
orfs <- orfs[grep("orf", orfs$ID),]
orfs <- orfs[order(orfs$group),]
ggplot(orfs, aes(x=group, y=Y, fill=r2_test)) + geom_tile() +
        geom_text(aes(label=round(r2_test, 2)), size=2.5, color="white") +
        theme_minimal(base_size=8) + ylab("Conditions") + xlab("Test Set")
ggsave("Scripts/Data_Vis/Test_comparisons_ORFs.pdf", width=10, height=9)

## ORFs copy number
cnvs <- cnvs[order(cnvs$group),]
ggplot(cnvs, aes(x=group, y=Y, fill=r2_test)) + geom_tile() +
        geom_text(aes(label=round(r2_test, 2)), size=2.5, color="white") +
        theme_minimal(base_size=8) + ylab("Conditions") + xlab("Test Set")
ggsave("Scripts/Data_Vis/Test_comparisons_CNVs.pdf", width=10, height=9)



################################################################################
# BASELINE (ALL FEATURES) MULTI-OUTPUT RF PREDICTION PERFORMANCES
################################################################################
## SNPs



## ORFs presence/absence



## ORFs copy number


################################################################################
# BASELINE (ALL FEATURES) MULTI-TRAIT BAYESIAN MODEL PREDICTION PERFORMANCES
################################################################################
# Generate subsets of environments to build multi-trait models with
pheno <- read.csv("Data/Peter_2018/pheno.csv", row.names=1) # fitness data
pCorEnvs <- data.frame(matrix(nrow=0, ncol=4)) # collect correlation statistics
colnames(pCorEnvs) <- c("Env1", "Env2", "PCC", "p.value")
i <- 0
for (env1 in colnames(pheno)){ # calculate env pair fitness based pearson's r
        for (env2 in colnames(pheno)){
                if (env1 != env2){
                        res <- cor.test(pheno[,env1], pheno[,env2])
                        pCorEnvs[i,] <- c(env1, env2, res$estimate, res$p.value)
                        i = i + 1
                }
        }
}
pCorEnvs$p.value <- as.numeric(pCorEnvs$p.value)
pCorEnvs$PCC <- as.numeric(pCorEnvs$PCC)
pCorEnvs <- pCorEnvs[order(pCorEnvs$p.value, decreasing=F),]
tmp <- as.data.frame(t(apply(pCorEnvs[1:2], 1, sort))) # sort env pair names
pCorEnvs$Env1 <- tmp$V1 # set sorted env pair columns to remove duplicates later
pCorEnvs$Env2 <- tmp$V2
pCorEnvs <- pCorEnvs[!duplicated(pCorEnvs[,1:2]),] # remove duplicate pairs
write.csv(pCorEnvs,"Scripts/Data_Vis/pheno_pairs_cor.csv", quote=F, row.names=F)

## SNPs

## ORFs presence/absence

## ORFs copy number





###############################################################################
# Get the genes, GO, and pwy info for benomyl top 10 genes
setwd("/mnt/home/seguraab/Shiu_Lab/Project/")
library(tidyverse)
go <- read.csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_go_all.tsv", sep="\t")
go_orf <- read.csv("Data/Peter_2018/ORFs_GO_genes.tsv", sep="\t")
master_pwy <- read.csv("/mnt/home/seguraab/Shiu_Lab/Co-function/Data/MetaCyc/All-pathways-S288c_descriptions.txt", sep="\t")
pwy <- read.csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_pwy_all.csv")
pwy_orf <- read.csv("Data/Peter_2018/ORFs_and_S288C_genes_pwy_all.csv")
snp_shap <- list.files("Scripts/Data_Vis/SNP_Figures/SHAP", "median", full.names=T)
orf_shap <- list.files("Scripts/Data_Vis/ORF_CNV_Figures/SHAP_orf_presence_absence", "median", full.names=T)
cnv_shap <- list.files("Scripts/Data_Vis/ORF_CNV_Figures/SHAP_orf_copy_number", "median", full.names=T)

myfunc <- function(path, snp='y'){
        df <- read.csv(path, sep="\t")
        if (snp=='y') {
                df2 <- left_join(df, go[,c('snp', 'gene', 'BP', 'CC', 'MF')],
                        by=c('X'='snp'))
                df3 <- df2[c("gene", "BP", "CC", "MF")] # subset (will add shap values later)
                df3 <- df3[!duplicated(df3),] # drop dupes before joining to pwy
                df3 <- merge(df3, pwy[,c('Accession.1', 'Pathways.of.gene')],
                        by.x="gene", by.y="Accession.1", all.x=T) # add pathway IDs
                df3 <- df3[!duplicated(df3),] # drop duplicates
                df3 <- merge(df2[,c("X", "X0", "gene")], df3, by.x="gene",
                        by.y="gene", all.x=T) # add snps & shap values
                df3 <- merge(df3, master_pwy[,c("Object.ID", "Pathways")],
                        by.x="Pathways.of.gene", by.y="Object.ID", all.x=T) # add pathway description
                df3 <- df3[!duplicated(df3),] # drop duplicates
                df3 <- df3[order(df3$X0, decreasing=T),] # sort by shap value
                path <- gsub(".tsv", "_go_pwy.txt", path)
                write.table(df3, path, quote=F, sep="\t", row.names=F)
        }
        if (snp=='n') {
                df2 <- left_join(df, go_orf[c('qacc', 'gene', 'BP')],
                        by=c('X'='qacc'))
                df3 <- df2[c("gene", "BP")] # subset (will add shap values later)
                df3 <- df3[!duplicated(df3),] # drop dupes before joining to pwy
                df3 <- merge(df3, pwy_orf[,c('Accession.1', 'Pathways.of.gene')],
                        by.x="gene", by.y="Accession.1", all.x=T) # add pathway IDs
                df3 <- df3[!duplicated(df3),] # drop duplicates
                df3 <- merge(df2[,c("X", "X0", "gene")], df3, by.x="gene",
                        by.y="gene", all.x=T) # add snps & shap values
                df3 <- merge(df3, master_pwy[,c("Object.ID", "Pathways")],
                        by.x="Pathways.of.gene", by.y="Object.ID", all.x=T) # add pathway description
                df3 <- df3[!duplicated(df3),] # drop duplicates
                df3 <- df3[order(df3$X0, decreasing=T),] # sort by shap value
                path <- gsub(".txt", "_go_pwy.tsv", path)
                write.table(df3, path, quote=F, sep="\t", row.names=F)
        }
}

for (i in 1:length(orf_shap)){
        myfunc(snp_shap[i], snp='y')
        myfunc(orf_shap[i], snp='n')
        myfunc(cnv_shap[i], snp='n')
}