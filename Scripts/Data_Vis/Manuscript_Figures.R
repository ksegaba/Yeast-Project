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
# line : RF feature selection (FS) curves
# line : RF after FS: Test set performance comparisons
# line : Heritability $ RF FS performance comparison
# line : GO Enrichment (RF after FS models)
# line : Pathway Enrichment (RF after FS models)
# line : Gene RF importance scores (after FS) comparisons across data types
# line : Fitness Factors impacting RF performance
# line : RF performance comparisons using different test sets (baseline using all features and randomized label)
# line : Multi-Output RF prediction performances (baseline using all features)

rm(list=ls())
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(RColorBrewer))
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
# RANDOM FOREST PREDICTION PERFORMANCE (figures in excel)
################################################################################
# SNP Models
path <- "/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results"
rf_base <- read.csv(file.path(path, "baseline/RESULTS_reg_baseline.txt"), sep="\t")
rf_base <- right_join(conds, rf_base, by=c("cond"="Y")) # add full condition names
write.table(rf_base, "Results/RESULTS_RF_SNPs_baseline.txt", sep="\t", row.names=F, quote=F)

# PC Models
rf_fs <- read.csv(file.path(path, "RESULTS_reg.txt"), sep="\t")
rf_pc <- rf_fs[grep("PCs", rf_fs$ID),]# population structure results
rf_pc <- right_join(conds, rf_pc, by=c("cond"="Y")) # add full condition names
write.table(rf_pc, "Results/RESULTS_RF_PCs.txt", sep="\t", row.names=F, quote=F)

# ORF Copy Number Models
path <- "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results"
rf_orf <- read.csv(file.path(path, "RESULTS_reg.txt"), sep="\t")
rf_orf <- rf_orf[grep("2022-0[4-6]-[0-2][5-9]",rf_orf$DateTime),] # extract correct runs by date
rf_cnv <- rf_orf[grep("_cno_baseline", rf_orf$ID),] # ORF copy number models
which(duplicated(rf_cnv$ID)) # YPDDANISO20 is duplicated (almost exact numbers due to the nature of the algorithm)
rf_cnv <- rf_cnv[!duplicated(rf_cnv$ID),] # remove the 2022-06-09 run; keep april date for consistency purposes
rf_cnv <- right_join(conds, rf_cnv, by=c("cond"="Y")) # add full condition names
write.table(rf_cnv, "Results/RESULTS_RF_CNVs_baseline.txt", sep="\t", row.names=F, quote=F)

# ORF Presence/Absence Models
rf_orf <- rf_orf[grep("_orf_baseline", rf_orf$ID),] # ORF presence/absence
rf_orf[which(duplicated(rf_orf$ID)),] # 8 traits duplicated
rf_orf <- rf_orf[!duplicated(rf_orf$ID),] # remove duplicates that are later by date
rf_orf <- right_join(conds, rf_orf, by=c("cond"="Y")) # add full condition names & remove duplicate rows
write.table(rf_orf, "Results/RESULTS_RF_ORFs_baseline.txt", sep="\t", row.names=F, quote=F)

################################################################################
# RANDOM FOREST FEATURE SELECTION CURVES
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

# SNP Models
path <- "/mnt/gs21/scratch/seguraab/yeast_project/SNP_yeast_RF_results/fs"
rf <- read.csv(file.path(path, "RESULTS_reg.txt"), sep="\t")
rf_fs <- rf[grep("rf_[0-9]+", rf$ID),] # SNP feature selection results
rf_fs <- right_join(conds, rf_fs, by=c("cond"="Y")) 
plot_fs(rf_fs, "_FS.pdf") # 2 to 40k features
rf_exp_fs <- rf_fs[grep("exp", rf_fs$ID),] 
plot_fs(rf_exp_fs, "_exp_FS.pdf") # 2 to 2048 features
out <- get_opt_feat(rf_fs) # optimal models
write.table(out, "Results/RESULTS_RF_SNPs_FS.txt", sep="\t", row.names=F, quote=F)

# ORF Copy Number Models
path <- "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results"
rf_orf <- read.csv(file.path(path, "RESULTS_reg.txt"), sep="\t")
rf_cnv_fs <- rf_orf[grep("_cno_[^b]", rf_orf$ID),]
plot_fs(rf_cnv_fs, "_cnv_FS.pdf")
out <- get_opt_feat(rf_cnv_fs) # optimal models
write.table(out, "Results/RESULTS_RF_CNVs_FS.txt", sep="\t", row.names=F, quote=F)

# ORF Presence/Absence Models
rf_orf_fs <- rf_orf[grep("_orf_[0-9]+", rf_orf$ID),]
plot_fs(rf_orf_fs, "_orf_FS.pdf")
out <- get_opt_feat(rf_orf_fs) # optimal models
write.table(out, "Results/RESULTS_RF_ORFs_FS.txt", sep="\t", row.names=F, quote=F)

################################################################################
# RANDOM FOREST AFTER FEATURE SELECTION: TEST SET PERFORMANCE COMPARISONS
################################################################################
pcs <- read.table("Results/RESULTS_RF_PCs.txt", sep="\t", header=T, row.names=3)
snp <- read.table("Results/RESULTS_RF_SNPs_FS.txt", sep="\t", header=T, row.names=2)
orf <- read.table("Results/RESULTS_RF_ORFs_FS.txt", sep="\t", header=T)
cnv <- read.table("Results/RESULTS_RF_CNVs_FS.txt", sep="\t", header=T)
orf <- right_join(conds, orf, by=c("cond"="Y"))
rownames(orf) <- orf$new_cond
cnv <- right_join(conds, cnv, by=c("cond"="Y"))
rownames(cnv) <- cnv$new_cond

# Bar plots with validation and PC performances for SNPs, ORFs, and CNVs
pcs <- pcs[order(rownames(pcs)),] # ensure all environments are in the same order
snp <- snp[rownames(pcs),]
orf <- orf[rownames(pcs),]
cnv <- cnv[rownames(pcs),]
rownames(pcs)==rownames(snp)
rownames(pcs)==rownames(orf)
rownames(pcs)==rownames(cnv)

ggplot() + geom_bar(aes(x=new_cond,y=r2_val), snp, fill="#90CE4F", stat="identity") +
        geom_errorbar(aes(x=new_cond, ymin=r2_val-r2_val_sd, 
                ymax=r2_val+r2_val_sd), snp, width=0.2) +
        geom_point(aes(x=rownames(pcs), y=r2_test), pcs, size=2.5, col="#702EA0") +
        geom_tile(aes(x=new_cond, y=r2_test, width=.9, height=.005), snp, fill="black") +
        theme_bw(10) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("Scripts/Data_Vis/snp_RF_performance_after_fs.pdf", width=8, height=4, units="in")

ggplot() + geom_bar(aes(x=new_cond,y=r2_val), orf, fill="#FFC000", stat="identity") +
        geom_errorbar(aes(x=new_cond, ymin=r2_val-r2_val_sd, 
                ymax=r2_val+r2_val_sd), orf, width=0.2) +
        geom_point(aes(x=rownames(pcs), y=r2_test), pcs, size=2.5, col="#702EA0") +
        geom_tile(aes(x=new_cond, y=r2_test, width=.9, height=.005), orf, fill="black") +
        theme_bw(10) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("Scripts/Data_Vis/orf_RF_performance_after_fs.pdf", width=8, height=4, units="in")

ggplot() + geom_bar(aes(x=new_cond,y=r2_val), cnv, fill="#5B9BD5", stat="identity") +
        geom_errorbar(aes(x=new_cond, ymin=r2_val-r2_val_sd, 
                ymax=r2_val+r2_val_sd), cnv, width=0.2) +
        geom_point(aes(x=rownames(pcs), y=r2_test), pcs, size=2.5, col="#702EA0") +
        geom_tile(aes(x=new_cond, y=r2_test, width=.9, height=.005), cnv, fill="black") +
        theme_bw(10) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("Scripts/Data_Vis/cnv_RF_performance_after_fs.pdf", width=8, height=4, units="in")

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
write.table(test, "Results/RESULTS_RF_FS_ALL.txt", sep="\t", quote=F, row.names=F)

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
h2 <- h2[order(h2$h2),] # sort in ascending order
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
# GENE IMPORTANCE SCORE (AFTER FS) COMPARISONS ACROSS DATA TYPES (what about across envs per data type? too much?)
################################################################################
## SNPs
traits <- c("YPACETATE", "YPD6AU", "YPD14", "YPD40", "YPD42", "YPDANISO10", 
        "YPDANISO20", "YPDANISO50", "YPDBENOMYL200", "YPDBENOMYL500", 
        "YPDCAFEIN40", "YPDCAFEIN50", "YPDCHX1", "YPDCHX05", "YPDCUSO410MM", 
        "YPDDMSO", "YPDETOH", "YPDFLUCONAZOLE", "YPDFORMAMIDE4", "YPDFORMAMIDE5", 
        "YPDHU", "YPDKCL2M", "YPDLICL250MM", "YPDMV", "YPDNACL1M", "YPDNACL15M", 
        "YPDNYSTATIN", "YPDSDS", "YPDSODIUMMETAARSENITE", "YPETHANOL", 
        "YPGALACTOSE", "YPGLYCEROL", "YPRIBOSE", "YPSORBITOL", "YPXYLOSE")
feats <- c(64, 3000, 512, 1024, 2000, 2000, 2000, 256, 1024, 4000, 3000, 2000, 
        512, 256, 1000, 512, 256, 1000, 128, 1024, 512, 512, 512, 1000, 3000, 
        1024, 32, 1000, 256, 1024, 512, 256, 512, 1000, 256)
map <- read.csv(
        "~/Shiu_Lab/Project/Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes.txt",
        header=F) # snp to gene map
colnames(map) <- c("snp", "chr", "pos", "gene")
# figure: how many snps per gene
counts <- map %>% group_by(gene) %>% tally()
ggplot(counts, aes(n)) + geom_histogram(bins=100) + xlim(c(0,100))
ggsave("Scripts/Data_Vis/SNPs_per_gene_counts.pdf")
table(counts$n) # most genes have about 40 or less snps

dir <- "/mnt/scratch/seguraab/yeast_project/SNP_yeast_RF_results/fs" # data directory
imp_df <- read.csv(file.path(dir,"YPACETATE_exp_rf_64_imp"), sep="\t") # first trait
colnames(imp_df) <- c("snp", "YPACETATE")
dim(imp_df)
imp_df <- imp_df[order(imp_df$YPACETATE, decreasing=T),]

combine <- function(imp_df, traits, feats, map, col, save){
        # function to combine the rest of the traits with imp_df
        for (i in 2:length(traits)){
                # df <- read.csv(list.files(
                #         dir, paste(traits[i], "[exp_rf]+", feats[i], "imp", sep="_"), #"orf"
                #         full.names=T, recursive=F), sep="\t")
                df <- read.csv(list.files(
                        dir, paste(traits[i], "cno", feats[i], "imp", sep="_"),
                        full.names=T, recursive=F), sep="\t")
                # select top 20 only and combine
                df <- df[order(df$mean_imp, decreasing=T),]
                colnames(df) <- c(col, traits[i])
                imp_df <- full_join(imp_df, df[1:20,], by=col)
                # somehow need to incorporate percentile information to make heatmap look better
        }
        print(dim(imp_df))
        
        # set NAs to zero, since they are not considered "top predictors"
        imp_df[is.na(imp_df)] <- 0

        # add gene info and aggregate by gene
        imp_df <- left_join(imp_df, map[c(col, "gene")], by=col)
        imp_df_gene <- imp_df[-1] %>% group_by(gene) %>% 
                summarise_all(median) %>% # list(median, mean, min, max, sd)
                as.data.frame()
        imp_df_gene <- imp_df_gene %>% group_by(gene) %>% 
                summarise_all(median) %>% # list(median, mean, min, max, sd)
                as.data.frame() # drop duplicate genes
        print(dim(imp_df_gene))
        imp_df_gene <- imp_df_gene[!is.na(imp_df_gene$gene),]
        print(dim(imp_df_gene))
        
        # write.csv(imp_df_gene, "Scripts/Data_Vis/SNPs_top20_all_env_imp.csv", quote=F, row.names=F)
        write.csv(imp_df_gene, "Scripts/Data_Vis/ORFs_pres_abs_top20_all_env_imp.csv", quote=F, row.names=F)
        write.csv(imp_df_gene, "Scripts/Data_Vis/ORFs_copy_num_top20_all_env_imp.csv", quote=F, row.names=F)
        
        rownames(imp_df_gene) <- imp_df_gene$gene
        imp_df_gene <- imp_df_gene[-1]
        
        
        # heatmap of importance scores
        palette_func <- colorRampPalette(c("red", "white"))
        palette <- rev(palette_func(2)) 
        col_fun <- colorRamp2(c(0, max(imp_df_gene)), palette)
        p <- Heatmap(as.matrix(imp_df_gene), 
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
                name="log10(q-values)",
                heatmap_legend_param=list(title="Feature Importance", 
                at=c(0, max(imp_df_gene))))
        pdf(save, height=unit(30, "in"), width=unit(10, "in"))
        draw(p)
        dev.off()
}
imp_df_orf <- combine(imp_df_orf[1:20,], traits, feats, "orf", "Scripts/Data_Vis/ORFs_top20_all_env_imp.pdf")

## ORF presence/absence
feats <- c(250, 1250, 250, 500, 250, 500, 500, 250, 750, 250, 500, 500, 250, 
        250, 500, 500, 250, 250, 250, 250, 250, 250, 750, 250, 250, 250, 750, 
        250, 250, 250, 500, 250, 250, 250, 250)
map <- read.csv("~/Shiu_Lab/Project/Data/Peter_2018/ORFs_GO_genes.tsv", sep="\t")
map <- map[c("qacc", "gene")]
map <- map[!duplicated(map),]
colnames(map) <- c("orf", "gene")
dir <- "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_RF_results" # data directory
imp_df_orf <- read.csv(file.path(dir,"YPACETATE_cno_250_imp"), sep="\t") # first trait
colnames(imp_df_orf) <- c("orf", "YPACETATE")
dim(imp_df_orf)
#imp_df_orf <- combine(imp_df_orf, traits, feats)

## ORF copy number
feats <- c(250, 500, 250, 250, 250, 500, 250, 250, 500, 250, 250, 250, 500, 
        250, 500, 500, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 
        250, 250, 250, 250, 250, 500, 250, 250)

## SNPs to ORF presence/absence
imp_snp <- read.csv("/mnt/scratch/seguraab/yeast_project/SNP_yeast_RF_results/fs/YPDCUSO410MM_rf_1000imp.csv")
imp_orf <- read.csv("Scripts/Data_Vis/ORFs_pres_abs_top20_all_env_imp.csv")
imp_cnv <- read.csv("Scripts/Data_Vis/ORFs_copy_num_top20_all_env_imp.csv")

# imp_snp <- read.csv("Scripts/Data_Vis/SNPs_top20_all_env_imp.csv", row.names=1)
# imp_orf <- read.csv("Scripts/Data_Vis/ORFs_pres_abs_top20_all_env_imp.csv", row.names=1)
# imp_cnv <- read.csv("Scripts/Data_Vis/ORFs_copy_num_top20_all_env_imp.csv", row.names=1)
palette_func <- colorRampPalette(c("red", "white"))
palette <- rev(palette_func(2)) 
col_fun <- colorRamp2(c(0, max(imp_orf)), palette)
p <- Heatmap(as.matrix(imp_orf), 
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
        name="log10(q-values)",
        heatmap_legend_param=list(title="Feature Importance", 
        at=c(0, max(imp_orf))))
pdf("Scripts/Data_Vis/SNPs_top20_all_env_imp.csv", height=unit(30, "in"), width=unit(10, "in"))
draw(p)
dev.off()
pdf("Scripts/Data_Vis/ORFs_pres_abs_top20_all_env_imp.csv", height=unit(30, "in"), width=unit(10, "in"))
draw(p)
dev.off()
pdf("Scripts/Data_Vis/ORFs_copy_num_top20_all_env_imp.csv", height=unit(30, "in"), width=unit(10, "in"))
draw(p)
dev.off()

toplot <- inner_join(imp_snp[c('gene', 'YPDCUSO410MM')], imp_orf[c('gene', 'YPDCUSO410MM')], by='gene')
dim(toplot)
ggplot(toplot, aes(x=YPDCUSO410MM.x, y=YPDCUSO410MM.y) ) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()
ggsave("Scripts/Data_Vis/SNPs_vs_ORF_pres_abs_YPDCUSO410MM_imp.pdf")

## SNPs to ORF copy number
toplot <- inner_join(imp_snp[c('gene', 'YPDCUSO410MM')], imp_cnv[c('gene', 'YPDCUSO410MM')], by='gene')
dim(toplot)
ggplot(toplot, aes(x=YPDCUSO410MM.x, y=YPDCUSO410MM.y) ) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()
ggsave("Scripts/Data_Vis/SNPs_vs_ORF_copy_num_YPDCUSO410MM_imp.pdf")

## ORF presence/absence to ORF copy number

################################################################################
# FITNESS FACTORS IMPACTING RF PERFORMANCE
################################################################################
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
## SNPs



## ORFs presence/absence



## ORFs copy number

################################################################################
# BASELINE (ALL FEATURES) MULTI-OUTPUT RF PREDICTION PERFORMANCES
################################################################################
## SNPs



## ORFs presence/absence



## ORFs copy number
