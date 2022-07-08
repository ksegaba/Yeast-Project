# Chapter 1 Figures
library(data.table)
library(gplots)
library(ggplot2)
library(viridis)
library(MASS)
library(dplyr)
library(tidyr)
library(matrixStats)

setwd("/mnt/home/seguraab/Shiu_Lab/Project")

# Isolate growth condition labels
cond <- c("YPACETATE", "YPD14", "YPD40", "YPD42", "YPD6AU", "YPDANISO10", 
        "YPDANISO20", "YPDANISO50", "YPDBENOMYL200", "YPDBENOMYL500", "YPDCAFEIN40", 
        "YPDCAFEIN50", "YPDCHX05", "YPDCHX1", "YPDCUSO410MM", "YPDDMSO", "YPDETOH", 
        "YPDFLUCONAZOLE", "YPDFORMAMIDE4", "YPDFORMAMIDE5", "YPDHU", "YPDKCL2M", 
        "YPDLICL250MM", "YPDMV", "YPDNACL15M", "YPDNACL1M", "YPDNYSTATIN", "YPDSDS", 
        "YPDSODIUMMETAARSENITE", "YPETHANOL", "YPGALACTOSE", "YPRIBOSE", "YPGLYCEROL", 
        "YPXYLOSE", "YPSORBITOL")
new_cond <- c("YP Acetate 2%", "YPD 14C", "YPD 40C", "YPD 42C", "YPD 6-Azauracile 600µg/ml",
        "YPD Anisomycin 10µg/ml", "YPD Anisomycin 20µg/ml", "YPD Anisomycin 50µg/ml",
        "YPD Benomyl 200µg/ml", "YPD Benomyl 500µg/ml", "YPD Caffeine 40mM", "YPD Caffeine 50mM",
        "YPD Cycloheximide 0.5µg/ml", "YPD Cycloheximide 1µg/ml", "YPD CuSO4 10mM", "YPD DMSO 6%",
        "YPD Ethanol 15%", "YPD Fluconazole 20µg/ml", "YPD Formamide 4%", "YPD Formamide 5%",
        "YPD Hydroxyurea 30mg/ml", "YPD KCL 2M", "YPD LiCl 250mM", "YPD Methylviologen 20mM",
        "YPD NaCl 1.5M", "YPD NaCl 1M", "YPD Nystatin 10µg/ml", "YPD SDS 0.2%", 
        "YPD Sodium metaarsenite 2.5mM", "YP Ethanol 2%", "YP Galactose 2%", "YP Ribose 2%",
        "YP Glycerol 2%", "YP Xylose 2%", "YP Sorbitol 2%")
conds <- as.data.frame(cbind(cond, new_cond))

#--------- LD Decay ---------#
ld <- fread("Data/Peter_2018/geno_LD_window_50.txt") # Window size 50, TASSEL5 output
ld_sub <- ld[which(ld$Dist_bp != "N/A" & ld$`R^2` != "NaN"), c("Dist_bp", "R^2", "DPrime", "pDiseq")] # remove rows with missing values
ld_sub$Dist_bp <- as.numeric(ld_sub$Dist_bp)
ld_sub$bin <- cut(ld_sub$Dist_bp, breaks=c(seq(1, 30976, by=50))) # bin by distance
ld_sub$bin <- as.character(lapply(strsplit(as.character(ld_sub$bin), split=","),head, n=1)) # rename bins
ld_sub$bin <- gsub("\\(", "", ld_sub$bin)

pdf("Scripts/Data_Vis/ld_w50.pdf")
ld_sub %>% group_by(bin) %>% summarise(med_R2=median(`R^2`), sd_R2=sd(`R^2`)) %>%
ggplot(aes(x=as.numeric(bin), y=med_R2)) + geom_point() + theme_bw() +
        xlab("Distance (bp)") + ylab(parse(text="R^2"))
dev.off()

#--------- KINSHIP ---------#
kin <- fread("Data/Peter_2018/geno_transposed.csv_ordered_kinship.txt")
kin <- as.matrix(kin, rownames=1, colnames=1)
colnames(kin) <- rownames(kin)
write.csv(kin, "Scripts/Data_Vis/kinship.csv", quote=F, row.names=T)

pdf("Scripts/Data_Vis/kinship.pdf")
hm <- heatmap.2(kin, 
                col = colorRampPalette(c("blue","white","red"))(21),
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

#--------- FITNESS CORRELATIONS ---------#
pCorEnvs <- read.csv("Data/Peter_2018/pheno_corr_envs.csv", header=T, row.names=1) # Correlation between conditions
colnames(pCorEnvs) <- rownames(pCorEnvs)
pCorIso <- read.csv("Data/Peter_2018/pheno_corr_isolates.csv", header=T, row.names=1) # Correlation between isolates across all conditions
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


pdf("Scripts/Data_Vis/pheno_corr_envs_v2.pdf")
hm3 <- heatmap.2(as.matrix(pCorEnvs),
                col = colorRampPalette(c("blue","white","red"))(21),
                #hclustfun = function(x) hclust(x, method="complete"),
                trace = "none",
                dendrogram = "none",
                notecex = 1,
                cex.axis = 1,
                cexRow = 0.6,
                cexCol = 0.6,
                density.info = "none",
                keysize = 1,
                key.title = "PCC",
                key.par = list(cex=1, mar=c(3,1,3,0)),
                margins = c(11,11))
dev.off()

#--------- ISOLATE FITNESS CORRELATION VS KINSHIP ---------#
k_bin <- reshape2::melt(kin) # reshape kinship matrix
#k_quant <- quantile(k_bin, seq(0,1,0.05)) # kinship quantiles
k_breaks <- seq(min(k_bin$value), max(k_bin$value), 0.5) # regular bins
#k_bin$bin <- cut(k_bin$value, breaks=k_quant) # quantile bins
k_bin$breaks <- cut(k_bin$value, breaks=k_breaks) # regular bins
#k_bin$bin <- as.character(lapply(strsplit(as.character(k_bin$bin),split=","),head,n=1)) # rename bins
k_bin$breaks <- as.character(lapply(strsplit(as.character(k_bin$breaks), split=","), head, n=1))
#k_bin$bin <- gsub("\\(", "", k_bin$bin) # remove
k_bin$breaks <- gsub("\\(", "", k_bin$breaks)
#k_bin$bin <- as.numeric(k_bin$bin) # convert to numeric
k_bin$breaks <- as.numeric(k_bin$breaks)

pCorIso <- as.matrix(pCorIso) # original, no values have been reset to 0
i_bin <- reshape2::melt(pCorIso) # reshape fitness across isolates correlation matrix
#i_quant <- quantile(i_bin$value, seq(0,1,0.05)) # pCor quantiles
i_breaks <- seq(min(i_bin$value), 1.01, 0.12)
#i_bin$bin <- cut(i_bin$value, breaks=i_quant)
i_bin$breaks <- cut(i_bin$value, breaks=i_breaks)
#i_bin$bin <- as.character(lapply(strsplit(as.character(i_bin$bin),split=","),head,n=1))
i_bin$breaks <- as.character(lapply(strsplit(as.character(i_bin$breaks), split=","), head, n=1))
#i_bin$bin <- gsub("\\(", "", i_bin$bin)
i_bin$breaks <- gsub("\\(", "", i_bin$breaks)
#i_bin$bin <- as.numeric(i_bin$bin)
i_bin$breaks <- as.numeric(i_bin$breaks)

k_i_bin <- merge(k_bin, i_bin, by=c("Var1", "Var2")) # merge dataframes
colnames(k_i_bin) <- c("Isolate 1", "Isolate 2", "Kinship", "Kinship bin", "pCor", "pCor bin")

k.i.cor <- cor(k_i_bin$Kinship, k_i_bin$pCor, method="spearman") # Spearman correlation
k_i_bin_count <- aggregate(cbind(count = `Isolate 2`) ~ `Kinship bin`, # bin counts
        data=k_i_bin, FUN=function(x){NROW(x)})
k_i_bin_median <- k_i_bin %>% group_by(`Kinship bin`) %>% summarise(median=median(pCor)) # kinship bin median pCor

get_density <- function(x, y, ...) { # density
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])}
k_i_bin$density <- get_density(k_i_bin$Kinship, k_i_bin$pCor, n=100)

# 5% and 95% quantiles for each kinship bin
#quants <- as.data.frame(do.call("rbind", tapply(k_i_bin$pCor, k_i_bin$`Kinship bin`, quantile, c(0.05, 0.95))))
#names(quants) <- c("x5","x95")
#quants$`Kinship bin` <- as.numeric(row.names(quants))

# Density scatter plot
ggplot(k_i_bin, aes(x=Kinship, y=pCor, color=density)) + geom_point(size=0.2, alpha=0.2) +
        geom_smooth(method="lm", formula=y~x, se=FALSE) +
        scale_color_viridis(direction=-1) + theme_bw(8)
ggsave("Scripts/Data_Vis/kinpCorIso_Density.pdf", width=5, height=4, device="pdf", useDingbats=FALSE)

ggplot(k_i_bin, aes(x=Kinship, y=pCor, color=density)) + geom_point(size=0.2, alpha=0.2) +
        scale_color_viridis(direction=-1) + theme_bw(base_size=12) + 
        theme(axis.text.x=element_text(color="black"), axis.text.y = element_text(color="black")) +
        geom_line(data=k_i_bin_median, aes(x=`Kinship bin`, y=median, group=1), color="black")
        #annotate(geom="text", x=2, y=-0.3, label=paste("Spearman's rank correlation = ", round(k.i.cor, digits=2),sep=""))
ggsave("Scripts/Data_Vis/kinpCorIso_Density_median.pdf", width=5, height=4, device="pdf", useDingbats=FALSE)

ggplot(k_i_bin, aes(x=`Kinship bin`, y=pCor, color=density)) + geom_point(size=0.2, alpha=0.2) +
        scale_color_viridis(direction=-1) + theme_bw(8) + 
        geom_line(data=k_i_bin_median, aes(x=`Kinship bin`, y=median, group=1), color="red")
ggsave("Scripts/Data_Vis/kinpCorIso_Density_median_bin.pdf", width=5, height=4, device="pdf", useDingbats=FALSE)

# Binned plot with quantiles
ggplot(data = k_i_bin, aes(x=`Kinship bin`, y=pCor, group=`Kinship bin`))  + geom_violin() +
        geom_boxplot(outlier.shape=NA) + theme_bw(10) +
        geom_line(data=quants, aes(x=`Kinship bin`, y=x95, group = 1), color="red") +
        geom_line(data=quants, aes(x=`Kinship bin`, y=x5, group = 1), color="blue")
ggsave("Scripts/Data_Vis/kinpCorIso_Distribution.pdf", width=4, height = 4, device="pdf", useDingbats=FALSE)

# Bin counts plot
ggplot(k_i_bin_count, aes(x=`Kinship bin`, y=log(count))) + theme_bw(10) +geom_bar(stat="identity", position= "dodge")
ggsave("Scripts/Data_Vis/kinpCorIso_BinCount.pdf", width=3.42, height = 2.42, device="pdf", useDingbats=FALSE)

#--------- ORF PRESENCE/ABSENCE & COPY NUMBER VARIATIONS ---------#
orf <- fread("Data/Peter_2018/ORFs_pres_abs.csv") # ORF presence/absence
orf <- as.matrix(orf, rownames=1, colnames=1)
oCor <- cor(t(orf)) # between isolates
oCor[oCor < 0.9] <- 0.9

cno <- fread("Data/Peter_2018/ORFs_no_NA.csv") # ORF copy number variants
cno <- as.matrix(cno, rownames=1, colnames=1)
cCor <- cor(t(cno))
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

#--------- HERITABILITY ---------#


#--------- FITNESS VARIANCE ACROSS CONDITIONS ---------#
y <- read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv", row.names=1)
iso_var <- as.matrix(rowVars(as.matrix(y)))
rownames(iso_var) <- rownames(y)
iso_var <- as.data.frame(iso_var)

y_var <- as.matrix(colVars(as.matrix(y)))
rownames(y_var) <- colnames(y)
y_var <- as.data.frame(y_var)
y_var$Envs <- rownames(y_var)
y_var <- merge(y_var, conds, by.x="Envs", by.y="cond")
y_var <- y_var[order(y_var$V1),]

ggplot(iso_var, aes(x=V1)) + geom_histogram(bins=40) + 
        #scale_x_continuous(breaks=seq(0,400,25)) + 
        labs(x="Isolate Fitness Variance", y="Counts")
ggsave("Scripts/Data_Vis/isolate_fitness_var_dist.pdf", device="pdf", useDingbats=FALSE)

ggplot(y_var, aes(x=V1, y=new_cond)) + geom_histogram(stat="identity") +
        labs(x="Fitness Variance", y="Environment") #+ coord_flip()
ggsave("Scripts/Data_Vis/envs_fitness_var_dist.pdf", device="pdf", useDingbats=FALSE)

#--------- SNP COUNTS PER GENE DISTRIBUTION ---------#
genes <- read.csv("Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes.txt", sep=",", header=FALSE)
genes2 <- genes[!(genes$V4=="intergenic"),]# remove intergenic
snp_n <- as.data.frame(table(genes2$V4)) # snp counts per gene

ggplot(snp_n, aes(x=Freq)) + geom_histogram(bins=40) + 
        scale_x_continuous(breaks=seq(0,400,25)) + 
        labs(x="Number of SNPs", y="Counts")
ggsave("Scripts/Data_Vis/snp_counts_dist_genes.pdf", device="pdf", useDingbats=FALSE)

#--------- GENOMIC PREDICTION ACCURACY ---------#
# Read in SNP-based model results
rrBLUP_PCs_cv = read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Results/rrBLUP_PCs_average_cv_R2.csv") # rrBLUP population structure validation
rrBLUP_PCs_test = read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Results/rrBLUP_PCs_average_test_R2.csv") # rrBLUP pop. struc. test
rrBLUP_cv = read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Results/rrBLUP_average_cv_R2.csv")
rrBLUP_test = read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Results/rrBLUP_average_test_R2.csv")
#GBLUP = read.csv() # Genomic BLUP (got deleted from scratch or something; need to rerun; or just bglr gblup and see if it agrees)
BL_cv = read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Results/BL_average_cv_R2.csv") # Bayes LASSO
BL_test = read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Results/BL_average_test_R2.csv")
BC = read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_BayesC/Results/RESULTS_BayesC.csv") # BayesC
bc <- BC[,c(1,5)] ; colnames(bc) <- c("Trait", "R2_test_mean") ; bc <- merge(conds, bc, by.x="cond", by.y="Trait")
RF = read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Results/RF_average_R2_.csv") # Random forest
rf <- RF[,c(2,3)]
XGB = read.delim("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_XGBoost/RESULTS_xgboost.txt", sep="\t") # XGBoost
xgb <- XGB[grep("2022-04-2[7-8]", XGB$Date),]
xgb_snp <- xgb[which(xgb$Data=="SNP"),]
xgb_snp <- xgb_snp[,c(4,8)] ; xgb_snp <- merge(conds, xgb_snp, by.x="cond", by.y="Trait")
xgb_snp <- xgb_snp %>% group_by(new_cond) %>% summarise_all(mean)
#GCN = read.csv() # GCN (will include later)

# Read in ORF/CNO-based model results
orf_rrblup = read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Results/rrBLUP_ORF_pres_abs_average_test_R2.csv")
cno_rrblup = read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Results/rrBLUP_CNO_average_test_R2.csv")
#orf_blup = read.csv() # this data got deleted, so I won't include it, instead I'll just show rrblup, xgb and rf
#cno_blup = read.csv()
orf_xgb = xgb[which(xgb$Data=="ORF"),]
cno_xgb = xgb[which(xgb$Data=="CNO"),]
orf_xgb <- merge(conds, orf_xgb, by.x="cond", by.y="Trait")
cno_xgb <- merge(conds, cno_xgb, by.x="cond", by.y="Trait")
orf_xgb <- orf_xgb[,c(2,9)]
cno_xgb <- cno_xgb[,c(2,9)]
ORF_RF = read.csv("/mnt/scratch/seguraab/yeast_project/ORF_yeast_RF_results/RESULTS_reg.txt", sep="\t")
ORF_RF <- ORF_RF[grep("2022-0[4-6]-[0-2][5-9]",ORF_RF$DateTime),]
orf_rf <- ORF_RF[grep("_orf_",ORF_RF$ID),]
cno_rf <- ORF_RF[grep("_cno_",ORF_RF$ID),] # need YPDSODIUMMETAARSENITE, YPDNACL1M, YPDANISO20, YPDDMSO
orf_rf <- orf_rf[, c(1,3,29)] ; orf_rf <- orf_rf %>% group_by(ID) %>% summarise_all(mean)
cno_rf <- cno_rf[,c(1,3,29)] ; cno_rf <- cno_rf %>% group_by(ID) %>% summarise_all(mean)
orf_rf$ID <- gsub("_orf_baseline", "", orf_rf$ID)
cno_rf$ID <- gsub("_cno_baseline", "", cno_rf$ID)
orf_rf <- merge(conds, orf_rf, by.x="cond", by.y="ID")
cno_rf <- merge(conds, cno_rf, by.x="cond", by.y="ID")

# Combine datasets
snps_cv <- merge(merge(merge(merge(merge(merge(merge( # not showing cv in lab meeting for now
                rrBLUP_PCs_cv, rrBLUP_cv, by="X"), 
                GBLUP, by="X"), 
                BL_cv, by="X"), 
                BC, by="X"), 
                RF, by="X"), 
                XGB, by="X"),
                GCN, by="X")
snps_cv <- snps_cv[,-c(5,6,8,9)]
colnames(snps_cv) <- c("Condition", "rrBLUP PCs", "rrBLUP", "GBLUP", "BL", "BC", "RF", "XGB", "GCN")

snps_test <-merge(merge(merge(merge(#merge(
                rrBLUP_PCs_test, rrBLUP_test, by="X"), 
                BL_test, by="X"), 
                bc, by.x="X", by.y="new_cond"), 
                rf, by.x="X", by.y="ID")#, 
                #xgb_snp, by.x="X", by.y="new_cond") # add gcn later
snps_test <- snps_test[,-5]
#snps_test <- snps_test[,-c(5,8)]
colnames(snps_test) <- c("Condition", "rrBLUP PCs", "rrBLUP", "Bayes LASSO", "BayesC", "Random Forest")#, "XGBoost")

orfs_test <- merge(merge(orf_rrblup, orf_rf, by.x="X", by.y="new_cond"),
                orf_xgb, by.x="X", by.y="new_cond")
orfs_test <- orfs_test[,-c(3,4)]
colnames(orfs_test) <- c("Condition", "rrBLUP", "Random Forest", "XGBoost")

cnos_test <- merge(merge(cno_rrblup, cno_rf, by.x="X", by.y="new_cond"),
                cno_xgb, by.x="X", by.y="new_cond")
cnos_test <- cnos_test[,-c(3,4)]
colnames(cnos_test) <- c("Condition", "rrBLUP", "Random Forest", "XGBoost")

# Heatmap (vertical)
snps_test2 <- melt(snps_test) ; colnames(snps_test2) <- c("Condition", "Model", "R-squared")
orfs_test2 <- melt(orfs_test) ; colnames(orfs_test2) <- c("Condition", "Model", "R-squared")
cnos_test2 <- melt(cnos_test) ; colnames(cnos_test2) <- c("Condition", "Model", "R-squared")
#snps_test <- snps_test[-35,]
rownames(snps_test) <- snps_test$Condition
snps_test <- snps_test[,-1]
rownames(orfs_test) <- orfs_test$Condition
orfs_test <- orfs_test[,-1]
rownames(cnos_test) <- cnos_test$Condition
cnos_test <- cnos_test[,-1]

combined <- merge(merge(snps_test, orfs_test, by="row.names"),cnos_test, by.x="Row.names", by.y="row.names")
combined <- combined[,c(1,6,8,11)]
colnames(combined) <- c("Conditions", "SNPs", "ORFs", "CNVs")
rownames(combined) <- combined$Conditions
combined <- combined[,-1]
write.csv(combined, "Scripts/Data_Vis/RF_performances_snps_orfs.csv", quote=F)

library(RColorBrewer)
Colors=rev(brewer.pal(6,"YlOrRd"))
Colors=colorRampPalette(Colors)(100)
par(mar=c(30,30,30,30)+1) #bottom,left,top,right
pdf("Scripts/Data_Vis/all_model_performances_snps2.pdf", height=400, width=400)
#heatmap(as.matrix(snps_test), Colv=NA, xlab="Models", ylab="Environments")
heatmap.2(as.matrix(snps_test),           # cell labeling
          cellnote=round(as.matrix(snps_test),2),
          dendrogram="none",
          Colv=F,
          cexRow=36, cexCol=36, # font size
          notecex=36.0, # cell font size
          notecol="black",
          na.color=par("bg"),
          breaks=c(0,0.2,0.4,0.6),
          col=c("white", "orange", "red"),
          trace="none",
          margins=c(250,500), # reduce label cutoff
          key.title="R-squared",
          keysize=1.4,
          key.par = list(cex=36), mar=c(2,0,2,0))
graphics.off()

pdf("Scripts/Data_Vis/all_model_performances_orfs.pdf", height=11, width=11)
#heatmap(as.matrix(snps_test), Colv=NA, xlab="Models", ylab="Environments")
heatmap.2(as.matrix(orfs_test),           # cell labeling
          cellnote=round(as.matrix(orfs_test),2),
          dendrogram="row",
          Colv=F,
          notecex=1.0,
          notecol="black",
          na.color=par("bg"),
          col="viridis",
          trace="none",
          key.title="R-squared")
dev.off()

pdf("Scripts/Data_Vis/all_model_performances_cnos.pdf", height=11, width=11)
#heatmap(as.matrix(snps_test), Colv=NA, xlab="Models", ylab="Environments")
heatmap.2(as.matrix(cnos_test),           # cell labeling
          cellnote=round(as.matrix(cnos_test),2),
          dendrogram="row",
          Colv=F,
          notecex=1.0,
          notecol="black",
          na.color=par("bg"),
          col="viridis",
          trace="none",
          key.title="R-squared")
dev.off()

# Random forest SNPs, ORFs, CNVs baseline model performances
combined <- reshape2::melt(combined)
med_snps <- median(combined[combined$variable=="SNPs","value"]) # median SNPs performance
med_orfs <- median(combined[combined$variable=="ORFs","value"]) # median ORFs performance
med_cnos <- median(combined[combined$variable=="CNVs","value"]) # median CNOs performance
ggplot(combined, aes(x=Conditions, y=value, fill=variable)) + geom_bar(stat="identity") + 
        theme_minimal(11) + facet_grid(variable ~ .) + ylab("Performance R-squared") +
        xlab("Environments") +
        geom_abline(aes(slope=0, intercept=med_snps, col="SNPs"), lty=2) +
        geom_abline(aes(slope=0, intercept=med_orfs, col="ORFs"), lty=2) +
        geom_abline(aes(slope=0, intercept=med_cnos, col="CNVs"), lty=2) + 
        theme(axis.text.x = element_text(angle=45, hjust=1, color="black", face="bold"),
         axis.text.y=element_text(color="black"), plot.margin = unit(c(0.5,0,0,0.5), "cm")) # top, right, bottom, left
ggsave("Scripts/Data_Vis/RF_performances_snps_orfs_bar.pdf", device="pdf", useDingbats=FALSE,
        width=18, height=22, units="cm")

par(mar=c(30,30,30,30)+1) #bottom,left,top,right
pdf("Scripts/Data_Vis/RF_performances_snps_orfs.pdf", height=400, width=400)
heatmap.2(as.matrix(combined), # the values are not matching the resulting image
          cellnote=round(as.matrix(snps_test),2),
          dendrogram="row",
          Colv=F,
          cexRow=36, cexCol=36, # font size
          notecex=36.0, # cell font size
          notecol="black",
          na.color=par("bg"),
          sepcolor="black",
          breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6),
          col=c("white","yellow","#ffae42","orange","#ff4500","red"),
          trace="none",
          margins=c(250,500), # reduce label cutoff
          key.title="R-squared",
          keysize=1.4,
          key.par = list(cex=36), mar=c(2,2,2,0))
graphics.off()

library(ComplexHeatmap)
library(circlize)
at = c(0,0.4,0.8)
col_fun_prop = colorRamp2(at, c("white","yellow","green"))
pdf("Scripts/Data_Vis/RF_performances_snps_orfs.pdf", height=30, width=30)
Heatmap(as.matrix(combined), cluster_columns=F, col=col_fun_prop,
        cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.2f", combined[i, j]), x, y, gp = gpar(fontsize = 36))
        }) + rowAnnotation(foo = anno_text(rownames(combined), gp = gpar(fontsize = 36,fontface="bold"))) +
        Legend(at=c(0,0.8), legend_gp=gpar(fill=col_fun_prop(at)), title="R-squared", legend_height = unit(4, "cm"),
        title_gp=gpar(fontsize="36"), labels_gp = gpar(fontsize = 36,fontface = "bold"))
graphics.off()

#--------- ACCURACY VS NARROW-SENSE HERITABILITY ---------#
accuracy <- read.csv("Scripts/Data_Vis/RF_performances_snps_orfs.csv") # random forest R-squared values
heritability <- read.csv("yeast_rrBLUP_results/SNPs_as_Features/Heritability_h2_H2_sommer.csv") # sommer narrow and broad sense heritability
heritability <- merge(conds, heritability, by.x="cond", by.y="Conditions")
data <- merge(heritability, accuracy, by.x="new_cond", by.y="Conditions")
data$new_cond <- factor(data$new_cond)
ggplot(data, aes(x=SNPs, y=h2, color=new_cond, shape=new_cond)) + 
        scale_shape_manual(values=1:nlevels(data$new_cond)/2) +
        geom_point(size=2)
ggsave("Scripts/Data_Vis/accuracy_vs_h2_sommer.pdf", width=7, height=4,
                device="pdf", useDingbats=FALSE)

#--------- FEATURE SELECTION ---------#
fs <- read.delim("/mnt/scratch/seguraab/yeast_project/yeast_rf_results/RESULTS_reg.txt", sep="\t")
fs <- fs[order(fs$ID),]
adj_r2 <- function(n, p, r2){
        # Adjusted R-squared
        # n is number of instances
        # p is number of features
        # r2 is model performance r-squared 
        return ( 1-(((1-r2)*(n-1))/(n-p-1)) )
}

# 500 to 50k features
for (i in 1:length(cond)){
        print(cond[i])
        print(new_cond[i])
        df <- fs[grep(cond[i],fs$ID),] # subset data corresponding to condition
        ggplot(tmp, aes(x=FeatureNum, y=r2_test)) + geom_line(color="red") + theme_bw() + 
                labs(title=new_cond[i],x="Features", y = "Performance") +
                geom_line(aes(x=FeatureNum, y=r2_val), color="black")
        ggsave(paste("Scripts/Data_Vis/",cond[i],"_FS.pdf", sep=""), width=7, height=4,
                device="pdf", useDingbats=FALSE)
        # adjusted r-squared 
        adj.r2_test <- adj_r2(125, df$FeatureNum, df$r2_test)
        adj.r2_val <- adj_r2(625, df$FeatureNum, df$r2_val)
        ggplot(df, aes(x=FeatureNum, y=adj.r2_test)) + geom_line(color="red") + theme_bw() + 
                labs(title=new_cond[i],x="Features", y = "Adjusted R-squared") +
                geom_line(aes(x=FeatureNum, y=adj.r2_val), color="black")
        ggsave(paste("Scripts/Data_Vis/",cond[i],"_FS_adj_r.pdf", sep=""), width=7, height=4,
                device="pdf", useDingbats=FALSE)
}

# 2 to 2064 features
sub <- fs[grep("_exp_", fs$ID),]
for (i in 1:length(cond)){
        print(cond[i])
        print(new_cond[i])
        sub_df <- sub[grep(cond[i],sub$ID),]
        ggplot(sub_df, aes(x=FeatureNum, y=r2_test)) + geom_line(color="red") + theme_bw() + 
                labs(title=new_cond[i],x="Features", y = "Performance") +
                geom_line(aes(x=FeatureNum, y=r2_val), color="black")
        ggsave(paste("Scripts/Data_Vis/",cond[i],"_exp_FS.pdf", sep=""), width=7, height=4,
                device="pdf", useDingbats=FALSE)
        # adjusted r-squared 
        adj.r2_test <- adj_r2(125, sub_df$FeatureNum, sub_df$r2_test)
        adj.r2_val <- adj_r2(625, sub_df$FeatureNum, sub_df$r2_val)
        ggplot(sub_df, aes(x=FeatureNum, y=adj.r2_test)) + geom_line(color="red") + theme_bw() + 
                labs(title=new_cond[i],x="Features", y = "Adjusted R-squared") +
                geom_line(aes(x=FeatureNum, y=adj.r2_val), color="black")
        ggsave(paste("Scripts/Data_Vis/",cond[i],"_exp_FS_adj_r.pdf", sep=""), width=7, height=4, 
                device="pdf", useDingbats=FALSE)
}

#--------- SNP IMPORTANCE VARIATION ACROSS CONDITIONS ---------#
# read importance score files

# rank SNPs



#--------- PATHWAY INFORMATION ----------#
pwys <- read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/biallelic_snps_diploid_S288C_genes_pwys_unique_descriptions.csv")
snps <- read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_pwys.csv", header=TRUE)
out <- merge(snps, pwys, by.x="V5", by.y="Pathway.ID")
write.csv(out, "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes_pwys_descriptions.csv", quote=F, row.names=F)

# ypdbenomyl500
sshap <- read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/SHAP_top_snps_genes_sorted_YPDBENOMYL500_training.txt", sep="\t")
sshap <- merge(sshap, out, by.x="X", by.y="V1")
write.csv(sshap, "Scripts/Data_Vis/00_shap_ypdbenomyl500_pathways.csv", quote=F, row.names=F)

oshap <- read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Data_Vis/SHAP_values_sorted_average_YPDBENOMYL200_orf_training.txt", sep="\t")
cshap

# ypdcafein40
sshap <- read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/SHAP_top_snps_genes_sorted_YPDCAFEIN40_training.txt", sep="\t")
sshap <- merge(sshap, out, by.x="X", by.y="V1")
write.csv(sshap, "Scripts/Data_Vis/00_shap_ypdcafein40_pathways.csv", quote=F, row.names=F)

# ypdcafein50
sshap <- read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_RF/SHAP/SHAP_top_snps_genes_sorted_YPDCAFEIN50_training.txt", sep="\t")
sshap <- merge(sshap, out, by.x="X", by.y="V1")
write.csv(sshap, "Scripts/Data_Vis/00_shap_ypdcafein50_pathways.csv", quote=F, row.names=F)



#------------ GENE COPY NUMBER DISTRIBUTION

cno <- read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs_no_NA.csv")
orf <- read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs_pres_abs.csv")

cno <- melt(cno)
orf <- melt(orf)

ggplot(cno, aes(x=value)) + geom_histogram(bins=40) + 
        labs(x="ORF Copy Number", y="Counts")
ggsave("Scripts/Data_Vis/orf_copy_number_dist.pdf", device="pdf", useDingbats=FALSE)

ggplot(orf, aes(x=value)) + geom_histogram(bins=40) + 
        labs(x="ORF presence/absence", y="Counts")
ggsave("Scripts/Data_Vis/orf_pres_abs_dist.pdf", device="pdf", useDingbats=FALSE)
