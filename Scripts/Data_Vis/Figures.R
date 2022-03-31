# Chapter 1 Figures
library(data.table)
library(gplots)
library(ggplot2)
library(viridis)
library(MASS)
library(dplyr)

setwd("/mnt/home/seguraab/Shiu_Lab/Project")

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

pdf("Scripts/Data_Vis/kinship.pdf")
hm <- heatmap.2(kin,trace="none", col = colorRampPalette(c("blue","white","red"))(21),
        dendrogram="both",cex.axis=0.4,notecex=0.4,cexRow=0.5,cexCol=0.4,Rowv = TRUE,
        Colv=TRUE,key.title="Kinship") #,lmat=lmat,lwid=lwid,lhei=lhei)
dev.off()

#--------- FITNESS CORRELATIONS ---------#
pCorEnvs <- read.csv("Data/Peter_2018/pheno_corr_envs.csv", header=T, row.names=1) # Correlation between conditions
pCorIso <- read.csv("Data/Peter_2018/pheno_corr_isolates.csv", header=T, row.names=1) # Correlation between isolates across all conditions

pdf("Scripts/Data_Vis/pheno_corr_isolates.pdf")
hm2 <- heatmap.2(as.matrix(pCorIso),trace="none",col = colorRampPalette(c("blue","white","red"))(21),
        dendrogram="both",cex.axis=0.4,notecex=0.4,cexRow=0.5,cexCol=0.4,key.title="PCC",
        Rowv=hm$rowDendrogram,Colv=hm$rowDendrogram)# ,lmat=lmat,lwid=lwid,lhei=lhei)
dev.off()

pdf("Scripts/Data_Vis/pheno_corr_envs.pdf")
hm3 <- heatmap.2(as.matrix(pCorEnvs),trace="none",col=colorRampPalette(c("blue","white","red"))(21),
        dendrogram="both",cex.axis=0.4,notecex=0.4,cexRow=0.5,cexCol=0.4,key.title="PCC")
dev.off()
#--------- ISOLATE FITNESS CORRELATION VS KINSHIP ---------#
k_bin <- reshape2::melt(kin) # reshape kinship matrix
k_quant <- quantile(k_bin, seq(0,1,0.05)) # kinship quantiles
k_bin$bin <- cut(k_bin$value, breaks=k_quant) # bins
k_bin$bin <- as.character(lapply(strsplit(as.character(k_bin$bin),split=","),head,n=1)) # rename bins
k_bin$bin <- gsub("\\(", "", k_bin$bin) # remove (
k_bin$bin <- as.numeric(k_bin$bin) # convert to numeric

pCorIso <- as.matrix(pCorIso)
i_bin <- reshape2::melt(pCorIso) # reshape fitness across isolates correlation matrix
i_quant <- quantile(i_bin$value, seq(0,1,0.05)) # pCor quantiles
i_bin$bin <- cut(i_bin$value, breaks=i_quant)
i_bin$bin <- as.character(lapply(strsplit(as.character(i_bin$bin),split=","),head,n=1))
i_bin$bin <- gsub("\\(", "", i_bin$bin)
i_bin$bin <- as.numeric(i_bin$bin)

k_i_bin <- merge(k_bin, i_bin, by=c("Var1", "Var2")) # merge dataframes
names(k_i_bin) <- c("Isolate 1", "Isolate 2", "Kinship", "Kinship bin", "pCor", "pCor bin")

k.i.cor <- cor(k_i_bin$Kinship, k_i_bin$pCor, method="spearman") # Spearman correlation
k_i_bin_count <- aggregate(cbind(count = `Isolate 2`) ~ `Kinship bin`, # bin counts
        data=k_i_bin, FUN=function(x){NROW(x)})
k_i_bin_median <- k_i_bin %>% group_by(`Kinship bin`) %>% summarise(median=median(pCor)) # bin pCor median

get_density <- function(x, y, ...) { # density
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])}
k_i_bin$density <- get_density(k_i_bin$Kinship, k_i_bin$pCor, n=100)

# 5% and 95% quantiles for each kinship bin
quants <- as.data.frame(do.call("rbind", tapply(k_i_bin$pCor, k_i_bin$`Kinship bin`, quantile, c(0.05, 0.95))))
names(quants) <- c("x5","x95")
quants$`Kinship bin` <- as.numeric(row.names(quants))

# Density scatter plot
ggplot(k_i_bin, aes(x=Kinship, y=pCor, color=density)) + geom_point(size=0.2, alpha=0.2) +
        geom_smooth(method="lm", formula=y~x, se=FALSE) +
        scale_color_viridis(direction=-1) + theme_bw(8)
ggsave("Scripts/Data_Vis/kinpCorIso_Density.pdf", width=5, height=4, device="pdf", useDingbats=FALSE)

ggplot(k_i_bin, aes(x=Kinship, y=pCor, color=density)) + geom_point(size=0.2, alpha=0.2) +
        scale_color_viridis(direction=-1) + theme_bw(8) +
        geom_line(data=k_i_bin_median, aes(x=`Kinship bin`, y=median, group=1), color="red") +
        annotate(geom="text", x=3, y=-0.3, label=paste("Spearman's rank cor = ", round(k.i.cor, digits=2),sep=""))
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

#--------- HERITABILITY ---------#


#--------- GENOMIC PREDICTION ACCURACY ---------#


#--------- ACCURACY VS HERITABILITY ---------#


#--------- FEATURE SELECTION ---------#
fs <- read.delim("/mnt/scratch/seguraab/yeast_project/yeast_rf_results/RESULTS_reg.txt", sep="\t")
fs <- fs[order(fs$ID),]
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
