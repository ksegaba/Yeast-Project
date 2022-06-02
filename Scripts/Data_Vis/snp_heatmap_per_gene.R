# Visualize SNP distribution heatmaps across diploid yeast isolates and across
# environments for each gene.
# Input:
# - genes mapped to SNPs file from match_SNPs_to_S288C.py
# - genotype file

library(data.table)
library(circlize)
library(ComplexHeatmap)

# Read in data
dir <- "/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/"
dir2 <- "/mnt/home/seguraab/Shiu_Lab/Project/Data/S288C_reference_genome_R64-2-1_20150113/"
df <- read.csv(paste(dir,"biallelic_snps_diploid_and_S288C_genes.txt", sep=""), header=F) # snps mapped to genes
geno <- fread(paste(dir, "geno.csv", sep="")) # snps genotype
biop <- read.csv(paste(dir2, "All_genes_and_pathways_in_S._cerevisiae_S288c.txt", sep=""), header=T, sep="\t") # gene pathway information

t_geno <- transpose(as.data.frame(geno)) # T=transpose data
colnames(t_geno) <- as.character(t_geno[1,]) # set column names as isolate IDs
t_geno <- t_geno[-1,] # drop first row containing isolate IDs
rownames(t_geno) <- colnames(geno)[2:length(geno)] # set row names as snp IDs

# Add pathway information
#meta <- merge(df, biop, by.x="V4", by.y="Genes", all.x=T)
#length(which(meta$V4=="intergenic")) # number of "intergenic" snps (5726)
#meta <- subset(meta, V4!="intergenic") # only genic snps (58730)
#meta <- meta[order(meta$V1),] # sort by snp name

# Exclude intergenic snps (5726)
t_geno <- subset(t_geno, df$V4 != "intergenic") # only genic snps (58730)
df <- subset(df, V4 != "intergenic")
rownames(df) <- df$V1 # set snp IDs as row names
df <- df[,4, drop=F] # drop extra columns

# subsets
geno_sub <- t_geno[1:200,]
df_sub <- df[1:200,, drop=F]

# Heatmap
colors = structure(c("blue", "white", "red"), names=c("-1", "0", "1"))
ann = HeatmapAnnotation(type = df$V4, annotation_name_side = "left" ) # gene annotations
pdf(paste(dir, "snp_genes_vs_isolates.pdf", sep=""), height=100, width = 50)
h = Heatmap(as.matrix(t_geno), name = "genotype", 
            column_title = "Yeast gene genotypes across isolates", 
            column_title_gp = gpar(fontsize = 36),
            col = colors, show_column_names = F, show_row_names = F,
            cluster_col = "hclust", show_column_dend = T) 
h2 = Heatmap(as.matrix(df), name = "genes", show_row_names = F, cluster_row = "hclust",
            width = 20, show_row_dend = T, row_dend_side = "left", row_dend_width = 20)
    #rowAnnotation(width = unit(1, "cm"))
draw(h+h2, main_heatmap = "genes")
dev.off()

