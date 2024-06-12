library(data.table)
library(BGLR)
setwd("/mnt/home/seguraab/Shiu_Lab/Project")

envs <- c("YPACETATE", "YPD14", "YPD40", "YPD42", "YPD6AU",
		   "YPDANISO10", "YPDANISO20", "YPDANISO50", "YPDBENOMYL200",
		   "YPDBENOMYL500", "YPDCAFEIN40", "YPDCAFEIN50", "YPDCHX05",
		   "YPDCHX1", "YPDCUSO410MM", "YPDDMSO", "YPDETOH",
		   "YPDFLUCONAZOLE", "YPDFORMAMIDE4", "YPDFORMAMIDE5", "YPDHU",
		   "YPDKCL2M", "YPDLICL250MM", "YPDMV", "YPDNACL15M", "YPXYLOSE",
		   "YPDNACL1M", "YPDNYSTATIN", "YPDSDS", "YPDSODIUMMETAARSENITE",
		   "YPETHANOL", "YPGALACTOSE", "YPRIBOSE", "YPGLYCEROL", "YPSORBITOL")

################################################################################
### TABLE S3 perhaps
################################################################################
############### Multi-trait GxE CV1: changes in env coefficients ###############
geno <- fread("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/ORFs_pres_abs.csv")
d <- "/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_GxE_results"
for (env in envs){
	files = list.files(path=paste(d, "/multi-trait/cv1/", sep=""), 
		pattern=paste("GxE_cv1_pav_5_envs_", env, "_.*.RDS", sep=""))

	# Get feature coefficients
	model <- files[1]
	mod <- readRDS(paste("/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_GxE_results/multi-trait/cv1/",
		model, sep="")) # multi-trait model
	coefs <- as.data.frame(mod$ETA[[1]]$beta) # feature coefficients
	rownames(coefs) <- colnames(geno)[2:ncol(geno)]
	colnames(coefs) <- strsplit(model, "_")[[1]][6:10]
	coefs <- reshape2::melt(coefs)
	for (model in files[2:length(files)]){ # coefficients from other model repetitions
		mod <- readRDS(paste("/mnt/gs21/scratch/seguraab/yeast_project/ORF_yeast_GxE_results/multi-trait/cv1/",
			model, sep="")) # multi-trait model
		coefs2 <- as.data.frame(mod$ETA[[1]]$beta) # feature coefficients
		colnames(coefs2) <- strsplit(model, "_")[[1]][6:10]
		coefs2$Feature <- colnames(geno)[2:ncol(geno)]
		coefs2 <- reshape2::melt(coefs2, id="Feature")
		coefs <- cbind(coefs, coefs2[,3])
	}

	# Calculate summary statistics
	colnames(coefs) <- c("Env", seq(1, ncol(coefs)-1, 1))
	coefs$mean <- rowMeans(coefs[,-1])
	coefs$sd <- genefilter::rowSds(coefs[,c(2:11)])
	coefs$var <- matrixStats::rowVars(as.matrix(coefs[,c(2:11)]))
	coefs$Feature <- coefs2$Feature
	
	# Sort by average feature importance and save
	coefs <- coefs[order(coefs$mean, decreasing=T),c("Feature", "Env", "mean", "sd", "var")]
	write.csv(coefs, sprintf("Scripts/Data_Vis/Section_4/Coefs_multi-trait_%s.csv", env), quote=F, row.names=F)
	
	# Create a plot for the top 5% features to show changes in feature importance
	tsh <- quantile(coefs$mean, 0.95)[[1]] # top 5%
	top5p <- coefs[coefs$mean >= tsh,c("Feature", "Env", "mean", "sd", "var")]
	head(tidyr::pivot_wider(top5p, id_cols="Feature", values_from="mean", names_from="Env"))
}