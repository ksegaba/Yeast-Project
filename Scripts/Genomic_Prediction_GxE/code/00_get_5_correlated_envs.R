# Get the top 4 correlated environments for each environment
setwd("/mnt/home/seguraab/Shiu_Lab/Project")
Y <- read.csv("Data/Peter_2018/pheno.csv", row.names=1) # phenotype matrix
Y_cor <- read.csv("Scripts/Data_Vis/Pheno_Figures/pheno_pairs_cor.csv") # phenotype matrix correlations

file.create("Scripts/Genomic_Prediction_GxE/data/5_env_fitness.txt")
file.create("Scripts/Genomic_Prediction_GxE/data/5_env_fitness_corr.txt")
write(paste("Target", "Env1", "Env2", "Env3", "Env4", "Target_PCC", "PCC_Env1",
            "PCC_Env2", "PCC_Env3", "PCC_Env4", "PCC_av", "PCC_sd", sep=","),
      file="Scripts/Genomic_Prediction_GxE/data/5_env_fitness_corr.txt", append=TRUE)
for (trait_name in colnames(Y)){
    index <- which(Y_cor$p.value < 0.05 & Y_cor$Env1==trait_name | Y_cor$Env2==trait_name)
    if (length(index) >= 4){
        y_cor <- Y_cor[index,] # subset correlated traits
        y_cor <- y_cor[order(y_cor$PCC, decreasing=T),] # sort
        sub <- y_cor[1:4,] # top 4 correlated traits
        envs <- unique(c(sub$Env1, sub$Env2)) # correlated traits
        envs <- envs[which(envs!=trait_name)] # names of correlated traits only
        write(paste(c(trait_name, envs), collapse=","),
              file="Scripts/Genomic_Prediction_GxE/data/5_env_fitness.txt",
              append=TRUE)
        pcc_list <- c(1)
        for (env in envs) pcc_list <- c(pcc_list, sub[sub$Env1==env | sub$Env2==env,]$PCC)
        line <- paste(c(paste(c(trait_name, envs), collapse=","),
                        paste(c(pcc_list, mean(sub$PCC), sd(sub$PCC)), collapse=",")),
                      collapse=",")
        write(line, file="Scripts/Genomic_Prediction_GxE/data/5_env_fitness_corr.txt",
              append=TRUE)
    }
}