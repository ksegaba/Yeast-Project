# Calculate sliding window variances
library(BGLR)
library(data.table)
library(tidyr)

X <- fread("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv") # genotype
X <- as.matrix(X, rownames=1, colnames=1)
X <- scale(X, center=TRUE, scale=FALSE)
y <- read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv", row.names=1) # phenotype

# Get SNP positions
snp_ids <- data.frame(SNPs=colnames(X))
snp_ids <- snp_ids %>% separate(SNPs, c("Chr", "bp"), remove=FALSE)
snp_ids$bp <- as.numeric(snp_ids$bp)
snp_ids <- snp_ids[order(snp_ids$bp),] # order ascending

setwd("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_BayesC/Results")

B <- readBinMat("YPDANISO20_bayesC_mapping_12000_ETA_1_b.bin")
colnames(B) <- colnames(X) 

########## Per Chromosome ##########
dir <- "/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_BayesC/Results_chr1"
####### Non-overlapping windows
#calc_wv <- function(chr="chromosome1", path=dir, title="Chromosome 1", wsize=10000){
#    setwd(path)
#    ch <- snp_ids[which(snp_ids$Chr == chr),] # SNPs in chromosome 1
#    startbp <- min(ch$bp)
#    endbp <- max(ch$bp)
#    nW <- ceiling(endbp/wsize) # number of windows (size 10kb)
#    res <- data.frame(Chr=1, Chr_start=startbp, Chr_end=endbp, 
#                    Window=1:nW, W_size=wsize, W_start=NA, W_end=NA,
#                    n_SNPs=NA)
#    yHats <- vector(mode="list", length=nrow(B))
#    wv <- c()
#    for (i in 1:nW){
#        sprintf("Window %i", i)
#        # window i start and end bp
#        wendbp <- startbp + wsize * i
#        res$W_start[i] <- wendbp - wsize
#        res$W_end[i] <- wendbp - 1

#        # get snps in window i
#        cols <- snp_ids[which(snp_ids$bp <= wendbp - 1 & snp_ids$bp >= wendbp - 10000),] # snps in window
#        X2 <- X[,which(colnames(X) %in% cols$SNPs)] # subset genotype matrix
#        res$n_SNPs[i] <- dim(X2)[2]

#        if (dim(X2)[2] != 0){
#            tmp <- c()
#            for (j in 1:nrow(B)){ # model iterations
#                # Window variances across SNPs per bayesc iteration
#                g <- X2 %*% B[j,colnames(X2)] # local genotypic value for each individual
#                tmp[j] <- var(g)
#                yHats[[j]] <- cbind(yHats[[j]],g)
#            }
#            pdf(paste(i,"_window_trace.pdf", sep=""))
#            plot(tmp, main=paste("Window ", i, sep=""))
#            dev.off()
#            wv <- cbind(wv, tmp)
#        }
#    }
#    pdf(paste(chr, "_manhattan_plot.pdf", sep=""))
#    plot(colMeans(wv), main=title)
#    dev.off()

#    # Variance of genotypic values across individuals in each iteration
#    vars <- c()
#    for (i in 1:nrow(B)){
#        vars[i] <- var(rowSums(yHats[[i]]))
#    }
#    pdf(paste(chr, "_manhattan_plot2.pdf", sep=""))
#    plot(vars, main=title)
#    dev.off()
#    return(wv, yHats)
#}

#calc_wv()


####### Sliding overlapping windows
chr="chromosome1"; path=dir; title="Chromosome 1"; wsize=1000
#calc_swv <- function(){
    setwd(path)
    ch <- snp_ids[which(snp_ids$Chr == chr),] # SNPs in chromosome 1
    startbp <- min(ch$bp)
    endbp <- max(ch$bp)
    nW <- length(seq(startbp, endbp, 500)) # number of windows (size 1kb, 100bp shift)
    res <- data.frame(Chr=1, Chr_start=startbp, Chr_end=endbp, Window=1:nW, W_size=wsize, 
                    W_start=seq(startbp, endbp, 500), W_end=NA, n_SNPs=NA)
    res$W_end <- res$W_start + wsize
    yHats <- vector(mode="list", length=nrow(B))
    wv <- c()
    for (i in 1:nW){
        sprintf(paste("Window ", i, sep=""))
        # get snps in window i
        cols <- snp_ids[which(snp_ids$bp <= res$W_end[i] & snp_ids$bp >= res$W_start[i]),] # snps in window
        X2 <- X[,which(colnames(X) %in% cols$SNPs)] # subset genotype matrix
        res$n_SNPs[i] <- dim(X2)[2]
        
        if (dim(X2)[2] != 0){
            tmp <- c()
            for (j in 1:nrow(B)){ # model iterations
                # Window variances across SNPs per bayesc iteration
                g <- X2 %*% B[j,colnames(X2)] # local genotypic value for each individual
                tmp[j] <- var(g)
                yHats[[j]] <- cbind(yHats[[j]],g)
            }
            #pdf(paste(i, wsize, "window_trace.pdf", sep="_"))
            #plot(tmp, main=paste("Window ", i, sep=""))
            #dev.off()
            wv <- cbind(wv, tmp)
        }
    }
    pdf(paste(chr, wsize, "manhattan_plot.pdf", sep="_"))
    plot(colMeans(wv), main=title)
    dev.off()

    # Variance of genotypic values across individuals in each iteration
    vars <- c()
    for (i in 1:nrow(B)){
        vars[i] <- var(rowSums(yHats[[i]]))
    }
    pdf(paste(chr, wsize, "manhattan_plot2.pdf", sep="_"))
    plot(vars, main=title)
    dev.off()
    #return(wv, yHats)
#}
write.csv(wv, "WV_bayesc_2kb_200slide.csv")
write.csv(res, "res_bayesc_2kb_200slide.csv")

#calc_swv(chr="chromosome1", path=dir, title="Chromosome 1", wsize=1000)


########## All Chromosomes together ##########
#setwd("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_BayesC/Results_1kb")
#n=1
#WV=data.frame(matrix(nrow=nrow(B), ncol=5*length(seq(1000, max(snp_ids$bp)+1000, 1000)))) # collect window variances
#columns=c()
#for (i in seq(1000, max(snp_ids$bp)+1000, 1000)){ 
#    cols <- snp_ids[which(snp_ids$bp <= i & snp_ids$bp > i-1000),] # snps to include
#    X2 <- X[,which(colnames(X) %in% cols$SNPs)] # subset genotype matrix
#    print(dim(X2))
#    columns[n] <- paste(i-1000, i, sep="-")

#    if (dim(X2)[2] != 0){
#        print(paste(i-1000, i, sep=" - "))
#        for (k in 1:nrow(B)){
#            # Window variances across SNPs per bayesc iteration
#            g <- X2 %*% B[k,colnames(X2)]
#            WV[k,n] <- var(g) 
#        }
        # trace plot of window variances
#        print("Plotting...")
#        pdf(paste(i-1000, i, "bayesC_YPDANISO20_wv_5kb.pdf", sep="_"))
#        plot(x=1:2000, y=WV[,n], main="Window variances (size 1000)")
#        dev.off()
#        print(mean(WV[,n]))
#        print(head(WV[,n]))
#        # what we want is a manhattan plot
#    }
#    n = n+1
#}

#write.csv(WV, file="WV_bayesc_5kb.csv", quote=FALSE)

#winvar <- read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_BayesC/Results/WV_bayesc_5kb.csv", row.names=1)
#pdf("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_BayesC/WV_bayesc_5kb_iter1.pdf")
#plot(winvar[1,], main="Iteration 1 WV (size 5kb)")
#dev.off()

#pdf("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_BayesC/Results/posterior_means.pdf")
#plot(colMeans(winvar))
#dev.off()