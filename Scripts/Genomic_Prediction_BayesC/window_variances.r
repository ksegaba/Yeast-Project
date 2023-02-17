# Calculate sliding window variances of snp effect sizes across iterations
library(BGLR)
library(data.table)
library(tidyr)
library(parallel)
library(stringr)


X <- fread("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv") # genotypes
X <- as.matrix(X, rownames=1, colnames=1)
X <- scale(X, center=TRUE, scale=FALSE)
y <- read.csv("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv", row.names=1) # phenotype


### BayesC
# ETA <- list(list(X=Xs, model="BayesC", count=110, probIn=1/100, saveEffects=TRUE))
# bayesc <- function(Y){ # Function to run BayesC on each trait
#     print(Y)
#     fit <- BGLR(y=as.matrix(y[Y]), ETA=ETA, nIter=12000, burnIn=2000, saveAt=paste(Y, "bayesC_mapping_12000_", sep="_"), verbose=FALSE)

#     # plot posterior probabilities
#     pdf(paste(Y, "bayesC_mapping_posterior.pdf", sep="_"))
#     plot(fit$ETA[[1]]$d, main=Y)
#     dev.off()

#     # save model
#     saveRDS(fit, paste(Y, "_bayesC_model_12000.rds", sep=""))
# }
# system.time(tmp <- mclapply(colnames(y), bayesc, mc.cores=35))


# Get SNP chromosomal positions
snp_ids <- data.frame(SNPs=colnames(X))
snp_ids <- snp_ids %>% separate(SNPs, c("Chr", "bp"), remove=FALSE)
snp_ids$bp <- as.numeric(snp_ids$bp)
snp_ids$Chr <- as.numeric(gsub("chromosome", "", snp_ids$Chr))
snp_ids <- snp_ids[order(snp_ids$Chr),] # order by chromosome

setwd("/mnt/home/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_BayesC")
B <- readBinMat("Results/YPDCAFEIN40_bayesC_mapping_12000_ETA_1_b.bin") # snp effects
colnames(B) <- colnames(X)

########## Per Chromosome ##########
####### Sliding overlapping windows
calc_swv <- function(chr){
    # create output folder
    path <- file.path(getwd(), paste("Results_YPDCAFEIN40_chr", chr, sep=""))
    dir.create(path)
    setwd(path)

    # subset snps within chromosome
    ch <- snp_ids[which(snp_ids$Chr == chr),]
    
    # get windows
    startbp <- min(ch$bp)
    endbp <- max(ch$bp)
    nW <- length(seq(startbp, endbp, 100)) # number of windows (size 1kb, 100bp shift)
    res <- data.frame(Chr=chr, Chr_start=startbp, Chr_end=endbp, Window=1:nW, W_size=1000,
                    W_start=seq(startbp, endbp, 100), W_end=NA, n_SNPs=NA, WV_mean=NA,
                    WV_sd=NA, WV_min=NA, WV_max=NA)
    res$W_end <- res$W_start + 1000
    
    # calculate window variances
    for (i in 1:nW){
        sprintf(paste("Window ", i, sep=""))
        
        # get snps in window i
        cols <- ch[which(ch$bp <= res$W_end[i] & ch$bp >= res$W_start[i]),] # snps in window
        res$n_SNPs[i] <- nrow(cols)
        #X2 <- X[,which(colnames(X) %in% cols$SNPs)] # subset genotype matrix
        
        # calculate genotypic values
        if (nrow(cols)>=2){ # at least two snps in window
            tmp <- c() # collect window variances in each iteration
            for (j in 1:nrow(B)){ # model iterations
                # Window variances across SNPs per model iteration
                g <- X2 %*% B[j,cols$SNPs] # local genotypic value for each individual
                tmp[j] <- var(g)
            }
            res$WV_mean[i] <- mean(unlist(tmp))
            res$WV_sd[i] <- sd(unlist(tmp))
            res$WV_min[i] <- min(unlist(tmp))
            res$WV_max[i] <- max(unlist(tmp))
            pdf(paste(i, res$W_start[i], res$W_end[i], "window_trace.pdf", sep="_"))
            plot(unlist(tmp), main=paste("Window ", i, sep=""), xlab="nIter", ylab="genomic variances")
            dev.off()
            #write.csv(unlist(tmp), "window_variances_per_iteration.csv")
        }
    }
    # mean window variances across all 2000 iterations
    pdf("wv_mean_manhattan_plot.pdf")
    plot(res$WV_mean, main=paste("Chromosome", chr), xlab="nW", ylab="mean WV across all 2000 iterations")
    dev.off()
    # number of snps per window
    pdf("n_snps_per_window_dist.pdf")
    hist(as.numeric(res$n_SNPs), breaks=10, xlab="number of SNPs")
    dev.off()
    # save results table
    write.csv(res, "res_1kb_500bp_shift.csv", quote=F, row.names=F)
    setwd("/mnt/ufs18/home-056/seguraab/Shiu_Lab/Project/Scripts/Genomic_Prediction_BayesC")
}
mclapply(1:12, calc_swv, mc.cores=12)
#lapply(1:12, calc_swv)


####### Calculate chromosome-level LD
X <- fread("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv") # genotypes
X <- as.matrix(X, rownames=1, colnames=1)
map <- read.csv('/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes.txt') # snp metadata
chrs <- unique(map$chr)
info <- data.frame(chr=1:12, n_snps=0) # snps per chromosome

# ensure map snps are in the same order as in cols of X
map <- map[order(match(map$snp,colnames(X))),]

# get snps in chromosome
se <- function(x) {return(sd(x)/sqrt(length(x)))} # standard error
calc_LD <- function(ch, tmp, half="whole"){
    R2=cor(X[,which(tmp)])^2 # snp genotype correlations
    D=as.matrix(dist(map$pos[which(tmp)])) # distance between snps
    dim(D)
    dim(R2)
    tmp2=col(D)>row(D) # get upper triangular matrix
    LD=data.frame(dist=D[tmp2],R2=R2[tmp2]) # linkage disequilibrium
    if (dim(LD)[1]==sum(tmp)*(sum(tmp)-1)/2) print("We're good!")

    tmp3=lowess(y=LD$R2,x=LD$dist) # locally weighted scatterplot smoothing
    range(LD$dist)
    range(LD$dist)/1e6
    LD$range=cut(LD$dist,breaks=seq(from=0,to=range(LD$dist)[2]+1,by=10000)) # distance bins
    table(LD$range)

    print("Generating plots...")
    R2_med <- tapply(FUN=median,X=LD$R2,INDEX=LD$range) # median R2
    dist_med <- tapply(FUN=median,X=LD$dist,INDEX=LD$range) # median distance
    pdf(paste(ch, half, "LD_median.pdf", sep="_")) # Plot LD
    plot(R2_med~dist_med,type='o')
    dev.off()

    N <- tapply(FUN=length,X=LD$dist,INDEX=LD$range) # bin counts
    R2_av <- tapply(FUN=mean,X=LD$R2,INDEX=LD$range) # mean R2
    dist_av <- tapply(FUN=mean,X=LD$dist,INDEX=LD$range) # mean distance
    pdf(paste(ch, half, "LD_mean.pdf", sep="_")) # Plot LD
    plot(R2_av~dist_av,type='o')
    dev.off()

    SD <- tapply(FUN=sd,X=LD$R2,INDEX=LD$range) # standard deviation
    summary(SD) ; summary(R2_av)
    pdf(paste(ch, half, "LD_mean_sd.pdf", sep="_"))
    plot(R2_av~dist_av,type='o')#,ylim=c(0, max(SD)))
    arrows(x0=dist_av, y0=R2_av-SD, x1=dist_av, y1=R2_av+SD, col="blue")
    dev.off()

    SE <- tapply(FUN=se,X=LD$R2,INDEX=LD$range) # standard deviation
    summary(SE) ; summary(R2_av)
    pdf(paste(ch, half, "LD_mean_se.pdf", sep="_"))
    plot(R2_av~dist_av,type='o')#,ylim=c(0, max(SD)))
    arrows(x0=dist_av, y0=R2_av-SE, x1=dist_av, y1=R2_av+SE, col="blue")
    dev.off()
}

apply_calc_LD <- function(ch){
    print(ch)
    tmp <- map$chr==ch
    index <- as.numeric(str_extract(ch, "[0-9]+$"))
    #info$n_snps[index] <- length(which(tmp)) # no. of snps

    calc_LD(ch, tmp) # whole chromosome
    
    # divide chromosome into 2 halves
    minimum <- min(map[map$chr==ch,]$pos)
    maximum <- max(map[map$chr==ch,]$pos)
    half <- round((maximum-minimum)/2, 0)
    tmpa <- (map$chr==ch & map$pos <= half) # first half of chromosome
    tmpb <- (map$chr==ch & map$pos > half) # second half of chromosome
    calc_LD(ch, tmpa, "a")
    calc_LD(ch, tmpb, "b")

    return(c(ch, index, length(which(tmp))))
}
info <- mclapply(chrs, apply_calc_LD, mc.cores=12)
info <- unlist(info)
info <- data.frame(matrix(info, nrow=2, byrow=T)
write.csv(t(info), "num_snps_per_chr.csv", quote=F)

# Singular value decomposition of genotype matrix
X <- scale(X,center=TRUE,scale=FALSE)
SVD <- svd(X,nu=3)
dim(SVD$u)
dim(X)
pdf("SVD_geno.pdf")
plot(SVD$u[,1:2])
dev.off()


########## Chromosome-sized windows ##########
B <- readBinMat("Results/YPDCAFEIN40_bayesC_mapping_12000_ETA_1_b.bin") # snp effects
colnames(B) <- colnames(X)
map <- read.csv('/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/biallelic_snps_diploid_and_S288C_genes.txt') # snp metadata
chrs <- unique(map$chr)
calc_chrwv <- function(ch){
    # create output folder
    path <- file.path(getwd(), paste("Results_YPDCAFEIN40_", ch, sep=""))
    if (dir.exists(path)==FALSE) dir.create(path)

    # subset snps within chromosome
    chrom <- map[which(map$chr == ch),]
    
    # calculate window variances
    X2 <- X[,which(colnames(X) %in% chrom$snp)] # subset genotype matrix
    tmp <- c() # collect window variances in each model iteration
    for (j in 1:nrow(B)){ # model iterations
        # Window variances across SNPs per model iteration
        g <- X2 %*% B[j,chrom$snp] # local genotypic value for each individual
        tmp[j] <- var(g) # genomic variances across individuals
    }

    # generate genomic variance plot
    pdf(file.path(path, paste(ch, "window_trace.pdf", sep="_")))
    plot(unlist(tmp), main=ch, xlab="nIter", ylab="genomic variances")
    abline(a=mean(unlist(tmp)), b=0, col="blue") # mean window variance across the 2000 iterations
    dev.off()
    #write.csv(unlist(tmp), file.path(path, paste(ch,"as_window_var_per_iter_YPDCAFEIN40.csv", sep="_")), quote=F, row.names=F)

    # out <- c(min(chrom$pos), max(ch$pos), mean(unlist(tmp)), sd(unlist(tmp)), min(unlist(tmp)), max(unlist(tmp)))
    print(c(min(chrom$pos), max(chrom$pos), mean(unlist(tmp)), sd(unlist(tmp)), min(unlist(tmp)), max(unlist(tmp))))
    return(tmp) # out
}
res <- lapply(chrs, calc_chrwv)
res <- unlist(res)
res <- matrix(res, nrow=12, byrow=T)
rownames(res) <- chrs
colnames(res) <- 1:2000
#colnames(res) <- c('Chr', 'Chr_start', 'Chr_end', 'WV_mean', 'WV_sd', 'WV_min', 'WV_max') # for 'out'
write.csv(t(res), "chr_as_windows_geno_var.csv", quote=F)

# 12x12 covariance matrix of genotypic value across model iterations
res <- read.csv("chr_as_windows_geno_var.csv", header=T, row.names=1)
res.cov <- cov(res)
summary(res.cov)
A <- matrix(0, nrow=12, ncol=12)
A <- (A + res.cov)/2000 # posterior means
write.csv(A, "chr_as_windows_post_means.csv", quote=F)

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