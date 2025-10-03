# Linkage Disequilibrium plot using output from TASSEL5
# Source: https://github.com/mohsinali1990/My_scripts/blob/main/LD%20decay%20Plot%20from%20TASSEL%20LDoutput.R
library(data.table)
library(ggplot2)
library(dplyr)

# remove missing values from R^2 and Distance variables
setwd("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018")
ld <- fread('geno_LD_window_50.txt')
# ld <- read.delim('geno_LD_window_50.txt', stringsAsFactors=FALSE, header=TRUE, sep="\t")
# ld100 <- read.delim('geno_LD_window_100.txt', stringsAsFactors=FALSE, header=TRUE, sep="\t")
ld <- ld[ld$R^2 != "NaN",]
# ld100 <- ld100[ld100$R^2 != "NaN",]
ld$Dist_bp <- as.numeric(ld$Dist_bp)
# ld100$Dist_bp <- as.numeric(ld100$Dist_bp)
ld <- ld[ld$Dist_bp != "NAN",]
# ld100 <- ld[ld100$Dist_bp != "NAN",]

pdf("~/Shiu_Lab/Project/Scripts/Data_Vis/Section_3/geno_LD_decay.pdf")
plot(x=ld$Dist_bp, y=ld$`R^2`, pch=".", cex=2, xlab="Distance window size 50", ylab="R^2")
dev.off()

pdf("geno_LD_decay_100.pdf")
plot(x=ld100$Dist_bp, y=ld100$`R^2`, pch=".", cex=2, xlab="Distance window size 50", ylab="R^2")
dev.off()

# Bin by distance (50 bp bins)
ld_sub <- ld[, c("Dist_bp", "R^2", "DPrime", "pDiseq")] # remove rows with missing values
ld_sub$bin <- cut(ld_sub$Dist_bp, breaks=c(seq(1, 40543, by=50))) # bin by distance
ld_sub$bin <- as.character(lapply(strsplit(as.character(ld_sub$bin), split=","),head, n=1)) # rename bins
ld_sub$bin <- gsub("\\(", "", ld_sub$bin)

pdf("~/Shiu_Lab/Project/Scripts/Data_Vis/Section_3/geno_LD_decay_binned_median.pdf")
ld_sub %>% group_by(bin) %>% summarise(med_R2=median(`R^2`)) %>%
ggplot(aes(x=as.numeric(bin), y=med_R2)) + geom_point() + theme_bw() +
        xlab("Distance (50 bp bins)") + ylab(parse(text="MedianR^2"))
dev.off()

ld_sub %>% group_by(bin) %>% summarise(med_R2=median(`R^2`),
        mean_R2=mean(`R^2`), sd_R2=sd(`R^2`)) %>% write.table(
        "geno_LD_window_50_binned.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# plot the mean_R2 values and the bars with no whiskers from the sd_R2 values
ld_sub <- fread("geno_LD_window_50_binned.txt")
pdf("~/Shiu_Lab/Project/Scripts/Data_Vis/Section_3/geno_LD_decay_binned_mean.pdf")
ld_sub %>%
ggplot(aes(x=as.numeric(bin), y=mean_R2)) + geom_point() + theme_bw() +
        xlab("Distance (50 bp bins)") + ylab(parse(text="MeanR^2"))
dev.off()

pdf("~/Shiu_Lab/Project/Scripts/Data_Vis/Section_3/geno_LD_decay_binned_mean_sd.pdf")
ld_sub %>%
ggplot(aes(x=as.numeric(bin), y=mean_R2)) + geom_point(color="blue") + 
        geom_errorbar(aes(ymin=mean_R2-sd_R2, ymax=mean_R2+sd_R2), width=0.05, alpha=0.3) +
        theme_bw() + xlab("Distance (50 bp bins)") + ylab(parse(text="MeanR^2"))
dev.off()

# Fit a non linear model using the arbitrary C value
# N is the number of the genotypes that have the SNP site
file <- ld[,c(1,2,7,8,13:17)]
Cstart <- c(C=0.1) # C values range from 0.5 to 2, start with 0.1

modelC <- nls(`R^2` ~ ( (10+C*Dist_bp)/( (2+C*Dist_bp) * (11+C*Dist_bp) ) ) * 
            ( 1+( (3+C*Dist_bp) * (12+12*C*Dist_bp+(C*Dist_bp)^2) ) / ( 2*N*(2+C*Dist_bp) * (11+C*Dist_bp) ) ), 
            data=file, start=Cstart, control=nls.control(maxiter=100))

# extract rho, the recombination parameter in 4Nr
rho <- summary(modelC)$parameters[1]

# feed in the new value of rho to obtain LD values adjusted for their Dist_bpances along the chromosome/genome
newrsq <- ( (10+rho*file$Dist_bp) / ( (2+rho*file$Dist_bp) * (11+rho*file$Dist_bp) ) ) *
        ( 1 + ( (3+rho * file$Dist_bp) * (12+12*rho*file$Dist_bp + (rho*file$Dist_bp)^2) ) / 
        (2*file$N*(2+rho*file$Dist_bp) * (11+rho*file$Dist_bp) ) )

newfile <- data.frame(file$Dist_bp, newrsq)

maxld <- max(newfile$newrsq,na.rm=TRUE) #using max LD value from adjusted data
halfdecay = maxld*0.5
halfdecayDist_bp <- newfile$file.Dist_bp[which.min(abs(newfile$newrsq-halfdecay))]
newfile <- newfile[order(newfile$file.Dist_bp),]

# plotting the values
pdf("~/Shiu_Lab/Project/Scripts/Data_Vis/Section_3/LD_decay_window_50.pdf", height=5, width = 5)
png("~/Shiu_Lab/Project/Scripts/Data_Vis/Section_3/LD_decay_window_50.png", type='cairo')
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 
plot(file$Dist_bp, file$rsq, pch=".", cex=2, xlab="Distance (bp)", ylab=expression(LD ~ (r^2)), col="grey")
lines(newfile$file.Dist_bp, newfile$newrsq, col="red", lwd=2)
abline(h=0.1, col="blue") # if you need to add horizontal line
abline(v=halfdecayDist_bp, col="green")
mtext(round(halfdecayDist_bp,2), side=1, line=0.05, at=halfdecayDist_bp, cex=0.75, col="green")
dev.off()