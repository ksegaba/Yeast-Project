# Linkage Disequilibrium plot using output from TASSEL5
# Source: https://github.com/mohsinali1990/My_scripts/blob/main/LD%20decay%20Plot%20from%20TASSEL%20LDoutput.R
library(data.table)

# remove missing values from R^2 and Distance variables
setwd("/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018")
ld <- read.delim('geno_LD_window_50.txt', stringsAsFactors=FALSE, header=TRUE, sep="\t")
ld100 <- read.delim('geno_LD_window_100.txt', stringsAsFactors=FALSE, header=TRUE, sep="\t")
ld <- ld[ld$R^2 != "NaN",]
ld100 <- ld100[ld100$R^2 != "NaN",]
ld$Dist_bp <- as.numeric(ld$Dist_bp)
ld100$Dist_bp <- as.numeric(ld100$Dist_bp)
ld <- ld[ld$Dist_bp != "NAN",]
ld100 <- ld[ld100$Dist_bp != "NAN",]

pdf("geno_LD_decay.pdf")
plot(x=ld$Dist_bp, y=ld$`R^2`, pch=".", cex=2, xlab="Distance window size 50", ylab="R^2")
dev.off()

pdf("geno_LD_decay_100.pdf")
plot(x=ld100$Dist_bp, y=ld100$`R^2`, pch=".", cex=2, xlab="Distance window size 50", ylab="R^2")
dev.off()

file <- ld[,c(1,2,7,8,13:17)]

# C values range from 0.5 to 2, start with 0.1
Cstart <- c(C=0.1)

# fit a non linear model using the arbitrary C value
# N is the number of the genotypes that have the SNP site
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
pdf("LD_decay_window_50.pdf", height=5, width = 5)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 
plot(file$Dist_bp, file$rsq, pch=".", cex=2, xlab="Distance (bp)", ylab=expression(LD ~ (r^2)), col="grey")
lines(newfile$file.Dist_bp, newfile$newrsq, col="red", lwd=2)
abline(h=0.1, col="blue") # if you need to add horizental line
abline(v=halfdecayDist_bp, col="green")
mtext(round(halfdecayDist_bp,2), side=1, line=0.05, at=halfdecayDist_bp, cex=0.75, col="green")
dev.off()