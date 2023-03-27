
library(data.table)
library(BGData)
library(BGLR)
library(ggplot2)

X=fread('/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv',data.table=FALSE)
<<<<<<< HEAD
# X=fread('/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno_pruned.csv',data.table=FALSE)
=======

>>>>>>> 2f27eb9783697f60426388411650f4fdb22e190b
ID=X[,1]
X=as.matrix(X[,-1])
rownames(X)=ID

MAP=fread('../data/feat_map_for_GDLC.csv',data.table=FALSE)
tmp=matrix(byrow=TRUE,ncol=2,data=unlist(strsplit(MAP[,1],split='_')))
getChr=function(x){
  n=length(unlist(strsplit(tmp,split='')))
  return(substr(x,start=11,stop=11+n-1))
}
MAP=MAP[,1,drop=FALSE]
MAP$chr=as.integer(getChr(tmp[,1]))
MAP$pos=as.numeric(tmp[,2])

<<<<<<< HEAD
# subset MAP (keep only SNPs not pruned)
MAP <- MAP[MAP$ID %in% colnames(X),]

=======
>>>>>>> 2f27eb9783697f60426388411650f4fdb22e190b
Y=read.csv('/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv')
stopifnot(all(Y$ID==rownames(X)))
X=X+1

VAR_COMP=data.frame(trait=colnames(Y)[-1],h2=NA,vG=NA,vE=NA,vQTL=NA)

models=list.files('../output/',recursive=T,full.names=TRUE,pattern='.RData')
<<<<<<< HEAD
# models=list.files('../output/',recursive=T,full.names=TRUE,pattern='pruned.RData')
=======
>>>>>>> 2f27eb9783697f60426388411650f4fdb22e190b

setwd('../output')

MAP$Mbp=0
for(i in unique(MAP$chr)){
  MAP$Mbp[MAP$chr==i]=max(MAP$Mbp)+MAP$pos[MAP$chr==i]/1e6+ifelse(i>1,.1,0)
}
traits=colnames(Y)[-1]

for(i in 1:length(models)){
<<<<<<< HEAD
       load(models[i])
       trait=strsplit(models[i], "/")[[1]][4]
       print(trait)
       setwd(trait)
       B=readBinMat('ETA_2_b.bin')!=0
       # B=readBinMat('prunedETA_2_b.bin')!=0
       pdf('posterior_probs.pdf')
       # pdf('posterior_probs_pruned.pdf')
       plot(colMeans(B))
       dev.off()

       VAR_COMP[VAR_COMP$trait==trait,]$h2=1-fm$varE 

       LFDR=1-colMeans(B) # 1 - posterior probability = local FDR
       DS=segments(LFDR,chr=MAP$chr,bp=MAP$Mbp,threshold=0.97,gap=.05)
       
       DS$setPIP=NA # posterior inclusion probability
=======

       load(models[i])
       trait=strsplit(models[i], "/")[[1]][4]
       setwd(trait)
       B=readBinMat('ETA_2_b.bin')!=0
       VAR_COMP[VAR_COMP$trait==trait,]$h2=1-fm$varE 

       LFDR=1-colMeans(B)
       DS=segments(LFDR,chr=MAP$chr,bp=MAP$Mbp,threshold=0.97,gap=.05)
       
       DS$setPIP=NA
>>>>>>> 2f27eb9783697f60426388411650f4fdb22e190b
       for(i in 1:nrow(DS)){
         DS$setPIP[i]=mean(apply(FUN=any,X=B[,DS$start[i]:DS$end[i]]!=0,MARGIN=1))
       }
       
       isQTL=fm$ETA[[2]]$d>.03
       
       DF=data.frame(Mbp=MAP$Mbp,PIP=colMeans(B),chr=MAP$chr)
       
<<<<<<< HEAD
      p=ggplot(DF,aes(x=Mbp,y=PIP,group=chr, col=chr))+
        geom_point(size=.7)+ #color='skyblue',
=======
      p=ggplot(DF,aes(x=Mbp,y=PIP,group=chr))+
        geom_point(color='skyblue',size=.7)+
>>>>>>> 2f27eb9783697f60426388411650f4fdb22e190b
        theme(legend.position="none")+
        ylim(c(0,1))+
        geom_hline(yintercept=.02,linetype='dashed')
      
      for(j in 1:nrow(DS)){
        p=p+annotate("rect", xmin = MAP$Mbp[DS$start[j]], xmax =MAP$Mbp[DS$end[j]], ymin =0, ymax = DS$setPIP[j],
                      alpha = .1,fill = "blue")
      }
       ggsave(p,file=paste0('manhattan_',trait,'_.tiff'),device='tiff')
<<<<<<< HEAD
       # ggsave(p,file=paste0('manhattan_',trait,'_pruned.tiff'),device='tiff')
       write.csv(DS,file=paste0('DS_',trait,'_.csv'))
       # write.csv(DS,file=paste0('DS_',trait,'_pruned.csv'))
     setwd('../')
}

# consider ld pruning, or do ld blocks and pick a couple of snps per block (can be done in plink)
# can do by chromosome to make sure blocks don't have snps that are too far away
# change figure so that it colors by chromosome
=======
       write.csv(DS,file=paste0('DS_',trait,'_.csv'))
     setwd('../')
}
>>>>>>> 2f27eb9783697f60426388411650f4fdb22e190b
