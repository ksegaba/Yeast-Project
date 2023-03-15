#!/usr/bin/env Rscript
#SBATCH --job-name=CS
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --array=1-35
#SBATCH --output=../logs/%x_%A_%a
#SBATCH --constraint=[intel18|intel16]

system("module load GCC/9.3.0  OpenMPI/4.0.3  R/4.0.3")

jobID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))


library(data.table)
library(BGData)
library(BGLR)
# X=fread('/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno.csv',data.table=FALSE)
X=fread('/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/geno_pruned.csv',data.table=FALSE)

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

Y=read.csv('/mnt/home/seguraab/Shiu_Lab/Project/Data/Peter_2018/pheno.csv')
stopifnot(all(Y$ID==rownames(X)))
X=X+1
G=getG(X,nCores=1,chunkSize=3000,verbose=TRUE,center=TRUE,scale=TRUE,scaleG=TRUE)

traits=colnames(Y)[-1]
trait=traits[jobID]
y=scale(Y[,trait])

setwd('../output')
dir.create(trait)
setwd(trait)

# This model includes a GRM ()
fm=BGLR(y=y,ETA=list(
  list(K=G,model='RKHS'),
  list(X=X,model='BayesC',saveEffects=TRUE,probIn=1/100,counts=1e10)),
  nIter=32000,burnIn=2000,saveAt='pruned')
  # nIter=32000,burnIn=2000)
# save(fm,file=paste0('fm_',trait,'.RData'))
save(fm,file=paste0('fm_',trait,'_pruned.RData'))

unlink('*.dat')
quit(save='no')

