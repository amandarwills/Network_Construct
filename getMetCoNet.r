#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if(length(args)==0)
{
	cat("cmdline e.g. -> ./getMetCoNet.r neg_mcap_tp1_dam.txt sample_info.txt T1\n")
	stop("Need arguments!", call.=FALSE)
}

#Load libraries
suppressWarnings(suppressMessages(library(DGCA, quietly = TRUE)))
suppressWarnings(suppressMessages(library(WGCNA, quietly = TRUE)))
suppressWarnings(suppressMessages(library(matrixStats, quietly = TRUE)))

#Parameters
maxPval=0.05

#Read data
metabolites=read.table(args[1],h=T,row.names=1,sep="\t",check.names=F)
samplesInfo=read.table(args[2],h=T,row.names=1,sep="\t")
timepoint=args[3]
A=samplesInfo[which(samplesInfo$TP==timepoint & samplesInfo$COND=="A"),]
H=samplesInfo[which(samplesInfo$TP==timepoint & samplesInfo$COND=="H"),]
A.norm=sweep(metabolites[,rownames(A)],2,as.numeric(A$WEIGHT),FUN='/')
H.norm=sweep(metabolites[,rownames(H)],2,as.numeric(H$WEIGHT),FUN='/')
data=cbind(A.norm,H.norm)
dataA=t(data)
corAMB=matCorr(dataA,corrType="pearson")
pairscorAMB=data.frame(n1=rownames(corAMB)[row(corAMB)],n2=colnames(corAMB)[col(corAMB)],cor=c(corAMB))
corAMBrow=rownames(corAMB)
corAMBcol=colnames(corAMB)
nsampleA=matNSamp(dataA)
corAMBpval=matCorSig(corAMB,nsampleA)
colnames(corAMBpval)=corAMBcol
rownames(corAMBpval)=corAMBrow
pairsPval=data.frame(n1=rownames(corAMBpval)[row(corAMBpval)],n2=colnames(corAMBpval)[col(corAMBpval)],pval=c(corAMBpval))
corAMBpvalVec=as.vector(corAMBpval)
corAMBAdjPval=adjustPVals(corAMBpvalVec,adjust="BH")
corAMBAdjPval=as.numeric(format.pval(corAMBAdjPval,digits=2,nsmall=3))
dim(corAMBAdjPval)=dim(corAMBpval)
colnames(corAMBAdjPval)=corAMBcol
rownames(corAMBAdjPval)=corAMBrow
pairsAdjPval=data.frame(n1=rownames(corAMBAdjPval)[row(corAMBAdjPval)],n2=colnames(corAMBAdjPval)[col(corAMBAdjPval)],adjPval=c(corAMBAdjPval))
corAMBval=cbind(pairscorAMB,pval=pairsPval$pval,adjPval=pairsAdjPval$adjPval)
corAMBvalFinal=corAMBval[complete.cases(corAMBval),]
corAMBvalFinalFiltered=corAMBvalFinal[corAMBvalFinal$adjPval <= maxPval,]
output=paste(args[1],".cor.txt",sep = "")
write.table(corAMBvalFinalFiltered,file=output,quote=FALSE,sep="\t",row.names=FALSE)
tree<-hclust(as.dist(1-corAMB),method="average")
module_labels <- cutreeDynamicTree(dendro=tree, minModuleSize=10,deepSplit=TRUE)
m <- as.data.frame(module_labels)
rownames(m)=rownames(corAMB)
output=paste(args[1],".modules.txt",sep = "")
write.table(m,file=output,quote=FALSE,sep="\t",row.names=T,col.names=F)
