#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if(length(args)==0)
{
	cat("cmdline e.g.-> ./getGeneCoNet.r combined_salmon.isoform.tpm.matrix.renamed TP1.regulated_FC1_P005.txt TP1 outfilename\n")
	stop("Need arguments!", call.=FALSE)
}
  
#Load libraries
suppressMessages(suppressWarnings(library(DGCA, quietly = TRUE)))
suppressMessages(suppressWarnings(library(WGCNA, quietly = TRUE)))
suppressMessages(suppressWarnings(library(matrixStats, quietly = TRUE)))

#Parameters
maxPval=0.05

#Read data
tpm.data=as.matrix(read.table(args[1],h=T,row.names=1,sep="\t"))
deg.data=read.table(args[2],h=T,row.names=1,sep="\t")

#sigG=as.vector(read.table(args[2])[,1])
deg.ids=rownames(deg.data)
timepoint=args[3]

#Remove non numerical data and reformat the matrix
indx <- grepl(timepoint, colnames(tpm.data))
TP=which(indx==T)
tpm.data=tpm.data[deg.ids,TP]
C=colnames(tpm.data)
R=rownames(tpm.data)
dims=dim(tpm.data)
tpm.data=as.numeric(tpm.data)
dim(tpm.data)=dims
colnames(tpm.data)=C
rownames(tpm.data)=R

##Start computing co-expression networks##
t.data=t(tpm.data)
cor.data=matCorr(t.data,corrType="pearson")
pairs.cor.data=data.frame(n1=rownames(cor.data)[row(cor.data)],n2=colnames(cor.data)[col(cor.data)],cor=c(cor.data))
cor.data.row=rownames(cor.data)
cor.data.col=colnames(cor.data)
nsample.data=matNSamp(t.data)
cor.data.pval=matCorSig(cor.data,nsample.data)
colnames(cor.data.pval)=cor.data.col
rownames(cor.data.pval)=cor.data.row
pairsPval=data.frame(n1=rownames(cor.data.pval)[row(cor.data.pval)],n2=colnames(cor.data.pval)[col(cor.data.pval)],pval=c(cor.data.pval))
cor.data.pval.vec=as.vector(cor.data.pval)
cor.data.adjPval=adjustPVals(cor.data.pval.vec,adjust="BH")
cor.data.adjPval=as.numeric(format.pval(cor.data.adjPval,digits=2,nsmall=3))
dim(cor.data.adjPval)=dim(cor.data.pval)
colnames(cor.data.adjPval)=cor.data.col
rownames(cor.data.adjPval)=cor.data.row
pairsAdjPval=data.frame(n1=rownames(cor.data.adjPval)[row(cor.data.adjPval)],n2=colnames(cor.data.adjPval)[col(cor.data.adjPval)],adjPval=c(cor.data.adjPval))
#
cor.data.val=cbind(pairs.cor.data,pval=pairsPval$pval,adjPval=pairsAdjPval$adjPval)
cor.data.val.final=cor.data.val[complete.cases(cor.data.val),]
cor.data.val.final.filtered=cor.data.val.final[cor.data.val.final$adjPval <= maxPval,]
#
output=paste(arg[4],".corData.txt",sep = "")
write.table(cor.data.val.final.filtered,file=output,quote=FALSE,sep="\t",row.names=FALSE)
#
gene_tree<-hclust(as.dist(1-cor.data),method="average")
module_labels <- cutreeDynamicTree(dendro=gene_tree, minModuleSize=10,deepSplit=TRUE)
m <- as.data.frame(module_labels)
rownames(m)=rownames(cor.data)
output=paste(arg[4],".modules.txt",sep = "")
write.table(m,file=output,quote=FALSE,sep="\t",row.names=T,col.names=F)
#END
