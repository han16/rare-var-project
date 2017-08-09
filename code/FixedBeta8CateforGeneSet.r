######################
# This code computes BF of every gene with fixed parameter within each variant category.
######################
rm(list=ls())
library(RSQLite)
library(dplyr)
library(knitr)

intergrand=function(aa, var.case, var.contr, bar.gamma, sig, N1, N0)
{
  ff=dbinom(var.case, sum(var.case, var.contr), aa*N1/(aa*N1+N0))*dgamma(aa, bar.gamma*sig, sig)
  return(ff)
}
# calculate the bayes factor of a single variant via integration
BF.gene.inte=function(var.case, var.contr, bar.gamma, sig, N1, N0)
{
  marglik0.CC <- dbinom(var.case, sum(var.case, var.contr), N1/(N1+N0))    # Under H0: gamma=1

  marglik1.CC <- integrate(intergrand, var.case, var.contr, bar.gamma, sig, N1, N0, low=0, upper=100, stop.on.error=F)$value # Under H1: gamma~gamma(gamma.mean*sigma, sigma)
  BF.var <- marglik1.CC/marglik0.CC

  return(BF.var)
}
fixed.beta.func=function(new.data, N1, N0, gamma.mean, sigma, beta, num.group) # new.data has one column specifying its group index
{
  #######################
  beta.k=beta
  full.info.genevar=list()
  gene.list=new.data$Gene; unique.gene=unique(gene.list) # find the gene list
  num.gene=length(unique.gene)
  BF.gene=numeric()
  ########################
  # calculate the Bayes factor for variant (i,j) and gene i as initials.
  for (i in 1:num.gene)
  {
    cat(i, "th gene of ", "\t", num.gene, "\t", "is running", "\n")
    var.index.list=which(gene.list==unique.gene[i])
    indi.gene=new.data[var.index.list,]   # note var.index.list matches new.data
    bb=1; var.BF=numeric()
    if (length(var.index.list)>0) # calculate Bayes factor for variant (i,j)
      for (j in 1:length(var.index.list))
      {
        if (new.data$group.index[var.index.list[j]]<=2)
          var.BF[j]=BF.gene.inte(new.data$No.case[var.index.list[j]], new.data$No.contr[var.index.list[j]], bar.gamma=6, sig=sigma, N1, N0)
        if (new.data$group.index[var.index.list[j]]>2)
          var.BF[j]=BF.gene.inte(new.data$No.case[var.index.list[j]], new.data$No.contr[var.index.list[j]], bar.gamma=gamma.mean, sig=sigma, N1, N0)
        bb=bb*((1-beta.k[new.data$group.index[var.index.list[j]]])+beta.k[new.data$group.index[var.index.list[j]]]*var.BF[j])
      }
    full.info.genevar[[i]]=cbind(indi.gene, var.BF)
    BF.gene[i]=bb
  }

  return(result=list(BayesFactor=data.frame(Gene=unique.gene, BF=BF.gene), Full.info=full.info.genevar))
}

#All.Anno.Data=read.table("D:\\ResearchWork\\StatisticalGenetics\\Rare-variant-project\\AnnotatedTrans.txt", header=T)
All.Anno.Data=read.table("C:\\han\\ResearchWork\\StatGene\\AutismData\\AnnotatedTrans.txt", header=T)
#All.Anno.Data=read.table("C:\\han\\ResearchWork\\StatGene\\SCZData\\AnnotatedSCZ.txt", header=T)
Anno.Data=All.Anno.Data[which(All.Anno.Data$ExacAF<0.05 & All.Anno.Data$Annotation!="synonymous SNV"),] # use AF cutoff and exclude synonumous SNV
#N1=5072/2; N0=5086/2
N1=4315; N0=4315
var.data=data.frame(ID=Anno.Data$ID, No.case=Anno.Data$No.case, No.contr=Anno.Data$No.contr)
#gene.set=read.csv("C:\\han\\ResearchWork\\StatGene\\AutismData\\lessconstraint1000genes.csv", header=T)
#gene.set=as.character(read.csv("D:\\ResearchWork\\StatisticalGenetics\\Rare-variant-project\\rare-var-project\\data\\GeneSet\\Samocha_2014NG_contraintgene.csv", header=T)$gene)
gene.set=as.character(read.table("C:\\han\\ResearchWork\\StatGene\\AutismData\\contraint.gene.txt", header=T)[[1]])
gene.FDR=read.table("C:\\han\\ResearchWork\\StatGene\\170726_to_Shengtong_gene_FDR_based_on_Sanders_Neuron_denovo_coding_SNV", header=T)  # use 1-FDR as prior


vart.set=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% gene.set)])
cand.data=Anno.Data[which(Anno.Data$ID %in% vart.set),]
par.evid=list()
LoF.def=c("stopgain", "frameshift substitution", "splicing", "stoploss")
par.evid[[1]]=which(cand.data$Annotation %in% LoF.def==T & cand.data$ExacAF<0.05 & cand.data$ExacAF>=0.01 )
par.evid[[2]]=which(cand.data$Annotation %in% LoF.def==T &  cand.data$ExacAF<0.01)
par.evid[[3]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))>=0.957 & cand.data$ExacAF>=0.01 & cand.data$ExacAF<0.05)
par.evid[[4]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))>=0.957 & cand.data$ExacAF>=0.001 & cand.data$ExacAF<0.01)
par.evid[[5]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))>=0.957 & cand.data$ExacAF<0.001)
par.evid[[6]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))<0.957 & cand.data$ExacAF>=0.01 & cand.data$ExacAF<0.05)
par.evid[[7]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))<0.957 & cand.data$ExacAF>=0.001 & cand.data$ExacAF<0.01)
par.evid[[8]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))<0.957 & cand.data$ExacAF<0.001)
group.index=rep(NA, nrow(cand.data))
for (i in 1:length(par.evid))
  group.index[par.evid[[i]]]=i
group.index[is.na(group.index)]=6  # note here
gene.data=data.frame(ID=cand.data$ID, Gene=cand.data$Gene, No.case=cand.data$No.case, No.contr=cand.data$No.contr, group.index=group.index)
num.group=8; beta.init=c(0.5, 0.9, 0.1, 0.2, 0.4, 0.05, 0.1, 0.2)
gene.result=fixed.beta.func(gene.data, N1, N0, gamma.mean=3, sigma=2, beta.init, num.group)
gene.BF=gene.result[[1]]
gene.fdr=numeric()
for (i in 1:length(gene.BF$Gene))
 if (length(which(gene.FDR$genename==as.character(gene.BF$Gene)[i]))>0)
 gene.fdr[i]=gene.FDR$FDR[which(gene.FDR$genename==as.character(gene.BF$Gene)[i])]

gene.comb=cbind(gene.BF, gene.fdr)
BF.with.prior=gene.comb$gene.fdr+(1-gene.comb$gene.fdr)*gene.comb$BF
gene.res=cbind(gene.comb, BF.with.prior)
gene.res.order=gene.res[order(gene.res$BF.with.prior, decreasing=T),]
top20gene=as.character(gene.res.order$Gene)[1:20]
################ rank of top genes in other gene sets
denovo=read.table("C:\\Users\\han\\Dropbox\\StatisticalGenetics\\TADA_SNV_CNV_combined_Feb7.txt", header=T)
order.denovo=denovo[order(denovo$qvalue.combined, decreasing=F),]
denovo.rank=numeric()
for (i in 1:length(top20gene))
 denovo.rank[i]=which(order.denovo$RefSeqName==top20gene[i])

ASD.prior=read.table("C:\\Users\\han\\Dropbox\\StatisticalGenetics\\ASD_Priors.txt", header=T)
order.ASD.prior=ASD.prior[order(ASD.prior$post, decreasing=T),] 
ASD.prior.rank=numeric()
for (i in 1:length(top20gene))
if (length(which(order.ASD.prior$Gene==top20gene[i]))>0) 
 ASD.prior.rank[i]=which(order.ASD.prior$Gene==top20gene[i])







