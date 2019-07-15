##### suppose the current folder is C:\Shengtong\Research\rare-var\rare-var-project\code
library(RSQLite)
library(dplyr)
library(knitr)
library(kableExtra)
library(RColorBrewer)
library(gplots)
library(tidyverse)
library(gridExtra)
library(ggpubr)
set.seed(123)
source("MIRAGE_burden_only_variant_level.R")

test.func=function(evid, Data, N1, N0) # given evid, and sample size, perform the burden analysis
{
  evid.data=Data[Data$ID %in% evid,]
  count=c(sum(evid.data$No.case), sum(evid.data$No.contr))
  Time=c(N1, N0)
  pois.test=poisson.test(count, Time, r=1, alternative="greater")
  return(result=list(odds.ratio=pois.test$estimate, p.value=pois.test$p.value, rate.case=count[1]/N1, rate.contr=count[2]/N0))
}


### read into the data

All.Anno.Data=as_tibble(read.table("../../AnnotatedTrans.txt", header=T))
#All.Anno.Data=read.table("../../AnnotatedTrans.txt", header=T)
#All.Anno.Data=read.table("C:\\Shengtong\\Research\\rare-var\\AnnotatedTrans.txt", header=T)


##### pre-process the data
N1=4315; N0=4315
All.Anno.Data[All.Anno.Data =="."] <- NA
All.Anno.Data$ExacAF[is.na(All.Anno.Data$ExacAF)]=0 # set AF of NA to zero
Anno.Data=All.Anno.Data[which(All.Anno.Data$ExacAF<0.05 & All.Anno.Data$Annotation!="synonymous SNV"),] # use AF cutoff and exclude synonymous SNV
Anno.Data=All.Anno.Data[which(All.Anno.Data$Annotation!="synonymous SNV"),] # use AF cutoff and exclude synonumous SNV
var.data=data.frame(ID=Anno.Data$ID, No.case=Anno.Data$No.case, No.contr=Anno.Data$No.contr)



### different gene set

GeneDB=src_sqlite(path="C:\\Shengtong\\Research\\rare-var\\gene.list.db", create=F)
gene_cate1=data.frame(collect(tbl(GeneDB, "SFARI_HighConf")))
gene_cate2=data.frame(collect(tbl(GeneDB, "SFARI_StrongCand")))
gene_cate3=data.frame(collect(tbl(GeneDB, "SFARI_cate3_gene")))
gene_cate4=data.frame(collect(tbl(GeneDB, "SFARI_cate4_gene")))
gene_cate5=data.frame(collect(tbl(GeneDB, "SFARI_cate5_gene")))
gene_cate6=data.frame(collect(tbl(GeneDB, "SFARI_cate6_gene")))
gene_cateS=data.frame(collect(tbl(GeneDB, "SFARI_cateS_gene")))
IDGene=data.frame(collect(tbl(GeneDB, "Pinto14AJHG_IDgene")))
TADAGene=data.frame(collect(tbl(GeneDB, "TADAGenelist")))
Qlessthan5percentgene=TADAGene$TadaName[TADAGene$qvalue.combined<0.05]
Qlessthan20percentgene=TADAGene$TadaName[TADAGene$qvalue.combined<0.2]
Qlessthan30percentgene=TADAGene$TadaName[TADAGene$qvalue.combined<0.3]
Qlessthan40percentgene=TADAGene$TadaName[TADAGene$qvalue.combined<0.4]
Qlessthan50percentgene=TADAGene$TadaName[TADAGene$qvalue.combined<0.5]
Qlargerthan90percentgene=TADAGene$TadaName[TADAGene$qvalue.combined>0.9]
purcell.genelist=data.frame(collect(tbl(GeneDB, "Purcell2014_genelist"))) ## PSD gene, SCZdenovo gene
ASD.gene=data.frame(collect(tbl(GeneDB, "AutismKB_gene")))
constraint.gene=data.frame(collect(tbl(GeneDB, "Samocha_2014NG_constraintgene")))$gene
RVIS.Allgene=data.frame(collect(tbl(GeneDB, "RVIS_gene")))
RVIS.gene=RVIS.Allgene$GeneID[RVIS.Allgene$RVIS.percentile<5] # top 5% gene
haploinsuff.gene=data.frame(collect(tbl(GeneDB, "Petrovski_plosgen_haploinsuff_gene")))
purcell.genelist=data.frame(collect(tbl(GeneDB, "Purcell2014_genelist"))) ## PSD gene, SCZdenovo gene
ASD.gene=data.frame(collect(tbl(GeneDB, "AutismKB_gene")))
constraint.gene=data.frame(collect(tbl(GeneDB, "Samocha_2014NG_constraintgene")))$gene
RVIS.Allgene=data.frame(collect(tbl(GeneDB, "RVIS_gene")))
RVIS.gene=RVIS.Allgene$GeneID[RVIS.Allgene$RVIS.percentile<5] # top 5% gene
haploinsuff.gene=data.frame(collect(tbl(GeneDB, "Petrovski_plosgen_haploinsuff_gene")))
gene.set=c("IDgene","High", "Mod", "PSD", "FMRP", "AutismKB", "constraintgene", "RVIS", "Haploinsuff", "SCZgene")
max.gene=length(gene.set)
LoF.def=c("stopgain", "frameshift substitution", "splicing", "stoploss")
LoF.var=as.character(Anno.Data$ID[which(Anno.Data$Annotation %in% LoF.def==T)])
gene.summy=matrix(nrow=max.gene, ncol=4)
gene.evid=list(); var.index=1


gene.evid[[1]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% IDGene$GeneID )])
high.conf=union(union(gene_cate1$GeneID, gene_cate2$GeneID),Qlessthan5percentgene)
mod.conf=setdiff(union(union(gene_cate3$GeneID, gene_cateS$GeneID), Qlessthan30percentgene),Qlessthan5percentgene)
gene.evid[[2]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% high.conf )])
gene.evid[[3]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% mod.conf )])
psd.gene=purcell.genelist$Gene_symbol[purcell.genelist$PSD=="Y"]
FMRP.gene=purcell.genelist$Gene_symbol[purcell.genelist$FMRP.target=="Y"]
gene.evid[[4]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% psd.gene )])
gene.evid[[5]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% FMRP.gene )])
gene.evid[[6]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% ASD.gene$GeneID )])
gene.evid[[7]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% constraint.gene )])
gene.evid[[8]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% RVIS.gene )])
gene.evid[[9]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% haploinsuff.gene$GeneID )])
gene.evid[[10]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% purcell.genelist$Gene_symbol )])






## Exon level analysis

# All subsequence analsyis foucs on variant with AF < 5% and exluding synonymous mutautions. **

#  burden analysis-exon level
exon.cutoff=c(0.9, 0.8, 0.7, 0.6, 0.5)
exon.fea=c(paste("ExonTop", exon.cutoff*100, "%",sep=""), "CriticalExon")
exon.summ=matrix(nrow=(1+length(exon.cutoff)), ncol=4)
evid.exon=list()
colnames(exon.summ)=c("OR", "p.value", "rate.ca", "rate.co")
rownames(exon.summ)=c(paste("Top", exon.cutoff*100, "%", sep=""), "critical exon")
for (i in 1:length(exon.cutoff))
{
  threshold=quantile(Anno.Data$Exon.Exac.Cons, prob=exon.cutoff[i], na.rm=T)
  evid.exon[[i]]=Anno.Data$ID[which(Anno.Data$Exon.Exac.Cons>threshold)]
  # cat(i, "is running", "\n")
  pois.test=test.func(evid.exon[[i]], var.data, N1, N0)
  exon.summ[i,]=c(pois.test$odds.ratio, pois.test$p.value, pois.test$rate.case, pois.test$rate.contr)
}
threshold=quantile(Anno.Data$Exon.Exac.Cons, prob=0.75, na.rm=T)
evid.exon[[1+length(exon.cutoff)]]=Anno.Data$ID[which( Anno.Data$Exon.Exac.Cons>threshold)]
exon.thrshd=quantile(Anno.Data$Exon.Brain.Exp, prob=0.75, na.rm=T)
expre.exon=Anno.Data$ID[which( Anno.Data$Exon.Brain.Exp>exon.thrshd)]
critical.exon=union(evid.exon[[1+length(exon.cutoff)]],expre.exon)
pois.test=test.func(critical.exon, var.data, N1, N0)
exon.summ[(1+length(exon.cutoff)),]=c(pois.test$odds.ratio, pois.test$p.value, pois.test$rate.case, pois.test$rate.contr)






## Combined feature analysis



#################  record burden and mirage for each cate split from combined feature of exon+gene set ###########
Burden.cate=list()
MIRAGE.pvalue=list()
MIRAGE.para.est=list()
MIRAGE.cate.index=list()
nn=0
for (i in 1:length(exon.fea))
{ # i=1


    for (j in 1:length(gene.set))
    { #  j=1
      cat(i, "th exon.fea  ", length(exon.fea),  j, "th gene of ", length(gene.set),  " is running", "\n")
        comb.evid=intersect(evid.exon[[i]], gene.evid[[j]]) # Exon interacts with gene set
        ###################### MIRAGE burden
        cand.data=Anno.Data[which(Anno.Data$ID %in% comb.evid),]
        gene.data=eight.partition(cand.data)   # partition by 8 categories
        overlap.data=gene.data[gene.data$ID %in% comb.evid,]
        order.overlap.data=overlap.data[order(overlap.data$group.index, decreasing=F),]

        psbl.index=unique(order.overlap.data$group.index); actu.num.group=length(psbl.index)  # re-order the group index
        delta.init=runif(1); beta.init=runif(actu.num.group)
        order.overlap.data$original.group.index=order.overlap.data$group.index

        if (nrow(order.overlap.data)>0)
        {
          burden.matrix=matrix(nrow=actu.num.group, ncol=4)  # this is burden for every category split from combined feature
          colnames(burden.matrix)=c("OR", "p.value", "rate.case", "rate.contr")
          rownames(burden.matrix)=paste(exon.fea[i], gene.set[j], "cate", psbl.index, sep="_")
          for (jj in 1:actu.num.group)
          {
            order.overlap.data$group.index[order.overlap.data$group.index==psbl.index[jj]]=jj # re-index the group labels
            pois.test=test.func(order.overlap.data[order.overlap.data$original.group.index==psbl.index[jj],]$ID, var.data, N1, N0)
            burden.matrix[jj,]=c(pois.test$odds.ratio, pois.test$p.value, pois.test$rate.case, pois.test$rate.contr)
          }

          para.est=multi.group.func.for.variant(order.overlap.data, N1, N0, gamma.mean=3, sigma=2, delta=0.2, beta.init, actu.num.group)
          nn=nn+1
          MIRAGE.pvalue[[nn]]=para.est$cate.pvalue
          MIRAGE.para.est[[nn]]=para.est$beta.est
          MIRAGE.cate.index[[nn]]=psbl.index
          Burden.cate[[nn]]=burden.matrix
        }

    } ################## end of j

}  ############# end of i

save(MIRAGE.cate.index, MIRAGE.para.est, MIRAGE.pvalue, Burden.cate, file="..\\output\\CombinedFeature\\Burden.MIRAGE.for.Combined.exon.geneset.partition.RData")
