rm(list=ls())
library(RSQLite)
library(dplyr)
library(knitr)
library(kableExtra)
library(RColorBrewer)
library(gplots)
library(tidyverse)
library(gridExtra)
library(ggpubr)
source("MIRAGE_burden_only_variant_level.R")
set.seed(123)


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}



wrap.it <- function(x, len)
{
  sapply(x, function(y) paste(strwrap(y, len),
                              collapse = "\n"),
         USE.NAMES = FALSE)
}


# Call this function with a list or vector
wrap.labels <- function(x, len)
{
  if (is.list(x))
  {
    lapply(x, wrap.it, len)
  } else {
    wrap.it(x, len)
  }
}

test.func=function(evid, Data, N1, N0) # given evid, and sample size, perform the burden analysis
{
  evid.data=Data[Data$ID %in% evid,]
  count=c(sum(evid.data$No.case), sum(evid.data$No.contr))
  Time=c(N1, N0)
  pois.test=poisson.test(count, Time, r=1, alternative="greater")
  return(result=list(odds.ratio=pois.test$estimate, p.value=pois.test$p.value, rate.case=count[1]/N1, rate.contr=count[2]/N0))
}


### Annotated data quality. 
All.Anno.Data=as_tibble(read.table("../../AnnotatedTrans.txt", header=T))
N1=N0=4315
All.Anno.Data[All.Anno.Data =="."] <- NA
All.Anno.Data$ExacAF[is.na(All.Anno.Data$ExacAF)]=0 # set AF of NA to zero 
Anno.Data=All.Anno.Data[which(All.Anno.Data$ExacAF<0.05 & All.Anno.Data$Annotation!="synonymous SNV"),] # use AF cutoff and exclude synonymous SNV
Anno.Data=All.Anno.Data[which(All.Anno.Data$Annotation!="synonymous SNV"),] # use AF cutoff and exclude synonumous SNV

var.data=data.frame(ID=Anno.Data$ID, No.case=Anno.Data$No.case, No.contr=Anno.Data$No.contr)



### Variant rate at different AF cutoffs
ExacAF.cutoff=c(1e-2, 1e-3, 1e-4, 1e-5)
AF.summary=matrix(nrow=length(ExacAF.cutoff)+2, ncol=5)
rownames(AF.summary)=c("NULL", paste("MAF<", ExacAF.cutoff, sep=""), "private"); rownames(AF.summary)[1]=""
colnames(AF.summary)=c("No.rows", "No.var.ca", "No.var.co",  "rate.ca", "rate.co")
AF.summary[1,]=c(nrow(All.Anno.Data),sum(All.Anno.Data$No.case),sum(All.Anno.Data$No.contr), sum(All.Anno.Data$No.case)/N1, sum(All.Anno.Data$No.contr)/N0)

for (i in 1:length(ExacAF.cutoff))
{
  #cat(i, "is running", "\n") 
  select.data=All.Anno.Data[which(All.Anno.Data$ExacAF<ExacAF.cutoff[i]),]  
  AF.summary[i+1,]=c(nrow(select.data),sum(select.data$No.case),sum(select.data$No.contr), sum(select.data$No.case)/N1, sum(select.data$No.contr)/N0)
}
select.data=All.Anno.Data[c(which(is.na(All.Anno.Data$ExacAF)==T), which(All.Anno.Data$ExacAF==0)),]
AF.summary[(length(ExacAF.cutoff)+2),]=c(nrow(select.data),sum(select.data$No.case),sum(select.data$No.contr), sum(select.data$No.case)/N1, sum(select.data$No.contr)/N0)



## Variant level analysis

### single variant 

var.fea=c("Prby Damaging", "Psbl Damaging", "SIFT<0.05", "CADD top10%", "BrainExp top10%", "Consensus", "LoF", "Private"); max.vart=length(var.fea)
var.evid=list()
var.evid[[1]]=as.character(Anno.Data$ID[which(as.numeric(as.character(Anno.Data$Polyphen2.HDIV.score))>=0.957 )]) # probably damaging >=0.957
var.evid[[2]]=as.character(Anno.Data$ID[which(as.numeric(as.character(Anno.Data$Polyphen2.HDIV.score))<=0.956 & as.numeric(as.character(Anno.Data$Polyphen2.HDIV.score))>=0.453)]) # possibly damaging
var.evid[[3]]=as.character(Anno.Data$ID[which(as.numeric(as.character(Anno.Data$SIFT.score))<0.05 )]) # deleterious SIFT<0.05
CADD.cutoff=quantile(as.numeric(as.character(Anno.Data$CADD.raw)), prob=0.9, na.rm=TRUE)
var.evid[[4]]=as.character(Anno.Data$ID[which(as.numeric(as.character(Anno.Data$CADD.raw))>CADD.cutoff)]) # CADD top 10% 
BEE.cutoff=quantile(as.numeric(Anno.Data$Exon.Brain.Exp), prob=0.9, na.rm=TRUE)
var.evid[[5]]=as.character(Anno.Data$ID[which(Anno.Data$Exon.Brain.Exp>BEE.cutoff)]) # Brain expression exon top 10%
var.evid[[6]]=union(union(var.evid[[1]], var.evid[[3]]), var.evid[[4]]) # consensus
LoF.def=c("stopgain", "frameshift substitution", "splicing", "stoploss")
var.evid[[7]]=as.character(Anno.Data$ID[which(Anno.Data$Annotation %in% LoF.def==T)]) 
var.evid[[8]]=union(as.character(All.Anno.Data$ID[which(All.Anno.Data$ExacAF==0)]), as.character(All.Anno.Data$ID[which(is.na(All.Anno.Data$ExacAF)==T)]))  
summy=matrix(nrow=max.vart, ncol=4)
colnames(summy)=c("OR", "p.value", "rate.ca", "rate.co")
rownames(summy)=var.fea
for (vart in 1:max.vart)
{
  # cat(vart, "is running", "\n")
  pois.test=test.func(var.evid[[vart]], var.data, N1, N0)
  summy[vart,]=c(pois.test$odds.ratio, pois.test$p.value, pois.test$rate.case, pois.test$rate.contr)
}




### different gene set 
GeneDB=src_sqlite(path="C:\\Shengtong\\Research\\rare-var\\gene.list.db", create=F)
#GeneDB=src_sqlite(path="../../gene.list.db", create=F)
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
gene.set=c("ID gene","High", "Mod", "PSD", "FMRP", "AutismKB", "constraint gene", "RVIS", "Haploinsuff", "SCZ gene", "Olfac.gene")
gene.fea=c("cate1", "cate2", "cate3", "cate4", "cate5", "cate6", "cateS", "TADAq<5%", "TADAq<20%", "TADAq<30%", "TADAq<40%", "TADAq<50%", "TADAq>90%", gene.set); max.gene=length(gene.fea)
LoF.def=c("stopgain", "frameshift substitution", "splicing", "stoploss")
LoF.var=as.character(Anno.Data$ID[which(Anno.Data$Annotation %in% LoF.def==T)])
gene.summy=matrix(nrow=max.gene, ncol=4)
gene.evid=list(); var.index=1
gene.evid[[1]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% gene_cate1$GeneID )])
gene.evid[[2]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% gene_cate2$GeneID  )])
gene.evid[[3]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% gene_cate3$GeneID )])
gene.evid[[4]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% gene_cate4$GeneID )])
gene.evid[[5]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% gene_cate5$GeneID )])
gene.evid[[6]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% gene_cate6$GeneID )])
gene.evid[[7]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% gene_cateS$GeneID )])
gene.evid[[8]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% Qlessthan5percentgene)])
gene.evid[[9]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% Qlessthan20percentgene )])
gene.evid[[10]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% Qlessthan30percentgene )])
gene.evid[[11]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% Qlessthan40percentgene )])
gene.evid[[12]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% Qlessthan50percentgene )])
gene.evid[[13]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% Qlargerthan90percentgene  )])
gene.evid[[14]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% IDGene$GeneID )])
high.conf=union(union(gene_cate1$GeneID, gene_cate2$GeneID),Qlessthan5percentgene)
mod.conf=setdiff(union(union(gene_cate3$GeneID, gene_cateS$GeneID), Qlessthan30percentgene),Qlessthan5percentgene)
gene.evid[[15]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% high.conf )])
gene.evid[[16]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% mod.conf )])
psd.gene=purcell.genelist$Gene_symbol[purcell.genelist$PSD=="Y"]
FMRP.gene=purcell.genelist$Gene_symbol[purcell.genelist$FMRP.target=="Y"]
gene.evid[[17]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% psd.gene )])
gene.evid[[18]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% FMRP.gene )])
gene.evid[[19]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% ASD.gene$GeneID )])
gene.evid[[20]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% constraint.gene )])
gene.evid[[21]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% RVIS.gene )])
gene.evid[[22]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% haploinsuff.gene$GeneID )])
gene.evid[[23]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% purcell.genelist$Gene_symbol )])
olfac.gene=data.frame(collect(tbl(GeneDB, "Olfac_gene")))
gene.evid[[24]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% olfac.gene$GeneID )])
#sczgene=as.character(read.table("D:\\ResearchWork\\StatisticalGenetics\\NumericAnalysis\\RealData\\SCZData\\GeneList\\SCZ.67gene.q0.3.txt", header=T)[[1]])
#gene.evid[[25]]=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% sczgene & Anno.Data$Annotation %in% LoF.def==T)])

colnames(gene.summy)=c("OR", "p.value", "rate.ca", "rate.co")
#rownames(gene.summy)=c(gene.fea, "67SCZriskgene")
rownames(gene.summy)=gene.fea
for (gene in 1:(max.gene))
  if (length(gene.evid[[gene]])>0)
  {
    # cat(gene, "is running", "\n")
    pois.test=test.func(gene.evid[[gene]], var.data, N1, N0)
    gene.summy[gene,]=c(pois.test$odds.ratio, pois.test$p.value, pois.test$rate.case, pois.test$rate.contr)
  }






## Exon level analysis

exon.cutoff=c(0.9, 0.8, 0.7, 0.6, 0.5,0)
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
evid.exon[[1+length(exon.cutoff)]]=critical.exon
pois.test=test.func(critical.exon, var.data, N1, N0)
exon.summ[(1+length(exon.cutoff)),]=c(pois.test$odds.ratio, pois.test$p.value, pois.test$rate.case, pois.test$rate.contr)




## Combined feature analysis


###  combination of gene set+AF+exon

comb.exon.cutoff=exon.cutoff;
#sub.com=c(exon.fea, "LoF", "Prby Damaging", "CriticalExon", "Consensus")
sub.com=c(exon.fea, var.fea)
comb.fea=paste(rep(sub.com, each=length(gene.set)), gene.set)
AF.cutoff=c("MAF<1e-2", "MAF<1e-3", "MAF<1e-4", "Private"); colname=c("OR", "p.value", "rate.ca", "rate.co")
comb.summ=matrix(nrow=length(comb.fea), ncol=4*length(AF.cutoff))
rownames(comb.summ)=comb.fea; colnames(comb.summ)=paste(rep(AF.cutoff, each=length(colname)), colname, sep="~")
comb.summ.MIRAGE=comb.summ
#################  record burden and mirage for each cate split from combined feature ###########
Burden.cate=list()
MIRAGE.pvalue=list()
MIRAGE.para.est=list()
MIRAGE.cate.index=list()
nn=0
all.evid=list(); ii=0; all.evid.fea=paste(rep(comb.fea, each=length(AF.cutoff)), AF.cutoff)
signal1.evid=list(); jj1=0
signal2.evid=list(); jj2=0
signal3.evid=list(); jj3=0
signal4.evid=list(); jj4=0
for (i in 1:length(sub.com))
{# i=1 
  cat(i, "th subcom of ", length(sub.com), " is running", "\n")
  
  if (i<=length(comb.exon.cutoff))    # this is for exon 
    for (j in 1:length(gene.set))
    { # j=1
      for (k in 1:(length(AF.cutoff)-1)) 
      { # k=1
        comb.evid=intersect(intersect(evid.exon[[i]], gene.evid[[13+j]]), var.evid[[k]]) 
        pois.test=test.func(comb.evid, var.data, N1, N0)
        comb.summ[(i-1)*length(gene.set)+j,((k-1)*4+1):(k*4)]=c(pois.test$odds.ratio, pois.test$p.value, pois.test$rate.case, pois.test$rate.contr)
        ii=ii+1
        all.evid[[ii]]=comb.evid
        ###################### MIRAGE burden
        cand.data=Anno.Data[which(Anno.Data$ID %in% comb.evid),]
        gene.data=eight.partition(cand.data)
        overlap.data=gene.data[gene.data$ID %in% comb.evid,]
        order.overlap.data=overlap.data[order(overlap.data$group.index, decreasing=F),]
        psbl.index=unique(order.overlap.data$group.index); actu.num.group=length(psbl.index)
        delta.init=runif(1); beta.init=runif(actu.num.group)
        order.overlap.data$original.group.index=order.overlap.data$group.index
        
        if (nrow(order.overlap.data)>0)
        {  
          burden.matrix=matrix(nrow=actu.num.group, ncol=4)  # this is burden for every category split from combined feature
          colnames(burden.matrix)=c("OR", "p.value", "rate.case", "rate.contr")
          rownames(burden.matrix)=paste(sub.com[i], gene.set[j], AF.cutoff[k], "cate", psbl.index, sep="_")
          
          for (jj in 1:actu.num.group)
          {    
            order.overlap.data$group.index[order.overlap.data$group.index==psbl.index[jj]]=jj # re-index the group labels
            pois.test=test.func(order.overlap.data[order.overlap.data$original.group.index==psbl.index[jj],]$ID, var.data, N1, N0)
            burden.matrix[jj,]=c(pois.test$odds.ratio, pois.test$p.value, pois.test$rate.case, pois.test$rate.contr)
          }
          
          para.est=multi.group.func.for.variant(order.overlap.data, N1, N0, gamma.mean=3, sigma=2, delta=0.2, beta.init, actu.num.group)
          comb.summ.MIRAGE[(i-1)*length(gene.set)+j,((k-1)*4+1):(k*4)]=c(pois.test$odds.ratio, para.est$pvalue[length(para.est$pvalue)], pois.test$rate.case, pois.test$rate.contr)
          nn=nn+1
          MIRAGE.pvalue[[nn]]=para.est$cate.pvalue
          MIRAGE.para.est[[nn]]=para.est$beta.est
          MIRAGE.cate.index[[nn]]=psbl.index
          Burden.cate[[nn]]=burden.matrix
        }
        if (nrow(order.overlap.data)==0)  
          comb.summ.MIRAGE[(i-1)*length(gene.set)+j,((k-1)*4+1):(k*4)]=NA
        
        
        ######################
        
        if (pois.test$p.value<0.05)
        {
          if (k==1)
          {
            jj1=jj1+1
            signal1.evid[[jj1]]=comb.evid
          }
          if (k==2)
          {
            jj2=jj2+1
            signal2.evid[[jj2]]=comb.evid
          }
          if (k==3)
          {
            jj3=jj3+1
            signal3.evid[[jj3]]=comb.evid
          }
        }
      }  
      k=length(AF.cutoff)
      comb.evid=intersect(intersect(evid.exon[[i]], gene.evid[[13+j]]), var.evid[[8]]) 
      pois.test=test.func(comb.evid, var.data, N1, N0)
      comb.summ[(i-1)*length(gene.set)+j,((k-1)*4+1):(k*4)]=c(pois.test$odds.ratio, pois.test$p.value, pois.test$rate.case, pois.test$rate.contr)
      ii=ii+1
      all.evid[[ii]]=comb.evid 
      
      ###################### MIRAGE burden
      cand.data=Anno.Data[which(Anno.Data$ID %in% comb.evid),]
      gene.data=eight.partition(cand.data)
      overlap.data=gene.data[gene.data$ID %in% comb.evid,]
      order.overlap.data=overlap.data[order(overlap.data$group.index, decreasing=F),]
      psbl.index=unique(order.overlap.data$group.index); actu.num.group=length(psbl.index)
      delta.init=runif(1); beta.init=runif(actu.num.group)
      order.overlap.data$original.group.index=order.overlap.data$group.index
      
      if (nrow(order.overlap.data)>0)
      {  
        burden.matrix=matrix(nrow=actu.num.group, ncol=4)  # this is burden for every category split from combined feature
        colnames(burden.matrix)=c("OR", "p.value", "rate.case", "rate.contr")
        rownames(burden.matrix)=paste(sub.com[i], gene.set[j], AF.cutoff[k], "cate", psbl.index, sep="_")
        
        for (jj in 1:actu.num.group)
        {    
          order.overlap.data$group.index[order.overlap.data$group.index==psbl.index[jj]]=jj # re-index the group labels
          pois.test=test.func(order.overlap.data[order.overlap.data$original.group.index==psbl.index[jj],]$ID, var.data, N1, N0)
          burden.matrix[jj,]=c(pois.test$odds.ratio, pois.test$p.value, pois.test$rate.case, pois.test$rate.contr)
        }
        #  delta=runif(1)
        para.est=multi.group.func.for.variant(order.overlap.data, N1, N0, gamma.mean=3, sigma=2, delta=0.2, beta.init, actu.num.group)
        comb.summ.MIRAGE[(i-1)*length(gene.set)+j,((k-1)*4+1):(k*4)]=c(pois.test$odds.ratio, para.est$pvalue[length(para.est$pvalue)], pois.test$rate.case, pois.test$rate.contr)
        nn=nn+1
        MIRAGE.pvalue[[nn]]=para.est$cate.pvalue
        MIRAGE.para.est[[nn]]=para.est$beta.est
        MIRAGE.cate.index[[nn]]=psbl.index
        Burden.cate[[nn]]=burden.matrix
      }
      if (nrow(order.overlap.data)==0)  
        comb.summ.MIRAGE[(i-1)*length(gene.set)+j,((k-1)*4+1):(k*4)]=NA 
      
      
      
      if (pois.test$p.value<0.05)
      {
        jj4=jj4+1
        signal4.evid[[jj4]]=comb.evid
      }
    }
  
  if (sub.com[i] %in% var.fea==T)   # this is for variant feature 
    for (j in 1:length(gene.set))
    {    
      for (k in 1:(length(AF.cutoff)-1))  
      {
        comb.evid=intersect(intersect(var.evid[[which(sub.com[i]==var.fea)]], gene.evid[[13+j]]), var.evid[[k]]) 
        pois.test=test.func(comb.evid, var.data, N1, N0)
        comb.summ[(i-1)*length(gene.set)+j,((k-1)*4+1):(k*4)]=c(pois.test$odds.ratio, pois.test$p.value, pois.test$rate.case, pois.test$rate.contr)
        ii=ii+1
        all.evid[[ii]]=comb.evid
        
        ###################### MIRAGE burden
        cand.data=Anno.Data[which(Anno.Data$ID %in% comb.evid),]
        gene.data=eight.partition(cand.data)
        overlap.data=gene.data[gene.data$ID %in% comb.evid,]
        order.overlap.data=overlap.data[order(overlap.data$group.index, decreasing=F),]
        psbl.index=unique(order.overlap.data$group.index); actu.num.group=length(psbl.index)
        delta.init=runif(1); beta.init=runif(actu.num.group)
        order.overlap.data$original.group.index=order.overlap.data$group.index
        
        if (nrow(order.overlap.data)>0)
        {     
          burden.matrix=matrix(nrow=actu.num.group, ncol=4)  # this is burden for every category split from combined feature
          colnames(burden.matrix)=c("OR", "p.value", "rate.case", "rate.contr")
          rownames(burden.matrix)=paste(sub.com[i], gene.set[j], AF.cutoff[k], "cate", psbl.index, sep="_")
          
          for (jj in 1:actu.num.group)
          {    
            order.overlap.data$group.index[order.overlap.data$group.index==psbl.index[jj]]=jj # re-index the group labels
            pois.test=test.func(order.overlap.data[order.overlap.data$original.group.index==psbl.index[jj],]$ID, var.data, N1, N0)
            burden.matrix[jj,]=c(pois.test$odds.ratio, pois.test$p.value, pois.test$rate.case, pois.test$rate.contr)
          }
          #  delta=runif(1)
          para.est=multi.group.func.for.variant(order.overlap.data, N1, N0, gamma.mean=3, sigma=2, delta=0.2, beta.init, actu.num.group)
          comb.summ.MIRAGE[(i-1)*length(gene.set)+j,((k-1)*4+1):(k*4)]=c(pois.test$odds.ratio, para.est$pvalue[length(para.est$pvalue)], pois.test$rate.case, pois.test$rate.contr)
          nn=nn+1
          MIRAGE.pvalue[[nn]]=para.est$cate.pvalue
          MIRAGE.para.est[[nn]]=para.est$beta.est
          MIRAGE.cate.index[[nn]]=psbl.index
          Burden.cate[[nn]]=burden.matrix
        }
        if (nrow(order.overlap.data)==0)  
          comb.summ.MIRAGE[(i-1)*length(gene.set)+j,((k-1)*4+1):(k*4)]=NA
        
        
        if (pois.test$p.value<0.05)
        {
          if (k==1)
          {
            jj1=jj1+1
            signal1.evid[[jj1]]=comb.evid
          }
          if (k==2)
          {
            jj2=jj2+1
            signal2.evid[[jj2]]=comb.evid
          }
          if (k==3)
          {
            jj3=jj3+1
            signal3.evid[[jj3]]=comb.evid
          }
        }
      } 
      k=length(AF.cutoff)
      comb.evid=intersect(intersect(var.evid[[which(sub.com[i]==var.fea)]], gene.evid[[13+j]]), var.evid[[8]]) 
      pois.test=test.func(comb.evid, var.data, N1, N0)
      comb.summ[(i-1)*length(gene.set)+j,((k-1)*4+1):(k*4)]=c(pois.test$odds.ratio, pois.test$p.value, pois.test$rate.case, pois.test$rate.contr)
      ii=ii+1
      all.evid[[ii]]=comb.evid 
      
      ###################### MIRAGE burden
      cand.data=Anno.Data[which(Anno.Data$ID %in% comb.evid),]
      gene.data=eight.partition(cand.data)
      overlap.data=gene.data[gene.data$ID %in% comb.evid,]
      order.overlap.data=overlap.data[order(overlap.data$group.index, decreasing=F),]
      psbl.index=unique(order.overlap.data$group.index); actu.num.group=length(psbl.index)
      delta.init=runif(1); beta.init=runif(actu.num.group)
      order.overlap.data$original.group.index=order.overlap.data$group.index
      
      if (nrow(order.overlap.data)>0)
      {  
        burden.matrix=matrix(nrow=actu.num.group, ncol=4)  # this is burden for every category split from combined feature
        colnames(burden.matrix)=c("OR", "p.value", "rate.case", "rate.contr")
        rownames(burden.matrix)=paste(sub.com[i], gene.set[j], AF.cutoff[k], "cate", psbl.index, sep="_")
        
        for (jj in 1:actu.num.group)
        {    
          order.overlap.data$group.index[order.overlap.data$group.index==psbl.index[jj]]=jj # re-index the group labels
          pois.test=test.func(order.overlap.data[order.overlap.data$original.group.index==psbl.index[jj],]$ID, var.data, N1, N0)
          burden.matrix[jj,]=c(pois.test$odds.ratio, pois.test$p.value, pois.test$rate.case, pois.test$rate.contr)
        }
        
        para.est=multi.group.func.for.variant(order.overlap.data, N1, N0, gamma.mean=3, sigma=2, delta=0.2, beta.init, actu.num.group)
        comb.summ.MIRAGE[(i-1)*length(gene.set)+j,((k-1)*4+1):(k*4)]=c(pois.test$odds.ratio, para.est$pvalue[length(para.est$pvalue)], pois.test$rate.case, pois.test$rate.contr)
        nn=nn+1
        MIRAGE.pvalue[[nn]]=para.est$cate.pvalue
        MIRAGE.para.est[[nn]]=para.est$beta.est
        MIRAGE.cate.index[[nn]]=psbl.index
        Burden.cate[[nn]]=burden.matrix
      }
      if (nrow(order.overlap.data)==0)  
        comb.summ.MIRAGE[(i-1)*length(gene.set)+j,((k-1)*4+1):(k*4)]=NA
      
      
      
      
      if (pois.test$p.value<0.05)
      {
        jj4=jj4+1
        signal4.evid[[jj4]]=comb.evid
      }
    }
  
  if (sub.com[i]=="CriticalExon")  # this is for critical exon
    for (j in 1:length(gene.set))
    {    
      for (k in 1:(length(AF.cutoff)-1))  
      {
        comb.evid=intersect(intersect(critical.exon, gene.evid[[13+j]]), var.evid[[k]]) 
        pois.test=test.func(comb.evid, var.data, N1, N0)
        comb.summ[(i-1)*length(gene.set)+j,((k-1)*4+1):(k*4)]=c(pois.test$odds.ratio, pois.test$p.value, pois.test$rate.case, pois.test$rate.contr)
        ii=ii+1
        all.evid[[ii]]=comb.evid
        
        ###################### MIRAGE burden
        cand.data=Anno.Data[which(Anno.Data$ID %in% comb.evid),]
        gene.data=eight.partition(cand.data)
        overlap.data=gene.data[gene.data$ID %in% comb.evid,]
        order.overlap.data=overlap.data[order(overlap.data$group.index, decreasing=F),]
        psbl.index=unique(order.overlap.data$group.index); actu.num.group=length(psbl.index)
        delta.init=runif(1); beta.init=runif(actu.num.group)
        order.overlap.data$original.group.index=order.overlap.data$group.index
        
        if (nrow(order.overlap.data)>0)
        {  
          burden.matrix=matrix(nrow=actu.num.group, ncol=4)  # this is burden for every category split from combined feature
          colnames(burden.matrix)=c("OR", "p.value", "rate.case", "rate.contr")
          rownames(burden.matrix)=paste(sub.com[i], gene.set[j], AF.cutoff[k], "cate", psbl.index, sep="_")
          
          for (jj in 1:actu.num.group)
          {    
            order.overlap.data$group.index[order.overlap.data$group.index==psbl.index[jj]]=jj # re-index the group labels
            pois.test=test.func(order.overlap.data[order.overlap.data$original.group.index==psbl.index[jj],]$ID, var.data, N1, N0)
            burden.matrix[jj,]=c(pois.test$odds.ratio, pois.test$p.value, pois.test$rate.case, pois.test$rate.contr)
          }
          #  delta=runif(1)
          para.est=multi.group.func.for.variant(order.overlap.data, N1, N0, gamma.mean=3, sigma=2, delta=0.2, beta.init, actu.num.group)
          comb.summ.MIRAGE[(i-1)*length(gene.set)+j,((k-1)*4+1):(k*4)]=c(pois.test$odds.ratio, para.est$pvalue[length(para.est$pvalue)], pois.test$rate.case, pois.test$rate.contr)
          nn=nn+1
          MIRAGE.pvalue[[nn]]=para.est$cate.pvalue
          MIRAGE.para.est[[nn]]=para.est$beta.est
          MIRAGE.cate.index[[nn]]=psbl.index
          Burden.cate[[nn]]=burden.matrix
        }
        if (nrow(order.overlap.data)==0)  
          comb.summ.MIRAGE[(i-1)*length(gene.set)+j,((k-1)*4+1):(k*4)]=NA
        
        
        
        
        if (pois.test$p.value<0.05)
        {
          if (k==1)
          {
            jj1=jj1+1
            signal1.evid[[jj1]]=comb.evid
          }
          if (k==2)
          {
            jj2=jj2+1
            signal2.evid[[jj2]]=comb.evid
          }
          if (k==3)
          {
            jj3=jj3+1
            signal3.evid[[jj3]]=comb.evid
          }
        }
      } 
      k=length(AF.cutoff)
      comb.evid=intersect(intersect(critical.exon, gene.evid[[13+j]]), var.evid[[8]]) 
      pois.test=test.func(comb.evid, var.data, N1, N0)
      comb.summ[(i-1)*length(gene.set)+j,((k-1)*4+1):(k*4)]=c(pois.test$odds.ratio, pois.test$p.value, pois.test$rate.case, pois.test$rate.contr)
      ii=ii+1
      all.evid[[ii]]=comb.evid 
      
      ###################### MIRAGE burden
      cand.data=Anno.Data[which(Anno.Data$ID %in% comb.evid),]
      gene.data=eight.partition(cand.data)
      overlap.data=gene.data[gene.data$ID %in% comb.evid,]
      order.overlap.data=overlap.data[order(overlap.data$group.index, decreasing=F),]
      psbl.index=unique(order.overlap.data$group.index); actu.num.group=length(psbl.index)
      delta.init=runif(1); beta.init=runif(actu.num.group)
      order.overlap.data$original.group.index=order.overlap.data$group.index
      
      if (nrow(order.overlap.data)>0)
      {     
        burden.matrix=matrix(nrow=actu.num.group, ncol=4)  # this is burden for every category split from combined feature
        colnames(burden.matrix)=c("OR", "p.value", "rate.case", "rate.contr")
        rownames(burden.matrix)=paste(sub.com[i], gene.set[j], AF.cutoff[k], "cate", psbl.index, sep="_")
        
        for (jj in 1:actu.num.group)
        {    
          order.overlap.data$group.index[order.overlap.data$group.index==psbl.index[jj]]=jj # re-index the group labels
          pois.test=test.func(order.overlap.data[order.overlap.data$original.group.index==psbl.index[jj],]$ID, var.data, N1, N0)
          burden.matrix[jj,]=c(pois.test$odds.ratio, pois.test$p.value, pois.test$rate.case, pois.test$rate.contr)
        }
        #  delta=runif(1)
        
        para.est=multi.group.func.for.variant(order.overlap.data, N1, N0, gamma.mean=3, sigma=2, delta=0.2, beta.init, actu.num.group)
        comb.summ.MIRAGE[(i-1)*length(gene.set)+j,((k-1)*4+1):(k*4)]=c(pois.test$odds.ratio, para.est$pvalue[length(para.est$pvalue)], pois.test$rate.case, pois.test$rate.contr)
        nn=nn+1
        MIRAGE.pvalue[[nn]]=para.est$cate.pvalue
        MIRAGE.para.est[[nn]]=para.est$beta.est
        MIRAGE.cate.index[[nn]]=psbl.index
        Burden.cate[[nn]]=burden.matrix
      }
      if (nrow(order.overlap.data)==0)  
        comb.summ.MIRAGE[(i-1)*length(gene.set)+j,((k-1)*4+1):(k*4)]=NA 
      
      
      
      if (pois.test$p.value<0.05)
      {
        jj4=jj4+1
        signal4.evid[[jj4]]=comb.evid
      }
    }
  
}  
#save(comb.summ, comb.summ.MIRAGE, MIRAGE.cate.index, MIRAGE.para.est, MIRAGE.pvalue, Burden.cate, file="..\\output\\CombinedFeature\\Burden.MIRAGE.for.Combined.exon.geneset.AF.partition.old.version.RData")
