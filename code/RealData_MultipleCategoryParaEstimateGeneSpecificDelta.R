#############################
rm(list=ls())
library(RSQLite)
library(dplyr)
library(knitr)
set.seed(123)
###################################
intergrand=function(aa, var.case, var.contr, bar.gamma, sig, N1, N0)
{
  ff=dbinom(var.case, sum(var.case, var.contr), aa*N1/(aa*N1+N0))*dgamma(aa, bar.gamma*sig, sig)
  return(ff)
}
# calculate the bayes factor of a single variant via integration
BF.var.inte=function(var.case, var.contr, bar.gamma, sig, N1, N0)
{
  marglik0.CC <- dbinom(var.case, sum(var.case, var.contr), N1/(N1+N0))    # Under H0: gamma=1

  marglik1.CC <- integrate(intergrand, var.case, var.contr, bar.gamma, sig, N1, N0, low=0, upper=100, stop.on.error=F)$value # Under H1: gamma~gamma(gamma.mean*sigma, sigma)
  BF.var <- marglik1.CC/marglik0.CC

  return(BF.var)
}
######################################################
multi.group.func=function(new.data, N1, N0, gamma.mean, sigma,beta.init, num.group) # new.data has one column specifying its group index
{
  ########################
  max.iter=1e4
  stop.cond=0; iter=1  # parameter settings
  thrshd=1e-5
  beta.k=matrix(nrow=max.iter, ncol=num.group)
  beta.k[1,]=beta.init
  full.info.genevar=list()
  gene.list=new.data$Gene; unique.gene=unique(gene.list) # find the gene list
  num.gene=length(unique.gene)
  gene.prior=numeric()
  BF.gene=matrix(nrow=max.iter, ncol=num.gene)
  ########################
  # calculate the Bayes factor for variant (i,j) and gene i as initials.
  for (i in 1:num.gene)
  {
    cat(i, "th gene of ", "\t", num.gene, "\t", "is running", "\n")
    var.index.list=which(gene.list==unique.gene[i])
    indi.gene=new.data[var.index.list,]
    bb=1; var.BF=numeric()
    gene.prior[i]=0
    if (length(var.index.list)>0) # calculate Bayes factor for variant (i,j)
     {
      for (j in 1:length(var.index.list))
      {

        if (new.data$group.index[var.index.list[j]]<=5)
          var.BF[j]=BF.var.inte(new.data$No.case[var.index.list[j]], new.data$No.contr[var.index.list[j]], bar.gamma=6, sig=sigma, N1, N0)
        if (new.data$group.index[var.index.list[j]]>5)
          var.BF[j]=BF.var.inte(new.data$No.case[var.index.list[j]], new.data$No.contr[var.index.list[j]], bar.gamma=gamma.mean, sig=sigma, N1, N0)
        bb=bb*((1-beta.k[1, new.data$group.index[var.index.list[j]]])+beta.k[1, new.data$group.index[var.index.list[j]]]*var.BF[j])
      }
      gene.prior[i]=new.data$gene.prior[var.index.list[1]]
     }
    full.info.genevar[[i]]=cbind(indi.gene, var.BF)
    BF.gene[1, i]=bb
    
  }
  ########################## EM algorithm
  ########################
  while (stop.cond==0)
  {
    iter=iter+1
    ############## EM algorithm: E step
    EUiZij=list() # expectation for variant (i,j), every gene may have varying number of variant
    EUi=numeric() # expectation for gene i.

    total.Zij=matrix(nrow=num.gene, ncol=num.group); total.Ui=matrix(nrow=num.gene, ncol=num.group)  # used to estimate beta

    for (i in 1:num.gene)
    {
      info.single.gene=full.info.genevar[[i]] # this is a small matrix for that single gene. each row is one variant
      bb=1
      UiZij=numeric()
      if (nrow(info.single.gene)>0)
        for (j in 1:nrow(info.single.gene))
        {

            numer=info.single.gene$var.BF[j]*beta.k[(iter-1), info.single.gene$group.index[j]]*gene.prior[i]
            denom=(gene.prior[i]+(1-gene.prior[i])/BF.gene[(iter-1),i])*(beta.k[(iter-1), info.single.gene$group.index[j]]*info.single.gene$var.BF[j]
                                                               +(1-beta.k[(iter-1), info.single.gene$group.index[j]]))
          UiZij[j]=numer/denom
          bb=bb*((1-beta.k[(iter-1), info.single.gene$group.index[j]])+beta.k[(iter-1), info.single.gene$group.index[j]]*info.single.gene$var.BF[j])

        }
      EUiZij[[i]]=UiZij
      BF.gene[iter,i]=bb
      EUi[i]=gene.prior[i]*bb/(gene.prior[i]*bb+1-gene.prior[i])
      ######################
      # Note here each gene may have multiple annotation groups
      tt=EUiZij[[i]]
      tt[is.na(tt)]=0
      for (g in 1:num.group)
      {
        total.Zij[i, g]=sum(tt[which(info.single.gene$group.index==g)])
        total.Ui[i, g]=sum(sum(tt[which(info.single.gene$group.index==g)]>0)*EUi[i])
      }

    }  # end of i
    ############## EM algorithm: M step
    for (g in 1:num.group)
    {
      if (sum(total.Ui[,g])!=0)
        beta.k[iter, g]=sum(total.Zij[,g])/sum(total.Ui[,g])
      if (sum(total.Ui[,g])==0)
        beta.k[iter, g]=0
    }
    ################
    if (num.group>1)
      diff=sum(abs(beta.k[iter,]-beta.k[(iter-1),]))
    if (diff<thrshd || iter>(max.iter-1))
      stop.cond=1
    cat(iter, "th iteration is running", "\n")
  } # end of iter
  ##############################
  if (iter<max.iter)
  {
    if (num.group>1)
      beta.k=beta.k[complete.cases(beta.k),]
    if (num.group==1)
      beta.k=beta.k[complete.cases(beta.k)]
  }
  ################## calculate the likelihood ratio test statistics and p value
  # beta.k[(iter-1), -7]=0
  lkhd=rep(1,num.gene); total.lkhd=0
  teststat=numeric(); pvalue=numeric()
  for (i in 1:num.gene)
  {
    #    i=33
    data=full.info.genevar[[i]]
    if (nrow(data)>0)
      for (j in 1:nrow(data))
        lkhd[i]=lkhd[i]*((1-beta.k[(iter-1), data$group.index[j]])+beta.k[(iter-1), data$group.index[j]]*data$var.BF[j])

        teststat[i]=2*log((1-gene.prior[i])+gene.prior[i]*lkhd[i]); # this is the test statistics of one gene
        total.lkhd=total.lkhd+log((1-gene.prior[i])+gene.prior[i]*lkhd[i])

      pvalue[i]=pchisq(teststat[i], 2, lower.tail=F)
  } # end of i
  teststat[num.gene+1]=2*total.lkhd
  pvalue[num.gene+1]=pchisq(teststat[num.gene+1], 2, lower.tail=F)
  ##################
  cate.lkhd=rep(1,num.group); cate.stat=numeric()
  cate.pvalue=numeric(num.group); sum.lkhd=0
  for (g in 1:num.group)
  { # g=2
    total.lkhd=0; lkhd.gene=rep(1, num.gene)
    for (i in 1:num.gene)
    {
      data=full.info.genevar[[i]]
      if (nrow(data)>0)
        for (j in 1:nrow(data))
          if (data$group.index[j]==g)
          {
            lkhd.gene[i]=lkhd.gene[i]*((1-beta.k[(iter-1), g])+beta.k[(iter-1), g]*data$var.BF[j])
            cate.lkhd[g]=cate.lkhd[g]*((1-beta.k[(iter-1), g])+beta.k[(iter-1), g]*data$var.BF[j])
          }

        total.lkhd=total.lkhd+log((1-gene.prior[i])+gene.prior[i]*lkhd.gene[i])
    } # end of i
    cate.stat[g]=2*total.lkhd
    cate.pvalue[g]=pchisq(2*total.lkhd, 1, lower.tail=F)
  } # end of g
  sum.lkhd=sum(cate.stat)
  ##############################################
  return(result=list(gene.prior=gene.prior, beta.est=beta.k[(iter-1),], beta.stat=cate.stat, cate.pvalue=cate.pvalue, BF.gene=data.frame(Gene=unique.gene, BF=BF.gene[(iter-1),]), full.info=full.info.genevar, Eui=EUi))
}
#####################################
fixed.beta.func=function(new.data, N1, N0, gamma.mean, sigma, beta) # new.data has one column specifying its group index
{
  #######################
  beta.k=beta
  full.info.genevar=list()
  gene.list=new.data$Gene; unique.gene=unique(gene.list) # find the gene list
  num.gene=length(unique.gene)
  BF.gene=numeric()
  ########################
  for (i in 1:num.gene)
  {
    cat(i, "th gene of ", "\t", num.gene, "\t", "is running", "\n")
    var.index.list=which(gene.list==unique.gene[i])
    indi.gene=new.data[var.index.list,]   # note var.index.list matches new.data
    bb=1; var.BF=numeric()
    if (length(var.index.list)>0) # calculate Bayes factor for variant (i,j)
      for (j in 1:length(var.index.list))
      {
        if (new.data$group.index[var.index.list[j]]<=5)
          var.BF[j]=BF.var.inte(new.data$No.case[var.index.list[j]], new.data$No.contr[var.index.list[j]], bar.gamma=6, sig=sigma, N1, N0)
        if (new.data$group.index[var.index.list[j]]>5)
          var.BF[j]=BF.var.inte(new.data$No.case[var.index.list[j]], new.data$No.contr[var.index.list[j]], bar.gamma=gamma.mean, sig=sigma, N1, N0)
        bb=bb*((1-beta.k[new.data$group.index[var.index.list[j]]])+beta.k[new.data$group.index[var.index.list[j]]]*var.BF[j])
      }
    full.info.genevar[[i]]=cbind(indi.gene, var.BF)
    BF.gene[i]=bb
  }

  return(result=list(BayesFactor=data.frame(Gene=unique.gene, BF=BF.gene), full.info=full.info.genevar))
}

#########################################
fifteen.partition=function(cand.data) # given gene data and annotations, do variant partitions
{
  par.evid=list()
  LoF.def=c("stopgain", "frameshift substitution", "splicing", "stoploss")
  par.evid[[1]]=which(cand.data$Annotation %in% LoF.def==T & cand.data$ExacAF<0.05 & cand.data$ExacAF>=0.01 )
  par.evid[[2]]=which(cand.data$Annotation %in% LoF.def==T &  cand.data$ExacAF<0.01 & cand.data$ExacAF>=0.001)
  par.evid[[3]]=which(cand.data$Annotation %in% LoF.def==T &  cand.data$ExacAF<0.001 & cand.data$ExacAF>=0.0001)
  par.evid[[4]]=which(cand.data$Annotation %in% LoF.def==T &  cand.data$ExacAF<0.0001 & cand.data$ExacAF>0)
  par.evid[[5]]=which(cand.data$Annotation %in% LoF.def==T &  cand.data$ExacAF==0)


  par.evid[[6]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))>=0.957 & cand.data$ExacAF>=0.01 & cand.data$ExacAF<0.05)
  par.evid[[7]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))>=0.957 & cand.data$ExacAF>=0.001 & cand.data$ExacAF<0.01)
  par.evid[[8]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))>=0.957 & cand.data$ExacAF>=0.0001 & cand.data$ExacAF<0.001)
  par.evid[[9]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))>=0.957 & cand.data$ExacAF>0 & cand.data$ExacAF<0.0001)
  par.evid[[10]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))>=0.957 & cand.data$ExacAF==0 )


  par.evid[[11]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))<0.957 & cand.data$ExacAF>=0.01 & cand.data$ExacAF<0.05)
  par.evid[[12]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))<0.957 & cand.data$ExacAF>=0.001 & cand.data$ExacAF<0.01)
  par.evid[[13]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))<0.957 & cand.data$ExacAF>=0.0001 & cand.data$ExacAF<0.001)
  par.evid[[14]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))<0.957 & cand.data$ExacAF>0 & cand.data$ExacAF<0.0001)
  par.evid[[15]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))<0.957 & cand.data$ExacAF==0)

  group.index=rep(NA, nrow(cand.data))
  for (i in 1:length(par.evid))
    group.index[par.evid[[i]]]=i
  gene.data=data.frame(ID=cand.data$ID, Gene=cand.data$Gene, No.case=cand.data$No.case, No.contr=cand.data$No.contr, group.index=group.index)
  gene.data=gene.data[complete.cases(gene.data),]
  return(gene.data)
}

eight.partition=function(cand.data) # given gene data and annotations, do variant partitions
{
  par.evid=list()
  LoF.def=c("stopgain", "frameshift substitution", "splicing", "stoploss")
  par.evid[[1]]=which(cand.data$Annotation %in% LoF.def==T & cand.data$ExacAF<0.05 & cand.data$ExacAF>=0.01 )
  par.evid[[2]]=which(cand.data$Annotation %in% LoF.def==T &  cand.data$ExacAF<0.01)

  par.evid[[3]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))>=0.957 & cand.data$ExacAF>=0.01 & cand.data$ExacAF<0.05)
  par.evid[[4]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))>=0.957 & cand.data$ExacAF>=0.001 & cand.data$ExacAF<0.01)
  par.evid[[5]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))>=0.957 &  cand.data$ExacAF<0.001)

  par.evid[[6]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))<0.957 & cand.data$ExacAF>=0.01 & cand.data$ExacAF<0.05)
  par.evid[[7]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))<0.957 & cand.data$ExacAF>=0.001 & cand.data$ExacAF<0.01)
  par.evid[[8]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))<0.957 & cand.data$ExacAF<0.001)

  group.index=rep(NA, nrow(cand.data))
  for (i in 1:length(par.evid))
    group.index[par.evid[[i]]]=i
  gene.data=data.frame(ID=cand.data$ID, Gene=cand.data$Gene, No.case=cand.data$No.case, No.contr=cand.data$No.contr, group.index=group.index)
  gene.data=gene.data[complete.cases(gene.data),]
  return(gene.data)
}
#########################################
#All.Anno.Data=read.table("D:\\ResearchWork\\StatisticalGenetics\\Rare-variant-project\\AnnotatedTrans.txt", header=T)
All.Anno.Data=read.table("C:\\han\\ResearchWork\\StatGene\\AutismData\\AnnotatedTrans.txt", header=T)
N1=4315; N0=4315
All.Anno.Data[All.Anno.Data =="."] <- NA
All.Anno.Data$ExacAF[is.na(All.Anno.Data$ExacAF)]=0 # set AF of NA to zero
Anno.Data=All.Anno.Data[which(All.Anno.Data$ExacAF<0.05 & All.Anno.Data$Annotation!="synonymous SNV"),] # use AF cutoff and exclude synonumous SNV
var.data=data.frame(ID=Anno.Data$ID, No.case=Anno.Data$No.case, No.contr=Anno.Data$No.contr)
########################################
#gene.set=as.character(read.csv("D:\\ResearchWork\\StatisticalGenetics\\Rare-variant-project\\rare-var-project\\data\\GeneSet\\Samocha_2014NG_contraintgene.csv", header=T)$gene)
#gene.set=as.character(read.csv("C:\\Users\\han\\Dropbox\\StatisticalGenetics\\Samocha_2014NG_contraintgene.csv", header=T)$gene)
gene.set=as.character(unique(All.Anno.Data$Gene))
#gene.set=as.character(read.table("C:\\Users\\han\\Dropbox\\StatisticalGenetics\\RVISGene\\RVIS.quantilelessthan50.gene.txt", header=T)[[1]])
gene.FDR=read.table("C:\\han\\ResearchWork\\StatGene\\170726_to_Shengtong_gene_FDR_based_on_Sanders_Neuron_denovo_coding_SNV.txt", header=T)  # use 1-FDR as prior
#gene.FDR=read.table("D:\\ResearchWork\\StatisticalGenetics\\Rare-variant-project\\170726_to_Shengtong_gene_FDR_based_on_Sanders_Neuron_denovo_coding_SNV.txt", header=T)  # use 1-FDR as prior
vart.set=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% gene.set)])
cand.data=Anno.Data[which(Anno.Data$ID %in% vart.set),]
gene.data=eight.partition(cand.data)

gene.prior=numeric()                           # add gene priors in the last column
for (i in 1:nrow(gene.data))
 if (length(which(gene.FDR$genename==as.character(gene.data$Gene)[i])))
   gene.prior[i]=1-gene.FDR$FDR[which(gene.FDR$genename==as.character(gene.data$Gene)[i])]
gene.prior[is.na(gene.prior)]=0

new.data=cbind(gene.data, gene.prior)

{ #i=1
  overlap.id=new.data$ID
  if (length(overlap.id)>0)
  {
    overlap.data=new.data[new.data$ID %in% overlap.id,]
    order.overlap.data=overlap.data[order(overlap.data$group.index, decreasing=F),]
    psbl.index=unique(order.overlap.data$group.index); actu.num.group=length(psbl.index)
    beta.init=runif(actu.num.group)
    for (j in 1:actu.num.group)
      order.overlap.data$group.index[order.overlap.data$group.index==psbl.index[j]]=j # re-index the group labels
    para.est=multi.group.func(order.overlap.data, N1, N0, gamma.mean=3, sigma=2, beta.init, actu.num.group)
    cutoff=1
#    for (j in 1:actu.num.group)
#    {
#      if (para.est$beta.stat[j]<1)
#        para.est$beta.est[j]=0
#      pASD.cons.para.est[(1+psbl.index[j]),((i-1)*2+1)]=para.est$beta.est[j]
#      pASD.cons.para.est[(1+psbl.index[j]),((i-1)*2+2)]=para.est$beta.stat[j]

#    }
    ########################## after category selection, enrichment of top 20 genes with denovo genes
#    beta.est=pASD.cons.para.est[2:16,((i-1)*2+1)]
#    beta.stat=pASD.cons.para.est[2:16,((i-1)*2+2)]
#    beta.stat[is.na(beta.stat)]=0
#    beta.stat[which(beta.stat<1)]=0
#    beta.est[which(beta.stat==0)]=0
#    BF.calculation=fixed.beta.func(order.overlap.data, N1, N0, gamma.mean=3, sigma=2, beta=beta.est)
#    BF.gene=BF.calculation[[1]]
#    order.BF.gene=BF.gene[order(BF.gene$BF, decreasing=T),]
#    topgene=as.character(order.BF.gene$Gene[1:20]); totalgene=order.BF.gene$Gene
#    count=matrix(nrow=2, ncol=2)
#    count[1,1]=length(unique(intersect(top1000denovogene, topgene))); count[1,2]=20
#    count[2,1]=length(unique(intersect(top1000denovogene, totalgene))); count[2,2]=length(totalgene)
#    burden.summ[i, 1]=fisher.test(count)$estimate; burden.summ[i, 2]=fisher.test(count)$p.value
#    burden.summ[i, 3]=count[1,1]; burden.summ[i, 4]=count[2,1]
    ##########################
  }
#  if (length(overlap.id)==0)
#    pASD.cons.para.est[,i]=NA
}
##########################

BF.gene=para.est$BF.gene
order.BF.gene=BF.gene[order(BF.gene$BF, decreasing=T),]
top20gene=as.character(order.BF.gene$Gene)[1:100]

denovo=read.table("C:\\Users\\han\\Dropbox\\StatisticalGenetics\\TADA_SNV_CNV_combined_Feb7.txt", header=T)
order.denovo=denovo[order(denovo$qvalue.combined, decreasing=F),]
order.denovo=denovo[order(denovo$BF.SNV.dn, decreasing=T),]
top1000.denovo.gene=as.character(order.denovo$RefSeqName[1:1000])
top20.overlap=length(intersect(top20gene, top1000.denovo.gene))
all.overlap=length(intersect(gene.set, top1000.denovo.gene))
top20.overlap
all.overlap

ASD.prior=read.table("C:\\Users\\han\\Dropbox\\StatisticalGenetics\\ASD_Priors.txt", header=T)
order.ASD.prior=ASD.prior[order(ASD.prior$post, decreasing=T),]
top1000.ASD.prior.gene=as.character(order.ASD.prior$Gene[1:1000])
top20.overlap=length(intersect(top20gene, top1000.ASD.prior.gene))
all.overlap=length(intersect(gene.set, top1000.ASD.prior.gene))
top20.overlap
all.overlap