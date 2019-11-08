#############################
rm(list=ls())
setwd("C:\\Shengtong\\Research\\rare-var\\rare-var-project\\")
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
multi.group.func=function(new.data, N1, N0, gamma.mean, sigma, delta, beta.init, num.group) # new.data has one column specifying its group index
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
  BF.gene=matrix(nrow=max.iter, ncol=num.gene)
  LoF.BF.gene=matrix(nrow=max.iter, ncol=num.gene)
  nonLoF.BF.gene=matrix(nrow=max.iter, ncol=num.gene)
  delta.est=numeric(); delta.est[1]=delta

  ########################
  # calculate the Bayes factor for variant (i,j) and gene i as initials.
  for (i in 1:num.gene)
  {
    cat(i, "th gene of ", "\t", num.gene, "\t", "is running", "\n")
    var.index.list=which(gene.list==unique.gene[i])
    indi.gene=new.data[var.index.list,]
    bb=1; var.BF=numeric()
    bb.LoF=1; bb.nonLoF=1
    if (length(var.index.list)>0) # calculate Bayes factor for variant (i,j)
      for (j in 1:length(var.index.list))
      {

        if (new.data$group.index[var.index.list[j]]<=5)
          var.BF[j]=BF.var.inte(new.data$No.case[var.index.list[j]], new.data$No.contr[var.index.list[j]], bar.gamma=6, sig=sigma, N1, N0)
        if (new.data$group.index[var.index.list[j]]>5)
          var.BF[j]=BF.var.inte(new.data$No.case[var.index.list[j]], new.data$No.contr[var.index.list[j]], bar.gamma=gamma.mean, sig=sigma, N1, N0)
        bb=bb*((1-beta.k[1, new.data$group.index[var.index.list[j]]])+beta.k[1, new.data$group.index[var.index.list[j]]]*var.BF[j])
        ################## split BF of LoF  and non LoF
        if (new.data$group.index[var.index.list[j]]<=2)
          bb.LoF=bb.LoF*((1-beta.k[1, new.data$group.index[var.index.list[j]]])+beta.k[1, new.data$group.index[var.index.list[j]]]*var.BF[j])
        if (new.data$group.index[var.index.list[j]]>2)
          bb.nonLoF=bb.nonLoF*((1-beta.k[1, new.data$group.index[var.index.list[j]]])+beta.k[1, new.data$group.index[var.index.list[j]]]*var.BF[j])

      }
    full.info.genevar[[i]]=cbind(indi.gene, var.BF)
    BF.gene[1, i]=bb
    LoF.BF.gene[1,i]=bb.LoF
    nonLoF.BF.gene[1,i]=bb.nonLoF
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
      bb=1; bb.LoF=1; bb.nonLoF=1
      UiZij=numeric()
      if (nrow(info.single.gene)>0)
        for (j in 1:nrow(info.single.gene))
        {

            numer=info.single.gene$var.BF[j]*beta.k[(iter-1), info.single.gene$group.index[j]]*delta.est[iter-1]
            denom=(delta.est[iter-1]+(1-delta.est[iter-1])/BF.gene[(iter-1),i])*(beta.k[(iter-1), info.single.gene$group.index[j]]*info.single.gene$var.BF[j]
                                                               +(1-beta.k[(iter-1), info.single.gene$group.index[j]]))
          UiZij[j]=numer/denom
          bb=bb*((1-beta.k[(iter-1), info.single.gene$group.index[j]])+beta.k[(iter-1), info.single.gene$group.index[j]]*info.single.gene$var.BF[j])

          ########################### split into LoF and non-LoF two parts
          if (info.single.gene$group.index[j]<=2)
            bb.LoF=bb.LoF*((1-beta.k[(iter-1), info.single.gene$group.index[j]])+beta.k[(iter-1), info.single.gene$group.index[j]]*info.single.gene$var.BF[j])
          if (info.single.gene$group.index[j]>2)
            bb.nonLoF=bb.nonLoF*((1-beta.k[(iter-1), info.single.gene$group.index[j]])+beta.k[(iter-1), info.single.gene$group.index[j]]*info.single.gene$var.BF[j])

        }
      EUiZij[[i]]=UiZij
      BF.gene[iter,i]=bb
      LoF.BF.gene[iter,i]=bb.LoF
      nonLoF.BF.gene[iter,i]=bb.nonLoF
      EUi[i]=delta.est[iter-1]*bb/(delta.est[iter-1]*bb+1-delta.est[iter-1])
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
      delta.est[iter]=sum(EUi)/num.gene
    #  delta.est[iter]=delta
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
 #   cat(iter, "th iteration is running", "\n")
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

        teststat[i]=2*log((1-delta.est[iter-1])+delta.est[iter-1]*lkhd[i]); # this is the test statistics of one gene
        total.lkhd=total.lkhd+log((1-delta.est[iter-1])+delta.est[iter-1]*lkhd[i])

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

        total.lkhd=total.lkhd+log((1-delta.est[iter-1])+delta.est[iter-1]*lkhd.gene[i])
    } # end of i
    cate.stat[g]=2*total.lkhd
    cate.pvalue[g]=pchisq(2*total.lkhd, 1, lower.tail=F)
  } # end of g
  sum.lkhd=sum(cate.stat)
  ##############################################
  return(result=list(delta.est=delta.est[iter-1], delta.pvalue=pvalue[length(pvalue)], beta.est=beta.k[(iter-1),], beta.stat=cate.stat, beta.pvalue=cate.pvalue, BF.gene=data.frame(Gene=unique.gene, BF=BF.gene[(iter-1),], LoF.BF=LoF.BF.gene[(iter-1),], nonLoF.BF=nonLoF.BF.gene[(iter-1),]), full.info=full.info.genevar, Eui=EUi))
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
#All.Anno.Data=read.table("C:\\han\\ResearchWork\\StatGene\\AutismData\\AnnotatedTrans.txt", header=T)
All.Anno.Data=read.table("..\\AnnotatedTrans.txt", header=T)
#All.Anno.Data=read.table("C:\\han\\ResearchWork\\StatGene\\SCZData\\AnnotatedSCZ.txt", header=T)
N1=4315; N0=4315
#N1=2536; N0=2543
All.Anno.Data[All.Anno.Data =="."] <- NA
All.Anno.Data$ExacAF[is.na(All.Anno.Data$ExacAF)]=0 # set AF of NA to zero
Anno.Data=All.Anno.Data[which(All.Anno.Data$ExacAF<0.05 & All.Anno.Data$Annotation!="synonymous SNV"),] # use AF cutoff and exclude synonumous SNV
var.data=data.frame(ID=Anno.Data$ID, No.case=Anno.Data$No.case, No.contr=Anno.Data$No.contr)
CADD.cutoff=quantile(as.numeric(as.character(All.Anno.Data$CADD.raw)), prob=0.9, na.rm=TRUE)
########################################
################ whole genome
#gene.set=as.character(unique(All.Anno.Data$Gene))
gene.set=as.character(read.csv("data\\GeneSet\\Samocha_2014NG_contraintgene.csv", header=T)$gene)
#gene.set=as.character(read.csv("C:\\Users\\han\\Dropbox\\StatisticalGenetics\\Samocha_2014NG_contraintgene.csv", header=T)$gene)
vart.set=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% gene.set)])
cand.data=Anno.Data[which(Anno.Data$ID %in% vart.set),]
gene.data=eight.partition(cand.data)

  overlap.id=gene.data$ID
  if (length(overlap.id)>0)
  {
    overlap.data=gene.data[gene.data$ID %in% overlap.id,]
    order.overlap.data=overlap.data[order(overlap.data$group.index, decreasing=F),]
    psbl.index=unique(order.overlap.data$group.index); actu.num.group=length(psbl.index)
    delta.init=runif(1); beta.init=runif(actu.num.group)
    for (j in 1:actu.num.group)
      order.overlap.data$group.index[order.overlap.data$group.index==psbl.index[j]]=j # re-index the group labels
  #  delta=runif(1)

    ###########################################
    permutate.BF.gene=list()
    num.permu=100
    for (ii in 1:num.permu)
    {
       cat(ii, "th permutation of ", num.permu, "is running", "\n")
    ###################### permutate the data: re-distribute the variants in cases and controls  ##############
    permutate.order.overlap.data=order.overlap.data

   # for (i in 1:10)
     for ( i in 1:nrow(order.overlap.data))
     {
       total.count=order.overlap.data[i,]$No.case+order.overlap.data[i,]$No.contr
       if (total.count>0)
         {
          permutate.order.overlap.data[i,]$No.case=sum(rbinom(total.count, 1, 0.5))
          permutate.order.overlap.data[i,]$No.contr=total.count-permutate.order.overlap.data[i,]$No.case
       }
       if (total.count==0)
         permutate.order.overlap.data[i,]=order.overlap.data[i,]
     }

    ###########################################################################################################
    #para.est=multi.group.func(order.overlap.data, N1, N0, gamma.mean=3, sigma=2, delta=0.2, beta.init, actu.num.group)
     para.est=multi.group.func(permutate.order.overlap.data, N1, N0, gamma.mean=3, sigma=2, delta=0.2, beta.init, actu.num.group)
     permutate.BF.gene[[ii]]=para.est$BF.gene
    } # end of ii

#    cutoff=1
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
##########################
############# use other methods to find risk genes #################
#cand.data=gene.data[gene.data$group.index<6,]  # removing non-damaging variants  why?????????

gene_var_data=as_tibble(cand.data %>% select(ID, Gene, No.case, No.contr))
gene=levels(factor(gene_var_data$Gene))

SKATO.pvalue=numeric()
CMC.pvalue=numeric()
ASUM.pvalue=numeric()
Fisher.pvalue=numeric()
Fisher.adj.pvalue=numeric()
CALPHA.pvalue=numeric()
for (i in 1:length(gene))
#for (i in 1:1000)
{
 # i=1
  cat(i, "is running", "\n")
  single_gene_var=gene_var_data %>% filter(Gene==gene[i])
  var_in_gene=levels(factor((single_gene_var%>% select(ID))$ID))
  No_var=length(var_in_gene)
  geno=matrix(0,nrow=(N1+N0), ncol=No_var)
  for (j in 1:No_var)
  {
    if (single_gene_var[j,3]>0)
    {
      case_var_pos=sample(seq(1,N1),as.numeric(single_gene_var[j,3]), replace=F)
      geno[case_var_pos,j]=1
    }
    if (single_gene_var[j,4]>0)
    {
     contr_var_pos=N1+sample(seq(1,N0),as.numeric(single_gene_var[j,4]), replace=F)
     geno[contr_var_pos,j]=1
    }

  }  # end of j

  pheno_geno=list(pheno=c(rep(1,N1), rep(0, N0)), geno=geno)

  ############ SKATO
  library(SKAT)
  library(AssotesteR)
  detach("package:AssotesteR", unload=TRUE)
  obj<-SKAT_Null_Model(pheno_geno$pheno~ 1, out_type="D")  # without covariates
  if (ncol(pheno_geno$geno)>1)
   SKATO.pvalue[i]=SKAT(pheno_geno$geno, obj, method="SKATO")$p.value
   if (ncol(pheno_geno$geno)==1)
     SKATO.pvalue[i]=NA

  library(AssotesteR)
  ############ CALPHA
  if (ncol(pheno_geno$geno)>1)
  {
    c.alpha.test=CALPHA(pheno_geno$pheno, pheno_geno$geno,  perm=NULL)
    CALPHA.pvalue[i]=c.alpha.test$asym.pval
  }
  if (ncol(pheno_geno$geno)==1)
    CALPHA.pvalue[i]=NA

  ############ CMC
  if (ncol(pheno_geno$geno)>1)
  {
    cmc.test=CMC(pheno_geno$pheno, pheno_geno$geno, maf=10^(-4)*c(0.05, 0.2, 0.5), perm=100)
    CMC.pvalue[i]=cmc.test$asym.pval
  }
  if (ncol(pheno_geno$geno)==1)
    CMC.pvalue[i]=NA

  ############  ASUM test
  if (ncol(pheno_geno$geno)>1)
  {
    ASUM.test=ASUM(pheno_geno$pheno, pheno_geno$geno, perm=100)
    if (ASUM.test$asum.stat>0)
      ASUM.pvalue[i]=ASUM.test$perm.pval
    if (ASUM.test$asum.stat<0)
      ASUM.pvalue[i]=1

  }
  if (ncol(pheno_geno$geno)==1)
    ASUM.pvalue[i]=NA


  ########### Burden test
  cont_matrix=matrix(nrow=2, ncol=2)
  cont_matrix[1,]=c(sum(single_gene_var$No.case), sum(single_gene_var$No.contr))
  cont_matrix[2,]=c(N1, N0)
  Fisher.pvalue[i]=fisher.test(cont_matrix)$p.value

  ########## Burden test for every category
  single_gene_cate=gene.data[gene.data$Gene==gene[i],]
  var_cate=unique(single_gene_cate$group.index)
  cate_pvalue=numeric()
  for (j in 1:length(var_cate))
  {
    cate_data=single_gene_cate[single_gene_cate$group.index==var_cate[j],]
    cate_pvalue[j]=fisher.test(matrix(c(sum(cate_data$No.case), N1, sum(cate_data$No.contr), N0), nrow=2))$p.value
  }
  Fisher.adj.pvalue[i]=min(p.adjust(cate_pvalue, method = "bonferroni", n = length(cate_pvalue)))

}  # end of j


all_pvalue=data.frame(Gene=gene,SKATO.pvalue=SKATO.pvalue, CMC.pvalue=CMC.pvalue, ASUM.pvalue=ASUM.pvalue, Fisher.pvalue=Fisher.pvalue, Fisher.adj.pvalue=Fisher.adj.pvalue,
                      CALPHA.pvalue=CALPHA.pvalue)
#all_pvalue[is.na(all_pvalue)] <- 1
all_pvalue_adj=all_pvalue
for (i in 1:ncol(all_pvalue))
  all_pvalue_adj[,i]=p.adjust(all_pvalue[,i], method="BH", n=length(all_pvalue[,i]))
fdr.level=0.05
cbind(cons_post, all_pvalue_adj)

