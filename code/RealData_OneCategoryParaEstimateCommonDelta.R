#############################
# this code estimates delta and beta for one category only for real data
#############################
#rm(list=ls())
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
single.group.func=function(new.data, N1, N0, gamma.mean, sigma, delta, beta.init) # new.data has one column specifying its group index
{
  #gamma.mean=3; sigma=2; delta=0.3; beta.init=0.5; num.group=1
  ########################
  max.iter=1e4
  stop.cond=0; iter=1  # parameter settings
  thrshd=1e-5
  beta.k=numeric()
  delta.est=numeric()
  delta.est[1]=delta
  beta.k[1]=beta.init
  full.info.genevar=list()
  gene.list=new.data$Gene; unique.gene=unique(gene.list) # find the gene list
  num.gene=length(unique.gene)
  BF.gene=matrix(nrow=max.iter, ncol=num.gene)
  ########################
  # calculate the Bayes factor for variant (i,j) and gene i as initials.
  for (i in 1:num.gene)
  {
    cat(i, "th gene of ", "\t", num.gene, "\t", "is running", "\n")
    var.index.list=which(gene.list==unique.gene[i])
    indi.gene=new.data[var.index.list,]
    bb=1; var.BF=numeric()
    if (length(var.index.list)>0) # calculate Bayes factor for variant (i,j)
      for (j in 1:length(var.index.list))
      {
        
        if (new.data$group.index[var.index.list[j]]<=2)
          var.BF[j]=BF.var.inte(new.data$No.case[var.index.list[j]], new.data$No.contr[var.index.list[j]], bar.gamma=6, sig=sigma, N1, N0)
        if (new.data$group.index[var.index.list[j]]>2)
          var.BF[j]=BF.var.inte(new.data$No.case[var.index.list[j]], new.data$No.contr[var.index.list[j]], bar.gamma=gamma.mean, sig=sigma, N1, N0)
        bb=bb*((1-beta.k[1])+beta.k[1]*var.BF[j])
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
    
    total.Zij=numeric(num.gene); total.Ui=numeric(num.gene)  # used to estimate beta
    
    for (i in 1:num.gene)
    {
      info.single.gene=full.info.genevar[[i]] # this is a small matrix for that single gene. each row is one variant
      bb=1
      UiZij=numeric()
      if (nrow(info.single.gene)>0)
        for (j in 1:nrow(info.single.gene))
        {
          
          numer=info.single.gene$var.BF[j]*beta.k[(iter-1)]*delta.est[iter-1]
          denom=(delta.est[iter-1]+(1-delta.est[iter-1])/BF.gene[(iter-1),i])*(beta.k[(iter-1)]*info.single.gene$var.BF[j]+(1-beta.k[(iter-1)]))
          
          UiZij[j]=numer/denom
          bb=bb*((1-beta.k[(iter-1)])+beta.k[(iter-1)]*info.single.gene$var.BF[j])
          
        }
      EUiZij[[i]]=UiZij
      BF.gene[iter,i]=bb
      EUi[i]=delta.est[iter-1]*bb/(delta.est[iter-1]*bb+1-delta.est[iter-1])
      ######################
      # Note here each gene may have multiple annotation groups
      tt=EUiZij[[i]]
      tt[is.na(tt)]=0
      
      total.Zij[i]=sum(tt)
      total.Ui[i]=sum(sum(tt>0)*EUi[i])
      
      
    }  # end of i
    ############## EM algorithm: M step
    delta.est[iter]=sum(EUi)/num.gene
    #delta.est[iter]=delta
    
    if (sum(total.Ui)!=0)
      beta.k[iter]=sum(total.Zij)/sum(total.Ui)
    if (sum(total.Ui)==0)
      beta.k[iter]=0
    
    ################
    diff=sum(abs(beta.k[iter]-beta.k[(iter-1)]))
    if (diff<thrshd || iter>(max.iter-1))
      stop.cond=1
    cat(iter, "th iteration is running", "\n")
  } # end of iter
  ##############################
  if (iter<max.iter)
    beta.k=beta.k[complete.cases(beta.k)]
  ################## calculate the likelihood ratio test statistics and p value
  # beta.k[(iter-1), -7]=0
  #delta=delta.est
  lkhd=rep(1,num.gene); total.lkhd=0
  teststat=numeric(); pvalue=numeric()
  for (i in 1:num.gene)
  {
    #    i=33
    data=full.info.genevar[[i]]
    if (nrow(data)>0)
      for (j in 1:nrow(data))
        lkhd[i]=lkhd[i]*((1-beta.k[(iter-1)])+beta.k[(iter-1)]*data$var.BF[j])
      
      teststat[i]=2*log((1-delta.est[iter-1])+delta.est[iter-1]*lkhd[i]); # this is the test statistics of one gene
      total.lkhd=total.lkhd+log((1-delta.est[iter-1])+delta.est[iter-1]*lkhd[i])
      
      pvalue[i]=pchisq(teststat[i], 2, lower.tail=F)
  } # end of i
  teststat[num.gene+1]=2*total.lkhd
  pvalue[num.gene+1]=pchisq(teststat[num.gene+1], 2, lower.tail=F)
  ##################
  cate.lkhd=0; cate.stat=numeric()
  cate.pvalue=1
  
  total.lkhd=0; lkhd.gene=rep(1, num.gene)
  for (i in 1:num.gene)
  {
    data=full.info.genevar[[i]]
    if (nrow(data)>0)
      for (j in 1:nrow(data))
      {
        lkhd.gene[i]=lkhd.gene[i]*((1-beta.k[(iter-1)])+beta.k[(iter-1)]*data$var.BF[j])
        cate.lkhd=cate.lkhd*((1-beta.k[(iter-1)])+beta.k[(iter-1)]*data$var.BF[j])
      }
    
    total.lkhd=total.lkhd+log((1-delta.est[iter-1])+delta.est[iter-1]*lkhd.gene[i])
  } # end of i
  cate.stat=2*total.lkhd
  cate.pvalue=pchisq(2*total.lkhd, 1, lower.tail=F)
  ##############################################
  # BF.total=prod(1-delta.est[iter-1]+delta.est[iter-1]*BF.gene[(iter-1),which(unique.gene %in% RiskGene)])*prod(1-delta.est[iter-1]+delta.est[iter-1]*BF.gene[(iter-1),which(unique.gene %in% ControlGene)])
  return(result=list(delta.est=delta.est[iter-1], beta.est=beta.k[(iter-1)], cate.stat=cate.stat, cate.pvalue=cate.pvalue, BF.gene=BF.gene[(iter-1),], full.info=full.info.genevar, Eui=EUi))
}
#########################################
All.Anno.Data=read.table("D:\\ResearchWork\\StatisticalGenetics\\Rare-variant-project\\AnnotatedTrans.txt", header=T)
N1=4315; N0=4315
All.Anno.Data[All.Anno.Data =="."] <- NA
All.Anno.Data$ExacAF[is.na(All.Anno.Data$ExacAF)]=0 # set AF of NA to zero
Anno.Data=All.Anno.Data[which(All.Anno.Data$ExacAF<0.05 & All.Anno.Data$Annotation!="synonymous SNV"),] # use AF cutoff and exclude synonumous SNV
var.data=data.frame(ID=Anno.Data$ID, No.case=Anno.Data$No.case, No.contr=Anno.Data$No.contr)
########################################
#gene.set=as.character(read.table("D:\\ResearchWork\\StatisticalGenetics\\Rare-variant-project\\rare-var-project\\data\\GeneSet\\TADAq.lessthan5per.gene.txt", header=T)[[1]])
#gene.set=as.character(read.csv("D:\\ResearchWork\\StatisticalGenetics\\Rare-variant-project\\rare-var-project\\data\\GeneSet\\Samocha_2014NG_contraintgene.csv", header=T)$gene)
#LoF.def=c("stopgain", "frameshift substitution", "splicing", "stoploss")
#cand.data=Anno.Data[which(Anno.Data$Gene %in% gene.set),]
#LoF.evid=which(cand.data$Annotation %in% LoF.def==T & cand.data$ExacAF<0.01)
#group.index=rep(0, nrow(cand.data))
#group.index[LoF.evid]=1
#gene.data=data.frame(ID=cand.data$ID, Gene=cand.data$Gene, No.case=cand.data$No.case, No.contr=cand.data$No.contr, group.index=group.index)
#new.new.data=gene.data[gene.data$group.index==1, ]
ID=signal1.evid[[1]]
gene.data=All.Anno.Data[All.Anno.Data$ID %in% ID, ]  
new.new.data=data.frame(ID=ID, Gene=gene.data$Gene, No.case=gene.data$No.case, No.contr=gene.data$No.contr, group.index=rep(1, length(ID)))
delta.range=seq(0.1, 1, by=0.1)
para.summary=matrix(nrow=length(delta.range), ncol=2)
#for (i in 1:length(delta.range))
{   i=1
 para.est=single.group.func(new.new.data, N1, N0, gamma.mean=3, sigma=2, delta=delta.range[i], beta.init=0.5)
 para.summary[i,1]=para.est$beta.est
 para.summary[i,2]=para.est$cate.stat
 }
 