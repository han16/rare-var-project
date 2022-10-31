
```{r, echo=F, warning=F, message=F}
rm(list=ls())
set.seed(123)
library(knitr)
library(RSQLite)
library(dplyr)
library(knitr)
library(kableExtra)
library(RColorBrewer)
library(gplots)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(DT)
#library("devtools")
#install.packages("Rtools")
#install_github('xinhe-lab/mirage')
library(mirage)
```
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
        
        if (new.data$group.index[var.index.list[j]]<=6)
          var.BF[j]=BF.var.inte(new.data$No.case[var.index.list[j]], new.data$No.contr[var.index.list[j]], bar.gamma=6, sig=sigma, N1, N0)
        if (new.data$group.index[var.index.list[j]]>6)
          var.BF[j]=BF.var.inte(new.data$No.case[var.index.list[j]], new.data$No.contr[var.index.list[j]], bar.gamma=gamma.mean, sig=sigma, N1, N0)
        bb=bb*((1-beta.k[1, new.data$group.index[var.index.list[j]]])+beta.k[1, new.data$group.index[var.index.list[j]]]*var.BF[j])
        ################## split BF of LoF  and non LoF
      #  if (new.data$group.index[var.index.list[j]]<=2)
      #    bb.LoF=bb.LoF*((1-beta.k[1, new.data$group.index[var.index.list[j]]])+beta.k[1, new.data$group.index[var.index.list[j]]]*var.BF[j])
      #  if (new.data$group.index[var.index.list[j]]>2)
      #    bb.nonLoF=bb.nonLoF*((1-beta.k[1, new.data$group.index[var.index.list[j]]])+beta.k[1, new.data$group.index[var.index.list[j]]]*var.BF[j])
        
      }
    full.info.genevar[[i]]=cbind(indi.gene, var.BF)
    bb=ifelse(bb==Inf, 3*10^300, bb) # set the upper limit when overflow
    BF.gene[1, i]=bb
 #   LoF.BF.gene[1,i]=bb.LoF
  #  nonLoF.BF.gene[1,i]=bb.nonLoF
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
      bb=1; bb.LoF=1; #bb.nonLoF=1
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
        #  if (info.single.gene$group.index[j]<=2)
        #    bb.LoF=bb.LoF*((1-beta.k[(iter-1), info.single.gene$group.index[j]])+beta.k[(iter-1), info.single.gene$group.index[j]]*info.single.gene$var.BF[j])
        #  if (info.single.gene$group.index[j]>2)
        #    bb.nonLoF=bb.nonLoF*((1-beta.k[(iter-1), info.single.gene$group.index[j]])+beta.k[(iter-1), info.single.gene$group.index[j]]*info.single.gene$var.BF[j])
          
        }
      EUiZij[[i]]=UiZij
      bb=ifelse(bb==Inf, 3*10^300, bb) # set the upper limit when overflow
      BF.gene[iter,i]=bb
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
   # delta.est[iter]=sum(EUi)/num.gene
      delta.est[iter]=delta # use fixed proportion of risk genes 
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
       cat("distance between iterations is ", diff, "\n")
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
  return(result=list(total.lkhd=total.lkhd, delta.est=delta.est[iter-1], delta.pvalue=pvalue[length(pvalue)], beta.est=beta.k[(iter-1),], beta.stat=cate.stat, beta.pvalue=cate.pvalue, BF.gene=data.frame(Gene=unique.gene, BF=BF.gene[(iter-1),]), full.info=full.info.genevar, Eui=EUi))
}

#########################################

```{r, echo=F, cache=T, warning=F, message=F}
### read into new sample 
path="C:\\Shengtong\\research\\rare-var\\"
ASC_new_sample3=as_tibble(read.table(paste(path, "ASC\\ASC_v17_raw_TDT_by_parent_2021-05-10.txt", sep=""), sep='\t',  header=T, fill=T)) # this is the correct way that Kyle read the data
ASC_new_sample3=ASC_new_sample3 %>% filter(isSyn=="FALSE") # filter syn variants first 
```

```{r, echo=F, warning=F, message=F}

##### count transmitted variants and untransmitted in proband by summing over male, female, dad, mom and ind 
######### use ASC_new_sample3 rather than ASC_new_sample2
Transmitted_proband=ASC_new_sample3%>% select(starts_with("t_proband"))%>%  mutate(Transmitted_proband = select(., t_proband_male_dad:t_proband_female_ind) %>% rowSums())%>% select(Transmitted_proband) # extract transmitted variants for proband 

Untransmitted_proband=ASC_new_sample3%>% select(starts_with("u_proband"))%>%  mutate(Untransmitted_proband = select(., u_proband_male_dad:u_proband_female_ind) %>% rowSums())%>% select(Untransmitted_proband) # extract untransmitted variants for proband 

Transmitted_sibling=ASC_new_sample3%>% select(starts_with("t_sibling"))%>%  mutate(Transmitted_sibling = select(., t_sibling_male_dad:t_sibling_female_ind) %>% rowSums())%>% select(Transmitted_sibling) # extract transmitted variants for sibling 

Untransmitted_sibling=ASC_new_sample3%>% select(starts_with("u_sibling"))%>%  mutate(Untransmitted_sibling = select(., u_sibling_male_dad:u_sibling_female_ind) %>% rowSums())%>% select(Untransmitted_sibling) # extract untransmitted variants for sibling

##### add new columns of total transmitted and non-transmitted variants 
#ASC_new_sample=ASC_new_sample3%>%mutate(Transmitted_proband=Transmitted_proband$Transmitted_proband, Untransmitted_proband=Untransmitted_proband$Untransmitted_proband, Transmitted_sibling=Transmitted_sibling$Transmitted_sibling, Untransmitted_sibling=Untransmitted_sibling$Untransmitted_sibling) %>% filter(isIndel=="FALSE" & VQSLOD<9) # exclude indel & VQSLOD filter 

################################################

N1=N0=7291
ASC_new_sample=ASC_new_sample3%>%mutate(Transmitted_proband=Transmitted_proband$Transmitted_proband, Untransmitted_proband=Untransmitted_proband$Untransmitted_proband, Transmitted_sibling=Transmitted_sibling$Transmitted_sibling, Untransmitted_sibling=Untransmitted_sibling$Untransmitted_sibling) %>% filter(isIndel=="FALSE" & QD >= 10 & AS_SOR <= 3 & AS_ReadPosRankSum >= -0.8)%>%filter(!str_detect(Consequence, 'intron_variant'))# exclude indel & apply other filters for SNP's 
```

```{r, echo=F, message=F, warning=F, eval=T}
# add one column of variant group index 
ASC_new_sample_GI2=ASC_new_sample %>% filter(!str_detect(Consequence, 'UTR')& Gnomad_non_neuro_AF<0.05)%>% add_column(Group_Index=NA)

############# LOF variants 
ASC_new_sample_GI2$Group_Index[which(ASC_new_sample_GI2$pLI>=0.995 & ASC_new_sample_GI2$isPTV==T & ASC_new_sample_GI2$Gnomad_non_neuro_AF<0.01)]=1
ASC_new_sample_GI2$Group_Index[which(ASC_new_sample_GI2$pLI>=0.995 & ASC_new_sample_GI2$isPTV==T & ASC_new_sample_GI2$Gnomad_non_neuro_AF>=0.01)]=2

ASC_new_sample_GI2$Group_Index[which(ASC_new_sample_GI2$pLI>=0.5 & ASC_new_sample_GI2$pLI<0.995 & ASC_new_sample_GI2$isPTV==T & ASC_new_sample_GI2$Gnomad_non_neuro_AF<0.01)]=3
ASC_new_sample_GI2$Group_Index[which(ASC_new_sample_GI2$pLI>=0.5 & ASC_new_sample_GI2$pLI<0.995 & ASC_new_sample_GI2$isPTV==T & ASC_new_sample_GI2$Gnomad_non_neuro_AF>=0.01)]=4

ASC_new_sample_GI2$Group_Index[which(ASC_new_sample_GI2$pLI<0.5 & ASC_new_sample_GI2$isPTV==T & ASC_new_sample_GI2$Gnomad_non_neuro_AF<0.01)]=5
ASC_new_sample_GI2$Group_Index[which(ASC_new_sample_GI2$pLI<0.5 & ASC_new_sample_GI2$isPTV==T & ASC_new_sample_GI2$Gnomad_non_neuro_AF>=0.01)]=6
########## missense variants 
ASC_new_sample_GI2$Group_Index[which(ASC_new_sample_GI2$MPC>=2 & ASC_new_sample_GI2$Gnomad_non_neuro_AF<0.001)]=7
ASC_new_sample_GI2$Group_Index[which(ASC_new_sample_GI2$MPC>=2 & ASC_new_sample_GI2$Gnomad_non_neuro_AF>=0.001 & ASC_new_sample_GI2$Gnomad_non_neuro_AF<0.01)]=8
ASC_new_sample_GI2$Group_Index[which(ASC_new_sample_GI2$MPC>=2 & ASC_new_sample_GI2$Gnomad_non_neuro_AF>=0.01 )]=9
ASC_new_sample_GI2$Group_Index[which(ASC_new_sample_GI2$MPC>=1 & ASC_new_sample_GI2$MPC<2 & ASC_new_sample_GI2$Gnomad_non_neuro_AF<0.001)]=10
ASC_new_sample_GI2$Group_Index[which(ASC_new_sample_GI2$MPC>=1 & ASC_new_sample_GI2$MPC<2 & ASC_new_sample_GI2$Gnomad_non_neuro_AF>=0.001 & ASC_new_sample_GI2$Gnomad_non_neuro_AF<0.01)]=11
ASC_new_sample_GI2$Group_Index[which(ASC_new_sample_GI2$MPC>=1 & ASC_new_sample_GI2$MPC<2 & ASC_new_sample_GI2$Gnomad_non_neuro_AF>=0.01)]=12

ASC_new_sample_GI2$Group_Index[which(ASC_new_sample_GI2$MPC<1 & ASC_new_sample_GI2$Gnomad_non_neuro_AF<0.001)]=13
ASC_new_sample_GI2$Group_Index[which(ASC_new_sample_GI2$MPC<1 & ASC_new_sample_GI2$Gnomad_non_neuro_AF>=0.001 & ASC_new_sample_GI2$Gnomad_non_neuro_AF<0.01)]=14
ASC_new_sample_GI2$Group_Index[which(ASC_new_sample_GI2$MPC<1 & ASC_new_sample_GI2$Gnomad_non_neuro_AF>=0.01)]=15
num.family=7291

actu.num.group=15
new.data=ASC_new_sample_GI2 %>% filter(Gnomad_non_neuro_AF<0.05) %>% select(Variant, Gene, Transmitted_proband, Untransmitted_proband, Group_Index) %>% drop_na() # drop all variants neither LoF nor missense variants 
```

delta=0.1  # use fixed proportion of risk genes 
beta.init=rep(0.1, 15)
colnames(new.data)[3:5]=c("No.case", "No.contrl", "group.index")
para.est=multi.group.func(new.data, N1, N0, gamma.mean=3, sigma=2, delta=0.3, beta.init, actu.num.group)
