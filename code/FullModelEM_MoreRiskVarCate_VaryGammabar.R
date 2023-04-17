# Ui; latent binary variant for gene , Ui~Ber(delta) Ui=1: risk gene Ui=0: non risk gene
# Zij: latent binary variant for variant, Zij~Ber(pi), Zij=1, variant j of gene i is causal variant, Zij=0 otherwise  
# m: number of variant 
# N0: number of controls; N1: number of cases
# num.gene: number of genes in the sample 
##################################
rm(list=ls())
#set.seed(123)
library(SKAT)
library(ACAT)
start.time=date()
cat("program starts at", date(),"\n\n")
ptm <- proc.time()
intergrand=function(aa, indi.var, pheno, bar.gamma, sig)
{
  ff=dbinom(sum(indi.var[pheno==1]), sum(indi.var), aa*N1/(aa*N1+N0))*dgamma(aa, bar.gamma*sig, sig)
  return(ff)
}
# calculate the bayes factor of a single variant via integration
BF.gene.inte=function(geno.var, pheno, bar.gamma, sig)
{
  marglik0.CC <- dbinom(sum(geno.var[pheno==1]), sum(geno.var), N1/(N1+N0))    # Under H0: gamma=1
  
  marglik1.CC <- integrate(intergrand, geno.var, pheno, bar.gamma, sig, low=0, upper=100, stop.on.error=F)$value # Under H1: gamma~gamma(gamma.mean*sigma, sigma) 
  BF.var <- marglik1.CC/marglik0.CC
  
  return(BF.var)
}
##################################
gene.simu=function(N0, N1, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, pi, model, num.group, split.ratio)
{

  pheno=c(rep(0,N0), rep(1,N1))
  q <- rbeta(m, alpha0, beta0)
  gamma <- rep(1,m)
  Zij=rep(0,m)
  if (model==1)
  {
    for (k in 1:num.group)
    {
      from=split.ratio[k]*m+1; to=split.ratio[k+1]*m 
      Zij[from:to]=rbinom((to-from+1), 1, pi[k])
      Zij_1_index=which(Zij[from:to]==1)
      gamma[from+Zij_1_index]=rgamma(length(Zij_1_index), gamma.mean[k]*sigma, sigma)
      #gamma[Zij[from:to]==1] <- rgamma(sum(Zij[from:to]==1), gamma.mean[k]*sigma, sigma)
    }
    q[Zij==1] <- rbeta(sum(Zij==1), alpha, beta)
    #gamma[Zij==1] <- rgamma(sum(Zij==1), gamma.mean*sigma, sigma)
  }
  
  x <- array(0, dim=c(length(pheno), m))
  for (j in 1:m)
  {
    x[pheno==0,j] <- sample(c(1,0), N0, replace=TRUE, prob=c(2*q[j],1-2*q[j]))
    x[pheno==1,j] <- sample(c(1,0), N1, replace=TRUE, prob=c(2*q[j]*gamma[j],1-2*q[j]*gamma[j]))
  }
  
  var.index=which(colSums(x!=0)>0) # get the variant index
  x=x[,!apply(x==0,2,all)]  # filter variants with 0 counts in both cases and controls since these kind of variants are noninformative 
  x=as.matrix(x)
  return (list(geno=x,pheno=pheno,q=q,gamma=gamma, Zij=Zij, var.index=var.index)) 
  
}
##################################
num.gene=1000
m=100
N0=3000; N1=3000
delta=0.1
alpha0 <- 0.1
beta0 <- 1000
alpha <- 0.1
beta <- 2000
gamma.mean <- c(3,3,5) # three variant categories, each has its own gamma.mean 
sigma <- 1
num.group=3
split.ratio=c(0, 0.6, 0.9, 1)
pi=numeric(num.group)
pi[1]=0.05; pi[2]=0.2; pi[3]=0.5
max.run=5
MAC_threshold=0
all.pi=matrix(nrow=max.run, ncol=(num.group+1))
all.teststat=matrix(nrow=max.run, ncol=num.group)
actu.pi=matrix(nrow=max.run, ncol=num.group)
pvalue.fish=matrix(nrow=max.run, ncol=num.group)


all.pvalue=numeric()
all.conti.table=list()
odds.ratio=numeric()
log.lkhd=numeric() # used to check the concavity of likelihood function of beta 
TP=matrix(nrow=max.run, ncol=3)

MIRAGE.pvalue=list()
MIRAGE.BF=list()
Fisher.pvalue=list()
SKATO.pvalue=list()
Gene.Risk.Status=list()
CMC.pvalue=list()
ASUM.pvalue=list()
Fisher.separate.pvalue=list()
Regres.pvalue=list()
ACAT.pvalue=list()
########################  generate data 
for (run in 1:max.run) {
actu.no.var=matrix(nrow=num.gene, ncol=num.group)
all.data=list()
causal.var.index=matrix(nrow=num.gene, ncol=m)
BF.var=matrix(0,nrow=num.gene, ncol=m);
Ui=rbinom(num.gene, 1, delta)
stop.cond=0; iter=1
thrshd=1e-5
max.iter=100000
beta.k=matrix(nrow=max.iter, ncol=num.group)
beta.k[1, 1]=0.7; beta.k[1, 2]=0.2; beta.k[1,3]=0.2 # initial settings 
delta.est=numeric(); delta.est[1]=0.8
BF.gene=matrix(nrow=max.iter, ncol=num.gene)
final.causal.var=matrix(nrow=num.gene, ncol=num.group)
conti.matx=list()
skat.pvalue=numeric()
acat_gene_pvalue=numeric()
########################
for (i in 1:num.gene)
{
  cat(i, "th gene of ", "\t", num.gene, "\t", "is running", "\t", run, "th run", "\n")
  data=gene.simu(N0, N1, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, pi, Ui[i], num.group, split.ratio)
  all.data[[i]]=data
  causal.var.index[i,]=data$Zij
  conti.gene.wide=list()
  acat_var_group_pvalue=numeric()
  for (k in 1:num.group)
  {
    from=split.ratio[k]*m+1; to=split.ratio[k+1]*m 
    xx=which(data$Zij==1); xx=xx[xx>=from & xx<=to] # xx; causal variant  before filtering
    yy=data$var.index[data$var.index>=from & data$var.index<=to] # yy: final variant index after filtering
    actu.no.var[i,k]=length(yy); zz=intersect(xx, yy) # zz; causal variant after filtering
    final.causal.var[i,k]=length(zz) 
    
    #################### fisher exact test 
    conti.table=matrix(nrow=2,ncol=2)  # used for fisher exact test 
    conti.table[1,1]=sum(data$geno[(data$pheno==1),which(data$var.index %in% yy)]) # number of variants in cases 
    conti.table[1,2]=sum(data$pheno==1)*length(yy)-conti.table[1,1]                
    conti.table[2,1]=sum(data$geno[(data$pheno==0),which(data$var.index %in% yy)])  # number of variants in control 
    conti.table[2,2]=sum(data$pheno==0)*length(yy)-conti.table[2,1]
    conti.gene.wide[[k]]=conti.table  
    ##################### ACAT test 
    acat_var_pvalue=numeric()
    if (length(yy)>1)
      yy_less_than10=yy[which(colSums(data$geno[,which(data$var.index %in% yy)])<MAC_threshold)] # find var index that has MAC<10 and put them together 
    if (length(yy)==1)
      yy_less_than10=yy[which(sum(data$geno[,which(data$var.index %in% yy)])<MAC_threshold)] # find var index that has MAC<10 and put them together 
    
    
    if (length(yy_less_than10)==length(yy)) #case1:  all single variants has MAC<10, then pool them together to do binomial test
      { 
      acat_var_pvalue=fisher.test(conti.table)$p.value # pool all variants together, which is conti.table for burden test
      acat_var_group_pvalue[k]=acat_var_pvalue
    }
    if (length(yy_less_than10)==0) # case2: no variants has MAC<10
    {
      ## do burden test for all the rest one by one 
      yy_more_than10=yy
      
      for (t in 1:length(yy_more_than10))
      {
        conti.table=matrix(nrow=2,ncol=2)  # used for fisher exact test 
        conti.table[1,1]=sum(data$geno[(data$pheno==1),which(data$var.index ==yy_more_than10[t])]) # number of variants in cases 
        conti.table[1,2]=sum(data$pheno==1)-conti.table[1,1]                
        conti.table[2,1]=sum(data$geno[(data$pheno==0),which(data$var.index ==yy_more_than10[t])])  # number of variants in control 
        conti.table[2,2]=sum(data$pheno==0)-conti.table[2,1]
        acat_var_pvalue[t]=fisher.test(conti.table)$p.value 
      }
      if (max(acat_var_pvalue)>=1)
        acat_var_group_pvalue[k]=1-1/(length(acat_var_pvalue))
      if (max(acat_var_pvalue)<1)
        acat_var_group_pvalue[k]=ACAT(acat_var_pvalue)
      
    }
    
    if ( length(yy_less_than10)>0 & length(yy_less_than10)<length(yy)) # case3: some not all variants that have MAC<10
    {
      ## combine all variants with MAC<10 into one variant, then do burden test 
      conti.table=matrix(nrow=2,ncol=2)  # used for fisher exact test 
      conti.table[1,1]=sum(data$geno[(data$pheno==1),which(data$var.index %in% yy_less_than10)]) # number of variants in cases 
      conti.table[1,2]=sum(data$pheno==1)*length(yy_less_than10)-conti.table[1,1]                
      conti.table[2,1]=sum(data$geno[(data$pheno==0),which(data$var.index %in% yy_less_than10)])  # number of variants in control 
      conti.table[2,2]=sum(data$pheno==0)*length(yy_less_than10)-conti.table[2,1]
      acat_var_pvalue[1]=fisher.test(conti.table)$p.value # pool all variants together, which is conti.table for burden test
      ## do burden test for all the rest one by one 
      yy_more_than10=yy[which(colSums(data$geno[,which(data$var.index %in% yy)])>=MAC_threshold)]
      for (t in 1:length(yy_more_than10))
      {
        conti.table=matrix(nrow=2,ncol=2)  # used for fisher exact test 
        conti.table[1,1]=sum(data$geno[(data$pheno==1),which(data$var.index ==yy_more_than10[t])]) # number of variants in cases 
        conti.table[1,2]=sum(data$pheno==1)-conti.table[1,1]                
        conti.table[2,1]=sum(data$geno[(data$pheno==0),which(data$var.index ==yy_more_than10[t])])  # number of variants in control 
        conti.table[2,2]=sum(data$pheno==0)-conti.table[2,1]
        acat_var_pvalue[1+t]=fisher.test(conti.table)$p.value 
      }
      if (max(acat_var_pvalue)>=1)
        acat_var_group_pvalue[k]=1-1/(length(acat_var_pvalue))
      if (max(acat_var_pvalue)<1)
        acat_var_group_pvalue[k]=ACAT(acat_var_pvalue)
      
    } # end of  if (length(yy_less_than10)<length(yy)) # exist variants that have MAC>=10
    
    
  }  # end of  for (k in 1:num.group)
  ################ ACAT; combine p value from group pvalues 
  acat_var_group_pvalue=acat_var_group_pvalue[!is.na(acat_var_group_pvalue)] # remove NA which is caused by no variants in the group after filtering
  if (max(acat_var_group_pvalue)>=1)
    acat_gene_pvalue[i]=1-1/length(acat_var_group_pvalue)
  if (max(acat_var_group_pvalue)<1)
    acat_gene_pvalue[i]=ACAT(acat_var_group_pvalue)
  ###############
  
  conti.matx[[i]]=conti.gene.wide
  
  #############  SKAT
  obj<-SKAT_Null_Model(all.data[[i]]$pheno~ 1, out_type="D")  # without covariates 
  skat.pvalue[i]=SKAT(all.data[[i]]$geno, obj, method="SKATO")$p.value
  
  
  bb=1
  if (ncol(data$geno)>0) # calculate Bayes factor for variant (i,j)
    for (j in 1:ncol(data$geno)) 
    {
      orig.var.index=data$var.index[j]
     # BF.var[i,orig.var.index]=BF.gene.inte(data$geno[,j], data$pheno, bar.gamma=gamma.mean, sig=sigma)
     for (k in 1:num.group)
     {
       from=split.ratio[k]*m+1; to=split.ratio[k+1]*m 
       if (orig.var.index>=from & orig.var.index<=to)
          gp.index=k
     }
      BF.var[i,orig.var.index]=BF.gene.inte(data$geno[,j], data$pheno, bar.gamma=gamma.mean[gp.index], sig=sigma) # use cate-specific gamma.mean 
      bb=bb*((1-beta.k[1, gp.index])+beta.k[1, gp.index]*BF.var[i,orig.var.index])
    }
  BF.gene[1, i]=bb
}  ### end of i 

################
for (k in 1:num.group)
 actu.pi[run,k]=sum(final.causal.var[Ui==1,k])/sum(actu.no.var[Ui==1,k])
########################## EM algorithm ############################

########################
while (stop.cond==0)
{
  iter=iter+1
  ############## EM algorithm: E step 
  EUiZij=matrix(0,nrow=num.gene, ncol=m)
  EUi=numeric()
  for (i in 1:num.gene)
  { 
    data=all.data[[i]]
    bb=1
    if (ncol(data$geno)>0)
      for (j in 1:ncol(data$geno))  
      {
        orig.var.index=data$var.index[j]
        for (k in 1:num.group)
        {
          from=split.ratio[k]*m+1; to=split.ratio[k+1]*m 
          if (orig.var.index>=from & orig.var.index<=to)
            gp.index=k
        }
        numer=BF.var[i,orig.var.index]*beta.k[(iter-1), gp.index]*delta.est[iter-1]
        denom=(delta.est[iter-1]+(1-delta.est[iter-1])/BF.gene[(iter-1),i])*(beta.k[(iter-1), gp.index]*BF.var[i,orig.var.index]+(1-beta.k[(iter-1), gp.index]))
        EUiZij[i, orig.var.index]=numer/denom
        bb=bb*((1-beta.k[(iter-1), gp.index])+beta.k[(iter-1), gp.index]*BF.var[i,orig.var.index])
      }
    BF.gene[iter,i]=bb
    EUi[i]=delta.est[iter-1]*bb/(delta.est[iter-1]*bb+1-delta.est[iter-1])
  }
  ############## EM algorithm: M step
  delta.est[iter]=sum(EUi)/num.gene
  #delta.est[iter]=delta
  for (k in 1:num.group) 
  {
    from=split.ratio[k]*m+1; to=split.ratio[k+1]*m 
    beta.k[iter,k]=sum(EUiZij[,(from:to)])/sum(actu.no.var[,k]*EUi) 
  }
  
  diff=sum(abs(beta.k[iter,]-beta.k[(iter-1),]))+abs(delta.est[iter]-delta.est[iter-1])
  if (diff<thrshd || iter>(max.iter-1))  # stop condition 
    stop.cond=1
  cat(iter, "th iteration is running", "\t", run, "th run", "\n")
} # end of iter
if (iter<max.iter)
  beta.k=beta.k[complete.cases(beta.k),]
################## calculate the likelihood ratio test statistics 
  lkhd=matrix(1,nrow=num.gene, ncol=num.group); total.lkhd=rep(0, num.group)
  for (i in 1:num.gene)
  {
#    i=33
    data=all.data[[i]]
    
    for (j in 1:ncol(data$geno))
    {
      orig.var.index=data$var.index[j]
      for (k in 1:num.group) # find the variant group for variant (i,j)
      {
        from=split.ratio[k]*m+1; to=split.ratio[k+1]*m 
        if (orig.var.index>=from & orig.var.index<=to)
          gp.index=k
      }
      lkhd[i,gp.index]=lkhd[i, gp.index]*((1-beta.k[(iter-1), gp.index])+beta.k[(iter-1), gp.index]*BF.var[i,orig.var.index])
    } # end of j 
    total.lkhd=total.lkhd+log((1-delta.est[iter-1])+delta.est[iter-1]*lkhd[i,])
    } # end of i
2*total.lkhd
all.teststat[run,]=2*total.lkhd
##################
# calculate pvalues by fisher exact test 
for (k in 1:num.group)
{
  contigen=matrix(0, nrow=2, ncol=2)
  for (i in 1:num.gene)
    contigen=contigen+conti.matx[[i]][[k]]
  
  pvalue.fish[run,k]=fisher.test(contigen)$p.value
}
##################
############# gene level test ###############
mirage.test.stat=2*rowSums(lkhd)   # mirage 
mirage.pvalue=pchisq(mirage.test.stat, 4, lower.tail=F)  # 4 is the number of unknown parameters 
fisher.odds.ratio=numeric(); fisher.pvalue=numeric()
fisher.separate.odds.ratio=matrix(nrow=num.gene, ncol=3); fisher.separate.pvalue=matrix(nrow=num.gene, ncol=3)
cmc.pvalue=numeric()
asum.pvalue=numeric()
regres.pvalue=numeric()
library(AssotesteR)
for (i in 1:num.gene)
{
  ###############    fisher 
  fish.test=fisher.test(conti.matx[[i]][[1]]+conti.matx[[i]][[2]]+conti.matx[[i]][[3]])
  fisher.odds.ratio[i]=fish.test$estimate
  fisher.pvalue[i]=fish.test$p.value
  
  for (ii in 1:3)
  {
  fish.separate.test=fisher.test(conti.matx[[i]][[ii]])
  fisher.separate.odds.ratio[i,ii]=fish.separate.test$estimate
  fisher.separate.pvalue[i,ii]=fish.separate.test$p.value
  }
  ############ regress on group specific variants by logit regression
  no.var=c(conti.matx[[i]][[1]][1,1], conti.matx[[i]][[2]][1,1], conti.matx[[i]][[3]][1,1])
  risk.status=Ui[i]
  gene.data=data.frame(risk.status=risk.status, no.var=no.var)
  full.model <- glm(risk.status ~no.var,family=binomial(link='logit'),data=gene.data)
  null.model<-glm(risk.status ~1,family=binomial(link='logit'),data=gene.data)
  regres.pvalue[i]=pchisq(deviance(null.model)-deviance(full.model),
                       df.residual(null.model)-df.residual(full.model),
                       lower.tail=FALSE)
  ############# CMC test 
  cmc.test=CMC(all.data[[i]]$pheno, all.data[[i]]$geno, maf=10^(-4)*c(0.05, 0.2, 0.5), perm=100)
  cmc.pvalue[i]=cmc.test$asym.pval
  
  ############  ASUM test 
  asum.pvalue[i]=ASUM(all.data[[i]]$pheno, all.data[[i]]$geno, perm=100)$perm.pval
  
  #######################
}
detach("package:AssotesteR", unload=TRUE)
#####################

par(mfrow=c(1,3))
if (num.group==1)
{
 plot(beta.k, ylim=c(0, 1), main=expression(paste (beta)), xlab="Iteration index", ylab="")
 abline(h=pi[1], col="red")
}
if (num.group>1)
  plot(beta.k[,1], ylim=c(0, 1), main=expression(paste (beta[1])), xlab="Iteration index", ylab="")
#abline(h=pi[1], col="red")
abline(h=actu.pi[run,1], col="red")
plot(beta.k[,2], ylim=c(0, 1), main=expression(paste (beta[2])), xlab="Iteration index", ylab="")
#abline(h=pi[2], col="red")
abline(h=actu.pi[run,2], col="red")
#plot(beta.k[,3], ylim=c(0, 1), main=expression(paste (beta[3])), xlab="Iteration index", ylab="")
#abline(h=pi[3], col="red")
plot(delta.est, ylim=c(0, 1), main=expression(delta), xlab="Iteration index", ylab="")
abline(h=delta, col="red")
#beta.k
#pi
all.pi[run,]=c(beta.k[(iter-1),], delta.est[iter-1])

TP[run,]=c(sum(mirage.pvalue<0.05), sum(fisher.pvalue<0.05), sum(skat.pvalue<0.05))/num.gene

###################### store statistics for further analysis
Gene.Risk.Status[[run]]=Ui
MIRAGE.pvalue[[run]]=mirage.pvalue
MIRAGE.BF[[run]]=BF.gene[iter,]
Fisher.pvalue[[run]]=fisher.pvalue
SKATO.pvalue[[run]]=skat.pvalue
CMC.pvalue[[run]]=cmc.pvalue
ASUM.pvalue[[run]]=asum.pvalue
Fisher.separate.pvalue[[run]]=fisher.separate.pvalue
Regres.pvalue[[run]]=regres.pvalue
ACAT.pvalue[[run]]=acat_gene_pvalue

} # end of run

#save(Gene.Risk.Status, MIRAGE.pvalue, MIRAGE.BF, Fisher.pvalue, SKATO.pvalue, CMC.pvalue, 
#    ASUM.pvalue, Fisher.separate.pvalue, ACAT.pvalue, 
#     file="C:\\Shengtong\\Research\\rare-var\\RareVariant\\rare-var-project\\output\\Var_specific_bargamma_MixedGene_MixtureVariant\\7methods_delta0.1_MAC_threshold0.gammamean335_5rep1.RData")
######################
end.time=date()
cat("program ends at", date(),"\n\n")
time.spent<-proc.time()-ptm
time.spent

