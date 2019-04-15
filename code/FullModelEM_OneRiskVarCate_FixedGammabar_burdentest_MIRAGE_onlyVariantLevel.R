# Ui; latent binary variant for gene , Ui~Ber(delta) Ui=1: risk gene Ui=0: non risk gene
# Zij: latent binary variant for variant, Zij~Ber(pi), Zij=1, variant j of gene i is causal variant, Zij=0 otherwise  
# m: number of variant 
# N0: number of controls; N1: number of cases
# num.gene: number of genes in the sample 
##################################
rm(list=ls())
set.seed(123)
library(SKAT)

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
  pheno=c(rep(0,N0), rep(1,N1))  # filtering step
  q <- rbeta(m, alpha0, beta0)
  x1x0=matrix(nrow=(N0+N1), ncol=m)
  for (i in 1:m)  # generate the data with certain variants full of zeros                
    x1x0[,i]=rpois((N0+N1), q[i])
  filter.x1x0=x1x0[,!apply(x1x0==0,2,all)]  # filter variants with 0 counts in both cases and controls since these variants are noninformative 
  mm=ncol(filter.x1x0)
  ##################
  gamma <- rep(1,mm)
  Zij=rep(0,mm)
  qq=rbeta(mm, alpha0, beta0)
  if (model==1)
  {
    for (k in 1:num.group)
    {
      from=split.ratio[k]*mm+1; to=split.ratio[k+1]*mm 
      Zij[from:to]=rbinom((to-from+1), 1, pi[k])
    }
    qq[Zij==1] <- rbeta(sum(Zij==1), alpha, beta)
    gamma[Zij==1] <- rgamma(sum(Zij==1), gamma.mean*sigma, sigma)
  }
  
  x <- array(0, dim=c(length(pheno), mm))
  for (j in 1:mm)
  {
    x1=rbinom(1, sum(filter.x1x0[,j]), N1*gamma[j]/(gamma[j]*N1+N0))
    #x1=sum(filter.x1x0[,j])-x0
    x0=sum(filter.x1x0[,j])-x1
    x[sample.int(N0,x0),j]=1
    x[N0+sample.int(N1,x1),j]=1
  }
  
  var.index=which(colSums(x1x0!=0)>0) # get the variant index in the original data before filtering
  x=as.matrix(x)
  return (list(geno=x,pheno=pheno,q=q,gamma=gamma, Zij=Zij, var.index=var.index)) 
  
}
##################################
num.gene=1
m=100
N0=3000; N1=3000
delta=0
alpha0 <- 0.1
beta0 <- 1000
alpha <- 0.1
beta <- 2000
gamma.mean <- 4
sigma <- 1
num.group=1
split.ratio=c(0, 1)
pi=0
max.run=100
all.pi=matrix(nrow=max.run, ncol=1)
all.teststat=numeric()
pvalue.fish=numeric()
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
########################  generate data 
for (run in 1:max.run) {
  actu.no.var=numeric()
  all.data=list()
  BF.var=numeric()
  Ui=rbinom(num.gene, 1, delta)
  stop.cond=0; iter=1
  thrshd=1e-5
  max.iter=5000
  beta.k=numeric(); beta.k[1]=runif(1)
  delta.est=numeric(); delta.est[1]=delta
  BF.gene=matrix(nrow=max.iter, ncol=num.gene)
  conti.matx=list()
  ########################
  data=gene.simu(N0, N1, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, pi, Ui, num.group, split.ratio)

    # construct  the contigency table for each gene
    xx=which(data$Zij==1)    # xx; causal variant 
    yy=data$var.index   # yy: final variant index 
    actu.no.var=length(yy)
    
    
    conti.table=matrix(nrow=2,ncol=2)  # used for fisher exact test 
    conti.table[1,1]=sum(data$geno[(data$pheno==1),which(yy%in%data$var.index)]) # var.counts in cases
    conti.table[1,2]=sum(data$pheno==1)*length(yy)-conti.table[1,1]
    conti.table[2,1]=sum(data$geno[(data$pheno==0),which(yy%in%data$var.index)]) # var.counts in controls
    conti.table[2,2]=sum(data$pheno==0)*length(yy)-conti.table[2,1]
    
    conti.matx=conti.table
    
    bb=1
    if (ncol(data$geno)>0) # calculate Bayes factor for variant j
      for (j in 1:ncol(data$geno)) 
      {
        orig.var.index=data$var.index[j]
        BF.var[orig.var.index]=BF.gene.inte(data$geno[,j], data$pheno, bar.gamma=gamma.mean, sig=sigma)
        
        bb=bb*((1-beta.k[1])+beta.k[1]*BF.var[orig.var.index])
      }

  ################
  ########################## EM algorithm
  ########################
  while (stop.cond==0)
  {
    iter=iter+1

    EZj=numeric() # expectation for variant j 
    
      bb=1
      if (ncol(data$geno)>0)
        for (j in 1:ncol(data$geno))  
        {
          orig.var.index=data$var.index[j]
          numer=BF.var[orig.var.index]*beta.k[iter-1]
          denom=beta.k[iter-1]*BF.var[orig.var.index]+(1-beta.k[iter-1])
          EZj[j]=numer/denom
          #EUiZij[i, orig.var.index]=numer/denom
          bb=bb*((1-beta.k[iter-1])+beta.k[iter-1]*BF.var[orig.var.index])
        }
    ############## EM algorithm: M step
    beta.k[iter]=sum(EZj)/(actu.no.var)
    ####################
    diff=sum(abs(beta.k[iter]-beta.k[iter-1]))
    if (diff<thrshd || iter>(max.iter-1))
      stop.cond=1
    cat(iter, "th iteration is running", "\t", run, "th run", "\n")
  } # end of iter
  if (iter<max.iter)
    beta.k=beta.k[complete.cases(beta.k)]
  ################## calculate the likelihood ratio test statistics 
  lkhd=1; total.lkhd=0
    for (j in 1:ncol(data$geno))
    {
      orig.var.index=data$var.index[j]
      lkhd=lkhd*((1-beta.k[iter-1])+beta.k[iter-1]*BF.var[orig.var.index])
    } # end of j 

  ##################
  ############# gene level test ###############
  mirage.test.stat=2*log(lkhd)   # mirage 
  mirage.pvalue=pchisq(mirage.test.stat, 1, lower.tail=F)  # degrees of freedom depends on the number of unknown parameters 
  fisher.odds.ratio=numeric(); fisher.pvalue=numeric()
  skat.pvalue=numeric()
  
    fish.test=fisher.test(conti.matx)
    fisher.odds.ratio=fish.test$estimate
    fisher.pvalue=fish.test$p.value
    #############  SKAT
    obj<-SKAT_Null_Model(data$pheno~ 1, out_type="D")  # without covariates 
    skat.pvalue=SKAT(data$geno, obj, method="SKATO")$p.value
    
  
  
  
  #################
  plot(beta.k, ylim=c(0, 1), main=expression(paste(eta)), xlab="Iteration index", ylab="")
  abline(h=pi, col="red")
  
  all.pi[run]=beta.k[iter-1]
  
  ###################### store statistics for further analysis
  Gene.Risk.Status[[run]]=Ui
  MIRAGE.pvalue[[run]]=mirage.pvalue
  Fisher.pvalue[[run]]=fisher.pvalue
  SKATO.pvalue[[run]]=skat.pvalue
  
} # end of run
######################

summary=tibble(eta_est=all.pi, MIRAGE.pvalue=MIRAGE.pvalue, Fisher.pvalue=Fisher.pvalue, SKATO.pvalue=SKATO.pvalue)

end.time=date()
cat("program ends at", date(),"\n\n")
time.spent<-proc.time()-ptm
time.spent


