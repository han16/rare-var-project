# Ui; latent binary variant for gene , Ui~Ber(delta) Ui=1: risk gene Ui=0: non risk gene
# Zij: latent binary variant for variant, Zij~Ber(pi), Zij=1, variant j of gene i is causal variant, Zij=0 otherwise
# m: number of variant
# N0: number of controls; N1: number of cases
# num.gene: number of genes in the sample
##################################
rm(list=ls())
set.seed(123)
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
    }
    q[Zij==1] <- rbeta(sum(Zij==1), alpha, beta)
    gamma[Zij==1] <- rgamma(sum(Zij==1), gamma.mean*sigma, sigma)
  }

  x <- array(0, dim=c(length(pheno), m))
  for (j in 1:m)
  {
    x[pheno==0,j] <- sample(c(1,0), N0, replace=TRUE, prob=c(2*q[j],1-2*q[j]))
    x[pheno==1,j] <- sample(c(1,0), N1, replace=TRUE, prob=c(2*q[j]*gamma[j],1-2*q[j]*gamma[j]))
  }

  return (list(geno=x,pheno=pheno,q=q,gamma=gamma, Zij=Zij))

}
##################################
num.gene=500
m=100
N0=3000; N1=3000
delta=1
alpha0 <- 0.1
beta0 <- 1000
alpha <- 0.1
beta <- 2000
gamma.mean <- 6
sigma <- 1
num.group=1
split.ratio=c(0, 1)
############ 

#pi=tau.true
pi=0.1
max.run=100
all.pi=matrix(nrow=max.run, ncol=(num.group+1))
all.teststat=numeric()
pvalue.fish=numeric()
all.conti.table=list()
odds.ratio=numeric()
log.lkhd=numeric() # used to check the concavity of likelihood function of beta
########################  generate data
for (run in 1:max.run) {
  all.data=list()
  causal.var.index=matrix(nrow=num.gene, ncol=m)
  BF.var=matrix(0,nrow=num.gene, ncol=m);
  Ui=rbinom(num.gene, 1, delta)
  stop.cond=0; iter=1
  thrshd=1e-5
  max.iter=1000
  beta.k=numeric(); beta.k[1]=0.7
  delta.est=numeric(); delta.est[1]=0.8
  BF.gene=matrix(nrow=max.iter, ncol=num.gene)
  conti.matx=list()
  ########################
  for (i in 1:num.gene)
  {
    cat(i, "th gene of ", "\t", num.gene, "\t", "is running", "\t", run, "th run of ", max.run, "\n")
    data=gene.simu(N0, N1, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, pi, Ui[i], num.group, split.ratio)
    all.data[[i]]=data

      # construct  the contigency table for each gene
    #  xx=which(data$Zij==1)    # xx; causal variant  before filtering
    #  yy=data$var.index   # yy: final variant index after filtering
    #  final.causal.var[i]=length(zz)

    #  conti.table=matrix(nrow=2,ncol=2)  # used for fisher exact test
    #  conti.table[1,1]=sum(data$geno[(data$pheno==1),which(yy%in%data$var.index)])
    #  conti.table[1,2]=sum(data$pheno==1)*length(yy)-conti.table[1,1]
    #  conti.table[2,1]=sum(data$geno[(data$pheno==0),which(yy%in%data$var.index)])
    #  conti.table[2,2]=sum(data$pheno==0)*length(yy)-conti.table[2,1]

    #   conti.matx[[i]]=conti.table

    bb=1
      for (j in 1:ncol(data$geno))
      {
        BF.var[i,j]=BF.gene.inte(data$geno[,j], data$pheno, bar.gamma=gamma.mean, sig=sigma)

        bb=bb*((1-beta.k[1])+beta.k[1]*BF.var[i,j])
      }
    BF.gene[1, i]=bb   # bayes facor for that gene
  }
  ########################## EM algorithm
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
        for (j in 1:m)
        {
          numer=BF.var[i,j]*beta.k[iter-1]*delta.est[iter-1]
          denom=(delta.est[iter-1]+(1-delta.est[iter-1])/BF.gene[(iter-1),i])*(beta.k[iter-1]*BF.var[i,j]+(1-beta.k[iter-1]))
          EUiZij[i, j]=numer/denom
          bb=bb*((1-beta.k[iter-1])+beta.k[iter-1]*BF.var[i,j])
        }
      BF.gene[iter,i]=bb
      EUi[i]=delta.est[iter-1]*bb/(delta.est[iter-1]*bb+1-delta.est[iter-1])
    }
    ############## EM algorithm: M step
    #delta.est[iter]=sum(EUi)/num.gene
    delta.est[iter]=1
    beta.k[iter]=sum(EUiZij)/sum(m*EUi)

    ####################
    diff=sum(abs(beta.k[iter]-beta.k[iter-1]))+abs(delta.est[iter]-delta.est[iter-1])
    if (diff<thrshd || iter>(max.iter-1))
      stop.cond=1
    cat(iter, "th iteration is running", "\t", run, "th run", "\n")
  } # end of iter
  if (iter<max.iter)
    beta.k=beta.k[complete.cases(beta.k)]
  ################## calculate the likelihood ratio test statistics
  lkhd=rep(1,num.gene); total.lkhd=0
  for (i in 1:num.gene)
  {
    #    i=33
    data=all.data[[i]]
    for (j in 1:ncol(data$geno))
      lkhd[i]=lkhd[i]*((1-beta.k[iter-1])+beta.k[iter-1]*BF.var[i,j])

    total.lkhd=total.lkhd+log((1-delta.est[iter-1])+delta.est[iter-1]*lkhd[i])
  } # end of i
  2*total.lkhd
  all.teststat[run]=2*total.lkhd
  ##################
  # calculate pvalues by fisher exact test
  #  contigen=matrix(0, nrow=2, ncol=2)
  #  for (i in 1:num.gene)
  #    contigen=contigen+conti.matx[[i]]
  #
  #  pvalue.fish[run]=fisher.test(contigen)$p.value
  #  odds.ratio[run]=fisher.test(contigen)$estimate
  #  all.conti.table[[run]]=contigen
  ##################

  par(mfrow=c(1,2))
  plot(beta.k, ylim=c(0, 1), main=expression(paste (beta)), xlab="Iteration index", ylab="")
  abline(h=pi, col="red")
  plot(delta.est, ylim=c(0, 1), main=expression(delta), xlab="Iteration index", ylab="")
  abline(h=delta, col="red")

  all.pi[run,]=c(beta.k[iter-1], delta.est[iter-1])
} # end of run
######################
end.time=date()
cat("program ends at", date(),"\n\n")
time.spent<-proc.time()-ptm
time.spent
