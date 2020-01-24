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
num.gene=1  # number of genes-gene level
m=1000
N0=N1=5000
delta=1
alpha0 <- 0.1
beta0 <- 1000
alpha <- 0.1
beta <- 2000
gamma.mean <- 5
sigma <- 1
num.group=1
split.ratio=c(0, 1)
pi=0.1
max.run=1
all.pi=matrix(nrow=max.run, ncol=(num.group+1))
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
  BF.var=matrix(0,nrow=num.gene, ncol=m);
  Ui=rbinom(num.gene, 1, delta)
  stop.cond=0; iter=1
  thrshd=1e-5
  max.iter=5000
  beta.k=numeric(); beta.k[1]=runif(1)
  delta.est=numeric(); delta.est[1]=delta
  BF.gene=matrix(nrow=max.iter, ncol=num.gene)
  conti.matx=list()
  ########################
  for (i in 1:num.gene)
  {
    cat(i, "th gene of ", "\t", num.gene, "\t", "is running", "\t", run, "th run of ", max.run, "\n")
    data=gene.simu(N0, N1, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, pi, Ui[i], num.group, split.ratio)
    all.data[[i]]=data

    # construct  the contigency table for each gene
    xx=which(data$Zij==1)    # xx; causal variant
    yy=data$var.index   # yy: final variant index
    actu.no.var[i]=length(yy)


    conti.table=matrix(nrow=2,ncol=2)  # used for fisher exact test
    conti.table[1,1]=sum(data$geno[(data$pheno==1),which(yy%in%data$var.index)]) # var.counts in cases
    conti.table[1,2]=sum(data$pheno==1)*length(yy)-conti.table[1,1]
    conti.table[2,1]=sum(data$geno[(data$pheno==0),which(yy%in%data$var.index)]) # var.counts in controls
    conti.table[2,2]=sum(data$pheno==0)*length(yy)-conti.table[2,1]

    conti.matx[[i]]=conti.table

    bb=1
    if (ncol(data$geno)>0) # calculate Bayes factor for variant (i,j)
      for (j in 1:ncol(data$geno))
      {
        orig.var.index=data$var.index[j]
        BF.var[i,orig.var.index]=BF.gene.inte(data$geno[,j], data$pheno, bar.gamma=gamma.mean, sig=sigma)

        bb=bb*((1-beta.k[1])+beta.k[1]*BF.var[i,orig.var.index])
      }
    BF.gene[1, i]=bb   # bayes facor for that gene
  }  # end of i
  ################
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
      if (ncol(data$geno)>0)
        for (j in 1:ncol(data$geno))
        {
          orig.var.index=data$var.index[j]
          numer=BF.var[i,orig.var.index]*beta.k[iter-1]*delta.est[iter-1]
          denom=(delta.est[iter-1]+(1-delta.est[iter-1])/BF.gene[(iter-1),i])*(beta.k[iter-1]*BF.var[i,orig.var.index]+(1-beta.k[iter-1]))
          EUiZij[i, orig.var.index]=numer/denom
          bb=bb*((1-beta.k[iter-1])+beta.k[iter-1]*BF.var[i,orig.var.index])
        }
      BF.gene[iter,i]=bb
      EUi[i]=delta.est[iter-1]*bb/(delta.est[iter-1]*bb+1-delta.est[iter-1])
    }
    ############## EM algorithm: M step
    #delta.est[iter]=sum(EUi)/num.gene
    delta.est[iter]=delta
    beta.k[iter]=sum(EUiZij)/sum(actu.no.var*EUi)

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
    {
      orig.var.index=data$var.index[j]
      lkhd[i]=lkhd[i]*((1-beta.k[iter-1])+beta.k[iter-1]*BF.var[i,orig.var.index])
    } # end of j
    total.lkhd=total.lkhd+log((1-delta.est[iter-1])+delta.est[iter-1]*lkhd[i])
  } # end of i
  2*total.lkhd
  all.teststat[run]=2*total.lkhd
  all.pvalue[run]=pchisq(2*total.lkhd, 1, lower.tail=F)
  ##################
  # calculate pvalues by fisher exact test
  contigen=matrix(0, nrow=2, ncol=2)
  for (i in 1:num.gene)
    contigen=contigen+conti.matx[[i]]

  fish.test=fisher.test(contigen)
  pvalue.fish[run]=fish.test$p.value
  odds.ratio[run]=fish.test$estimate
  all.conti.table[[run]]=contigen
  ##################
  ############# gene level test ###############
  mirage.test.stat=2*lkhd   # mirage
  mirage.pvalue=pchisq(mirage.test.stat, 1, lower.tail=F)  # degrees of freedom depends on the number of unknown parameters
  fisher.odds.ratio=numeric(); fisher.pvalue=numeric()
  skat.pvalue=numeric()
  for (i in 1:num.gene)
  {
    ###############    fisher
    fish.test=fisher.test(conti.matx[[i]])
    fisher.odds.ratio[i]=fish.test$estimate
    fisher.pvalue[i]=fish.test$p.value
    #############  SKAT
    obj<-SKAT_Null_Model(all.data[[i]]$pheno~ 1, out_type="D")  # without covariates
    skat.pvalue[i]=SKAT(all.data[[i]]$geno, obj, method="SKATO")$p.value

  }



  #################

  par(mfrow=c(1,2))
  plot(beta.k, ylim=c(0, 1), main=expression(paste(eta)), xlab="Iteration index", ylab="")
  abline(h=pi, col="red")
  plot(delta.est, ylim=c(0, 1), main=expression(delta), xlab="Iteration index", ylab="")
  abline(h=delta, col="red")

  all.pi[run,]=c(beta.k[iter-1], delta.est[iter-1])

  TP[run,]=c(sum(mirage.pvalue<0.05), sum(fisher.pvalue<0.05), sum(skat.pvalue<0.05))/num.gene
  ###################### store statistics for further analysis
  Gene.Risk.Status[[run]]=Ui
  MIRAGE.pvalue[[run]]=mirage.pvalue
  MIRAGE.BF[[run]]=BF.gene[iter,]
  Fisher.pvalue[[run]]=fisher.pvalue
  SKATO.pvalue[[run]]=skat.pvalue

} # end of run
######################
row_sub= apply(TP, 1, function(row) all(row !=0 )) # remove rows with zero
TP_effect=TP[row_sub,]
colMeans(TP_effect)


end.time=date()
cat("program ends at", date(),"\n\n")
time.spent<-proc.time()-ptm
time.spent


