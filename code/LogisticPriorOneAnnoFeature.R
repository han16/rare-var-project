rm(list=ls())
# Ui; latent binary variant for gene , Ui~Ber(delta) Ui=1: risk gene Ui=0: non risk gene
# Zij: latent binary variant for variant, Zij~Ber(pi), Zij=1, variant j of gene i is causal variant, Zij=0 otherwise
# m: number of variant
# N0: number of controls; N1: number of cases
# num.gene: number of genes in the sample
library(DEoptim)
#################################
gene.simu=function(N0, N1, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, pi, model, num.group, split.ratio)
{
  pheno=c(rep(0,N0), rep(1,N1))  # filtering step
  q <- rbeta(m, alpha0, beta0)
  ############ generate counts to store in matrix, but may take large memory size  
  #x1x0=matrix(nrow=(N0+N1), ncol=m)
  #for (i in 1:m)  # generate the data with certain variants full of zeros
  #  x1x0[,i]=rpois((N0+N1), q[i])
  #filter.x1x0=x1x0[,!apply(x1x0==0,2,all)]  # filter variants with 0 counts in both cases and controls since these variants are noninformative
  #mm=ncol(filter.x1x0)
  #var.index=which(colSums(x1x0!=0)>0) # get the variant index in the original data before filtering
  #############
  #filter.x1x0=rep(NA, (N0+N1))
  #var.index=numeric(); ii=1
  #for (i in 1:m)
  #  {
  #   x1x0=rpois((N0+N1), q[i])
  #   if (sum(x1x0)>0)
  #   {
  #     ii=ii+1
  #     var.index[ii]=i
  #   }
  #   filter.x1x0=cbind(filter.x1x0, x1x0)
  #  }
  ################## only save the counts 
  total.count=numeric()
  for (i in 1:m)  # generate the data with certain variants full of zeros
    total.count[i]=sum(rpois((N0+N1), q[i]))
  var.index=which(total.count>0) # get the variant index in the original data before filtering
  mm=length(var.index)
  filter.count=total.count[var.index]
  
  gamma <- rep(1,mm)
  Zij=rep(0,mm)
  qq=rbeta(mm, alpha0, beta0)
  if (model==1)
  {
    #  for (k in 1:num.group) # when one gene has multiple annotation groups
    #  {
    #    from=split.ratio[k]*mm+1; to=split.ratio[k+1]*mm
    #    Zij[from:to]=rbinom((to-from+1), 1, pi[k])
    #  }
    for (i in 1:mm)
      Zij[i]=rbinom(1,1,pi[var.index[mm]])
    
    qq[Zij==1] <- rbeta(sum(Zij==1), alpha, beta)
    gamma[Zij==1] <- rgamma(sum(Zij==1), gamma.mean*sigma, sigma)
  }
  
 # x <- array(0, dim=c(length(pheno), mm)) # store individual genotype 
  x <- array(0, dim=c(2, mm))             # store the total counts in case and controls
  for (j in 1:mm)  # re-distribute generated total counts x1+x0 into cases and controls by conditional distribution
  {
    #if (Zij[j]==1)
    #  filter.count[j]=sum(rpois((N0+N1), qq[j]))
    #x1=rbinom(1, sum(filter.x1x0[,j]), N1*gamma[j]/(gamma[j]*N1+N0))
    x1=rbinom(1, filter.count[j], N1*gamma[j]/(gamma[j]*N1+N0))
    #x1=sum(filter.x1x0[,j])-x0
    x0=filter.count[j]-x1
    #x[sample.int(N0,x0),j]=1
    #x[N0+sample.int(N1,x1),j]=1
    x[1,j]=x0
    x[2,j]=x1
  }
  x=as.matrix(x)
  return (list(geno=x,pheno=pheno,q=q,gamma=gamma, Zij=Zij, var.index=var.index))
  
}
################################################
intergrand=function(aa, var.case, var.contr, bar.gamma, sig, N1, N0)
{
  ff=dbinom(var.case, sum(var.case, var.contr), aa*N1/(aa*N1+N0))*dgamma(aa, bar.gamma*sig, sig)
  return(ff)
}
# calculate the bayes factor of a single variant via integration
BF.var.inte=function(var.case, var.contr, bar.gamma, sig, N1, N0)
{
  marglik0.CC <- dbinom(var.case, sum(var.case, var.contr), N1/(N1+N0), log=T)    # Under H0: gamma=1
  
  marglik1.CC <- integrate(intergrand, var.case, var.contr, bar.gamma, sig, N1, N0, low=0, upper=100, stop.on.error=F)$value # Under H1: gamma~gamma(gamma.mean*sigma, sigma) 
  #BF.var=log(marglik1.CC)/marglik0.CC
  BF.var <- exp(min(log(marglik1.CC)-marglik0.CC, 709)) # avoid overflow 
  
  
  return(BF.var)
}
################################################
objec.func=function(beta.est) # this is the log likelihood function to be maximized
{
  lkhd=0
  for (i in 1:num.var)
    #lkhd=lkhd+log(1+BF.var[i]*exp(Ajk.effect[i]*beta.est))-log(1+exp(Ajk.effect[i]*beta.est))
     lkhd=lkhd+min(exp(709), log(1+BF.var[i]*exp(Ajk.effect[i]*beta.est))-log(1+exp(Ajk.effect[i]*beta.est))) # avoid overflow
  lkhd=lkhd-lambda/2*sum(beta.est^2)
  
  return(lkhd)  
}

objec.func.min=function(beta.est) # this is the log likelihood function to be maximized
{
  lkhd=0
  for (i in 1:num.var)
    lkhd=lkhd-log(1+BF.var[i]*exp(Ajk.effect[i]*beta.est))+log(1+exp(Ajk.effect[i]*beta.est))
  
  lkhd=lkhd+lambda/2*sum(beta.est^2)
  
  return(lkhd)  
}

deriv.objec.func=function(beta.est) # derivative function of parameters 
{
  par.beta=numeric(anno.num)
  for (k in 1:anno.num) # k: annotattion index
  {
    par.beta.k=0
    for (j in 1:num.var) # j: variant index
    {
      delta.j=exp(Ajk.effect[j]*beta.est)
      par.beta.k=par.beta.k+Ajk.effect[j]*delta.j*(BF.var[j]/(1+delta.j*BF.var[j])-1/(1+delta.j))
    }  
    par.beta[k]=par.beta.k-lambda*beta.est[k] # k th element of first derivative
  }
  
  return(par.beta)
}
################################################
num.gene=1000
m=200
N0=10000; N1=10000
delta=1
alpha0 <- 0.1
beta0 <- 1000
alpha <- 0.1
beta <- 2000
gamma.mean <- 10
sigma <- 1
split.ratio=c(0,1)
num.group=length(split.ratio)-1
anno.num=1
############  generate feature covariate 
Ajk=rep(0, (num.gene*m))
Ajk[sample(length(Ajk), (1*length(Ajk)))]=1
############## generate beta 
beta.true=numeric(anno.num)
beta.true[1]=2
############# compute tau: risk of every variant being risk variant
tau.true=exp(Ajk*beta.true)/(1+exp(Ajk*beta.true))
######################
max.rep=10
all.beta=numeric(max.rep)
ratio=seq(0, 1, by=0.05)
#ratio=0
for (rep in 1:max.rep)
{  
  
  cat(rep, "th replicate is running", "\n")  
  Ui=rbinom(num.gene, 1, delta)
  all.data=list(); var.orig.index=numeric(); var.orig.index[1]=NA
  for (i in 1:num.gene)
  {
    data=gene.simu(N0, N1, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, pi=tau.true[((i-1)*m+1):(i*m)], Ui[i], num.group, split.ratio)
    all.data[[i]]=data
    var.orig.index=c(var.orig.index, (i-1)*m+data$var.index)
  }
  var.orig.index=var.orig.index[-1]
  Ajk.effect=Ajk[var.orig.index]
  ################### calculate bayes factor for every variant 
  k=0; BF.var=numeric(); var.count=matrix(nrow=(m*num.gene), ncol=2)
  for (i in 1:length(all.data))
  {
    for (j in 1:ncol(all.data[[i]]$geno))
    {
      k=k+1
      BF.var[k]=BF.var.inte(all.data[[i]]$geno[2,j],all.data[[i]]$geno[1,j], bar.gamma=6, sig=sigma, N1, N0)
      var.count[k,]=c(all.data[[i]]$geno[2,j],all.data[[i]]$geno[1,j] )
    }
  } # end of i
 var.count=var.count[-((k+1):nrow(var.count)),]
 ###############################
 ################################
  num.var=length(BF.var)
  theta.est=list()
  #theta.est[[1]]=runif(anno.num, -1,1)
  #lambda.range=seq(0, 5000, by=100)
  lambda.range=0
  #beta.est=theta.est[[1]]
  error=numeric()
  beta.fin.est=numeric(length(lambda.range))
  #beta.est=beta.true
###################################
#  lambda=0
#  beta.range=seq(0, 1, by=0.01)
#  log.lkhd=numeric(length(beta.range))
#  for (i in 1:length(beta.range))
#    log.lkhd[i]=objec.func(beta.range[i])
  
################################### use GLM 
  Y=c(NA)
  for (i in 1:length(all.data))
    Y=c(Y, all.data[[i]]$Zij)
  Response=Y[-1]
  fit <- glm(Response~0+Ajk.effect,family=binomial())
 all.beta[rep]=fit$coefficients 
  
###################################  
  
#  for (k in 1:length(lambda.range))
#  {
#    cat(k, "is running", "\n")  
#    lambda=lambda.range[k] 
#    beta.est=runif(anno.num, 1,3)
#    theta.est=optim(beta.est, objec.func, deriv.objec.func, method="BFGS", control=list(fnscale=-1))
#    beta.fin.est[k]=theta.est$par
#    theta.est=DEoptim(fn=objec.func.min, lower=rep(-1, anno.num), upper=rep(200, anno.num), control=list(NP=100, itermax=100,trace=FALSE))
#    beta.fin.est[k,]=theta.est$optim$bestmem
#    error[k]=sum((beta.fin.est[k]-beta.true)^2)
    
#  }
#  plot(beta.true, ylim=c(-1,1), type="o")
#  lines(beta.fin.est[which.min(error)], ylim=c(-1,1), type="o", col=2)
  #beta.fin.est
#  all.beta[rep]=beta.fin.est[which.min(error)]
} # end of rep 