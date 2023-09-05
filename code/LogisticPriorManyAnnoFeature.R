# Ui; latent binary variant for gene , Ui~Ber(delta) Ui=1: risk gene Ui=0: non risk gene
# Zij: latent binary variant for variant, Zij~Ber(pi), Zij=1, variant j of gene i is causal variant, Zij=0 otherwise
# m: number of variants in each gene. 
# N0: number of controls; N1: number of cases
# num.gene: number of genes in the sample
# q: allele frequency q~Beta(alpha0, beta0)
library(DEoptim)
#################################
rm(list=ls())
gene.simu=function(N0, N1, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, pi, model, num.group, split.ratio)
{
  pheno=c(rep(0,N0), rep(1,N1))  # filtering step
  q <- rbeta(m, alpha0, beta0)
  #x1x0=matrix(nrow=(N0+N1), ncol=m)
  #for (i in 1:m)  # generate the data with certain variants full of zeros
  #  x1x0[,i]=rpois((N0+N1), q[i])
  #filter.x1x0=x1x0[,!apply(x1x0==0,2,all)]  # filter variants with 0 counts in both cases and controls since these variants are non-informative
  #mm=ncol(filter.x1x0)
  #var.index=which(colSums(x1x0!=0)>0) # get the variant index in the original data before filtering
  ################## only save the counts 
  total.count=numeric()
  for (i in 1:m)  # generate the data with certain variants full of zeros that will be filtered out later 
    total.count[i]=sum(rpois((N0+N1), q[i]))
  var.index=which(total.count>0) # get the variant index in the original data before filtering
  mm=length(var.index) # number of effective variants with at least one count 
  filter.count=total.count[var.index]
  ##################
  #mm=m
  gamma <- rep(1,mm)
  Zij=rep(0,mm)
  qq=rbeta(mm, alpha0, beta0)
  if (model==1)
  {
    #  for (k in 1:num.group)
    #  {
    #    from=split.ratio[k]*mm+1; to=split.ratio[k+1]*mm
    #    Zij[from:to]=rbinom((to-from+1), 1, pi[k])
    #  }
    for (i in 1:mm)
      Zij[i]=rbinom(1,1,pi[var.index[mm]])
    
    qq[Zij==1] <- rbeta(sum(Zij==1), alpha, beta)
    gamma[Zij==1] <- rgamma(sum(Zij==1), gamma.mean*sigma, sigma)
  }
  
  x <- array(0, dim=c(2, mm)) 
  for (j in 1:mm)
  {
    #x1=rbinom(1, sum(filter.x1x0[,j]), N1*gamma[j]/(gamma[j]*N1+N0))
    #x1=sum(filter.x1x0[,j])-x0
    #x0=sum(filter.x1x0[,j])-x1
    #x[sample.int(N0,x0),j]=1
    #x[N0+sample.int(N1,x1),j]=1
    
    x1=rbinom(1, filter.count[j], N1*gamma[j]/(gamma[j]*N1+N0))
    x0=filter.count[j]-x1
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
  marglik0.CC <- dbinom(var.case, sum(var.case, var.contr), N1/(N1+N0))    # Under H0: gamma=1
  
  marglik1.CC <- integrate(intergrand, var.case, var.contr, bar.gamma, sig, N1, N0, low=0, upper=100, stop.on.error=F)$value # Under H1: gamma~gamma(gamma.mean*sigma, sigma) 
  BF.var <- marglik1.CC/marglik0.CC
  
  return(BF.var)
}
################################################
objec.func=function(beta.est) # this is the log likelihood function to be maximized
{
  lkhd=0
  for (i in 1:num.var)
    lkhd=lkhd+log(1+BF.var[i]*exp(Ajk.effect[i,]%*%beta.est))-log(1+exp(Ajk.effect[i,]%*%beta.est))
  
  lkhd=lkhd-lambda/2*sum(beta.est^2)
  
  return(lkhd)  
}

objec.func.min=function(beta.est) # this is the log likelihood function to be maximized
{
  lkhd=0
  for (i in 1:num.var)
    lkhd=lkhd-log(1+BF.var[i]*exp(Ajk.effect[i,]%*%beta.est))+log(1+exp(Ajk.effect[i,]%*%beta.est))
  
  lkhd=lkhd+lambda/2*sum(beta.est^2)
  
  return(lkhd)  
}

deriv.objec.func=function(beta.est)
{
  par.beta=numeric(anno.num)
  for (k in 1:anno.num) # k: annotation index
  {
    par.beta.k=0
    for (j in 1:num.var) # j: variant index
    {
      delta.j=exp(Ajk.effect[j,]%*%beta.est)
      
      par.beta.k=par.beta.k+Ajk.effect[j,k]*delta.j*(BF.var[j]/(1+delta.j*BF.var[j])-1/(1+delta.j))
    }  
    par.beta[k]=par.beta.k-lambda*beta.est[k] # k th element of first derivative
  }
  
  return(par.beta)
}
################################################
num.gene=1000  # number of genes 
m=200  # number of variants 
N0=500000; N1=500000
delta=1
alpha0 <- 0.1
beta0 <- 1000
alpha <- 0.1
beta <- 2000
gamma.mean <- 10
sigma <- 1
split.ratio=c(0,1)
num.group=length(split.ratio)-1
anno.num=1 # number of annotation groups 
############  generate feature covariates 
Ajk=matrix(0, nrow=(num.gene*m), ncol=anno.num) # Ajk are indicator matrix with elements of being 0, or 1
for (i in 1:ncol(Ajk))
  Ajk[sample(nrow(Ajk), (0.5*nrow(Ajk))),i]=1 # half of variants across all genes are "1", i.e. risk variants 
############## generate beta 
beta.true=numeric(anno.num) # beta coefficients for all annotation groups 
beta.true[1]=1
#beta.true[2]=2
#beta.true[3]=3.5
#beta.true[4:anno.num]=runif(anno.num-4+1)
############# compute eta: risk of every variant being a  risk variant
eta.true=exp(Ajk%*%beta.true)/(1+exp(Ajk%*%beta.true)) 
# eta.true is variant specific, could be very big number 
######################
max.rep=10
all.beta=matrix(nrow=max.rep, ncol=anno.num)
for (rep in 1:max.rep)
{  
  ############## generate the sample 
  cat(rep, "th replicate is running", "\n")  
  Ui=rbinom(num.gene, 1, delta)
  all.data=list(); var.orig.index=numeric(); var.orig.index[1]=NA
  for (i in 1:num.gene)
  {
    data=gene.simu(N0, N1, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, pi=eta.true[((i-1)*m+1):(i*m)], Ui[i], num.group, split.ratio)
    all.data[[i]]=data
    var.orig.index=c(var.orig.index, (i-1)*m+data$var.index)
  }
  var.orig.index=var.orig.index[-1]
  Ajk.effect=Ajk[var.orig.index, ]
  ################### calculate Bayes factor for every variant 
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
  num.var=length(BF.var)
  theta.est=list()
  lambda.range=0
  error=numeric()
  beta.fin.est=matrix(nrow=length(lambda.range), ncol=anno.num)
  ####################### use glm to estimate beta
  Y=c(NA)
  for (i in 1:length(all.data))
    Y=c(Y, all.data[[i]]$Zij)
  Response=Y[-1]
  #fit <- glm(Response~0+Ajk.effect[,1]+Ajk.effect[,2],family=binomial())
  fit <- glm(Response~0+Ajk.effect,family=binomial())
  #fit <- glm(Response~0+Ajk.effect[,1]+Ajk.effect[,2]+Ajk.effect[,3]+Ajk.effect[,4]+Ajk.effect[,5]+Ajk.effect[,6]+Ajk.effect[,7]
  #           +Ajk.effect[,8]+Ajk.effect[,9]+Ajk.effect[,10],family=binomial())
  all.beta[rep,]=fit$coefficients 
  #######################  use optim to estimate beta 
 # for (k in 1:length(lambda.range))
  {
#    cat(k, "is running", "\n")  
#    lambda=lambda.range[k] 
#    beta.est=runif(anno.num, 0,3)
#    theta.est=optim(beta.est, objec.func, deriv.objec.func, method="BFGS", control=list(fnscale=-1))
#    beta.fin.est[k,]=theta.est$par
    #theta.est=DEoptim(fn=objec.func.min, lower=rep(-1, anno.num), upper=rep(5, anno.num), control=list(NP=100, itermax=100,trace=FALSE))
    #beta.fin.est[k,]=theta.est$optim$bestmem
#    error[k]=sum((beta.fin.est[k,]-beta.true)^2)
    
  }  # end of for (k in 1:length(lambda.range))
  
  #plot(beta.true, ylim=c(-1,1), type="o")
  #lines(beta.fin.est[which.min(error),], ylim=c(-1,1), type="o", col=2)
  #beta.fin.est
#  all.beta[rep,]=beta.fin.est[which.min(error),]
} # end of rep 
################################
#max.num=100
#count=matrix(nrow=max.num, ncol=2)
#count[,1]=seq(1,max.num); count[,2]=15*count[,1]
#BF=numeric()
#for (i in 1:nrow(count))
#  BF[i]=BF.var.inte(count[i,2], count[i,1], 6, 1, 100000, 100000)
#par(mfrow=c(1,2))
#plot(BF[1:50])
#plot(BF[51:100])