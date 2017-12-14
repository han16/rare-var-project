rm(list=ls())
set.seed(123)
#################################################################
# Simulation Functions
#################################################################

# Simulation of gene-level data
# N0, N: sample sizes of controls and cases
# m: number of variants
# model: 1 - causal gene (alternative); 0 - non-causal gene (null)
# Causal variants: q ~ Beta(alpha, beta)
# Non-causal variants: q ~ Beta(alpha0, beta0)
# Relative risk of causal variants: RR ~ Gamma(gamma.mean*sigma, sigma)
# pi: fraction of causal variants
simul.gene <- function(N0, N, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, pi=1, model=1) {
  # phenotype data
  y <- c(rep(0,N0), rep(1,N))

  # relative risk of each variants
  q <- rbeta(m, alpha0, beta0)
  gamma <- rep(1,m)
  if (model==1) {
    Z <- rbinom(m,1,pi)
    q[Z==1] <- rbeta(sum(Z==1), alpha, beta)
    gamma[Z==1] <- rgamma(sum(Z==1), gamma.mean*sigma, sigma)
  }

  # genotype data
  x <- array(0, dim=c(length(y), m))
  for (j in 1:m) {
    x[y==0,j] <- sample(c(1,0), N0, replace=TRUE, prob=c(2*q[j],1-2*q[j]))
    x[y==1,j] <- sample(c(1,0), N, replace=TRUE, prob=c(2*q[j]*gamma[j],1-2*q[j]*gamma[j]))
  }

  # TO DO: filter variants with (0,0) counts
  
  x=x[,!apply(x==0,2,all)]  # filter variants with 0 counts in both cases and controls since these kind of variants are noninformative
  x=as.matrix(x)

  return (list(x=x,y=y,q=q,gamma=gamma))
}

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
###################################
# the test statistic
burden.test.stat <- function(x, y, MAF.thr) {
  # weights of the variants based on MAF
  n <- length(y)
  MAF <- colSums(x, na.rm=T) / n
  weights <- (MAF < MAF.thr)

  # test statistic: genetic burden in cases
  x.case <- as.matrix(x[y==1,])
  x.var <- colSums(x.case, na.rm=T)
  T <- sum(x.var*weights)

  return (list(weights=weights, stat=T))
}

# Burden test with fixed threshold
# x: genotype matrix (n by m)
# y: phenotype
# MAF.thr: MAF threshold, default 1% (MAF from both cases and controls)
# nperm: number of permutations
burden.test <- function(x, y, MAF.thr=0.01, nperm=1000) {
  # observe test statistic
  T.obs <- burden.test.stat(x,y,MAF.thr)$stat

  # computing the p-value by permutation
  T.perm <- length(nperm)
  for (i in 1:nperm) {
    y.perm <- sample(y, length(y))
    T.perm[i] <- burden.test.stat(x,y.perm,MAF.thr)$stat
  }
  pval <- sum(T.perm >= T.obs) / nperm

  return (list(stat=T.obs, p.value=pval))
}

# SKAT
SKAT.run <- function(x, y, SKATO=FALSE) {
  obj<-SKAT_Null_Model(y ~ 1, out_type="D")
  if (SKATO==FALSE)  pval <- SKAT(x, obj)$p.value
  else pval <- SKAT(x, obj, method="optimal")$p.value

  return (list(p.value=pval))
}
###################################
# simulation of data and evaluate the performance of the Bayesian method
simulation.BF <- function(N0, N, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, ngenes=500, pi=1,
alpha.model=0.5, beta.model=200, gamma.mean.model, sigma.model=1, pi.model=1) {
  BF.pos <- length(ngenes)
  BF.neg <- length(ngenes)
  for (i in 1:ngenes) {
    data <- simul.gene(N0, N, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, pi=pi, model=0)

    BF.neg.prod=1
    for (j in 1:ncol(data$x))
      BF.neg.prod=BF.neg.prod*BF.gene.inte(data$x[,j],data$y,gamma.mean.model,sigma.model)
    BF.neg[i]=BF.neg.prod
      
   # BF.neg[i] <- BF.gene.inte(data$x,data$y,alpha.model,beta.model,gamma.mean.model,sigma.model,pi=pi.model)

    data <- simul.gene(N0, N, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, pi=pi, model=1)
    BF.pos.prod=1
    for (j in 1:ncol(data$x))
      BF.pos.prod=BF.pos.prod*BF.gene.inte(data$x[,j],data$y,gamma.mean.model,sigma.model)
    BF.pos[i]=BF.pos.prod
    #BF.pos[i] <- BF.gene.inte(data$x,data$y,alpha.model,beta.model,gamma.mean.model,sigma.model,pi=pi.model)
  }

  return (list(neg=BF.neg, pos=BF.pos))
}

# simulation of data and evaluate the performance of the fixed threshold method
simulation.FT <- function(N0, N, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, ngenes=500, pi=1, MAF.thr=0.01,nperm=1000) {
  FT.pos <- length(ngenes)
  FT.neg <- length(ngenes)
  for (i in 1:ngenes) {
    data <- simul.gene(N0, N, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, pi=pi, model=0)
    FT.neg[i] <- burden.test(data$x, data$y, MAF.thr=MAF.thr,nperm=nperm)$p.value

    data <- simul.gene(N0, N, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, pi=pi,model=1)
    FT.pos[i] <- burden.test(data$x, data$y, MAF.thr=MAF.thr,nperm=nperm)$p.value
  }

  return (list(neg=-log(FT.neg), pos=-log(FT.pos)))
}

#BF.gene=function(x, y, gamma.mean, beta, S, pi.model){
# BF.var=numeric()
# for (i in 1:ncol(x))
#   {
#    marglik0.CC <- dbinom(sum(x[y==1,i]), sum(x[,i]), N/(N+N0))    # Under H0: gamma=1
#
#    gamma <- rgamma(S, (gamma.mean*beta), beta)                 # Under H1: gamma~gamma(gamma.mean*sigma, sigma)
#    marglik1.CC <- numeric(S)
#    for (j in 1:S)
#       marglik1.CC[j] <- dbinom(sum(x[y==1,i]), sum(x[,i]), gamma[j] * N/(gamma[j] * N+N0))
#    BF.var[i] <- mean(marglik1.CC) / marglik0.CC
#   }
#
# return(prod((1-pi.model)+pi.model*BF.var))
#
#}


#############################
# simulation of data and evaluate the performance of SKAT
library(SKAT)
simulation.SKAT <- function(N0, N, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, ngenes=500, pi=1) {
  SKAT.pos <- length(ngenes)
  SKAT.neg <- length(ngenes)
  for (i in 1:ngenes) {
    data <- simul.gene(N0, N, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, pi=pi, model=0)
    SKAT.neg[i] <- SKAT.run(data$x, data$y)$p.value

    data <- simul.gene(N0, N, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, pi=pi,model=1)
    SKAT.pos[i] <- SKAT.run(data$x, data$y)$p.value
  }

  return (list(neg=-log(SKAT.neg), pos=-log(SKAT.pos)))
}

#################################################################
# Simulation: Basic Model
#################################################################

# simulation/true parameters
N0 <- 3000
N <- 3000
N1=N
alpha0 <- 1
beta0 <- 1000
alpha <- 1
beta <- 2000
gamma.mean <- 3
sigma <- 1
ngenes <- 500
pi <- 0.1
m <- 100

# Bayesian method        # parameters in the model 
alpha.model <- 0.2 
beta.model <- 200
gamma.mean.model <- 3
sigma.model <- 1
pi.model <- pi
BF.distr <- simulation.BF(N0, N, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, ngenes=ngenes, pi=pi, 
alpha.model=alpha.model, beta.model=beta.model, gamma.mean.model=gamma.mean.model, sigma.model=sigma.model, pi.model=pi.model)
BF.pos <- BF.distr$pos
BF.neg <- BF.distr$neg
AUC.BF <- mean(sample(log(BF.pos), 1000, replace=TRUE) > sample(log(BF.neg),1000,replace=TRUE))
summary(BF.pos)
summary(BF.neg)
mean(BF.pos>2)
mean(BF.neg>2)
AUC.BF

# Fixed threshold method
MAF.thr <- 0.01
nperm <- 1000
FT.distr <- simulation.FT(N0, N, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, ngenes=ngenes, pi=pi, MAF.thr=MAF.thr,nperm=nperm)
FT.pos <- FT.distr$pos
FT.neg <- FT.distr$neg

# SKAT
SKAT.distr <- simulation.SKAT(N0, N, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, ngenes=ngenes, pi=pi)
SKAT.pos <- SKAT.distr$pos
SKAT.neg <- SKAT.distr$neg

# Bayesian results
#postscript(file="simulation/BF_distr.eps",height=6,width=6)
BF.distr <- data.frame(model=c(rep('Non-risk gene',ngenes), rep('Risk gene',ngenes)), BF=c(log(BF.neg), log(BF.pos)))
boxplot(BF~model, data=BF.distr, ylab="log(BF)")
dev.off()
summary(log(BF.neg))
summary(log(BF.pos))


library(ROCR)
plotROC <- function(scores.pos, scores.neg) {
  npos <- length(scores.pos)
  nneg <- length(scores.neg)
  pred <- prediction(c(scores.pos,scores.neg), c(rep(1,npos),rep(0,nneg)))
  perf <- performance(pred,"tpr", "fpr")
  return (perf)
}

# TPR at a certain alpha (FPR)
TPR <- function(scores.pos, scores.neg, alpha) {
  scores.neg <- sort(scores.neg, decreasing=TRUE)
  cutoff <- scores.neg[as.integer(length(scores.neg)*alpha)]
  TPR <- sum(scores.pos > cutoff)/length(scores.pos)
  return (TPR)
}



# Comparison of methods
#postscript(file="simulation/AUC_risk_only.eps",height=6,width=6)
plot(plotROC(BF.pos,BF.neg))
plot(plotROC(FT.pos,FT.neg),add=TRUE,col='red')
plot(plotROC(SKAT.pos,SKAT.neg),add=TRUE,col='blue')
legend("bottomright", c("Bayes","Burden", "SKAT"), col=c('black','red','blue'), lty=1)
dev.off()
AUC.FT <- mean(sample(FT.pos, 1000, replace=TRUE) > sample(FT.neg,1000,replace=TRUE))
AUC.SKAT <- mean(sample(SKAT.pos, 1000, replace=TRUE) > sample(SKAT.neg,1000,replace=TRUE))
AUC.BF
AUC.FT
AUC.SKAT
FPR <- 0.01
TPR(BF.pos, BF.neg, FPR)
TPR(FT.pos, FT.neg, FPR)
TPR(SKAT.pos, SKAT.neg, FPR)
save(list=ls(), file="Comparewithothers.RData")
# for (i in 1:ngenes) {
#   # negative genes
#   data <- simul.gene(N0, N, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, pi=pi,model=0)
#   BF.neg[i] <- BF.gene(data$x,data$y,alpha.model,beta.model,gamma.mean,sigma,pi=pi.model,S=1000)$BF
#   FT.neg[i] <- burden.test(data$x, data$y, MAF.thr=0.01,nperm=nperm)$p.value
#   SKAT.neg[i] <- SKAT.run(data$x, data$y)$p.value

#   # positive genes
#   data <- simul.gene(N0, N, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, pi=pi,model=1)
#   BF.pos[i] <- BF.gene(data$x,data$y,alpha.model,beta.model,gamma.mean,sigma,pi=pi.model,S=1000)$BF
#   FT.pos[i] <- burden.test(data$x, data$y, MAF.thr=0.01,nperm=nperm)$p.value
#   SKAT.pos[i] <- SKAT.run(data$x, data$y)$p.value
# }
#################################################################
# Simulation: Power Analysis as a Function of m and pi
#################################################################

