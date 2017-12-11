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

  return (list(x=x,y=y,q=q,gamma=gamma))
}

# simulation of data and evaluate the performance of the Bayesian method
simulation.BF <- function(N0, N, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, ngenes=500, pi=1, alpha.model=0.5, beta.model=200, gamma.mean.model=1, sigma.model=1, pi.model=1) {
  BF.pos <- length(ngenes)
  BF.neg <- length(ngenes)
  for (i in 1:ngenes) {
    data <- simul.gene(N0, N, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, pi=pi, model=0)
    BF.neg[i] <- BF.gene(data$x,data$y,alpha.model,beta.model,gamma.mean.model,sigma.model,pi=pi.model,S=1000)$BF

    data <- simul.gene(N0, N, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, pi=pi, model=1)
    BF.pos[i] <- BF.gene(data$x,data$y,alpha.model,beta.model,gamma.mean.model,sigma.model,pi=pi.model,S=1000)$BF
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

# simulation parameters
N0 <- 1000
N <- 1000
alpha0 <- 1
beta0 <- 1000
alpha <- 1
beta <- 2000
gamma.mean <- 2
sigma <- 4
ngenes <- 200
pi <- 0.2
m <- 5

# Bayesian method
alpha.model <- 0.2  # used in the model
beta.model <- 200
gamma.mean.model <- 2
sigma.model <- 5
pi.model <- pi
BF.distr <- simulation.BF(N0, N, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, ngenes=ngenes, pi=pi, alpha.model=alpha.model, beta.model=beta.model, gamma.mean.model=gamma.mean.model, sigma.model=sigma.model, pi.model=pi.model)
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
postscript(file="simulation/BF_distr.eps",height=6,width=6)
BF.distr <- data.frame(model=c(rep('Non-risk gene',ngenes), rep('Risk gene',ngenes)), BF=c(log(BF.neg), log(BF.pos)))
boxplot(BF~model, data=BF.distr, ylab="log(BF)")
dev.off()
summary(log(BF.neg))
summary(log(BF.pos))

# Comparison of methods
postscript(file="simulation/AUC_risk_only.eps",height=6,width=6)
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

# for (i in 1:ngenes) {
#   # negative genes
#   data <- simul.gene(N0, N, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, pi=pi,model=0)
#   BF.neg[i] <- BF.gene(data$x,data$y,alpha.model,beta.model,gamma.mean,sigma,pi=pi.model,S=1000)$BF
#   FT.neg[i] <- burden.test(data$x, data$y, MAF.thr=0.01,nperm=nperm)$p.value
#   SKAT.neg[i] <- SKAT.run(data$x, data$y)$p.value
#
#   # positive genes
#   data <- simul.gene(N0, N, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, pi=pi,model=1)
#   BF.pos[i] <- BF.gene(data$x,data$y,alpha.model,beta.model,gamma.mean,sigma,pi=pi.model,S=1000)$BF
#   FT.pos[i] <- burden.test(data$x, data$y, MAF.thr=0.01,nperm=nperm)$p.value
#   SKAT.pos[i] <- SKAT.run(data$x, data$y)$p.value
# }
#################################################################
# Simulation: Power Analysis as a Function of m and pi
#################################################################

# simulation parameters
N <- 1000
N0 <- 5000
alpha0 <- 1
beta0 <- 1000
alpha <- 1
beta <- 2000
gamma.mean <- 2
sigma <- 4
ngenes <- 500
pi <- 0.5 #seq(from=0.1,to=0.5,by=0.1)
m <- 5 #c(5,10,20)

# Bayesian method
alpha.model <- 0.2  # used in the model
beta.model <- 200
gamma.mean.model <- 2
sigma.model <- 4
AUC.BF <- array(0,dim=c(length(pi),length(m)))
for (i in 1:length(pi)) {
  for (j in 1:length(m)) {
    pi.model <- pi[i]
    BF.distr <- simulation.BF(N0, N, m[j], alpha0, beta0, alpha, beta, gamma.mean, sigma, ngenes=ngenes, pi=pi[i], alpha.model=alpha.model, beta.model=beta.model, gamma.mean.model=gamma.mean.model, sigma.model=sigma.model, pi.model=pi.model)
    BF.pos <- BF.distr$pos
    BF.neg <- BF.distr$neg
    AUC.BF[i,j] <- mean(sample(log(BF.pos), 1000, replace=TRUE) > sample(log(BF.neg),1000,replace=TRUE))
  }
}
AUC.BF

postscript(file="simulation/AUC_mis3_simulation.eps",width=6, height=6)
colors <- 2:4
pchs <- c(0:3)
ltys <- c(1:4)
lwd <- 1.5
plot(pi, AUC.BF[,1], type='b', col=colors[1], pch=pchs[1], lty=ltys[1], lwd=lwd,  xlab="Fraction of causal variants (pi)", ylab="AUC", ylim=c(0.5,0.7))
lines(pi, AUC.BF[,2], type='b', col=colors[2], pch=pchs[2], lty=ltys[2], lwd=lwd)
lines(pi, AUC.BF[,3], type='b', col=colors[3], pch=pchs[3], lty=ltys[3], lwd=lwd)
legend("topleft", c("m=5", "m=10", "m=20"), col=colors, lty=ltys, pch=pchs, lwd=lwd)
dev.off()

#################################################################
# Simulation: Power Analysis as a Function of Sample Size
#################################################################

# simulation parameters
N <- 1000
N0 <- 1000 * 1:10
alpha0 <- 1
beta0 <- 1000
alpha <- 1
beta <- 2000
gamma.mean <- 2
sigma <- 4
ngenes <- 500
pi <- 0.5 #seq(from=0.1,to=0.5,by=0.1)
m <- c(5,10,20)

# Bayesian method
alpha.model <- 0.2  # used in the model
beta.model <- 200
gamma.mean.model <- gamma.mean
sigma.model <- sigma
pi.model <- pi
AUC.BF <- array(0,dim=c(length(N0),length(m)))
for (i in 1:length(N0)) {
  for (j in 1:length(m)) {
    BF.distr <- simulation.BF(N0[i], N, m[j], alpha0, beta0, alpha, beta, gamma.mean, sigma, ngenes=ngenes, pi=pi, alpha.model=alpha.model, beta.model=beta.model, gamma.mean.model=gamma.mean.model, sigma.model=sigma.model, pi.model=pi.model)
    BF.pos <- BF.distr$pos
    BF.neg <- BF.distr$neg
    AUC.BF[i,j] <- mean(sample(log(BF.pos), 1000, replace=TRUE) > sample(log(BF.neg),1000,replace=TRUE))
  }
}
AUC.BF

postscript(file="simulation/AUC_mis3_simulation_controls.eps",width=6, height=6)
colors <- 2:4
pchs <- c(0:3)
ltys <- c(1:4)
lwd <- 1.5
plot(N0, AUC.BF[,1], type='b', col=colors[1], pch=pchs[1], lty=ltys[1], lwd=lwd,  xlab="Sample size of controls", ylab="AUC", ylim=c(0.6, 1.0))
lines(N0, AUC.BF[,2], type='b', col=colors[2], pch=pchs[2], lty=ltys[2], lwd=lwd)
lines(N0, AUC.BF[,3], type='b', col=colors[3], pch=pchs[3], lty=ltys[3], lwd=lwd)
legend("topleft", c("m=5", "m=10", "m=20"), col=colors, lty=ltys, pch=pchs, lwd=lwd)
dev.off()
#################################################################
# Systematic Simulation: Sensitivity Analysis
#################################################################

# Dominated by risk variants, or significant fraction of protective variants
N0 <- 2000
N <- 2000
m <- 10
alpha0 <- 1
beta0 <- 1000
alpha <- 1
beta <- 1000
ngenes <- 1000
gamma.mean <- 1.5
sigma <- 1
pi <- 1
alpha.model <- 0.5  # used in the model
beta.model <- 200
gamma.mean.model <- 1.5 #c(1.2,1.4,1.6,1.8,2.0)
sigma.model <- c(0.6,0.8,1.0,1.2,1.4)
pi.model <- 1
# npar <- length(gamma.mean.model)
npar <- length(sigma.model)
AUC.BF <- length(npar)
TPR.BF <- length(npar)
FPR <- 0.01
for (i in 1:npar) {
  BF.distr <- simulation.BF(N0, N, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, ngenes=ngenes, pi=pi, alpha.model=alpha.model, beta.model=beta.model, gamma.mean.model=gamma.mean.model, sigma.model=sigma.model[i], pi.model=pi.model)
  BF.pos <- BF.distr$pos
  BF.neg <- BF.distr$neg
  AUC.BF[i] <- mean(sample(log(BF.pos), 1000, replace=TRUE) > sample(log(BF.neg),1000,replace=TRUE))
  TPR.BF[i] <- TPR(BF.pos, BF.neg, FPR)
}
# report <- data.frame(gamma.mean=gamma.mean.model, AUC=AUC.BF, TPR=TPR.BF)
report <- data.frame(sigma=sigma.model, AUC=AUC.BF, TPR=TPR.BF)
write.csv(report, file="simulation/sensitivity_risk_only.csv")

# Mixture of risk and non-functional variants
N0 <- 2000
N <- 2000
m <- 20
alpha0 <- 1
beta0 <- 1000
alpha <- 1
beta <- 5000
ngenes <- 1500
gamma.mean <- 5
sigma <- 1
pi <- 0.2
alpha.model <- 0.5  # used in the model
beta.model <- 200
gamma.mean.model <- 4 #c(1.2,1.4,1.6,1.8,2.0)
sigma.model <- 1.3 #c(0.6,0.8,1.0,1.2,1.4)
pi.model <- 0.4 #c(0.1,0.2,0.3,0.4,0.5)
npar <- length(gamma.mean.model)
# npar <- length(sigma.model)
# npar <- length(pi.model)
AUC.BF <- length(npar)
TPR.BF <- length(npar)
FPR <- 0.01
for (i in 1:npar) {
  BF.distr <- simulation.BF(N0, N, m, alpha0, beta0, alpha, beta, gamma.mean, sigma, ngenes=ngenes, pi=pi, alpha.model=alpha.model, beta.model=beta.model, gamma.mean.model=gamma.mean.model[i], sigma.model=sigma.model, pi.model=pi.model)
  BF.pos <- BF.distr$pos
  BF.neg <- BF.distr$neg
  AUC.BF[i] <- mean(sample(log(BF.pos), 1000, replace=TRUE) > sample(log(BF.neg),1000,replace=TRUE))
  TPR.BF[i] <- TPR(BF.pos, BF.neg, FPR)
}
# report <- data.frame(gamma.mean=gamma.mean.model, AUC=AUC.BF, TPR=TPR.BF)
# report <- data.frame(sigma=sigma.model, AUC=AUC.BF, TPR=TPR.BF)
# report <- data.frame(pi=pi.model, AUC=AUC.BF, TPR=TPR.BF)
report <- data.frame(gamma.mean=gamma.mean.model, sigma=sigma.model, pi=pi.model, AUC=AUC.BF, TPR=TPR.BF)

# plot the results: low-risk
report <- read.table("simulation/report.txt",header=TRUE)
postscript("simulation/sensitivity_risk_protective.eps", width=10, height=6)
par(mfrow=c(1,2))
plot(report$gamma.mean, report$AUC, type='b',ylim=c(0.70,0.9),xlab='Relative risk', ylab='AUC')
abline(h=0.741,col='red',lty=2)
text(1.9,0.75,"Burden")
abline(h=0.844,col='blue',lty=2)
text(1.9,0.85,"SKAT")
plot(report$gamma.mean, report$TPR, type='b',ylim=c(0.20,0.45),xlab='Relative risk', ylab='True positive rate')
abline(h=0.23,col='red',lty=2)
text(1.9,0.24,"Burden")
abline(h=0.298,col='blue',lty=2)
text(1.9,0.31,"SKAT")
par(mfrow=c(1,1))
dev.off()

# plot the results: mixture
report <- read.table("simulation/report.txt",header=TRUE)
report.pi <- read.table("simulation/report_pi.txt",header=TRUE)
postscript("simulation/sensitivity_risk_mixture.eps", width=10, height=14)
par(mfrow=c(2,2))
plot(report$gamma.mean, report$AUC, type='b',ylim=c(0.65,0.8),xlab='Relative risk', ylab='AUC')
abline(h=0.671,col='red',lty=2)
text(6,0.665,"Burden")
abline(h=0.674,col='blue',lty=2)
text(6, 0.685,"SKAT")
plot(report$gamma.mean, report$TPR, type='b',ylim=c(0.05,0.2),xlab='Relative risk', ylab='True positive rate')
abline(h=0.093,col='red',lty=2)
text(6,0.10,"Burden")
abline(h=0.085,col='blue',lty=2)
text(6,0.08,"SKAT")
plot(report.pi$pi, report$AUC, type='b',ylim=c(0.65,0.8),xlab=expression(pi), ylab='AUC')
abline(h=0.671,col='red',lty=2)
text(0.4,0.665,"Burden")
abline(h=0.674,col='blue',lty=2)
text(0.4, 0.685,"SKAT")
plot(report.pi$pi, report$TPR, type='b',ylim=c(0.05,0.2),xlab=expression(pi), ylab='True positive rate')
abline(h=0.093,col='red',lty=2)
text(0.4,0.1,"Burden")
abline(h=0.085,col='blue',lty=2)
text(0.4,0.08,"SKAT")
par(mfrow=c(1,1))
dev.off()

