## main EM algorithm 

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
        
        if (new.data$group.index[var.index.list[j]]<=5)
          var.BF[j]=BF.var.inte(new.data$No.case[var.index.list[j]], new.data$No.contr[var.index.list[j]], bar.gamma=6, sig=sigma, N1, N0)
        if (new.data$group.index[var.index.list[j]]>5)
          var.BF[j]=BF.var.inte(new.data$No.case[var.index.list[j]], new.data$No.contr[var.index.list[j]], bar.gamma=gamma.mean, sig=sigma, N1, N0)
        bb=bb*((1-beta.k[1, new.data$group.index[var.index.list[j]]])+beta.k[1, new.data$group.index[var.index.list[j]]]*var.BF[j])
        ################## split BF of LoF  and non LoF
        if (new.data$group.index[var.index.list[j]]<=2)
          bb.LoF=bb.LoF*((1-beta.k[1, new.data$group.index[var.index.list[j]]])+beta.k[1, new.data$group.index[var.index.list[j]]]*var.BF[j])
        if (new.data$group.index[var.index.list[j]]>2)
          bb.nonLoF=bb.nonLoF*((1-beta.k[1, new.data$group.index[var.index.list[j]]])+beta.k[1, new.data$group.index[var.index.list[j]]]*var.BF[j])
        
      }
    full.info.genevar[[i]]=cbind(indi.gene, var.BF)
    BF.gene[1, i]=bb
    LoF.BF.gene[1,i]=bb.LoF
    nonLoF.BF.gene[1,i]=bb.nonLoF
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
      bb=1; bb.LoF=1; bb.nonLoF=1
      UiZij=numeric()
      if (nrow(info.single.gene)>0)
        for (j in 1:nrow(info.single.gene))
        {
          
          numer=info.single.gene$var.BF[j]*beta.k[(iter-1), info.single.gene$group.index[j]]*delta.est[iter-1]
          denom=(delta.est[iter-1]+(1-delta.est[iter-1])/BF.gene[(iter-1),i])*(beta.k[(iter-1), info.single.gene$group.index[j]]*info.single.gene$var.BF[j]
                                                                               +(1-beta.k[(iter-1), info.single.gene$group.index[j]]))
          UiZij[j]=numer/denom
          bb=bb*((1-beta.k[(iter-1), info.single.gene$group.index[j]])+beta.k[(iter-1), info.single.gene$group.index[j]]*info.single.gene$var.BF[j])
          ###########################
          if (info.single.gene$group.index[j]<=2)
            bb.LoF=bb.LoF*((1-beta.k[(iter-1), info.single.gene$group.index[j]])+beta.k[(iter-1), info.single.gene$group.index[j]]*info.single.gene$var.BF[j])
          if (info.single.gene$group.index[j]>2)
            bb.nonLoF=bb.nonLoF*((1-beta.k[(iter-1), info.single.gene$group.index[j]])+beta.k[(iter-1), info.single.gene$group.index[j]]*info.single.gene$var.BF[j])
          
        }
      EUiZij[[i]]=UiZij
      BF.gene[iter,i]=bb
      LoF.BF.gene[iter,i]=bb.LoF
      nonLoF.BF.gene[iter,i]=bb.nonLoF
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
    delta.est[iter]=sum(EUi)/num.gene
    #  delta.est[iter]=delta
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
    cat(iter, "th iteration is running", "\n")
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
  return(result=list(delta.est=delta.est[iter-1], delta.stat=sum.lkhd, beta.est=beta.k[(iter-1),], beta.stat=cate.stat, cate.pvalue=cate.pvalue, BF.gene=data.frame(Gene=unique.gene, BF=BF.gene[(iter-1),], LoF.BF=LoF.BF.gene[(iter-1),], nonLoF.BF=nonLoF.BF.gene[(iter-1),]), full.info=full.info.genevar, Eui=EUi))
}