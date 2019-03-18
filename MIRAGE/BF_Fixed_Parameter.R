# calculate Bayes factor with given parameter 
fixed.beta.func=function(new.data, N1, N0, gamma.mean, sigma, beta) # new.data has one column specifying its group index
{
  #######################
  beta.k=beta
  full.info.genevar=list()
  gene.list=new.data$Gene; unique.gene=unique(gene.list) # find the gene list
  num.gene=length(unique.gene)
  BF.gene=numeric()
  ########################
  for (i in 1:num.gene)
  {
    cat(i, "th gene of ", "\t", num.gene, "\t", "is running", "\n")
    var.index.list=which(gene.list==unique.gene[i])
    indi.gene=new.data[var.index.list,]   # note var.index.list matches new.data
    bb=1; var.BF=numeric()
    if (length(var.index.list)>0) # calculate Bayes factor for variant (i,j)
      for (j in 1:length(var.index.list))
      {
        if (new.data$group.index[var.index.list[j]]<=5)
          var.BF[j]=BF.var.inte(new.data$No.case[var.index.list[j]], new.data$No.contr[var.index.list[j]], bar.gamma=6, sig=sigma, N1, N0)
        if (new.data$group.index[var.index.list[j]]>5)
          var.BF[j]=BF.var.inte(new.data$No.case[var.index.list[j]], new.data$No.contr[var.index.list[j]], bar.gamma=gamma.mean, sig=sigma, N1, N0)
        bb=bb*((1-beta.k[new.data$group.index[var.index.list[j]]])+beta.k[new.data$group.index[var.index.list[j]]]*var.BF[j])
      }
    full.info.genevar[[i]]=cbind(indi.gene, var.BF)
    BF.gene[i]=bb
  }
  
  return(result=list(BayesFactor=data.frame(Gene=unique.gene, BF=BF.gene), full.info=full.info.genevar))
}