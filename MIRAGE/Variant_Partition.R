fifteen.partition=function(cand.data) # given gene data and annotations, do variant partitions
{
  par.evid=list()
  LoF.def=c("stopgain", "frameshift substitution", "splicing", "stoploss")
  par.evid[[1]]=which(cand.data$Annotation %in% LoF.def==T & cand.data$ExacAF<0.05 & cand.data$ExacAF>=0.01 )
  par.evid[[2]]=which(cand.data$Annotation %in% LoF.def==T &  cand.data$ExacAF<0.01 & cand.data$ExacAF>=0.001)
  par.evid[[3]]=which(cand.data$Annotation %in% LoF.def==T &  cand.data$ExacAF<0.001 & cand.data$ExacAF>=0.0001)
  par.evid[[4]]=which(cand.data$Annotation %in% LoF.def==T &  cand.data$ExacAF<0.0001 & cand.data$ExacAF>0)
  par.evid[[5]]=which(cand.data$Annotation %in% LoF.def==T &  cand.data$ExacAF==0)
  
  
  par.evid[[6]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))>=0.957 & cand.data$ExacAF>=0.01 & cand.data$ExacAF<0.05)
  par.evid[[7]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))>=0.957 & cand.data$ExacAF>=0.001 & cand.data$ExacAF<0.01)
  par.evid[[8]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))>=0.957 & cand.data$ExacAF>=0.0001 & cand.data$ExacAF<0.001)
  par.evid[[9]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))>=0.957 & cand.data$ExacAF>0 & cand.data$ExacAF<0.0001)
  par.evid[[10]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))>=0.957 & cand.data$ExacAF==0 )
  
  
  par.evid[[11]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))<0.957 & cand.data$ExacAF>=0.01 & cand.data$ExacAF<0.05)
  par.evid[[12]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))<0.957 & cand.data$ExacAF>=0.001 & cand.data$ExacAF<0.01)
  par.evid[[13]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))<0.957 & cand.data$ExacAF>=0.0001 & cand.data$ExacAF<0.001)
  par.evid[[14]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))<0.957 & cand.data$ExacAF>0 & cand.data$ExacAF<0.0001)
  par.evid[[15]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))<0.957 & cand.data$ExacAF==0)
  
  group.index=rep(NA, nrow(cand.data))
  for (i in 1:length(par.evid))
    group.index[par.evid[[i]]]=i
  gene.data=data.frame(ID=cand.data$ID, Gene=cand.data$Gene, No.case=cand.data$No.case, No.contr=cand.data$No.contr, group.index=group.index)
  gene.data=gene.data[complete.cases(gene.data),]
  return(gene.data)
}

eight.partition=function(cand.data) # given gene data and annotations, do variant partitions
{
  par.evid=list()
  LoF.def=c("stopgain", "frameshift substitution", "splicing", "stoploss")
  par.evid[[1]]=which(cand.data$Annotation %in% LoF.def==T & cand.data$ExacAF<0.05 & cand.data$ExacAF>=0.01 )
  par.evid[[2]]=which(cand.data$Annotation %in% LoF.def==T &  cand.data$ExacAF<0.01)
  
  par.evid[[3]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))>=0.957 & cand.data$ExacAF>=0.01 & cand.data$ExacAF<0.05)
  par.evid[[4]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))>=0.957 & cand.data$ExacAF>=0.001 & cand.data$ExacAF<0.01)
  par.evid[[5]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))>=0.957 &  cand.data$ExacAF<0.001)
  
  par.evid[[6]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))<0.957 & cand.data$ExacAF>=0.01 & cand.data$ExacAF<0.05)
  par.evid[[7]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))<0.957 & cand.data$ExacAF>=0.001 & cand.data$ExacAF<0.01)
  par.evid[[8]]=which(cand.data$Annotation %in% LoF.def==F & as.numeric(as.character(cand.data$Polyphen2.HDIV.score))<0.957 & cand.data$ExacAF<0.001)
  
  group.index=rep(NA, nrow(cand.data))
  for (i in 1:length(par.evid))
    group.index[par.evid[[i]]]=i
  gene.data=data.frame(ID=cand.data$ID, Gene=cand.data$Gene, No.case=cand.data$No.case, No.contr=cand.data$No.contr, group.index=group.index)
  gene.data=gene.data[complete.cases(gene.data),]
  return(gene.data)
}