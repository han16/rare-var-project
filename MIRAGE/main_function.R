main.function=function(evid, gene.set)  # input must be either evid or gene set, but not both
{
  if (length(evid)>1)
    cand.data=Anno.Data[which(Anno.Data$ID %in% evid),]

  if (length(gene.set)>1)
  {
    vart.set=as.character(Anno.Data$ID[which(Anno.Data$Gene %in% gene.set)])
    cand.data=Anno.Data[which(Anno.Data$ID %in% vart.set),]
   }
gene.data=eight.partition(cand.data)

{
  overlap.id=gene.data$ID
  if (length(overlap.id)>0)
  {
    overlap.data=gene.data[gene.data$ID %in% overlap.id,]
    order.overlap.data=overlap.data[order(overlap.data$group.index, decreasing=F),]
    psbl.index=unique(order.overlap.data$group.index); actu.num.group=length(psbl.index)
    delta.init=runif(1); beta.init=runif(actu.num.group)
    for (j in 1:actu.num.group)
      order.overlap.data$group.index[order.overlap.data$group.index==psbl.index[j]]=j # re-index the group labels
    #  delta=runif(1)
    para.est=multi.group.func(order.overlap.data, N1, N0, gamma.mean=3, sigma=2, delta=0.2, beta.init, actu.num.group)

  }
}
return(para.est)
}
