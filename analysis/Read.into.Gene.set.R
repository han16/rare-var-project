# read into gene set 
GeneDB=src_sqlite(path="C:\\Shengtong\\han_desktop\\Research\\rare-var\\gene.list.db", create=F)
#GeneDB=src_sqlite(path="../../gene.list.db", create=F)
gene_cate1=data.frame(collect(tbl(GeneDB, "SFARI_HighConf")))
gene_cate2=data.frame(collect(tbl(GeneDB, "SFARI_StrongCand")))
gene_cate3=data.frame(collect(tbl(GeneDB, "SFARI_cate3_gene")))
gene_cate4=data.frame(collect(tbl(GeneDB, "SFARI_cate4_gene")))
gene_cate5=data.frame(collect(tbl(GeneDB, "SFARI_cate5_gene")))
gene_cate6=data.frame(collect(tbl(GeneDB, "SFARI_cate6_gene")))
gene_cateS=data.frame(collect(tbl(GeneDB, "SFARI_cateS_gene")))
TADAGene=data.frame(collect(tbl(GeneDB, "TADAGenelist")))
Qlessthan5percentgene=TADAGene$TadaName[TADAGene$qvalue.combined<0.05]
Qlessthan20percentgene=TADAGene$TadaName[TADAGene$qvalue.combined<0.2]
Qlessthan30percentgene=TADAGene$TadaName[TADAGene$qvalue.combined<0.3]
Qlessthan40percentgene=TADAGene$TadaName[TADAGene$qvalue.combined<0.4]
Qlessthan50percentgene=TADAGene$TadaName[TADAGene$qvalue.combined<0.5]
Qlargerthan90percentgene=TADAGene$TadaName[TADAGene$qvalue.combined>0.9]
###### ID gene ######
ID_Gene=as.character(data.frame(collect(tbl(GeneDB, "Pinto14AJHG_IDgene")))$GeneID)
#### high confidence gene 
High_Confidence_Gene=union(union(gene_cate1$GeneID, gene_cate2$GeneID),Qlessthan5percentgene)
#### moderate confidence gene 
Moderate_Confidence_Gene=setdiff(union(union(gene_cate3$GeneID, gene_cateS$GeneID), Qlessthan30percentgene),Qlessthan5percentgene)
#### psd gene 
purcell.genelist=data.frame(collect(tbl(GeneDB, "Purcell2014_genelist"))) ## PSD gene, SCZdenovo gene 
PSD_Gene=purcell.genelist$Gene_symbol[purcell.genelist$PSD=="Y"]
#### FMRP gene 
FMRP_Gene=purcell.genelist$Gene_symbol[purcell.genelist$FMRP.target=="Y"]
#### Autism KB gene 
ASD_KB_Gene=as.character(data.frame(collect(tbl(GeneDB, "AutismKB_gene")))$GeneID)
#### constraint gene 
Constraint_Gene=data.frame(collect(tbl(GeneDB, "Samocha_2014NG_constraintgene")))$gene
### RVIS gene 
RVIS.Allgene=data.frame(collect(tbl(GeneDB, "RVIS_gene")))
RVIS_Gene=RVIS.Allgene$GeneID[RVIS.Allgene$RVIS.percentile<5] # top 5% gene
#### Haplo insufficient gene 
Haploinsufficient_Gene=as.character(data.frame(collect(tbl(GeneDB, "Petrovski_plosgen_haploinsuff_gene")))$GeneID)
#### SCZ data 
SCZ_Gene=purcell.genelist$Gene_symbol
Gene_Set=c("ID gene","High", "Mod", "PSD", "FMRP", "AutismKB", "constraint gene", "RVIS", "Haploinsuff", "SCZ gene")
Gene_Set_List=list()
Gene_Set_List[[1]]=ID_Gene
Gene_Set_List[[2]]=High_Confidence_Gene
Gene_Set_List[[3]]=Moderate_Confidence_Gene
Gene_Set_List[[4]]=PSD_Gene
Gene_Set_List[[5]]=FMRP_Gene
Gene_Set_List[[6]]=ASD_KB_Gene
Gene_Set_List[[7]]=Constraint_Gene
Gene_Set_List[[8]]=RVIS_Gene
Gene_Set_List[[9]]=Haploinsufficient_Gene
Gene_Set_List[[10]]=SCZ_Gene