#
setwd('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq')
load('rdas/interactionModel_allFeatures_DE_sva.rda')
load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda') #rpkm data
intJxn$EntrezID <- geneMap[intJxn[,'newGeneID'], c("EntrezID")]
load('/users/ssemick/Tox_Expr/Adult_Smoking/RNA_Seq/rdas/Adult/ER Annotation.rda')
ann$RowName <- rownames(ER)
intER$EntrezID <- as.character(ann[match(rownames(intER),ann$RowName) ,'Entrez'])
intER$EnsemblGeneID <- rownames(geneMap)[match(intER$EntrezID, geneMap$EntrezID)]

intGene$Type="Gene"
intExon$Type="Exon"
intJxn$Type="Jxn"
intER$Type="ER"

JxnSave <- plyr::rename(intJxn, c("seqnames" = "Chr", 
											  "start" = "Start", 
											  "end" = "End", 
											  "strand" ="Strand", 
											  "width" =  "Length",
											  "newGeneSymbol" = "Symbol",
											  "newGeneID" = "EnsemblGeneID"
											  ) )
ERSave <- plyr::rename(intER, c("seqnames" = "Chr", 
											  "start" = "Start", 
											  "end" = "End", 
											  "strand" ="Strand", 
											  "width" =  "Length",
											  "nearestSymbol" = "Symbol"
											  ) )
colnames(intExon)[1] = "EnsemblGeneID"
intGene$EnsemblGeneID = rownames(intGene)
											  
common_columns <- Reduce(intersect, lapply(list(intGene,intExon,JxnSave,ERSave),colnames))
All <- rbind(intGene[,common_columns], intExon[,common_columns], JxnSave[,common_columns], ERSave[,common_columns]  )

All[All$adj.P.Val<0.10,]

load('rdas/Genes_DE_sva.rda')
adultGene2 = adultGene
colnames(adultGene2) = paste0(colnames(adultGene2), "_Adult") 
fetalGene2 = fetalGene
colnames(fetalGene2) = paste0(colnames(fetalGene2), "_Prenatal") 
geneToBind = cbind(adultGene2[rownames(intGene),8:13], fetalGene2[rownames(intGene),8:13])

load('rdas/Exons_DE_sva.rda')
adultExon2 = adultExon
colnames(adultExon2) = paste0(colnames(adultExon2), "_Adult") 
fetalExon2 = fetalExon
colnames(fetalExon2) = paste0(colnames(fetalExon2), "_Prenatal") 
exonToBind = cbind(adultExon2[rownames(intExon),9:14], fetalExon2[rownames(intExon),9:14])

load('rdas/Jxns_DE_sva.rda')
adultJxn2 = adultJxn
colnames(adultJxn2) = paste0(colnames(adultJxn2), "_Adult") 
fetalJxn2 = fetalJxn
colnames(fetalJxn2) = paste0(colnames(fetalJxn2), "_Prenatal") 
jxnToBind = cbind(adultJxn2[rownames(intJxn),22:27], fetalJxn2[rownames(intJxn),22:27])

load('rdas/ERs_DE_sva.rda')
adultER2 = adultER
colnames(adultER2) = paste0(colnames(adultER2), "_Adult") 
fetalER2 = fetalER
colnames(fetalER2) = paste0(colnames(fetalER2), "_Prenatal") 
erToBind = cbind(adultER2[rownames(intER),18:23], fetalER2[rownames(intER),18:23])

featureToBind = rbind(geneToBind, exonToBind, jxnToBind, erToBind)
all_bound = cbind(All, featureToBind[rownames(All),])

all_bound$Feature_ID = rownames(all_bound)
save(all_bound, file='/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/allFeatures_bound.rda' )

all_bound_sig = all_bound[all_bound$adj.P.Val<0.10, ]
all_bound_sig = all_bound_sig[,c('Symbol','EntrezID','EnsemblGeneID','Type','P.Value','adj.P.Val','t','logFC','P.Value_Prenatal','adj.P.Val_Prenatal','t_Prenatal','logFC_Prenatal','P.Value_Adult','adj.P.Val_Adult','t_Adult','logFC_Adult')]
write.csv(all_bound_sig, '/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/Paper/FDR10_significant_interaction_results.csv',row.names=FALSE)
##########
all_bound_sig[all_bound_sig$P.Value_Prenatal<0.05,]
all_bound_sig[all_bound_sig$P.Value_Adult<0.05,]
table(Prenatal= all_bound_sig$P.Value_Prenatal<0.05, Adult= all_bound_sig$P.Value_Adult<0.05 )
fisher.test(table(Adult= all_bound_sig$P.Value_Adult<0.05,Prenatal= all_bound_sig$P.Value_Prenatal<0.05 ))
fisher.test(cbind(table(all_bound_sig$P.Value_Adult<0.05),table(all_bound_sig$P.Value_Prenatal<0.05) ))

chisq.test(table(sign(all_bound_sig$logFC) ) )

library(ggplot2)
theme_set(theme_bw(base_size=14) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"))
dat = all_bound_sig
dat$LFC_Sign = 		sign(dat$logFC)		 
ggplot(data=dat, aes(x =t_Prenatal, y = t_Adult ) ) + geom_point() +facet_grid(Type~LFC_Sign ) +geom_vline(xintercept=0) + geom_hline(yintercept=0)
ggplot(data=dat, aes(x =(t_Prenatal+t_Adult)/2, y = t_Prenatal- t_Adult) ) + geom_point() +facet_grid(Type~LFC_Sign ) +geom_vline(xintercept=0) + geom_hline(yintercept=0)
				 


