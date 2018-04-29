#
library(tidyr)
library(dplyr)
library(ggplot2)
theme_set(theme_bw(base_size=18) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="bottom"	)) 
### Checking enrichment with Geschwind autism DE
setwd('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq')
load('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/allFeatures_bound.rda')

### Test prenatal features 
TESTED_UNIV <- unique(all_bound$EnsemblGeneID)
PRENATAL_DE <- unique(all_bound[all_bound$adj.P.Val_Prenatal <0.1,'EnsemblGeneID'] )
PRENATAL_NOT_DE <- TESTED_UNIV[!TESTED_UNIV %in% PRENATAL_DE]
### Test interaction features 
INT_DE <- unique(all_bound[all_bound$adj.P.Val <0.1,'EnsemblGeneID'] )
INT_DE = INT_DE[!is.na(INT_DE)]
INT_NOT_DE <- TESTED_UNIV[!TESTED_UNIV %in% INT_DE]

### PUll in geschwind data
load('/dcl01/lieber/ajaffe/psychENCODE_Data/UCLA_R01MH094714/rdas/ASD_DE_gene.rda') #out
out = out[rownames(out)%in%TESTED_UNIV, ]

############## Gene set Test with limma ############
limma::geneSetTest(which(rownames(out)%in%INT_DE), out[,'Adj_ba41_42_22.tstat'] )

int_p = sapply(colnames(out)[grepl("tstat",colnames(out))], function(x) {
limma::geneSetTest(which(rownames(out)%in%INT_DE), abs(out[,x]) ) })

prenatal_p = sapply(colnames(out)[grepl("tstat",colnames(out))], function(x) {
limma::geneSetTest(which(rownames(out)%in%PRENATAL_DE), abs(out[,x]) ) })

out$interactionDE = ifelse(rownames(out)%in%INT_DE, "Yes", "No")
out$prenatalDE = ifelse(rownames(out)%in%PRENATAL_DE, "Yes", "No")
out$EnsemblGeneID = rownames(out)

dat = tidyr::gather(out[,c(grep("tstat",colnames(out),value=T), "interactionDE", "prenatalDE", "EnsemblGeneID" ) ], analysis, t_stat, -interactionDE, -prenatalDE, -EnsemblGeneID) %>% separate(analysis, c("Model", "Region"), sep="\\_",extra = "merge" )

ann_text_interaction <- data.frame(t_stat = paste0("p=",as.character(signif(int_p,3) ) ),
								   Model = factor( stringr::str_split_fixed(str = names(int_p), pattern = "\\_", n = 2)[,1] ),
								   Region = factor( stringr::str_split_fixed(str = names(int_p), pattern = "\\_", n = 2)[,2] ),
								   myIntDensity="Yes"
								   ) 
								   
ann_text_prenatal <- data.frame(t_stat = paste0("p=",as.character(signif(prenatal_p,3) ) ),
								   Model = factor( stringr::str_split_fixed(str = names(prenatal_p), pattern = "\\_", n = 2)[,1] ),
								   Region = factor( stringr::str_split_fixed(str = names(prenatal_p), pattern = "\\_", n = 2)[,2] ),
								   myIntDensity="Yes"
								   ) 
#
mean_t <- dplyr::group_by(dat,Model, Region, interactionDE) %>%
                 summarize(
                 mean= mean(abs(t_stat)) )

								   
myIntDensity = ggplot(subset(dat,Model=="Qual"),aes(x= abs(t_stat), fill=interactionDE) )+ 
			geom_density(alpha=0.6, position="identity") + 
			facet_wrap(~Region, labeller=labeller(Region=c(ba41_42_22.tstat="BA41, BA42, BA22", ba9.tstat="BA9", vermis.tstat="Vermis") )) + 
			labs(x="Absolute T-Statistic",
				 y="Density") + 
			scale_fill_manual(values=c("black", "red")) +
			scale_colour_manual(values=c("black", "red")) +
			geom_text(data=subset(ann_text_interaction,Model=="Qual"), aes(x=1.15, y=0.74, label=t_stat), colour="black",inherit.aes=FALSE) +
			theme(legend.position='none' ) 
			
myPrenatalDensity = ggplot(data=subset(dat,Model=="Qual"),aes(x=abs(t_stat), fill=prenatalDE) )+ geom_density(alpha=0.6, position="identity") + 
			facet_wrap(~Region, labeller=labeller(Region=c(ba41_42_22.tstat="BA41, BA42, BA22", ba9.tstat="BA9", vermis.tstat="Vermis") )) + 
			labs(x="Absolute T-Statistic",
				 y="Density") + 
			scale_fill_manual(values=c("black", "red")) +
			geom_text(data=subset(ann_text_prenatal,Model=="Qual"), aes(x=0.95, y=0.71, label=t_stat), colour="black",inherit.aes=FALSE) +
			theme(legend.position='none' )
			
#ggsave(myIntDensity,filename='/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/plots/qualModel_interaction_GeschwindAutism_t_statistic_density.pdf', height=11, width=8.5)
#ggsave(myPrenatalDensity,filename='/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/plots/qualModel_prenatal_GeschwindAutism_t_statistic_density.pdf', height=11,width=8.5)

ggsave(myIntDensity,filename='/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/plots/qualModel_interaction_GeschwindAutism_t_statistic_density_NGCposter.pdf', height=8.5, width=8.5)
ggsave(myPrenatalDensity,filename='/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/plots/qualModel_prenatal_GeschwindAutism_t_statistic_density_NGCposter.pdf', height=8.5,width=8.5)

### Random permutation
perm = sapply(1:5000,  function(i) {
fake_index = sample(nrow(out),size=sum(rownames(out)%in%INT_DE) )
dat = sapply(colnames(out)[grepl("tstat",colnames(out))], function(x) {
limma::geneSetTest(fake_index, out[,x] ) })
return(dat)
 })
perm=t(perm)

 pdf('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/plots/interaction_GeschwindAutism_DE_Enrichment_5000Permutation_Null_Histograms_Limma_geneSetTest.pdf')
 par(mfrow=c(3,2))
 hist(perm[,'Adj_ba41_42_22.tstat'], breaks=20)
 hist(perm[,'Qual_ba41_42_22.tstat'], breaks=20)
 hist(perm[,'Adj_ba9.tstat'], breaks=20)
 hist(perm[,'Qual_ba9.tstat'], breaks=20)
 hist(perm[,'Adj_vermis.tstat'], breaks=20)
 hist(perm[,'Qual_vermis.tstat'], breaks=20)
 dev.off()
 
####################################################
setList = list(
ba41_42_22 = rownames(out)[out$Adj_ba41_42_22.fdr <0.05],
Qual_ba41_42_22 = rownames(out)[out$Qual_ba41_42_22.fdr <0.05],
ba9 = rownames(out)[out$Adj_ba9.fdr <0.05],
Qual_ba9 = rownames(out)[out$Qual_ba9.fdr <0.05],
vermis = rownames(out)[out$Adj_vermis.fdr <0.05],
Qual_vermis = rownames(out)[out$Qual_vermis.fdr <0.05]
)


prenatalEnrich = sapply(setList, function(x) {
DE_OVERLAP = c( sum(PRENATAL_DE%in% x),sum(!PRENATAL_DE%in% x))
NOT_DE_OVERLAP= c(sum(PRENATAL_NOT_DE%in% x), sum(!PRENATAL_NOT_DE%in% x) )
enrich_table = cbind(DE_OVERLAP, NOT_DE_OVERLAP)
res = fisher.test(enrich_table)
dat=c(res$p.value, 
	  res$estimate, 
	  nGeneEnriched =enrich_table[1] , 
	  nGeneSet = enrich_table[1] +enrich_table[3],  
	  percentEnrich = (enrich_table[1] / enrich_table[2])*100, 
	  backgroundEnrich = (enrich_table[3] / enrich_table[4] ) *100,
	  geneSymbolsEnriched = paste( PRENATAL_DE[PRENATAL_DE%in% x] , collapse=";") )

names(dat) <- c("P.Value","Odds Ratio", "nGeneEnriched", "nGeneSet", "percentEnrich", "percentBackground","geneSymbolsEnriched")
return(dat)
})
write.csv(t(prenatalEnrich), file='/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/csvs/prenatal_GeschwindAutism_DE_Enrichment.csv' ,row.names=TRUE)

intEnrich = sapply(setList, function(x) {
DE_OVERLAP = c( sum(INT_DE%in% x),sum(!INT_DE%in% x))
NOT_DE_OVERLAP= c(sum(INT_NOT_DE%in% x), sum(!INT_NOT_DE%in% x) )
enrich_table = cbind(DE_OVERLAP, NOT_DE_OVERLAP)
res = fisher.test(enrich_table)
dat=c(res$p.value, 
	  res$estimate, 
	  nGeneEnriched =enrich_table[1] , 
	  nGeneSet = enrich_table[1] +enrich_table[3],  
	  percentEnrich = (enrich_table[1] / enrich_table[2])*100, 
	  backgroundEnrich = (enrich_table[3] / enrich_table[4] ) *100,
	  geneSymbolsEnriched = paste( INT_DE[INT_DE%in% x] , collapse=";") )

names(dat) <- c("P.Value","Odds Ratio", "nGeneEnriched", "nGeneSet", "percentEnrich", "percentBackground","geneSymbolsEnriched")
return(dat)
})
write.csv(t(intEnrich), file='/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/csvs/interaction_GeschwindAutism_DE_Enrichment.csv' ,row.names=TRUE)
intEnrich = read.csv('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/csvs/interaction_GeschwindAutism_DE_Enrichment.csv')

### Int enrich permutation
perm = sapply(1:10000,  function(i) {
INT_FAKE_DE = sample(TESTED_UNIV,size=length(INT_DE))
INT_FAKE_NOT_DE = TESTED_UNIV[!TESTED_UNIV %in% INT_FAKE_DE]
intEnrich_Null = sapply(setList, function(x) {
DE_OVERLAP = c( sum(INT_FAKE_DE%in% x),sum(!INT_FAKE_DE%in% x))
NOT_DE_OVERLAP= c(sum(INT_FAKE_NOT_DE%in% x), sum(!INT_FAKE_NOT_DE%in% x) )
enrich_table = cbind(DE_OVERLAP, NOT_DE_OVERLAP)
res = fisher.test(enrich_table)
dat=c(res$p.value, 
	  res$estimate, 
	  nGeneEnriched =enrich_table[1] , 
	  nGeneSet = enrich_table[1] +enrich_table[3],  
	  percentEnrich = (enrich_table[1] / enrich_table[2])*100, 
	  backgroundEnrich = (enrich_table[3] / enrich_table[4] ) *100 )

names(dat) <- c("P.Value","Odds Ratio", "nGeneEnriched", "nGeneSet", "percentEnrich", "percentBackground")
return(dat)
})
intEnrich_Null['P.Value',]
 })
 perm = t(perm)
 pdf('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/plots/interaction_GeschwindAutism_DE_EnrichmentPermutation_Null_Histograms.pdf')
 par(mfrow=c(2,2))
 hist(perm[,'ba41_42_22'], breaks=20)
 hist(perm[,'Qual_ba41_42_22'], breaks=20)
 hist(perm[,'ba9'], breaks=20)
 hist(perm[,'Qual_ba9'], breaks=20)
 dev.off()
 
 shared_genes = intersect(rownames(all_bound),rownames(out) )
 par(mfrow=c(3,2)) 
 plot ( all_bound[shared_genes,"t_Prenatal"], out[shared_genes, "Adj_ba9.tstat"] )
 plot ( all_bound[shared_genes,"t_Prenatal"], out[shared_genes, "Adj_ba41_42_22.tstat"] )
 plot ( all_bound[shared_genes,"t_Prenatal"], out[shared_genes, "Adj_vermis.tstat"] )
 plot ( all_bound[shared_genes,"t_Prenatal"], out[shared_genes, "Qual_ba9.tstat"] )
 plot ( all_bound[shared_genes,"t_Prenatal"], out[shared_genes, "Qual_ba41_42_22.tstat"] )
 plot ( all_bound[shared_genes,"t_Prenatal"], out[shared_genes, "Qual_vermis.tstat"] )
 
 par(mfrow=c(3,2)) 
 plot ( -log10(all_bound[shared_genes,"P.Value"]), -log10(out[shared_genes, "Adj_ba9.pval"])  )
 plot ( -log10(all_bound[shared_genes,"P.Value"]), -log10(out[shared_genes, "Adj_ba41_42_22.pval"]) )
 plot ( -log10(all_bound[shared_genes,"P.Value"]), -log10(out[shared_genes, "Adj_vermis.pval"]) )
 plot ( -log10(all_bound[shared_genes,"P.Value"]), -log10(out[shared_genes, "Qual_ba9.pval"]) )
 plot ( -log10(all_bound[shared_genes,"P.Value"]), -log10(out[shared_genes, "Qual_ba41_42_22.pval"]) )
 plot ( -log10(all_bound[shared_genes,"P.Value"]), -log10(out[shared_genes, "Qual_vermis.pval"]) )
 
  plot ( (out[shared_genes, "Adj_vermis.tstat"]), (out[shared_genes, "Qual_vermis.tstat"]) )