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
load('/dcl01/lieber/ajaffe/psychENCODE_Data/UCLA_R01MH094714/rdas/ASD_DE_gene.rda') #out

asd_col_interest = c(grep("pval", colnames(out),value=TRUE ), grep("tstat", colnames(out),value=TRUE ))
asd_col_interest=asd_col_interest[grepl("Qual",asd_col_interest)]
nic_col_interest=c("Symbol",grep("P\\.Value", colnames(all_bound),value=TRUE ),"t","t_Adult","t_Prenatal" )

dat = cbind(all_bound, out[match(all_bound$EnsemblGeneID,rownames(out) ),asd_col_interest] )
head(dat[dat[,'adj.P.Val_Prenatal'] < 0.10,c(nic_col_interest,asd_col_interest )]) #demonstration

nic_asd_stats = dat[,c(nic_col_interest,asd_col_interest )]
#Saving
save(nic_asd_stats,file= '/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/nicotine_asdGeschwind_stats.rda')
write.csv(nic_asd_stats,file= '/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/tables/nicotine_asdGeschwind_stats.csv',row.names=TRUE)
write.csv(nic_asd_stats,file= gzfile('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/tables/nicotine_asdGeschwind_stats.csv.gz'),row.names=TRUE)

#
dat$significant_de  = ifelse(dat$adj.P.Val_Prenatal<0.1 & dat$adj.P.Val<0.1, "Prenatal and Interaction",
						ifelse(dat$adj.P.Val_Prenatal<0.1, "Prenatal", 
							ifelse(dat$adj.P.Val<0.1, "Interaction", "Not DE") ) )
dat$prenatal_sig =ifelse(dat$adj.P.Val_Prenatal<0.1, TRUE, FALSE)
dat$interaction_sig =ifelse(dat$adj.P.Val<0.1, TRUE, FALSE)
							
							
nic_asd_prenatal_sig=dat[dat$adj.P.Val_Prenatal<0.1|dat$adj.P.Val<0.1,c(nic_col_interest,asd_col_interest, "prenatal_sig","interaction_sig" )]

save(nic_asd_prenatal_sig,file= '/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/nicotine_asdGeschwind_stats_sigPrenatal_or_Interaction.rda')
write.csv(nic_asd_prenatal_sig,file= '/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/tables/nicotine_asdGeschwind_stats_sigPrenatal_or_Interaction.csv',row.names=TRUE)
write.csv(nic_asd_prenatal_sig,file= gzfile('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/tables/nicotine_asdGeschwind_stats_sigPrenatal_or_Interaction.csv.gz'),row.names=TRUE)

##
nic_asd_prenatal_sig = read.csv('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/tables/nicotine_asdGeschwind_stats_sigPrenatal_or_Interaction.csv')
load('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/allFeatures_bound.rda')
write.csv(cbind(nic_asd_prenatal_sig, all_bound[as.character(nic_asd_prenatal_sig$X),c('EntrezID','EnsemblGeneID','Type')] ), file= '/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/tables/nicotine_asdGeschwind_stats_sigPrenatal_or_Interaction.csv',row.names=TRUE)

