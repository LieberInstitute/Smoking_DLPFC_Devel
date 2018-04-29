#Checking that gene length is not a source of bias for enrichment analylsis
load('rdas/interactionModel_allFeatures_DE_sva.rda')
intER[!duplicated(intER$Geneid),]

load('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/allFeatures_bound.rda')

sig_length=all_bound[all_bound$adj.P.Val<0.10,]
sig_length = sig_length[!duplicated(sig_length$EntrezID),]
sig_length = sig_length[!is.na(sig_length$EntrezID),]

notsig_length=all_bound[all_bound$adj.P.Val>0.10,]
notsig_length = notsig_length[!duplicated(notsig_length$EntrezID),]
notsig_length = notsig_length[!is.na(notsig_length$EntrezID),]

t.test(sig_length$Length,notsig_length$Length)