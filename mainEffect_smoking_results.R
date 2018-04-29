#Results Analysis
setwd('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq')
load('rdas/Genes_DE_sva.rda')
load('rdas/Exons_DE_sva.rda')
load('rdas/Jxns_DE_sva.rda')
load('rdas/ERs_DE_sva.rda')
load('rdas/interactionModel_allFeatures_DE_sva.rda')
#Primary Analysis

### Genes ###
 colnames(adultGene)[colnames(adultGene)%in%c('logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')] <-  paste0("Adult_",colnames(adultGene)[colnames(adultGene)%in%c('logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')] )
 colnames(fetalGene)[colnames(fetalGene)%in%c('logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')] <-  paste0("Prenatal_",colnames(fetalGene)[colnames(fetalGene)%in%c('logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')] )
 
 shared_info <- c(colnames(adultGene), colnames(fetalGene))[ duplicated(c(colnames(adultGene), colnames(fetalGene))) ]
 sigGene <- unique(c(rownames(adultGene[adultGene$Adult_adj.P.Val<0.10,]), rownames(fetalGene[fetalGene$Prenatal_adj.P.Val<0.10,]) ))
 bothGene<-cbind(adultGene[sigGene,],fetalGene[sigGene,!colnames(fetalGene)%in%shared_info])
 bothGene$Cohort_Sig <- ifelse(bothGene$Adult_adj.P.Val<0.10 & bothGene$Prenatal_adj.P.Val<0.10,"Both", ifelse(bothGene$Prenatal_adj.P.Val<0.10, "Prenatal", "Adult") )
 bothGene_Reduced <- bothGene[,c('Cohort_Sig','Symbol','Chr','Start','End','Prenatal_logFC', 'Adult_logFC','Prenatal_adj.P.Val','Adult_adj.P.Val','Prenatal_P.Value','Adult_P.Value')]

### Exons ###
 colnames(adultExon)[colnames(adultExon)%in%c('logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')] <-  paste0("Adult_",colnames(adultExon)[colnames(adultExon)%in%c('logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')] )
 colnames(fetalExon)[colnames(fetalExon)%in%c('logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')] <-  paste0("Prenatal_",colnames(fetalExon)[colnames(fetalExon)%in%c('logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')] )
 
 shared_info <- c(colnames(adultExon), colnames(fetalExon))[ duplicated(c(colnames(adultExon), colnames(fetalExon))) ]
 sigExon <- unique(c(rownames(adultExon[adultExon$Adult_adj.P.Val<0.10,]), rownames(fetalExon[fetalExon$Prenatal_adj.P.Val<0.10,]) ))
 #NOTHING SIG
 #bothExon<-cbind(adultExon[sigExon,],fetalExon[sigExon,!colnames(fetalExon)%in%shared_info]
 #bothExon$Cohort_Sig <- ifelse(bothExon$Adult_adj.P.Val<0.10 & bothExon$Prenatal_adj.P.Val<0.10,"Both", ifelse(bothExon$Prenatal_adj.P.Val<0.10, "Prenatal", "Adult") )
 #bothExon$Reduced <- bothExon[,c('Cohort_Sig','Symbol','Chr','Start','End','Prenatal_logFC', 'Adult_logFC','Prenatal_adj.P.Val','Adult_adj.P.Val')]

### Junctions ###
 colnames(adultJxn)[colnames(adultJxn)%in%c('logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')] <-  paste0("Adult_",colnames(adultJxn)[colnames(adultJxn)%in%c('logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')] )
 colnames(fetalJxn)[colnames(fetalJxn)%in%c('logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')] <-  paste0("Prenatal_",colnames(fetalJxn)[colnames(fetalJxn)%in%c('logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')] )
 
 shared_info <- c(colnames(adultJxn), colnames(fetalJxn))[ duplicated(c(colnames(adultJxn), colnames(fetalJxn))) ]
 sigJxn <- unique(c(rownames(adultJxn[adultJxn$Adult_adj.P.Val<0.10,]), rownames(fetalJxn[fetalJxn$Prenatal_adj.P.Val<0.10,]) ))
 bothJxn<-cbind(adultJxn[sigJxn,],fetalJxn[sigJxn,!colnames(fetalJxn)%in%shared_info])
 bothJxn$Cohort_Sig <- ifelse(bothJxn$Adult_adj.P.Val<0.10 & bothJxn$Prenatal_adj.P.Val<0.10,"Both", ifelse(bothJxn$Prenatal_adj.P.Val<0.10, "Prenatal", "Adult") )
 bothJxn_Reduced <- bothJxn[,c('Cohort_Sig','newGeneSymbol','seqnames','start','end','Prenatal_logFC', 'Adult_logFC','Prenatal_adj.P.Val','Adult_adj.P.Val','Prenatal_P.Value','Adult_P.Value')]

### ERs ###
 colnames(adultER)[colnames(adultER)%in%c('logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')] <-  paste0("Adult_",colnames(adultER)[colnames(adultER)%in%c('logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')] )
 colnames(fetalER)[colnames(fetalER)%in%c('logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')] <-  paste0("Prenatal_",colnames(fetalER)[colnames(fetalER)%in%c('logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')] )
 
 shared_info <- c(colnames(adultER), colnames(fetalER))[ duplicated(c(colnames(adultER), colnames(fetalER))) ]
 sigER <- unique(c(rownames(adultER[adultER$Adult_adj.P.Val<0.10,]), rownames(fetalER[fetalER$Prenatal_adj.P.Val<0.10,]) ))

 bothER<- cbind(adultER[sigER,,drop=FALSE],fetalER[sigER,!colnames(fetalER)%in%shared_info])
 bothER$Cohort_Sig <- ifelse(bothER$Adult_adj.P.Val<0.10 & bothER$Prenatal_adj.P.Val<0.10,"Both", ifelse(bothER$Prenatal_adj.P.Val<0.10, "Prenatal", "Adult") )
 bothER_Reduced <- bothER[,c('Cohort_Sig','nearestSymbol','seqnames','start','end','Prenatal_logFC', 'Adult_logFC','Prenatal_adj.P.Val','Adult_adj.P.Val', 'Prenatal_P.Value','Adult_P.Value')]
## Combining into a single table
bothER_Reduced<-plyr::rename(bothER_Reduced, c("nearestSymbol"="Symbol","seqnames" = "Chr", "start" = "Start","end"="End") )
bothJxn_Reduced<-plyr::rename(bothJxn_Reduced, c("newGeneSymbol"="Symbol","seqnames" = "Chr", "start" = "Start","end"="End") )

bothGene_Reduced$Feature <- "Gene" 
bothER_Reduced$Feature <-"Expressed Region"
bothJxn_Reduced$Feature <-"Junction"

allFeatures_Reduced <-rbind(bothGene_Reduced, bothER_Reduced, bothJxn_Reduced)

### Merging on interaction information
colnames(intGene) <- paste0("Interaction_", colnames(intGene))
colnames(intJxn) <- paste0("Interaction_", colnames(intJxn))
colnames(intER) <- paste0("Interaction_", colnames(intER))
intFeatures = rbind(intGene[,c('Interaction_P.Value', 'Interaction_adj.P.Val')], intJxn[,c('Interaction_P.Value', 'Interaction_adj.P.Val')], intER[,c('Interaction_P.Value', 'Interaction_adj.P.Val')])
allFeatures_Reduced = cbind(allFeatures_Reduced, intFeatures[rownames(allFeatures_Reduced),])
allFeatures_Reduced_ss = allFeatures_Reduced[, c("Cohort_Sig","Symbol", "Feature", "Prenatal_logFC", "Adult_logFC", "Prenatal_adj.P.Val" , "Adult_adj.P.Val", "Interaction_P.Value", "Interaction_adj.P.Val")]

### Saving results
write.csv(allFeatures_Reduced_ss, file='/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/Paper/Adult_Prenatal_Seperate_SigFeatures_Interaction_Table1.csv',row.names=FALSE)
median(abs(allFeatures_Reduced_ss$Prenatal_logFC[allFeatures_Reduced_ss$Cohort_Sig=="Prenatal"]))

## Supplemental Table for NRCAM junction results
NRCAM_junctions = bothJxn[bothJxn$Cohort_Sig=="Prenatal",]
NRCAM_junctions = cbind(NRCAM_junctions, intFeatures[rownames(NRCAM_junctions),])
NRCAM_junctions$ensemblTx = sapply(NRCAM_junctions$ensemblTx, paste, collapse=";")

write.csv(NRCAM_junctions, file='/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/Paper/NRCAM_junctions_supplemental_table.csv',row.names=FALSE)

#Ordinal Model
load('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/ordinal_sensitivity_model_results.rda')

#Schizophrenia replication
load('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/smoking_schizophrenia_model_results.rda')

adultJxn_SCZ[rownames(adultJxn[adultJxn$Adult_adj.P.Val<0.1,]),]
adultJxn_Ordinal[rownames(adultJxn[adultJxn$Adult_adj.P.Val<0.1,]),]

