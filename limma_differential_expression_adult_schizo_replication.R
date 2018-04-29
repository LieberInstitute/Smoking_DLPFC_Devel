#qsub -l bluejay -l mf=40G,h_vmem=60G,h_stack=256M -cwd -b y R CMD BATCH --no-save limma_differential_expression_adult_schizo_replication.R

#Limma Smoking Differential Expression Analyses
setwd('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq')
load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rawCounts_szControlDLPFC.rda') #load this file for the mapObjects 
library(LIBDpheno)
library(ggplot2)
library(reshape2)
library(limma)
library(edgeR)
library(sva)
#####################
#Add in LIMS information to data on cluster which isn't otherwise available
LIMS = toxicant[[1]]
pd.model = pd[pd$Dx =="Control",]
id <- match( brnum(pd.model$BrNum), LIMS$brnumerical  )
id <- id[!is.na(id)]
pd.model <- cbind( LIMS[id,], pd.model)
pd.model$source <- droplevels(pd.model$source)
pd.model <- pd.model[,colSums(sapply(pd.model, is.na))!=nrow(pd.model)]
pd.model <- pd.model[ ,!colnames(pd.model) %in% 
						colnames(pd.model[,sapply(pd.model, is.logical)])[!sapply( pd.model[,sapply(pd.model, is.logical)], any, na.rm=TRUE)] ] #drop logicals without any trues
#seeing smoking and toxicant status

pd.model[pd.model$BrNum %in% c('Br1813','Br1826','Br2047', "Br2049"),
		 c("nicotine_comments","other_drugs_comments","final_dx_comments") ]   
pd.model$modelGroup <- pd.model$cotinine| pd.model$BrNum %in% c('Br1813','Br1826','Br2047', "Br2049")
pd.model$modelGroup <- as.factor(pd.model$modelGroup)
pd.model$modelGroup <- plyr::revalue(pd.model$modelGroup, c("TRUE" = "Smoker", "FALSE" = "Non-Smoker"))
pd.model$modelGroup[which(pd.model$nicotine)] <- "Smoker"
pd.model <- pd.model[!is.na(pd.model$modelGroup),]
pd.model <- pd.model[pd.model$agedeath<0|pd.model$agedeath>16,]
pd.model$AgeGroup = factor(ifelse(pd.model$agedeath<0, "Fetal","Adult"))
 pd.model<- pd.model[-which(pd.model$AgeGroup=="Adult"&pd.model$smoking & pd.model$modelGroup =="Non-Smoker"),]

pd.model<- pd.model[!pd.model$BrNum %in%c("Br1179","Br1105","Br2267"),]#outliers on SNP MDS plot droppped ADULTS
pd.model<- pd.model[!pd.model$BrNum %in%c("Br1779","Br1794"),]#outliers on SNP MDS plot droppped PRENATAL

### Splitting into three groups for fetal alone, adult alone, adult-fetal interaction
pd.adult = pd.model[pd.model$AgeGroup=="Adult",]
pd.fetal = pd.model[pd.model$AgeGroup=="Fetal",]

### Filtering to features tested before
gIndex <- which(rowMeans(cpm(geneCounts[,pd.model$RNum])>1)>0.1)
eIndex <- which(rowMeans(cpm(exonCounts[,pd.model$RNum])>1)>0.1)
jIndex <- which(rowMeans(cpm(as.data.frame(jCounts[,pd.model$RNum]))>1)>0.1 & jMap$code != "Novel" )

#####
pd.scz = pd[pd$Dx =="Schizo",]
pd.scz = pd.scz[pd.scz$Age >16,]
id <- match( brnum(pd.scz$BrNum), LIMS$brnumerical  )
id <- id[!is.na(id)]
pd.scz <- cbind( LIMS[id,], pd.scz)
pd.scz$source <- droplevels(pd.scz$source)
pd.scz <- pd.scz[,colSums(sapply(pd.scz, is.na))!=nrow(pd.scz)]
pd.scz <- pd.scz[ ,!colnames(pd.scz) %in% 
						colnames(pd.scz[,sapply(pd.scz, is.logical)])[!sapply( pd.scz[,sapply(pd.scz, is.logical)], any, na.rm=TRUE)] ] #drop logicals without any trues
#seeing smoking and toxicant status
table(Smoking = pd.scz$smoking, Toxicant = (pd.scz$nicotine | pd.scz$cotinine), useNA = 'ifany' )

pd.scz$modelGroup <- pd.scz$cotinine| pd.scz$nicotine 
pd.scz <- pd.scz[!is.na(pd.scz$modelGroup),]
pd.scz$modelGroup <- as.factor(pd.scz$modelGroup)
pd.scz$modelGroup <- plyr::revalue(pd.scz$modelGroup, c("TRUE" = "Smoker", "FALSE" = "Non-Smoker"))
pd.scz = pd.scz[-which(pd.scz$smoking & pd.scz$modelGroup =="Non-Smoker"),] #dropping ambiguous "non-smokers"
pd.scz <-pd.scz[!pd.scz$BrNum%in%"Br1178",] #drop outlier on snpPC################ qSV Creation


geneCounts.filt = geneCounts[gIndex,pd.scz$RNum] #keeping only the exons above the rpkm threshold (threshold 0.01 here)
geneMap.filt = geneMap[gIndex,] #keeping only the gene map for genes above the rpkm threshold

# Adult
dge.scz <- DGEList(counts = geneCounts.filt[,pd.scz$RNum],
    genes = geneMap.filt)
dge.scz = 	calcNormFactors(dge.scz)

adult_mod0 <- model.matrix(~ Age+ RIN + 
											mitoRate + totalAssignedGene + 
											Sex + snpPC1
											, data=pd.scz)

adult_mod <- model.matrix(~ modelGroup	+ Age+
											RIN + 
											mitoRate +  totalAssignedGene + 
											Sex + snpPC1 
											, data=pd.scz)

vGene = voom(dge.scz,adult_mod,plot=FALSE)
sv.scz = sva(vGene$E,adult_mod,adult_mod0)
mod_gene_scz <- cbind(adult_mod,sv.scz$sv)
vGene = voom(dge.scz,mod_gene_scz,plot=FALSE)

fitGene = lmFit(vGene, mod_gene_scz)
eBGene = eBayes(fitGene)

adultGene_SCZ = topTable(eBGene,sort.by='p',coef="modelGroupSmoker",
    n=Inf)

### Exons

exonCounts.filt = exonCounts[eIndex,pd.scz$RNum] #keeping only the exons above the rpkm threshold (threshold 0.01 here)
exonMap.filt = exonMap[eIndex,] #keeping only the exon map for exons above the rpkm threshold

# Adult
dge.scz <- DGEList(counts = exonCounts.filt[,pd.scz$RNum],
    genes = exonMap.filt)
dge.scz = 	calcNormFactors(dge.scz)

adult_mod0 <- model.matrix(~ Age+ RIN + 
											mitoRate + totalAssignedGene + 
											Sex + snpPC1
									   		, data=pd.scz)

adult_mod <- model.matrix(~ modelGroup	+ Age+
											RIN + 
											mitoRate + totalAssignedGene + 
											Sex + snpPC1 
											, data=pd.scz)

vExon = voom(dge.scz,adult_mod,plot=FALSE)
sv.scz = sva(vExon$E,adult_mod,adult_mod0)
mod_exon_scz <- cbind(adult_mod,sv.scz$sv)
vExon = voom(dge.scz,mod_exon_scz,plot=FALSE)

fitExon = lmFit(vExon, mod_exon_scz)
eBExon = eBayes(fitExon)

adultExon_SCZ = topTable(eBExon,sort.by='p',coef="modelGroupSmoker",
    n=Inf)

### Junctions


jCounts.filt = as.matrix(jCounts[jIndex,pd.scz$RNum]) #keeping only the exons above the rpkm threshold (threshold 0.01 here)
jMap.filt = as.data.frame(jMap[jIndex,]) #keeping only the exon map for exons above the rpkm threshold

# Adult
dge.scz <- DGEList(counts = jCounts.filt[,pd.scz$RNum],
    genes = jMap.filt)
dge.scz = 	calcNormFactors(dge.scz)

adult_mod0 <- model.matrix(~ Age+ RIN + 
											mitoRate + totalAssignedGene +
											Sex + snpPC1
											, data=pd.scz)

adult_mod <- model.matrix(~ modelGroup	+ Age+
											RIN + 
											mitoRate + totalAssignedGene +
											Sex + snpPC1 
											, data=pd.scz)

vJxn = voom(dge.scz,adult_mod,plot=FALSE)
sv.scz = sva(vJxn$E,adult_mod,adult_mod0)
mod_jxn_scz <- cbind(adult_mod,sv.scz$sv)
vJxn = voom(dge.scz,mod_jxn_scz,plot=FALSE)

fitJxn = lmFit(vJxn, mod_jxn_scz)
eBJxn = eBayes(fitJxn)
adultJxn_SCZ = topTable(eBJxn,sort.by='p',coef="modelGroupSmoker",
    n=Inf)

### ERs
load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/processed_covMat.rda')
regionMat.model <- round(regionMat[, pd.scz$RNum])

# Adult
dge.scz <- DGEList(counts = regionMat.model[,pd.scz$RNum],
    genes = as.data.frame(regions))
dge.scz = 	calcNormFactors(dge.scz)

adult_mod0 <- model.matrix(~ Age+ RIN + 
											mitoRate + totalAssignedGene +
											Sex + snpPC1
											, data=pd.scz)

adult_mod <- model.matrix(~ modelGroup	+ Age+
											RIN + 
											mitoRate + totalAssignedGene + 
											Sex + snpPC1 
											, data=pd.scz)

vER = voom(dge.scz,adult_mod,plot=FALSE)
sv.scz = sva(vER$E,adult_mod,adult_mod0)
mod_er_scz <- cbind(adult_mod,sv.scz$sv)
vER = voom(dge.scz,mod_er_scz,plot=FALSE)

fitER = lmFit(vER, mod_er_scz)
eBER = eBayes(fitER)

adultER_SCZ = topTable(eBER,sort.by='p',coef="modelGroupSmoker",
    n=Inf)
#save(adultGene_SCZ, adultExon_SCZ, adultJxn_SCZ, adultER_SCZ,file='rdas/smoking_schizophrenia_model_results.rda')
save(mod_gene_scz, mod_exon_scz, mod_jxn_scz, mod_er_scz, file='rdas/smoking_schizophrenia_mod_DE.rda')		#saving model objects