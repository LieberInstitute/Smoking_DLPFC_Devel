#qsub -l bluejay -l mf=50G,h_vmem=70G,h_stack=256M -cwd -b y R CMD BATCH limma_differential_expression_analyses_sva_cluster.R

#Limma Smoking Differential Expression Analyses
setwd('/users/ssemick/Tox_Expr/Adult_Smoking/RNA_Seq')
load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rawCounts_szControlDLPFC.rda') #load this file for the mapObjects 
#load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda') #load rpkm data to filter junctions before DESeq2 pipeline--speed up computation!
library(LIBDpheno)
library(ggplot2)
library(reshape2)
library(limma)
library(edgeR)

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

### Adding qSVs
library(sva)
degradationMat <- read.csv(gzfile('/users/ajaffe/Lieber/Projects/RNAseq/Consortium/PhaseI/degradationMat_DLPFC_polyA_Phase1.csv.gz')) #degredation matrix for the rna-seq 
degradationMat.model = degradationMat[ , pd.model[, 'RNum'] ]
qsv.model = qsva(degradationMat.model)
#pd.model <- cbind(pd.model, qsv.model)


degradationMat.adult = degradationMat[ , pd.adult[, 'RNum'] ]
qsv.adult = qsva(degradationMat.adult)
#pd.adult <- cbind(pd.adult, qsv.adult)

degradationMat.fetal = degradationMat[ , pd.fetal[, 'RNum'] ]
qsv.fetal = qsva(degradationMat.fetal)
#pd.fetal <- cbind(pd.fetal, qsv.fetal)

### Genes

#gIndex = which(rowMeans(geneRpkm[,pd.model$RNum]) > .1) # Filtering genes by expression levels 	
gIndex <- which(rowMeans(cpm(geneCounts[,pd.model$RNum])>1)>0.1)
#gIndex <- which(rowSums(cpm(geneCounts[,pd.model$RNum])>1)>=30)

#geneRpkm.filt = geneRpkm[gIndex,pd.model$RNum] #keeping only the exons above the rpkm threshold (threshold 0.01 here)
geneCounts.filt = geneCounts[gIndex,pd.model$RNum] #keeping only the exons above the rpkm threshold (threshold 0.01 here)
geneMap.filt = geneMap[gIndex,] #keeping only the gene map for genes above the rpkm threshold

# Adult
dge.adult <- DGEList(counts = geneCounts.filt[,pd.adult$RNum],
    genes = geneMap.filt)
dge.adult = 	calcNormFactors(dge.adult)

adult_mod0 <- model.matrix(~ Age+ RIN + 
											mitoRate + totalAssignedGene + 
											Sex + snpPC1
											, data=pd.adult)

adult_mod <- model.matrix(~ modelGroup	+ Age+
											RIN + 
											mitoRate +  totalAssignedGene + 
											Sex + snpPC1 
											, data=pd.adult)

vGene = voom(dge.adult,adult_mod,plot=FALSE)
sv.adult = sva(vGene$E,adult_mod,adult_mod0)
mod_gene_adult <- cbind(adult_mod,sv.adult$sv)
vGene = voom(dge.adult,mod_gene_adult,plot=FALSE)

fitGene = lmFit(vGene, mod_gene_adult)
eBGene = eBayes(fitGene)

adultGene = topTable(eBGene,sort.by='p',coef="modelGroupSmoker",
    n=Inf)

# Fetal

#gIndex.fetal = which(rowMeans(geneRpkm[,pd.fetal$RNum]) > .01) # Filtering genes by expression levels 	
#gIndex.fetal <- which(rowSums(cpm(geneCounts[,pd.fetal$RNum])>1)>=3 )
#gIndex.fetal <- which(rowMeans(cpm(geneCounts[,pd.fetal$RNum])>1)>0)
#gIndex.fetal <- which(rowMeans(cpm(geneCounts[,pd.model$RNum]))>.01)

#geneRpkm.fetal = geneRpkm[gIndex.fetal,pd.fetal$RNum] #keeping only the exons above the rpkm threshold (threshold 0.01 here)
#geneCounts.fetal = geneCounts[gIndex.fetal,pd.fetal$RNum] #keeping only the exons above the rpkm threshold (threshold 0.01 here)
#geneMap.fetal = geneMap[gIndex.fetal,] #keeping only the gene map for genes above the rpkm threshold

dge.fetal <- DGEList(counts = geneCounts.filt[,pd.fetal$RNum],
    genes = geneMap.filt)
dge.fetal = calcNormFactors(dge.fetal)

fetal_mod0 <- model.matrix(~ Age+
											RIN + 
											mitoRate + totalAssignedGene +
											Sex + snpPC1
											, data=pd.fetal)

fetal_mod <- model.matrix(~ modelGroup	+ Age+
											RIN + 
											mitoRate +  totalAssignedGene +
											Sex + snpPC1 
											, data=pd.fetal)

vGene = voom(dge.fetal,fetal_mod,plot=FALSE)
sv.fetal = sva(vGene$E,fetal_mod,fetal_mod0)
mod_gene_fetal <- cbind(fetal_mod,sv.fetal$sv)
vGene = voom(dge.fetal,mod_gene_fetal,plot=FALSE)

fitGene = lmFit(vGene, mod_gene_fetal)
eBGene = eBayes(fitGene)
fetalGene = topTable(eBGene,sort.by='p',coef="modelGroupSmoker",
    n=Inf)

## Saving results	
#save(fetalGene, adultGene, file= "/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/Genes_DE_sva.rda")

### Exons

eIndex <- which(rowMeans(cpm(exonCounts[,pd.model$RNum])>1)>0.1)
#eIndex = which(rowMeans(exonRpkm[,pd.model$RNum]) > .1) # Filtering genes by expression levels 	
#gIndex <- which(rowMeans(cpm(geneCounts[,pd.model$RNum]))>.1)

#exonRpkm.filt = exonRpkm[eIndex,pd.model$RNum] #keeping only the exons above the rpkm threshold (threshold 0.01 here)
exonCounts.filt = exonCounts[eIndex,pd.model$RNum] #keeping only the exons above the rpkm threshold (threshold 0.01 here)
exonMap.filt = exonMap[eIndex,] #keeping only the exon map for exons above the rpkm threshold

# Adult
dge.adult <- DGEList(counts = exonCounts.filt[,pd.adult$RNum],
    genes = exonMap.filt)
dge.adult = 	calcNormFactors(dge.adult)

adult_mod0 <- model.matrix(~ Age+ RIN + 
											mitoRate + totalAssignedGene + 
											Sex + snpPC1
									   		, data=pd.adult)

adult_mod <- model.matrix(~ modelGroup	+ Age+
											RIN + 
											mitoRate + totalAssignedGene + 
											Sex + snpPC1 
											, data=pd.adult)

vExon = voom(dge.adult,adult_mod,plot=FALSE)
sv.adult = sva(vExon$E,adult_mod,adult_mod0)
mod_exon_adult <- cbind(adult_mod,sv.adult$sv)
vExon = voom(dge.adult,mod_exon_adult,plot=FALSE)

fitExon = lmFit(vExon, mod_exon_adult)
eBExon = eBayes(fitExon)

adultExon = topTable(eBExon,sort.by='p',coef="modelGroupSmoker",
    n=Inf)

# Fetal

#eIndex.fetal <- which(rowSums(cpm(exonCounts[,pd.fetal$RNum])>.5)>=3 )
#
##gIndex <- which(rowMeans(cpm(geneCounts[,pd.model$RNum]))>.1)
#
#exonRpkm.fetal = exonRpkm[eIndex.fetal,pd.fetal$RNum] #keeping only the exons above the rpkm threshold (threshold 0.01 here)
#exonCounts.fetal = exonCounts[eIndex.fetal,pd.fetal$RNum] #keeping only the exons above the rpkm threshold (threshold 0.01 here)
#exonMap.fetal = exonMap[eIndex.fetal,] #keeping only the gene map for genes above the rpkm threshold


dge.fetal <- DGEList(counts = exonCounts.filt[,pd.fetal$RNum],
    genes = exonMap.filt)
dge.fetal = calcNormFactors(dge.fetal)

fetal_mod0 <- model.matrix(~ Age+
											RIN + 
											mitoRate + totalAssignedGene +
											Sex + snpPC1
											, data=pd.fetal)

fetal_mod <- model.matrix(~ modelGroup	+ Age+
											RIN + 
											mitoRate + totalAssignedGene + 
											Sex + snpPC1 
											, data=pd.fetal)

vExon = voom(dge.fetal,fetal_mod,plot=FALSE)
sv.fetal = sva(vExon$E,fetal_mod,fetal_mod0)
mod_exon_fetal <- cbind(fetal_mod,sv.fetal$sv)
vExon = voom(dge.fetal,mod_exon_fetal,plot=FALSE)

fitExon = lmFit(vExon, mod_exon_fetal)
eBExon = eBayes(fitExon)

fetalExon = topTable(eBExon,sort.by='p',coef="modelGroupSmoker",
    n=Inf)
## Saving results	
#save(fetalExon, adultExon, file= "/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/Exons_DE_sva.rda")

### Junctions

jIndex <- which(rowMeans(cpm(as.data.frame(jCounts[,pd.model$RNum]))>1)>0.1 & jMap$code != "Novel")
#jIndex = which(rowMeans(jRpkm[,pd.model$RNum]) > .1) # Filtering genes by expression levels 	
#gIndex <- which(rowMeans(cpm(geneCounts[,pd.model$RNum]))>.1)

#jRpkm.filt = jRpkm[jIndex,pd.model$RNum] #keeping only the exons above the rpkm threshold (threshold 0.01 here)
jCounts.filt = as.matrix(jCounts[jIndex,pd.model$RNum]) #keeping only the exons above the rpkm threshold (threshold 0.01 here)
jMap.filt = as.data.frame(jMap[jIndex,]) #keeping only the exon map for exons above the rpkm threshold

# Adult
dge.adult <- DGEList(counts = jCounts.filt[,pd.adult$RNum],
    genes = jMap.filt)
dge.adult = 	calcNormFactors(dge.adult)

adult_mod0 <- model.matrix(~ Age+ RIN + 
											mitoRate + totalAssignedGene +
											Sex + snpPC1
											, data=pd.adult)

adult_mod <- model.matrix(~ modelGroup	+ Age+
											RIN + 
											mitoRate + totalAssignedGene +
											Sex + snpPC1 
											, data=pd.adult)

vJxn = voom(dge.adult,adult_mod,plot=FALSE)
sv.adult = sva(vJxn$E,adult_mod,adult_mod0)
mod_jxn_adult <- cbind(adult_mod,sv.adult$sv)
vJxn = voom(dge.adult,mod_jxn_adult,plot=FALSE)

fitJxn = lmFit(vJxn, mod_jxn_adult)
eBJxn = eBayes(fitJxn)
adultJxn = topTable(eBJxn,sort.by='p',coef="modelGroupSmoker",
    n=Inf)

# Fetal
dge.fetal <- DGEList(counts = jCounts.filt[ ,pd.fetal$RNum] ,
    genes = jMap.filt)
dge.fetal = calcNormFactors(dge.fetal)

fetal_mod0 <- model.matrix(~ Age+
											RIN + 
											mitoRate + totalAssignedGene +
											Sex + snpPC1
											, data=pd.fetal)

fetal_mod <- model.matrix(~ modelGroup	+ Age+
											RIN + 
											mitoRate + totalAssignedGene + 
											Sex + snpPC1 
											, data=pd.fetal)

vJxn = voom(dge.fetal,fetal_mod,plot=FALSE)
sv.fetal = sva(vJxn$E,fetal_mod,fetal_mod0)
mod_jxn_fetal <- cbind(fetal_mod,sv.fetal$sv)
vJxn = voom(dge.fetal,mod_jxn_fetal,plot=FALSE)

fitJxn = lmFit(vJxn, mod_jxn_fetal)
eBJxn = eBayes(fitJxn)

fetalJxn = topTable(eBJxn,sort.by='p',coef="modelGroupSmoker",
    n=Inf)
## Saving results	
#save(fetalJxn, adultJxn, file= "/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/Jxns_DE_sva.rda")

### ERs
load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/processed_covMat.rda')
regionMat.model <- round(regionMat[, pd.model$RNum])

# Adult
dge.adult <- DGEList(counts = regionMat.model[,pd.adult$RNum],
    genes = as.data.frame(regions))
dge.adult = 	calcNormFactors(dge.adult)

adult_mod0 <- model.matrix(~ Age+ RIN + 
											mitoRate + totalAssignedGene +
											Sex + snpPC1
											, data=pd.adult)

adult_mod <- model.matrix(~ modelGroup	+ Age+
											RIN + 
											mitoRate + totalAssignedGene + 
											Sex + snpPC1 
											, data=pd.adult)

vER = voom(dge.adult,adult_mod,plot=FALSE)
sv.adult = sva(vER$E,adult_mod,adult_mod0)
mod_er_adult <- cbind(adult_mod,sv.adult$sv)
vER = voom(dge.adult,mod_er_adult,plot=FALSE)

fitER = lmFit(vER, mod_er_adult)
eBER = eBayes(fitER)

adultER = topTable(eBER,sort.by='p',coef="modelGroupSmoker",
    n=Inf)

# Fetal
dge.fetal <- DGEList(counts = regionMat.model[,pd.fetal$RNum],
    genes = as.data.frame(regions))
dge.fetal = calcNormFactors(dge.fetal)
fetal_mod0 <- model.matrix(~ Age+
											RIN + 
											mitoRate + totalAssignedGene +
											Sex + snpPC1
											, data=pd.fetal)

fetal_mod <- model.matrix(~ modelGroup	+ Age+
											RIN + 
											mitoRate + totalAssignedGene + 
											Sex + snpPC1 
											, data=pd.fetal)

vER = voom(dge.fetal,fetal_mod,plot=FALSE)
sv.fetal = sva(vER$E,fetal_mod,fetal_mod0)
mod_er_fetal <- cbind(fetal_mod,sv.fetal$sv)
vER = voom(dge.fetal,mod_er_fetal,plot=FALSE)

fitER = lmFit(vER, mod_er_fetal)
eBER = eBayes(fitER)
fetalER = topTable(eBER,sort.by='p', coef="modelGroupSmoker", n=Inf)
## Saving results	
#save(fetalER, adultER, file= "/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/ERs_DE_sva.rda")
## Saving mod

save(mod_gene_adult, mod_gene_fetal, mod_exon_adult, mod_exon_fetal, mod_jxn_adult, mod_jxn_fetal, mod_er_adult, mod_er_fetal, file = "/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/modObjects_DE_sva.rda")
