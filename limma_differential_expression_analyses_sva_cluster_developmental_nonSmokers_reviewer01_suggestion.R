#qsub -l bluejay  -l mf=50G,h_vmem=70G,h_stack=256M -cwd -b y R CMD BATCH limma_differential_expression_analyses_sva_cluster_developmental_nonSmokers_reviewer01_suggestion.R

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
## Adding qSVs

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


###


### Genes

#gIndex = which(rowMeans(geneRpkm[,pd.model$RNum]) > .1) # Filtering genes by expression levels 	
gIndex <- which(rowMeans(cpm(geneCounts[,pd.model$RNum])>1)>0.1)
#gIndex <- which(rowSums(cpm(geneCounts[,pd.model$RNum])>1)>=30)


#geneRpkm.filt = geneRpkm[gIndex,pd.model$RNum] #keeping only the exons above the rpkm threshold (threshold 0.01 here)
geneCounts.filt = geneCounts[gIndex,pd.model[pd.model$modelGroup=='Non-Smoker','RNum'] ] #keeping only the exons above the rpkm threshold (threshold 0.01 here)
geneMap.filt = geneMap[gIndex,] #keeping only the gene map for genes above the rpkm threshold

# Interaction
dge = DGEList(counts = geneCounts.filt, 
    genes = geneMap.filt)
dge = calcNormFactors(dge)

mod0 <- model.matrix(~ RIN + mitoRate + totalAssignedGene +
											Sex + snpPC1
											, data=pd.model[pd.model$modelGroup=='Non-Smoker',])

mod <- model.matrix(~  AgeGroup+
							RIN + 
							mitoRate + totalAssignedGene+  
							Sex + snpPC1, data=pd.model[pd.model$modelGroup=='Non-Smoker',])

vGene = voom(dge,mod,plot=FALSE)
sv.int_gene = sva(vGene$E,mod,mod0)
mod_gene <- cbind(mod,sv.int_gene$sv)
vGene = voom(dge,mod_gene,plot=FALSE)
fitGene = lmFit(vGene, mod_gene)
eBGene = eBayes(fitGene)
mainGene = topTable(eBGene,sort.by='p',coef="AgeGroupFetal",
    n=Inf)

### Exons
eIndex <- which(rowMeans(cpm(exonCounts[,pd.model$RNum])>1)>0.1)
#eIndex = which(rowMeans(exonRpkm[,pd.model$RNum]) > .1) # Filtering genes by expression levels 	
#gIndex <- which(rowMeans(cpm(geneCounts[,pd.model$RNum]))>.1)

#exonRpkm.filt = exonRpkm[eIndex,pd.model$RNum] #keeping only the exons above the rpkm threshold (threshold 0.01 here)
exonCounts.filt = exonCounts[eIndex,pd.model[pd.model$modelGroup=='Non-Smoker','RNum']] #keeping only the exons above the rpkm threshold (threshold 0.01 here)
exonMap.filt = exonMap[eIndex,] #keeping only the exon map for exons above the rpkm threshold

# Interaction
dge = DGEList(counts = exonCounts.filt, 
    genes = exonMap.filt)
dge = calcNormFactors(dge)

vExon = voom(dge,mod,plot=FALSE)
sv.int_exon = sva(vExon$E,mod,mod0)
mod_exon <- cbind(mod,sv.int_exon$sv)
vExon = voom(dge, mod_exon,plot=FALSE)
fitExon = lmFit(vExon, mod_exon)

eBExon = eBayes(fitExon)

mainExon = topTable(eBExon,sort.by='p',coef="AgeGroupFetal",
    n=Inf)

### Junctions

jIndex <- which(rowMeans(cpm(as.data.frame(jCounts[,pd.model$RNum]))>1)>0.1 & jMap$code != "Novel")
#jIndex = which(rowMeans(jRpkm[,pd.model$RNum]) > .1) # Filtering genes by expression levels 	
#gIndex <- which(rowMeans(cpm(geneCounts[,pd.model$RNum]))>.1)

#jRpkm.filt = jRpkm[jIndex,pd.model$RNum] #keeping only the exons above the rpkm threshold (threshold 0.01 here)
jCounts.filt = as.matrix(as.data.frame(jCounts[jIndex,pd.model[pd.model$modelGroup=='Non-Smoker','RNum']])) #keeping only the exons above the rpkm threshold (threshold 0.01 here)
jMap.filt = as.data.frame(jMap[jIndex,]) #keeping only the exon map for exons above the rpkm threshold

# Interaction
dge = DGEList(counts = jCounts.filt, 
    genes = jMap.filt )
dge = calcNormFactors(dge)

vJxn = voom(dge,mod,plot=FALSE)
sv.int_jxn = sva(vJxn$E,mod,mod0)
mod_jxn <- cbind(mod,sv.int_jxn$sv)

vJxn = voom(dge,mod_jxn,plot=FALSE)
fitJxn = lmFit(vJxn, mod_jxn)
eBJxn = eBayes(fitJxn)
mainJxn = topTable(eBJxn,sort.by='p',coef="AgeGroupFetal",
    n=Inf)
	
### ERs
load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/processed_covMat.rda')
regionMat.model <- round(regionMat[, pd.model[pd.model$modelGroup=='Non-Smoker','RNum'] ])

# Interaction
dge = DGEList(counts = regionMat.model, 
    genes = as.data.frame(regions))
dge = calcNormFactors(dge)

vER = voom(dge,mod,plot=FALSE)
sv.int_er = sva(vER$E,mod,mod0)
mod_er <- cbind(mod,sv.int_er$sv)

vER = voom(dge,mod_er,plot=FALSE)
fitER = lmFit(vER, mod_er)
eBER = eBayes(fitER)

mainER = topTable(eBER,sort.by='p',coef="AgeGroupFetal",
    n=Inf)

## Saving results	
save(mainGene, mainExon, mainJxn,  mainER,file= "/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/developmentModel_allFeatures_DE_sva_nonSmokeExposed.rda")

## saving model
save(mod_gene, mod_exon, mod_jxn,mod_er, file="/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/developmentModel_mod_DE_sva_nonSmokeExposed.rda")