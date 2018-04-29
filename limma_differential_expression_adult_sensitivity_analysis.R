#qsub -l bluejay -l mf=40G,h_vmem=60G,h_stack=256M -cwd -b y R CMD BATCH --no-save limma_differential_expression_adult_sensitivity_analysis.R

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
jIndex <- which(rowMeans(cpm(as.data.frame(jCounts[,pd.model$RNum]))>1)>0.1 & jMap$code != "Novel")


### Reducing to sensitivty model
pd.adult$cotinine_blood_number <- NA
pd.adult$cotinine_blood_number <- ifelse(grepl('ng/mL', pd.adult$nicotine_comments), pd.adult$nicotine_comments, NA)
pd.adult$cotinine_blood_number <- gsub('.*cotinine', '', pd.adult$cotinine_blood_number, ignore.case = TRUE) #removing all text before ; (usually nicotine said first)
pd.adult$cotinine_blood_number <- gsub("[^0-9|.]", "", pd.adult$cotinine_blood_number)
pd.adult$cotinine_blood_number <- as.numeric(pd.adult$cotinine_blood_number)
#Assigning groups
pd.adult$ordinalModel <- ifelse( !is.na(pd.adult$cotinine_blood_number) & pd.adult$AgeGroup=="Adult", ifelse(pd.adult$cotinine_blood_number > 200, "Heavy Smoker", "Light Smoker"), NA) 
pd.adult$ordinalModel[!(pd.adult$nicotine | pd.adult$cotinine) & pd.adult$AgeGroup=="Adult"] <- "Non-Smoker"

pd.adult$ordinalModel <- factor(pd.adult$ordinalModel, levels = c("Non-Smoker", "Light Smoker", "Heavy Smoker"))
pd.adult <- pd.adult[!is.na(pd.adult$ordinalModel),]
#pd.model = pd.model[-which(pd.model$smoking & pd.model$ordinalModel =="Non-Smoker"),] #dropping ambiguous non-smokers
pd.adult$ordinalModel  <- as.numeric(pd.adult$ordinalModel)
## Adding qSVs

geneCounts.filt = geneCounts[gIndex,pd.adult$RNum] #keeping only the exons above the rpkm threshold (threshold 0.01 here)
geneMap.filt = geneMap[gIndex,] #keeping only the gene map for genes above the rpkm threshold

# Adult
dge.adult <- DGEList(counts = geneCounts.filt[,pd.adult$RNum],
    genes = geneMap.filt)
dge.adult = 	calcNormFactors(dge.adult)

adult_mod0 <- model.matrix(~ Age+ RIN + 
											mitoRate + totalAssignedGene + 
											Sex + snpPC1
											, data=pd.adult)

adult_mod <- model.matrix(~ ordinalModel	+ Age+
											RIN + 
											mitoRate +  totalAssignedGene + 
											Sex + snpPC1 
											, data=pd.adult)

vGene = voom(dge.adult,adult_mod,plot=FALSE)
sv.adult = sva(vGene$E,adult_mod,adult_mod0)
mod_gene_ord <- cbind(adult_mod,sv.adult$sv)
vGene = voom(dge.adult,mod_gene_ord,plot=FALSE)

fitGene = lmFit(vGene, mod_gene_ord)
eBGene = eBayes(fitGene)

adultGene_Ordinal = topTable(eBGene,sort.by='p',coef="ordinalModel",
    n=Inf)

### Exons

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

adult_mod <- model.matrix(~ ordinalModel	+ Age+
											RIN + 
											mitoRate + totalAssignedGene + 
											Sex + snpPC1 
											, data=pd.adult)

vExon = voom(dge.adult,adult_mod,plot=FALSE)
sv.adult = sva(vExon$E,adult_mod,adult_mod0)
mod_exon_ord <- cbind(adult_mod,sv.adult$sv)
vExon = voom(dge.adult,mod_exon_ord,plot=FALSE)

fitExon = lmFit(vExon, mod_exon_ord)
eBExon = eBayes(fitExon)

adultExon_Ordinal = topTable(eBExon,sort.by='p',coef="ordinalModel",
    n=Inf)
### Junctions


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

adult_mod <- model.matrix(~ ordinalModel	+ Age+
											RIN + 
											mitoRate + totalAssignedGene +
											Sex + snpPC1 
											, data=pd.adult)

vJxn = voom(dge.adult,adult_mod,plot=FALSE)
sv.adult = sva(vJxn$E,adult_mod,adult_mod0)
mod_jxn_ord <- cbind(adult_mod,sv.adult$sv)
vJxn = voom(dge.adult,mod_jxn_ord,plot=FALSE)

fitJxn = lmFit(vJxn, mod_jxn_ord)
eBJxn = eBayes(fitJxn)
adultJxn_Ordinal = topTable(eBJxn,sort.by='p',coef="ordinalModel",
    n=Inf)

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

adult_mod <- model.matrix(~ ordinalModel	+ Age+
											RIN + 
											mitoRate + totalAssignedGene + 
											Sex + snpPC1 
											, data=pd.adult)

vER = voom(dge.adult,adult_mod,plot=FALSE)
sv.adult = sva(vER$E,adult_mod,adult_mod0)
mod_er_ord <- cbind(adult_mod,sv.adult$sv)
vER = voom(dge.adult,mod_er_ord,plot=FALSE)

fitER = lmFit(vER, mod_er_ord)
eBER = eBayes(fitER)

adultER_Ordinal = topTable(eBER,sort.by='p',coef="ordinalModel",
    n=Inf)
#save(adultGene_Ordinal, adultExon_Ordinal, adultJxn_Ordinal, adultER_Ordinal,file='rdas/ordinal_sensitivity_model_results.rda')	
save(mod_gene_ord, mod_exon_ord, mod_jxn_ord, mod_er_ord, file='rdas/ordinal_sensitivity_mod_DE.rda')	
