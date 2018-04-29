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
geneCounts.filt = geneCounts[gIndex,pd.model$RNum] #keeping only the exons above the rpkm threshold (threshold 0.01 here)
geneMap.filt = geneMap[gIndex,] #keeping only the gene map for genes above the rpkm threshold

# Adult
dge.adult <- DGEList(counts = geneCounts.filt[,pd.adult$RNum],
    genes = geneMap.filt)
dge.adult = 	calcNormFactors(dge.adult)
adult_mod <- model.matrix(~ modelGroup	+ Age+
											RIN + 
											mitoRate +  totalAssignedGene + 
											Sex + snpPC1 
											, data=pd.adult)

vGene = voom(dge.adult,adult_mod,plot=FALSE)

#sv.adult = sva(vGene$E,adult_mod,adult_mod0)
#mod_gene_adult <- cbind(adult_mod,sv.adult$sv)
#vGene = voom(dge.adult,mod_gene_adult,plot=FALSE)

fitGene = lmFit(vGene, adult_mod)
eBGene = eBayes(fitGene)

adultGene_noSVA = topTable(eBGene,sort.by='p',coef="modelGroupSmoker",
    n=Inf)

# Fetal

dge.fetal <- DGEList(counts = geneCounts.filt[,pd.fetal$RNum],
    genes = geneMap.filt)
dge.fetal = calcNormFactors(dge.fetal)

fetal_mod <- model.matrix(~ modelGroup	+ Age+
											RIN + 
											mitoRate +  totalAssignedGene +
											Sex + snpPC1 
											, data=pd.fetal)

vGene = voom(dge.fetal,fetal_mod,plot=FALSE)

fitGene = lmFit(vGene, fetal_mod)
eBGene = eBayes(fitGene)
fetalGene_noSVA = topTable(eBGene,sort.by='p',coef="modelGroupSmoker",
    n=Inf)
## Saving results	

load("/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/Genes_DE_sva.rda")


pdf('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/plots/prenatal_secondary_de_analysis_Sva_effectSizes_vs_noSva_effectSizes.pdf',height=8,width=8)
#plot(fetalGene$logFC, fetalGene_noSVA[rownames(fetalGene),'logFC'])
cor.test(fetalGene$logFC, fetalGene_noSVA[rownames(fetalGene),'logFC'])

sigGene = rownames(fetalGene[fetalGene$adj.P.Val<0.1,])
par(mar=c(5,6,2,2))
plot(fetalGene[sigGene,'logFC'], fetalGene_noSVA[sigGene,'logFC'], pch = 21, bg="grey", main="Log2 Fold Change for prenatal genes called DE",
	xlab="Model with SVs", ylab="Model without SVs",
	xlim=c(-1.6,1.6), ylim = c(-1.6,1.6),
	cex.axis=2,cex.lab=2)
lines(x = c(-1.6,1.6), y = c(-1.6,1.6),col="red")

cor.test(fetalGene[sigGene,'logFC'], fetalGene_noSVA[sigGene,'logFC'])
dev.off()

par(mar=c(5,6,2,2))
plot(allRegion_mergedStats[,'ALL_subset_mainEffect_t'], 
	 Aging_controlOnly_mergedStats[rownames(allRegion_mergedStats),"CRB_Control_Aging_t"], pch = 21, bg="grey",
	xlab="AD Effect", ylab="Aging Effect (control only)",
	xlim=c(-8,8), ylim = c(-8,8),
	cex.axis=2,cex.lab=2)
abline(v=0,h=0, lty=2,col="blue",lwd=2)
lines(x = c(-8,8), y = c(-8,8),col="red")
dev.off()
