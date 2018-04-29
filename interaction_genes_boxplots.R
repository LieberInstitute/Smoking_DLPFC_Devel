## Plotting a few fetal genes of interest
load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda') #rpkm data
load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/processed_covMat.rda')
regionMat <- round(regionMat)

setwd('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq')
load('rdas/Genes_DE_sva.rda')
load('rdas/Exons_DE_sva.rda')
load('rdas/Jxns_DE_sva.rda')
load('rdas/ERs_DE_sva.rda')
load("/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/interactionModel_mod_DE_sva.rda")

### load svs
load('rdas/modObjects_DE_sva.rda')
#mod_gene_fetal

### 
library(ggplot2)
theme_set(theme_bw(base_size=18) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 legend.position="none"	,
				 plot.title = element_text(hjust = 0.5, face='italic') )) 

### Plot for adult controls: smoker vs. non-smoker MARCO, + prenatal smoke exposed vs. prenatal unexposed
#Add in LIMS information to data on cluster which isn't otherwise available
library(LIBDpheno)
LIMS = toxicant[[1]]
pd.model = pd[pd$Dx =="Control",]
id <- match( brnum(pd.model$BrNum), LIMS$brnumerical  )
id <- id[!is.na(id)]
pd.model <- cbind( LIMS[id,], pd.model)
pd.model$source <- droplevels(pd.model$source)
pd.model <- pd.model[,colSums(sapply(pd.model, is.na))!=nrow(pd.model)]
pd.model <- pd.model[ ,!colnames(pd.model) %in% 
						colnames(pd.model[,sapply(pd.model, is.logical)])[!sapply( pd.model[,sapply(pd.model, is.logical)], any, na.rm=TRUE)] ] #drop logicals without any trues seeing smoking and toxicant status

###
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

pd.adult = pd.model[pd.model$AgeGroup=="Adult",]
pd.fetal = pd.model[pd.model$AgeGroup=="Fetal",]

#### Plot for main effect and interaction
clean_geneRpkm = jaffelab::cleaningY( y=geneRpkm[,pd.model$RNum], mod=mod_gene_int, P=2 )
####
modIndex = grep("modelGroupSmoker:AgeGroupFetal", colnames(mod_gene_int))
mod_gene_int = mod_gene_int[ ,c(1,2,modIndex, 3:(modIndex-1),(modIndex+1):ncol(mod_gene_int)) ]
mod_exon_int = mod_exon_int[ ,c(1,2,modIndex, 3:(modIndex-1),(modIndex+1):ncol(mod_exon_int)) ]
mod_jxn_int = mod_jxn_int[ ,c(1,2,modIndex, 3:(modIndex-1),(modIndex+1):ncol(mod_jxn_int)) ]
mod_er_int = mod_er_int[ ,c(1,2,modIndex, 3:(modIndex-1),(modIndex+1):ncol(mod_er_int)) ]

"modelGroupSmoker:AgeGroupFetal"

dat_int = cbind(pd.model, t(jaffelab::cleaningY( log2(geneRpkm[rownames(intGene[intGene$adj.P.Val<0.10,]), pd.model$RNum] +1), mod =mod_gene_int, P=2)  ),
t(jaffelab::cleaningY( log2(exonRpkm[rownames(intExon[intExon$adj.P.Val<0.10,]), pd.model$RNum,drop=FALSE]+1),mod =mod_exon_int, P=4) ),
t(jaffelab::cleaningY( log2(jRpkm[rownames(intJxn[intJxn$adj.P.Val<0.10,]), pd.model$RNum,drop=FALSE]+1),mod =mod_jxn_int, P=4) ),
t(jaffelab::cleaningY( log2(regionMat[rownames(intER[intER$adj.P.Val<0.10,]), pd.model$RNum,drop=FALSE]+1),mod =mod_er_int, P=4) ) )

kk = match(colnames(dat_int), rownames(intGene) )
colnames(dat_int)[!is.na(kk)] <- intGene[kk[!is.na(kk)],'Symbol']
dat_int$modelGroup = plyr::mapvalues(dat_int$modelGroup, c("Non-Smoker", "Smoker"), c("Unexposed", "Exposed"))

dat_fetal_genes = dat_fetal[,!is.na(kk)]

#############  
## GABRA4
y_axis_title <- expression(log[2]~"("~italic("GABRA4")~" RPKM+1)")
GABRA4 = ggplot(data=dat_int,aes(x=AgeGroup,fill=modelGroup, y= GABRA4) ) +geom_boxplot(outlier.colour = NA, alpha = 0) + geom_point(alpha = 1, aes(colour = modelGroup), position = position_jitterdodge()) + 
		labs(x = "Age Group", 
			 y = "Adjusted RPKM",
			 title="GABRA4"
			 ) +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") +
		theme(axis.title.x=element_blank(),
		legend.background = element_rect(colour = "black"),
		legend.title = element_blank(),
		legend.position='none') 

## KCNN2
y_axis_title <- expression(log[2]~"("~italic("C1orf115")~" RPKM+1)")
C1orf115 = ggplot(data=dat_int,aes(x=AgeGroup,fill=modelGroup, y= `C1orf115`) ) +
		geom_boxplot(outlier.colour = NA, alpha = 0) + geom_point(alpha = 1, aes(colour = modelGroup), position = position_jitterdodge()) + 
		labs(x = "Age Group", 
			 y = "Adjusted RPKM",
			 title="C1orf115"
			 ) +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") +
		theme(axis.title.x=element_blank(),
		legend.background = element_rect(colour = "black"),
		legend.title = element_blank(),
		legend.position='none') 
		
## NRCAM
y_axis_title <- expression(log[2]~"("~italic("NRCAM")~" RPKM+1)")
NRCAM = ggplot(data=dat_int,aes(x=AgeGroup,fill=modelGroup, y= NRCAM) ) +
		geom_boxplot(outlier.colour = NA, alpha = 0) + geom_point(alpha = 1, aes(colour = modelGroup), position = position_jitterdodge()) + 
		labs(x = "Age Group", 
			 y = "Adjusted RPKM",
			 title="NRCAM"
			 ) +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") +
		theme(axis.title.x=element_blank(),
		legend.background = element_rect(colour = "black"),
		legend.title = element_blank(),
		legend.position='none') 

## CNTN4
y_axis_title <- expression(log[2]~"("~italic("CDH8")~" RPKM+1)")
CDH8 = ggplot(data=dat_int,aes(x=AgeGroup,fill=modelGroup, y= CDH8) ) +
		geom_boxplot(outlier.colour = NA, alpha = 0) + geom_point(alpha = 1, aes(colour = modelGroup), position = position_jitterdodge()) + 
		labs(x = "Age Group", 
			 y = "Adjusted RPKM",
			 title="CDH8"
			 ) +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") +
		theme(axis.title.x=element_blank(),
		legend.background = element_rect(colour = "black"),
		legend.title = element_blank(),
        legend.justification = c(1, 1), 
		legend.position='bottom') 

#########
library(gridExtra)
library(grid)
			   
interaction_boxplot <- arrangeGrob(GABRA4, C1orf115, NRCAM,CDH8, ncol=2,nrow=2) #generates g
ggsave(interaction_boxplot, file="/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/Paper/Interaction_boxplot_panel_big.pdf", height=8.5,width=11) 	 
