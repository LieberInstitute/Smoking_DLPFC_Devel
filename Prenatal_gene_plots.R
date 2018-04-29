#Plotting a few fetal genes of interest
load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda') #rpkm data
load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/processed_covMat.rda')
regionMat <- round(regionMat)

setwd('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq')
load('rdas/Genes_DE_sva.rda')
load('rdas/Exons_DE_sva.rda')
load('rdas/Jxns_DE_sva.rda')
load('rdas/ERs_DE_sva.rda')
load('rdas/interactionModel_allFeatures_DE_sva.rda')

##load svs
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

pd.adult = pd.model[pd.model$AgeGroup=="Adult",]
pd.fetal = pd.model[pd.model$AgeGroup=="Fetal",]

#### Plot for main effect and interaction
clean_geneRpkm = jaffelab::cleaningY(y= geneRpkm[,pd.fetal$RNum], mod =mod_gene_fetal, P=2 )

dat_fetal = cbind(pd.fetal, t(jaffelab::cleaningY( log2(geneRpkm[rownames(fetalGene[fetalGene$adj.P.Val<0.10,]), pd.fetal$RNum] +1), mod =mod_gene_fetal, P=2)  ),t(jRpkm[rownames(fetalJxn[fetalJxn$adj.P.Val<0.10,]), pd.fetal$RNum]  ), t(jaffelab::cleaningY( log2(regionMat[rownames(fetalER[fetalER$adj.P.Val<0.10,]), pd.fetal$RNum,drop=FALSE]+1),mod =mod_gene_fetal, P=2)    ) )

kk = match(colnames(dat_fetal), rownames(fetalGene) )
colnames(dat_fetal)[!is.na(kk)] <- fetalGene[kk[!is.na(kk)],'Symbol']
colnames(dat_fetal)[ncol(dat_fetal)] <- rownames(fetalER[fetalER$adj.P.Val<0.10,])
dat_fetal$modelGroup = plyr::mapvalues(dat_fetal$modelGroup, c("Non-Smoker", "Smoker"), c("Unexposed", "Exposed"))

dat_fetal_genes = dat_fetal[,!is.na(kk)]

#############  
library(pheatmap)
ann = dat_fetal[,c('modelGroup','Race','Sex','Age','RIN')]
anno_colors <- list(modelGroup = c('red','blue') )

pdf('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/Paper/Prenatal_heatmap_Figure1.pdf')
pheatmap(t(scale(as.matrix(dat_fetal_genes), center = TRUE, scale = TRUE) ),color=rev(gray.colors(100, start = 0, end = 1, gamma = 2.2, alpha = NULL)), annotation_col=ann, labels_col =dat_fetal$BrNum )
dev.off()
 colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
#############
## GABRA4

y_axis_title <- expression(log[2]~"("~italic("GABRA4")~" RPKM+1)")
GABRA4 = ggplot(data=dat_fetal,aes(x=modelGroup, y= GABRA4) ) +geom_boxplot(outlier.colour = NA, alpha = 0) + 
		geom_point(alpha = 1, aes(colour = modelGroup), position = position_jitterdodge()) + 
		labs(x = "Age Group", 
			 y = "Adjusted RPKM",
			 title="GABRA4"
			 ) +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + 			 
		geom_text(size=5, aes(x=-Inf,y=Inf,hjust= -.2,vjust=1.8,label=paste0("p=", signif(fetalGene[fetalGene$Symbol=="GABRA4","P.Value"],3) ) )) +
		theme(axis.title.x=element_blank(),
		legend.background = element_rect(colour = "black"),
		legend.justification = c(1, 1), 
		legend.title = element_blank() ) 

## KCNN2
y_axis_title <- expression(log[2]~"("~italic("KCNN2")~" RPKM+1)")
KCNN2 = ggplot(data=dat_fetal,aes(x=modelGroup, y= KCNN2 ) ) +geom_boxplot(outlier.colour = NA, alpha = 0) + 
		geom_point(alpha = 1, aes(colour = modelGroup), position = position_jitterdodge()) + 
		labs(x = "Age Group", 
			 y = "Adjusted RPKM",
			 title="KCNN2") +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + 			 
		geom_text(size=5, aes(x=Inf,y=Inf,hjust= 1.2,vjust=1.8,label=paste0("p=", signif(fetalGene[fetalGene$Symbol=="KCNN2","P.Value"],3) ) )) +
		theme(axis.title.x=element_blank(),
		legend.background = element_rect(colour = "black"),
		legend.justification = c(1, 1), 
		legend.title = element_blank() ) 
		
## NRCAM
y_axis_title <- expression(log[2]~"("~italic("NRCAM")~" RPKM+1)")
NRCAM = ggplot(data=dat_fetal,aes(x=modelGroup, y= log2(NRCAM+1) ) ) +geom_boxplot(outlier.colour = NA, alpha = 0) + 
		geom_point(alpha = 1, aes(colour = modelGroup), position = position_jitterdodge()) + 
		labs(x = "Age Group", 
			 y = "Adjusted RPKM",
			 title="NRCAM") +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + 			 
		geom_text(size=5, aes(x=Inf,y=Inf,hjust= 1.2,vjust=1.8,label=paste0("p=", signif(fetalGene[fetalGene$Symbol=="NRCAM","P.Value"],3) ) )) +
		theme(axis.title.x=element_blank(),
		legend.background = element_rect(colour = "black"),
		legend.justification = c(1, 1), 
		legend.title = element_blank() ) 


## CNTN4
y_axis_title <- expression(log[2]~"("~italic("CNTN4")~" RPKM+1)")
CNTN4 = ggplot(data=dat_fetal,aes(x=modelGroup, y= log2(CNTN4+1) ) ) +geom_boxplot(outlier.colour = NA, alpha = 0) + 
		geom_point(alpha = 1, aes(colour = modelGroup), position = position_jitterdodge()) + 
		labs(x = "Age Group", 
			 y = "Adjusted RPKM",
			 title="CNTN4"
			 ) +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + 			 
		geom_text(size=5, aes(x=-Inf,y=Inf,hjust= -.2,vjust=1.8,label=paste0("p=", signif(fetalGene[fetalGene$Symbol=="CNTN4","P.Value"],3) ) )) +
		theme(axis.title.x=element_blank(),
		legend.background = element_rect(colour = "black"),
		legend.justification = c(1, 1), 
		legend.title = element_blank() ) 

## EPHA8 
y_axis_title <- expression(log[2]~"("~italic("EPHA8")~" RPKM+1)")
EPHA8 = ggplot(data=dat_fetal,aes(x=modelGroup, y= EPHA8 ) ) +geom_boxplot(outlier.colour = NA, alpha = 0) + 
		geom_point(alpha = 1, aes(colour = modelGroup), position = position_jitterdodge()) + 
		labs(x = "Age Group", 
			 y ="Adjusted RPKM",
			 title ="EPHA8"
			 ) +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + 			 
		geom_text(size=5,aes(x=-Inf,y=Inf,hjust= -.2,vjust=1.8,label=paste0("p=", signif(fetalGene[fetalGene$Symbol=="EPHA8","P.Value"],3) ) )) +
		theme(axis.title.x=element_blank(),
		legend.background = element_rect(colour = "black"),
		legend.justification = c(1, 1), 
		legend.title = element_blank() ) 		

## PCDH10
y_axis_title <- expression(log[2]~"("~italic("PCDH10")~" RPKM+1)")
PCDH10 = ggplot(data=dat_fetal,aes(x=modelGroup, y= log2(er101585+1) ) ) +geom_boxplot(outlier.colour = NA, alpha = 0) + 
		geom_point(alpha = 1, aes(colour = modelGroup), position = position_jitterdodge()) + 
		labs(x = "Age Group", 
			 y = "Adjusted Coverage",
			 title = "PCDH10"
			 ) +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + 			 
		geom_text(size=5,aes(x=Inf,y=Inf,hjust= 1.2,vjust=1.8,label=paste0("p=", signif(fetalER[rownames(fetalER)=="er101585","P.Value"],3) ) )) +
		theme(axis.title.x=element_blank(),
		legend.background = element_rect(colour = "black"),
		legend.justification = c(1, 1), 
        legend.position = "none",
		legend.title = element_blank() ) 		
#########
library(gridExtra)
library(grid)
			   
Prenatal_boxplot <- arrangeGrob(GABRA4, KCNN2, NRCAM,CNTN4,EPHA8,PCDH10, ncol=3,nrow=2) #generates g
ggsave(Prenatal_boxplot, file="/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/Paper/Prenatal_boxplot_Figure1_big.pdf", height=8.5,width=11) 	 
