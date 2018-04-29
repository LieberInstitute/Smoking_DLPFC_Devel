library(tidyr)
library(dplyr)
library(ggplot2)
theme_set(theme_bw(base_size=18) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="bottom"	)) 
				 
###
load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda') #rpkm data
load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/processed_covMat.rda')
regionMat <- round(regionMat)

setwd('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq')
load('rdas/Genes_DE_sva.rda')
load('rdas/Exons_DE_sva.rda')
load('rdas/Jxns_DE_sva.rda')
load('rdas/ERs_DE_sva.rda')
load('rdas/interactionModel_allFeatures_DE_sva.rda')
load('rdas/modObjects_DE_sva.rda')

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
clean_geneRpkm = jaffelab::cleaningY(y= log2(geneRpkm[,pd.fetal$RNum]+1), mod =mod_gene_fetal, P=2 )

dat_fetal = cbind(pd.fetal, t(jaffelab::cleaningY( log2(geneRpkm[rownames(fetalGene[fetalGene$adj.P.Val<0.10,]), pd.fetal$RNum] +1), mod =mod_gene_fetal, P=2)  ),t(jRpkm[rownames(fetalJxn[fetalJxn$adj.P.Val<0.10,]), pd.fetal$RNum]  ), regionMat[rownames(fetalER[fetalER$adj.P.Val<0.10,]), pd.fetal$RNum]    )				 				 
				 
### Checking enrichment with Geschwind autism DE
setwd('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq')
load('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/allFeatures_bound.rda')

### Test prenatal features 
TESTED_UNIV <- unique(all_bound$EnsemblGeneID)
prenatal_de_asd <- unique(all_bound[all_bound$adj.P.Val_Prenatal <0.1 & all_bound$Type=="Gene" ,'EnsemblGeneID'] )

### Pull in geschwind data
load("/dcl01/lieber/ajaffe/psychENCODE_Data/UCLA_R01MH094714/rpkmCounts_asd.rda")
load("/dcl01/lieber/ajaffe/psychENCODE_Data/UCLA_R01MH094714/rdas/degradation_mat_UCLA_ASD_RiboZero_v2.rda")
## filter
gIndex = which(rowMeans(geneRpkm) > 0.1)
geneRpkm = geneRpkm[gIndex,]
geneMap = geneMap[gIndex,]


## stat list quality
rIndexes = split(1:nrow(pd), pd$Region)
pd$Diagnosis = factor(pd$Diagnosis, levels = c("CTL", "ASD"))

cleanedQualData = lapply(rIndexes, function(ii) {
	cat(".")
	mod = model.matrix(~Diagnosis + totalAssignedGene +
		Sequencing.Batch + Brain.Bank + RIN + Age + Sex, data =pd[ii,])
	degPca = prcomp(t(log2(degCovAdj[,ii]+1))) # do PCA
	k = sva::num.sv(log2(degCovAdj[,ii]+1), mod) # in sva package
	qSVs = degPca$x[,1:k] # identify quality surrogate variables
	mod2 = cbind(mod, qSVs)
	res=jaffelab::cleaningY(y=log2(geneRpkm[,ii]+1), mod=mod2, P=2)
return(res)
})
cleanedQualData_all = do.call("cbind", cleanedQualData) 
asd_cleaned = cleanedQualData_all[prenatal_de_asd, ]
asd_cleaned = t(asd_cleaned)
asd_cleaned = cbind(asd_cleaned, pd[match(rownames(asd_cleaned),pd$Sample.Name),c('Diagnosis','Region') ] )

### Plotting
load('/dcl01/lieber/ajaffe/psychENCODE_Data/UCLA_R01MH094714/rdas/ASD_DE_gene.rda') #out
library(gridExtra)
library(grid)	   
library(ggsignif)
library(ggpubr)
lay <- rbind(c(1,1,2,2,2) )


dat_fetal$modelGroup=plyr::mapvalues(dat_fetal$modelGroup,c("Non-Smoker", "Smoker"), c("Unexposed","Exposed") )
pdf('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/plots/all_prenatal_de_genes_asd_nicotine_boxplots.pdf')
for (i in prenatal_de_asd){
tmp_asd =  asd_cleaned[,c('Region','Diagnosis',i)]
colnames(tmp_asd)[3] = "AdjustedRPKM"

tmp_nic = dat_fetal[, c('modelGroup',i)]
colnames(tmp_nic)[2] = "AdjustedRPKM"

asd_plot = ggplot(data=tmp_asd, aes(x=Region , y = AdjustedRPKM ,fill=Diagnosis ) ) +geom_boxplot(outlier.colour = NA, alpha = 0, aes(fill=Diagnosis)) + 
		geom_point(alpha = 1, aes(colour = Diagnosis), position = position_jitterdodge()) + 
		labs(x = "Brain Region", 
			 y = "Adjusted Expression",
			 title = paste0("Autism") ) +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + 			 
		theme(axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		legend.background = element_rect(colour = "black"),
		legend.justification = c(1, 1), 
        legend.position = c(1, 1),
		legend.title = element_blank(),
		axis.text.x = element_text(angle = 45, hjust = 1) )			  

nic_plot = ggplot(data=tmp_nic, aes(x=modelGroup, y = AdjustedRPKM ) ) +geom_boxplot(outlier.colour = NA, alpha = 0) + 
		geom_jitter(alpha = 1, width=0.2) + 
		labs(x = "Brain Region", 
			 y = "Adjusted Expression",
			 title =  paste0("Smoking") ) +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + 			 
		theme(axis.title.x=element_blank(),
		legend.background = element_rect(colour = "black"),
		legend.justification = c(1, 1), 
        legend.position = c(1, 1),
		legend.title = element_blank(),
		axis.text.x = element_text(angle = 45, hjust = 1))
grid.arrange(nic_plot,asd_plot,layout_matrix = lay, top=textGrob(geneMap[match(i,rownames(geneMap)),'Symbol'],gp=gpar(fontface='italic',fontsize=25))) 				
}
dev.off()

myplots <- list()  # new empty list
for (i in prenatal_de_asd[1:9]) local ( {
i <- i

tmp_asd =  asd_cleaned[,c('Region','Diagnosis',i)]
colnames(tmp_asd)[3] = "AdjustedRPKM"

tmp_nic = dat_fetal[, c('modelGroup',i)]
colnames(tmp_nic)[2] = "AdjustedRPKM"

asd_plot = ggplot(data=tmp_asd, aes(x=Region , y = AdjustedRPKM ,fill=Diagnosis ) ) +geom_boxplot(outlier.colour = NA, alpha = 0, aes(fill=Diagnosis)) + 
		geom_point(alpha = 1, aes(colour = Diagnosis), position = position_jitterdodge()) + 
		labs(x = "Brain Region", 
			 y = "Adjusted Expression",
			 title = paste0("Autism") ) +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + 			 
		theme(axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		legend.background = element_rect(colour = "black"),
		legend.justification = c(1, 1), 
        legend.position = c(1, 1),
		legend.title = element_blank(),
		axis.text.x = element_text(angle = 45, hjust = 1) )			  

nic_plot = ggplot(data=tmp_nic, aes(x=modelGroup, y = AdjustedRPKM ) ) +geom_boxplot(outlier.colour = NA, alpha = 0) + 
		geom_jitter(alpha = 1, width=0.2) + 
		labs(x = "Brain Region", 
			 y = "Adjusted Expression",
			 title =  paste0("Smoking") ) +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + 			 
		theme(axis.title.x=element_blank(),
		legend.background = element_rect(colour = "black"),
		legend.justification = c(1, 1), 
        legend.position = c(1, 1),
		legend.title = element_blank(),
		axis.text.x = element_text(angle = 45, hjust = 1))

myplots[[i]] <<- (list(asd_plot,nic_plot))
})

print(myplots[[1]][[1]])



		   
pdf('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/plots/5_Prenatal_DE_ASD_Genes.pdf')
#ENSG00000109158
ENSG00000109158_asd = ggplot(data=asd_cleaned, aes(x=Region , y = ENSG00000109158,fill=Diagnosis ) ) +geom_boxplot(outlier.colour = NA, alpha = 0, aes(fill=Diagnosis)) + 
		geom_point(alpha = 1, aes(colour = Diagnosis), position = position_jitterdodge()) + 
		labs(x = "Brain Region", 
			 y = "Adjusted Expression",
			 title = "Autism") +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + 			 
		theme(axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		legend.background = element_rect(colour = "black"),
		legend.justification = c(1, 1), 
        legend.position = c(1, 1),
		legend.title = element_blank() ) 
ENSG00000109158_nic = ggplot(data=dat_fetal, aes(x=modelGroup, y = ENSG00000109158 ) ) +geom_boxplot(outlier.colour = NA, alpha = 0) + 
		geom_jitter(alpha = 1, width=0.2) + 
		labs(x = "Brain Region", 
			 y = "Adjusted Expression",
			 title = "Nicotine Exposure") +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + 			 
		theme(axis.title.x=element_blank(),
		legend.background = element_rect(colour = "black"),
		legend.justification = c(1, 1), 
        legend.position = c(1, 1),
		legend.title = element_blank() ) 
grid.arrange(ENSG00000109158_nic,ENSG00000109158_asd,layout_matrix = lay, top=textGrob("GABRA4",gp=gpar(fontface='italic',fontsize=25))) 
ENSG00000109158 <- arrangeGrob(ENSG00000109158_nic,ENSG00000109158_asd, ncol=2,nrow=1, layout_matrix = lay) #generates g

###ENSG00000144619
ENSG00000144619_asd = ggplot(data=asd_cleaned, aes(x=Region , y = ENSG00000144619,fill=Diagnosis ) ) +geom_boxplot(outlier.colour = NA, alpha = 0, aes(fill=Diagnosis)) + 
		geom_point(alpha = 1, aes(colour = Diagnosis), position = position_jitterdodge()) + 
		labs(x = "Brain Region", 
			 y = "Adjusted Expression",
			 title = "Autism") +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + 			 
		theme(axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		legend.background = element_rect(colour = "black"),
		legend.justification = c(1, 1), 
        legend.position = c(1, 1),
		legend.title = element_blank() ) 
ENSG00000144619_nic = ggplot(data=dat_fetal, aes(x=modelGroup, y = ENSG00000144619 ) ) +geom_boxplot(outlier.colour = NA, alpha = 0) + 
		geom_jitter(alpha = 1, width=0.2) + 
		labs(x = "Brain Region", 
			 y = "Adjusted Expression",
			 title = "Nicotine Exposure") +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + 			 
		theme(axis.title.x=element_blank(),
		legend.background = element_rect(colour = "black"),
		legend.justification = c(1, 1), 
        legend.position = c(1, 1),
		legend.title = element_blank() ) 
grid.arrange(ENSG00000144619_nic,ENSG00000144619_asd,layout_matrix = lay, top=textGrob("CNTN4",gp=gpar(fontface='italic',fontsize=25)) ) 
ENSG00000144619 <- arrangeGrob(ENSG00000144619_nic,ENSG00000144619_asd, ncol=2,nrow=1, layout_matrix = lay, top=textGrob("CNTN4",gp=gpar(fontface='italic',fontsize=25))) #generates g
	
###ENSG00000091129
ENSG00000091129_asd = ggplot(data=asd_cleaned, aes(x=Region , y = ENSG00000091129,fill=Diagnosis ) ) +geom_boxplot(outlier.colour = NA, alpha = 0, aes(fill=Diagnosis)) + 
		geom_point(alpha = 1, aes(colour = Diagnosis), position = position_jitterdodge()) + 
		labs(x = "Brain Region", 
			 y = "Adjusted Expression",
			 title = "Autism") +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + 			 
		theme(axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		legend.background = element_rect(colour = "black"),
		legend.justification = c(1, 1), 
        legend.position = c(1, 1),
		legend.title = element_blank() ) 
ENSG00000091129_nic = ggplot(data=dat_fetal, aes(x=modelGroup, y = ENSG00000091129 ) ) +geom_boxplot(outlier.colour = NA, alpha = 0) + 
		geom_jitter(alpha = 1, width=0.2) + 
		labs(x = "Brain Region", 
			 y = "Adjusted Expression",
			 title = "Nicotine Exposure") +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + 			 
		theme(axis.title.x=element_blank(),
		legend.background = element_rect(colour = "black"),
		legend.justification = c(1, 1), 
        legend.position = c(1, 1),
		legend.title = element_blank() ) 
grid.arrange(ENSG00000091129_nic,ENSG00000091129_asd,layout_matrix = lay, top=textGrob("NRCAM",gp=gpar(fontface='italic',fontsize=25))  ) 
ENSG00000091129 <- arrangeGrob(ENSG00000091129_nic,ENSG00000091129_asd, ncol=2,nrow=1, layout_matrix = lay) #generates g

###ENSG00000189108
ENSG00000189108_asd = ggplot(data=asd_cleaned, aes(x=Region , y = ENSG00000189108,fill=Diagnosis ) ) +geom_boxplot(outlier.colour = NA, alpha = 0, aes(fill=Diagnosis)) + 
		geom_point(alpha = 1, aes(colour = Diagnosis), position = position_jitterdodge()) + 
		labs(x = "Brain Region", 
			 y = "Adjusted Expression",
			 title = "Autism") +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + 			 
		theme(axis.title.x=element_blank(),
		axis.title.y=element_blank(),
		legend.background = element_rect(colour = "black"),
		legend.justification = c(1, 1), 
        legend.position = c(1, 1),
		legend.title = element_blank() )		
ENSG00000189108_nic = ggplot(data=dat_fetal, aes(x=modelGroup, y = ENSG00000189108 ) ) +geom_boxplot(outlier.colour = NA, alpha = 0) + 
		geom_jitter(alpha = 1, width=0.2) + 
		labs(x = "Brain Region", 
			 y = "Adjusted Expression",
			 title = "Nicotine Exposure") +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + 			 
		theme(axis.title.x=element_blank(),
		legend.background = element_rect(colour = "black"),
		legend.justification = c(1, 1), 
        legend.position = c(1, 1),
		legend.title = element_blank() ) 
grid.arrange(ENSG00000189108_nic,ENSG00000189108_asd,layout_matrix = lay, top=textGrob("IL1RAPL2",gp=gpar(fontface='italic',fontsize=25)) ) 
ENSG00000109158 <- arrangeGrob(ENSG00000189108_nic,ENSG00000189108_asd, ncol=2,nrow=1, layout_matrix = lay) #generates g
dev.off()


 