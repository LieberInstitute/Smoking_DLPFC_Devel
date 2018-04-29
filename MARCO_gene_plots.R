#Plotting MARCO: only adult gene
load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda') #rpkm data
setwd('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq')
load('rdas/Genes_DE_sva.rda')
load('rdas/Exons_DE_sva.rda')
load('rdas/Jxns_DE_sva.rda')
load('rdas/ERs_DE_sva.rda')
load('rdas/interactionModel_allFeatures_DE_sva.rda')
### 
load('rdas/modObjects_DE_sva.rda')
###
library(ggplot2)
theme_set(theme_bw(base_size=14) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"				 )) 

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


###
clean_geneRpkm = jaffelab::cleaningY(y= log2(geneRpkm[,pd.adult$RNum]+1), mod =mod_gene_adult, P=2 )

#### Plot for main effect and interaction
dat = cbind(pd.adult, MARCO = clean_geneRpkm[rownames(adultGene[adultGene$adj.P.Val<.1,]), pd.adult$RNum] )
a = ggplot(data=dat,aes(x=modelGroup, y= MARCO) ) +geom_boxplot(outlier.colour = NA, alpha = 0) + 
		geom_point(alpha = 1, position = position_jitter(w = 0.2, h = 0) ) + 
		labs(x = "Smoking Status", 
			 y = "Adjusted RPKM",
			 title="Controls"
			 ) +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + 			 
		geom_text(aes(x=Inf,y=Inf,hjust=1.1,vjust=1.8,label=paste0("p=", signif(adultGene[adultGene$adj.P.Val<.1,"P.Value"],3) ) )) +
		theme(axis.title.x=element_blank(),
		legend.background = element_rect(colour = "black"),
		legend.justification = c(1, 1), 
        legend.position = c(1, 1),
		legend.title = element_blank() ) 

### Plot for adult controls: ordinal model boxplots ####
load('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/ordinal_sensitivity_model_results.rda')
load('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/ordinal_sensitivity_mod_DE.rda')
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


clean_geneRpkm = jaffelab::cleaningY(y= log2(geneRpkm[,pd.adult$RNum]+1), mod =mod_gene_ord, P=2 )

dat2 = cbind(pd.adult, MARCO = clean_geneRpkm[rownames(adultGene_Ordinal[adultGene_Ordinal$adj.P.Val<.1,]), pd.adult$RNum] )
b = ggplot(data=dat2,aes(x=ordinalModel, y= MARCO ) ) +geom_boxplot(outlier.colour = NA, alpha = 0 ) + 
		geom_point(alpha = 1,, position = position_jitter(width=0.2,h=0) ) + 
		labs(x = "Smoker Category", 
			 y = "Adjusted RPKM",
			 title="Controls"
			 ) +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + 			 
		geom_text(aes(x=Inf,y=Inf,hjust=1.1,vjust=1.8,label=paste0("p=", signif(adultGene_Ordinal[adultGene_Ordinal$adj.P.Val<.1,"P.Value"],3) ) )) +
		theme(axis.title.x=element_blank() )


#pd.model = pd.model[-which(pd.model$smoking & pd.model$ordinalModel =="Non-Smoker"),] #dropping ambiguous non-smokers
### Plot for adult schizo: smoker vs. non-smoker
load('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/smoking_schizophrenia_model_results.rda')
load('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/smoking_schizophrenia_mod_DE.rda')

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
pd.scz <-pd.scz[!pd.scz$BrNum%in%"Br1178",]

clean_geneRpkm = jaffelab::cleaningY(y= log2(geneRpkm[,pd.scz$RNum]+1), mod =mod_gene_scz, P=2 )
dat3 = cbind(pd.scz, MARCO = clean_geneRpkm['ENSG00000019169', pd.scz$RNum] )
scz = ggplot(data=dat3,aes(x=modelGroup, y= MARCO ) ) +geom_boxplot(outlier.colour = NA, alpha = 0 ) + 
		geom_point(alpha = 1, position = position_jitter(width=0.2,h=0) ) + 
		labs(x = "Smoking Status", 
			 y = "Adjusted RPKM",
			 title="Schizophrenia"
			 ) +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + 		
		theme(axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
		geom_text(aes(x=Inf,y=Inf,hjust=1.1,vjust=1.8,label=paste0("p=", signif(adultGene_SCZ['ENSG00000019169',"P.Value"],3) ) ))
#########
library(gridExtra)
library(grid)
lay <- rbind(c(1,2),
             c(3,3) )			   
MARCO_plot <- arrangeGrob(a, scz, b, ncol=2,nrow=2, layout_matrix = lay) #generates g
ggsave(MARCO_plot, file="/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/Paper/MARCO_Figure1_1016.pdf", height=8.5,width=8.5) 	 