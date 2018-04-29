#Schizophrenia smoking replication
##################### 
#Setup 
setwd('/users/ssemick/Tox_Expr/Adult_Smoking/RNA_Seq')
load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rawCounts_szControlDLPFC.rda') #load this file for the mapObjects 
load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda') #load rpkm data to filter junctions before DESeq2 pipeline--speed up computation!
library(LIBDpheno)
library(ggplot2)
library(reshape2)
library(ggplot2)
theme_set(theme_bw(base_size=18) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"))

#####################
#Add in LIMS information to data on cluster which isn't otherwise available
LIMS = toxicant[[1]]
pd.model = pd[pd$Dx =="Control",]
pd.model = pd.model[pd.model$Age <0,]
id <- match( brnum(pd.model$BrNum), LIMS$brnumerical  )
id <- id[!is.na(id)]
pd.model <- cbind( LIMS[id,], pd.model)
pd.model$source <- droplevels(pd.model$source)
pd.model <- pd.model[,colSums(sapply(pd.model, is.na))!=nrow(pd.model)]
pd.model <- pd.model[ ,!colnames(pd.model) %in% 
						colnames(pd.model[,sapply(pd.model, is.logical)])[!sapply( pd.model[,sapply(pd.model, is.logical)], any, na.rm=TRUE)] ] #drop logicals without any trues
#seeing smoking and toxicant status
table(pd.model$cotinine| pd.model$BrNum %in% c('Br1813','Br1826','Br2047', "Br2049") )

pd.model[pd.model$BrNum %in% c('Br1813','Br1826','Br2047', "Br2049"),
		 c("nicotine_comments","other_drugs_comments","final_dx_comments") ]   
pd.model$modelGroup <- pd.model$cotinine| pd.model$BrNum %in% c('Br1813','Br1826','Br2047', "Br2049")
pd.model <- pd.model[!is.na(pd.model$modelGroup),]
pd.model$modelGroup <- as.factor(pd.model$modelGroup)
pd.model$modelGroup <- plyr::revalue(pd.model$modelGroup, c("TRUE" = "Smoker", "FALSE" = "Non-Smoker"))

pd.model <- pd.model[!pd.model$BrNum%in%c("Br1779","Br1794"),]
###############

######Cotinine Brain Number
pd.model[pd.model$modelGroup=="Smoker","nicotine_comments"]
pd.model$cotinine_brain_number <- NA
pd.model$cotinine_brain_number[grepl('Negative in brain',pd.model$nicotine_comments)] <- 0
pd.model$cotinine_brain_number[grepl('pg/mg', pd.model$nicotine_comments)] <- pd.model$nicotine_comments[grepl('pg/mg', pd.model$nicotine_comments)]
pd.model$cotinine_brain_number <- gsub('pg/mg.*', '', pd.model$cotinine_brain_number, ignore.case = TRUE) #removing all text after ;
pd.model$cotinine_brain_number <- gsub("[^0-9|.]", "", pd.model$cotinine_brain_number)
pd.model$cotinine_brain_number <- as.numeric(pd.model$cotinine_brain_number)
#pd.model$nicotine_brain_number[which(!pd.model$nicotine & !df$cotinine & df$primarydx=='Control'& grepl('CB', df$nicotine_comments))] <- 0
###############
setwd('plots')
setwd('Fetal')
######### Cotinine distribution
smoke_conc_brain <- ggplot (data = subset(pd.model,modelGroup=="Smoker"), aes(x = "", y = cotinine_brain_number))
smoke_conc_brain <- smoke_conc_brain + geom_boxplot(outlier.colour = NA, alpha = 0.1) + 
		geom_point(alpha = 1, position = position_jitter(width = 0.2)) +
		theme_bw(base_size = 18) +
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
			  legend.position="none") + 
			  scale_fill_grey() +scale_colour_grey()+
		labs(x = 'Prenatal Nicotine Exposed', y="Cotinine Concentration in Brain (pg/mg)")
pdf(file="Cotinine_Distributions_QC.pdf",height=8.5,width=11)
smoke_conc_brain
dev.off()
######### Age
pd.model$modelGroup<-plyr::revalue(pd.model$modelGroup, c("Non-Smoker"="Unexposed", "Smoker"="Exposed") )
t.test(Age~modelGroup, data = pd.model) #NS
age_plot <- ggplot (data = pd.model, aes(x = modelGroup, y = 40+52*Age))
age_plot <- age_plot + geom_boxplot(outlier.colour = NA, alpha = 0.1, aes(fill = modelGroup)) + 
		geom_point(alpha = 1, aes(colour = modelGroup), position = position_jitter(width = 0.2)) +
		theme_bw(base_size = 18) + 
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
			  legend.position="none") +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + labs(x="Exposure Status",y = "Gestational Age (Weeks)")
		   
#################### RIN
t.test(RIN~modelGroup, data = pd.model) #sig
rin_plot <- ggplot (data = pd.model, aes(x = modelGroup, y = RIN))
rin_plot <- rin_plot + geom_boxplot(outlier.colour = NA, alpha = 0.1, aes(fill = modelGroup)) + 
		geom_point(alpha = 1, aes(colour = modelGroup), position = position_jitter(width = 0.2)) +
		theme_bw(base_size = 18) + 
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
			  axis.title.x=element_blank(),
			  legend.position="none") +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + labs(x="Exposure Status", y = "RIN")
################## PMI
t.test(PMI~modelGroup, data = pd.model) #Ns
pmi_plot <- ggplot (data = pd.model, aes(x = modelGroup, y = PMI))
pmi_plot <- pmi_plot + geom_boxplot(outlier.colour = NA, alpha = 0.1, aes(fill = modelGroup)) + 
		geom_point(alpha = 1, aes(colour = modelGroup), position = position_jitter(width = 0.2)) +
		theme_bw(base_size = 18) + 
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
			  axis.title.x=element_blank(),
			  legend.position="none") +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + labs(x="Exposure Status", y="PMI")
################## mitoRate

t.test(mitoRate~modelGroup, data = pd.model) #ns
mito_plot <- ggplot (data = pd.model, aes(x = modelGroup, y = mitoRate))
mito_plot <- mito_plot + geom_boxplot(outlier.colour = NA, alpha = 0.1, aes(fill = modelGroup)) + 
		geom_point(alpha = 1, aes(colour = modelGroup), position = position_jitter(width = 0.2)) +
		theme_bw(base_size = 18) + 
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
			  axis.title.x=element_blank(),
			  legend.position="none") +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + labs(x="Exposure Status", y= "Mitochondrial Mapping Rate")
############### Gene Assignment Rate
t.test(totalAssignedGene~modelGroup, data = pd.model) #Ns
tag_plot <- ggplot (data = pd.model, aes(x = modelGroup, y = totalAssignedGene))
tag_plot <- tag_plot + geom_boxplot(outlier.colour = NA, alpha = 0.1, aes(fill = modelGroup)) + 
		geom_point(alpha = 1, aes(colour = modelGroup), position = position_jitter(width = 0.2)) +
		theme_bw(base_size = 18) + 
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
			  axis.title.x=element_blank(),
			  legend.position="none") +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + labs(x="Nicotine Exposure", y="Gene Assignment Rate")
############### mapping Rate
t.test(mappingRate~modelGroup, data = pd.model) #NS
map_plot <- ggplot (data = pd.model, aes(x = modelGroup, y = mappingRate))
map_plot <- map_plot + geom_boxplot(outlier.colour = NA, alpha = 0.1, aes(fill = modelGroup)) + 
		geom_point(alpha = 1, aes(colour = modelGroup), position = position_jitter(width = 0.2)) +
		theme_bw(base_size = 18) + 
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
			  axis.title.x=element_blank(),
			  legend.position="none") +
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + labs(x="Exposure Status", y = "Mapping Rate")		   

###################### Supplementary figure 1
library(gridExtra)
age_plot <- age_plot + theme_bw(base_size = 18) + theme(plot.title = element_text(hjust = 0.5),
														panel.grid.major = element_blank(),
														panel.grid.minor = element_blank(),
														axis.title.x=element_blank(),
														legend.position="none")
rin_plot <- rin_plot + theme_bw(base_size = 18)  + theme(plot.title = element_text(hjust = 0.5),
														panel.grid.major = element_blank(),
														panel.grid.minor = element_blank(),
														axis.title.x=element_blank(),
														legend.position="none")
pmi_plot <- pmi_plot + theme_bw(base_size = 18) + theme(plot.title = element_text(hjust = 0.5),
														panel.grid.major = element_blank(),
														panel.grid.minor = element_blank(),
														axis.title.x=element_blank(),
														legend.position="none")
mito_plot <- mito_plot + theme_bw(base_size = 18) + theme(plot.title = element_text(hjust = 0.5),
														panel.grid.major = element_blank(),
														panel.grid.minor = element_blank(),
														axis.title.x=element_blank(),
														legend.position="none")
tag_plot <- tag_plot + theme_bw(base_size = 18) + theme(plot.title = element_text(hjust = 0.5),
														panel.grid.major = element_blank(),
														panel.grid.minor = element_blank(),
														axis.title.x=element_blank(),
														legend.position="none")
map_plot <- map_plot + theme_bw(base_size = 18) + theme(plot.title = element_text(hjust = 0.5),
														panel.grid.major = element_blank(),
														panel.grid.minor = element_blank(),
														axis.title.x=element_blank(),
														legend.position="none")

Demo_EDA <- arrangeGrob(age_plot, rin_plot, pmi_plot, mito_plot, tag_plot, map_plot, ncol=3,nrow=2)
ggsave( Demo_EDA, file="Exploratory_Data_Analysis_Confounds_v_SmokingExposure.pdf",
           height=8.5,width=11)		
########################## Ancestry PCs
t.test(snpPC1~modelGroup, data = pd.model) #Ns .09
t.test(snpPC2~modelGroup, data = pd.model) #Ns
t.test(snpPC3~modelGroup, data = pd.model) #Ns
t.test(snpPC4~modelGroup, data = pd.model) #Ns
t.test(snpPC5~modelGroup, data = pd.model) #Ns
###############################

""
#> pd.model[,sapply(pd.model, is.numeric)]
#> names(pd.model[,sapply(pd.model, is.numeric)])
s

##############################
fisher.test(table(pd.model$Race, pd.model$modelGroup)) # .08
p <- ggplot(data = pd.model, aes(x = Race,fill =modelGroup))
p <- p + geom_bar(aes(y = (..count..)/sum(..count..)),position = "dodge") +
            theme_bw(base_size = 14) + 
  						  theme(plot.title = element_text(hjust = 0.5),
								panel.grid.major = element_blank(),
								panel.grid.minor = element_blank())+
            labs(x = "Race",
                 y = "Count",
                 title = "Race Distribution for Schizophrenics" +
    				scale_colour_brewer(palette="Set1") + 
    				scale_fill_brewer(palette="Set1") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust  = 1))
ggsave(p, file = "Fetal Smoking Race.png",
		   type = "cairo-png",
           height=5.5,width=9,dpi=600,units="in")				   
fisher.test(table(pd.model$Sex, pd.model$modelGroup)) # ns .0922

##################### Supplemental Toxicant Information

#################################################################
jIndex = which(jMap$code != "Novel")#Drop all novel junctions
jRpkm = jRpkm[jIndex,] #keeping only the non-novel junctions

gRpkm_rowmean <- as.data.frame(log2(rowMeans(geneRpkm[,pd.model$RNum])+1))
colnames(gRpkm_rowmean) <- "Rowmeans"
gRpkm_rowmean$type <- "Gene"

eRpkm_rowmean <- as.data.frame(log2(rowMeans(exonRpkm[,pd.model$RNum])+1))
colnames(eRpkm_rowmean) <- "Rowmeans"
eRpkm_rowmean$type <- "Exon"

jRpkm_rowmean <- as.data.frame(log2(rowMeans(jRpkm[,pd.model$RNum])+1))
colnames(jRpkm_rowmean) <- "Rowmeans"
jRpkm_rowmean$type <- "Junction"

model_rowmean = rbind(gRpkm_rowmean, eRpkm_rowmean, jRpkm_rowmean)
rpkm_density <- ggplot(data = model_rowmean, aes(x = Rowmeans, fill = type ) )
rpkm_density <- rpkm_density + stat_ecdf(geom = "step", alpha = .8) + 
							   geom_vline(data= data.frame(Thresholds=c("0.01","0.1","1.0"),vals = c(log2(0.01 +1),log2(0.1 +1),log2(1.0 +1))),aes(xintercept=vals, linetype=Thresholds, col=Thresholds), show.legend=TRUE) +
							   guides(fill=FALSE) +
							   facet_wrap( ~ type, scales = "free") + 
							   labs(x = "Log2(Mean Normalized Expression+1)",
									y = "Cumulative Density") +
								theme_bw(base_size = 18) + 
								theme(plot.title = element_text(hjust = 0.5),
								panel.grid.major = element_blank(),
								panel.grid.minor = element_blank(),
								legend.title=element_text(size=12),
								legend.text=element_text(size=10),
								legend.justification=c(1,0), 
								legend.position=c(.99,.05)) +
								scale_linetype_manual(values=c(2,1,3)) + 
								scale_color_manual(values = c("black","red","black") ) +
								coord_cartesian(xlim = c(0, 10)) + scale_y_continuous(limits = c(0,1),expand=c(0,0)) 
ggsave(rpkm_density, file = "Prenatal RPKM Density.pdf",
           height=8.5,width=11)										
###########################
#Pre-PCA Filtering by RPKM
geneRpkm.model <- geneRpkm[,pd.model$RNum]
gIndex=which(rowMeans(geneRpkm.model) > .1) 
geneRpkm.model_filt <- geneRpkm.model[gIndex,]

exonRpkm.model <- exonRpkm[,pd.model$RNum]
jRpkm.model <- jRpkm[,pd.model$RNum]
length(gIndex)
###########################
getPcaVars <- function(pca, digits = 3) {
    stopifnot(is(pca) == 'prcomp')
    signif(pca$sdev^2 / sum(pca$sdev^2) * 100, digits = digits)
}
pca = prcomp(t(log2(geneRpkm.model_filt+1)))
dat1 <- cbind(pca$x, pd.model)
pcaVars = getPcaVars(pca)
sum(pcaVars>1) #first 19 PCs contain more than 1 % of variance
sum(pcaVars[1:5]) #first 5 explain 64 percent of variance
library(ggplot2)
dat1$modelGroup <- plyr::revalue(dat1$modelGroup, c("Non-Smoker"="Unexposed","Smoker"="Exposed") )
PC1_2 <- ggplot(data = dat1, aes(x = PC1, y =PC2) )
PC1_2 <- PC1_2 + geom_point(aes(colour=modelGroup), alpha = 1) +
				 labs(x = paste0("PC1: ", getPcaVars(pca)[1], "% Var Expl"),
					  y = paste0("PC2: ", getPcaVars(pca)[2], "% Var Expl"),
					  colour = "Nicotine Exposure") +
				scale_colour_brewer(palette="Set1") + theme(legend.position='none',
		legend.background = element_rect(colour='black',fill=NA, size=.5, linetype="solid")) 


#ggsave(PC1_2,file = "Prenatal_PC1_2 Smoking Highlighted PFC.pdf",
#           height=8.5,width=11)									
						

PC2_3 <- ggplot(data = dat1, aes(x = PC2, y =PC3) )
PC2_3 <- PC2_3 + geom_point(aes(colour=modelGroup), alpha = 1.0) +
				 labs(x = paste0("PC2: ", getPcaVars(pca)[2], "% Var Expl"),
					  y = paste0("PC3: ", getPcaVars(pca)[3], "% Var Expl"),
					  colour = "Nicotine Exposure") +
				scale_colour_brewer(palette="Set1") + theme(legend.position='none',
		legend.background = element_rect(colour='black',fill=NA, size=.5, linetype="solid")) 
#ggsave(PC2_3,file = "Prenatal_PC2_3 Smoking Highlighted PFC.pdf",
#           height=8.5,width=11)									
		   
PC3_4 <- ggplot(data = dat1, aes(x = PC3, y =PC4) )
PC3_4 <- PC3_4 + geom_point(aes(colour=modelGroup), alpha = 1.0) +
				 labs(x = paste0("PC3: ", getPcaVars(pca)[3], "% Var Expl"),
					  y = paste0("PC4: ", getPcaVars(pca)[4], "% Var Expl"),
					  colour = "Nicotine Exposure") +
					  scale_colour_brewer(palette="Set1")  + theme(legend.position='none',
		legend.background = element_rect(colour='black',fill=NA, size=.5, linetype="solid")) 
#ggsave(PC3_4, file = "Prenatal_PC3_4 Smoking Highlighted PFC.pdf",
#           height=8.5,width=11)									
								

PC4_5 <- ggplot(data = dat1, aes(x = PC4, y =PC5) )
PC4_5 <- PC4_5 + geom_point(aes(colour=modelGroup), alpha = 1.0) +
				 labs(x = paste0("PC4: ", getPcaVars(pca)[4], "% Var Expl"),
					  y = paste0("PC5: ", getPcaVars(pca)[5], "% Var Expl"),
					  colour = "Smoking") +
					  scale_colour_brewer(palette="Set1") + theme(legend.position='bottom') 
#ggsave(PC4_5, file = "Prenatal_PC4_5 Smoking Highlighted PFC.pdf",
#           height=8.5,width=11)									

### 
library(gridExtra)
library(grid)			   
exprs_pca_boxplot <- arrangeGrob(PC1_2, PC2_3, PC3_4,PC4_5, ncol=2,nrow=2) #generates g
ggsave(exprs_pca_boxplot, file="/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/Paper/pcaExprs_boxplot_panel_fetal.pdf", height=11,width=8.5) 	 
##########################
#Check PC1
PC = "PC1"
t.test(dat1[,PC] ~ dat1$modelGroup)
t.test(dat1[,PC] ~ as.factor(dat1$Sex))
t.test(dat1[,PC] ~ as.factor(dat1$Race))
cor.test(dat1[,PC], dat1$Age) #sig
cor.test(dat1[,PC], dat1$RIN)
cor.test(dat1[,PC], dat1$mitoRate)
cor.test(dat1[,PC], dat1$mitoMapped)
cor.test(dat1[,PC], dat1$mappingRate)
cor.test(dat1[,PC], dat1$totalAssignedGene)
cor.test(dat1[,PC], dat1$totalMapped)
cor.test(dat1[,PC], dat1$PMI)

PC1_Age <- ggplot(data = dat1, aes(x=PC1, y = Age) )
PC1_Age <- PC1_Age + geom_point() + 
				 geom_smooth(method="lm", se =FALSE) + 
				 scale_colour_brewer(palette="Set1") +
				 theme_bw(base_size = 18) + 
				 labs(x = paste0("PC1: ", getPcaVars(pca)[1], "% Var Expl"),
					  y = "Age",
					  title = paste0("PC1 Associates with Age", "; P = ",
					  signif(cor.test(dat1[,PC], dat1$Age)$p.value,3) ))
ggsave(PC1_Age, file = "PC1 Associates with Age.png",
		   type = "cairo-png",
           height=5.5,width=9,dpi=600,units="in")									
####################################
#Check PC2
PC = "PC2"
t.test(dat1[,PC] ~ dat1$modelGroup)
t.test(dat1[,PC] ~ as.factor(dat1$Sex))
t.test(dat1[,PC] ~ as.factor(dat1$Race))
cor.test(dat1[,PC], dat1$Age)
cor.test(dat1[,PC], dat1$RIN)
cor.test(dat1[,PC], dat1$mitoRate) #
cor.test(dat1[,PC], dat1$mitoMapped)
cor.test(dat1[,PC], dat1$mappingRate)
cor.test(dat1[,PC], dat1$totalAssignedGene)
cor.test(dat1[,PC], dat1$totalMapped)
cor.test(dat1[,PC], dat1$PMI)

PC_mitoRate <- ggplot(data = dat1, aes(x=PC2, y = mitoRate) )
PC_mitoRate <- PC_mitoRate + geom_point() + 
				 scale_colour_brewer(palette="Set1") +
				 theme_bw(base_size = 18) + 
				 geom_smooth(method="lm", se =FALSE) + 
				 labs(x = paste0("PC2: ", getPcaVars(pca)[2], "% Var Expl"),
					  y = "Mito Mapping Rate",
					  title = paste0("PC2 Associates with mitoRate", "; P = ",
					  signif(cor.test(dat1[,PC], dat1$mitoRate)$p.value,3) ))
ggsave(PC_mitoRate, file = "PC2 Associates with mitoRate.png",
		   type = "cairo-png",
           height=5.5,width=9,dpi=600,units="in")									
####################################
#Check PC3
PC = "PC3"
t.test(dat1[,PC] ~ dat1$modelGroup)
t.test(dat1[,PC] ~ as.factor(dat1$Sex))
t.test(dat1[,PC] ~ as.factor(dat1$Race))
cor.test(dat1[,PC], dat1$Age)
cor.test(dat1[,PC], dat1$RIN)
cor.test(dat1[,PC], dat1$mitoRate)
cor.test(dat1[,PC], dat1$mitoMapped)
cor.test(dat1[,PC], dat1$mappingRate)
cor.test(dat1[,PC], dat1$totalAssignedGene)#
cor.test(dat1[,PC], dat1$totalMapped)
cor.test(dat1[,PC], dat1$PMI)

PC_tag <- ggplot(data = dat1, aes(x=PC3, y = totalAssignedGene) )
PC_tag <- PC_tag + geom_point() + 
				 geom_smooth(method="lm", se =FALSE) + 
				 scale_colour_brewer(palette="Set1") +
				 theme_bw(base_size = 18) + 
				 labs(x = paste0("PC3: ", getPcaVars(pca)[3], "% Var Expl"),
					  y = "Total Assigned Genes",
					  title = paste0("PC3 Associates with Total Assigned Genes", "; P = ",
					  signif(cor.test(dat1[,PC], dat1$totalAssignedGene)$p.value,3) ))
ggsave(PC_tag, file = "PC3 Associates with totalAssignedGene.png",
		   type = "cairo-png",
           height=5.5,width=9,dpi=600,units="in")									

####################################
#Check PC4
PC = "PC4"
t.test(dat1[,PC] ~ dat1$modelGroup)
t.test(dat1[,PC] ~ as.factor(dat1$Sex))
t.test(dat1[,PC] ~ as.factor(dat1$Race))
cor.test(dat1[,PC], dat1$Age)
cor.test(dat1[,PC], dat1$RIN)
cor.test(dat1[,PC], dat1$mitoRate)
cor.test(dat1[,PC], dat1$mitoMapped)
cor.test(dat1[,PC], dat1$mappingRate)
cor.test(dat1[,PC], dat1$totalAssignedGene)
cor.test(dat1[,PC], dat1$totalMapped)
cor.test(dat1[,PC], dat1$PMI)

PC_mapRate <- ggplot(data = dat1, aes(x=PC4, y = mappingRate) )
PC_mapRate <- PC_mapRate + geom_point() + 
				 geom_smooth(method="lm", se =FALSE) +
				 scale_colour_brewer(palette="Set1") +
				 theme_bw(base_size = 18) + 				 
				 labs(x = paste0("PC4: ", getPcaVars(pca)[4], "% Var Expl"),
					  y = "mappingRate",
					  title = paste0("PC4 Associates with mappingRate", "; P = ",
					  signif(cor.test(dat1[,PC], dat1$mappingRate)$p.value,3) ))
ggsave(PC_mapRate, file = "PC4 Associates with mappingRate.png",
		   type = "cairo-png",
           height=5.5,width=9,dpi=600,units="in")									

#Check PC5
PC = "PC5"
t.test(dat1[,PC] ~ dat1$modelGroup)
t.test(dat1[,PC] ~ as.factor(dat1$Sex))
t.test(dat1[,PC] ~ as.factor(dat1$Race))
cor.test(dat1[,PC], dat1$Age)
cor.test(dat1[,PC], dat1$RIN)
cor.test(dat1[,PC], dat1$mitoRate)
cor.test(dat1[,PC], dat1$mitoMapped)
cor.test(dat1[,PC], dat1$mappingRate)
cor.test(dat1[,PC], dat1$totalAssignedGene)
cor.test(dat1[,PC], dat1$totalMapped)
cor.test(dat1[,PC], dat1$PMI)

#doesn't strongly associate with anything

####################################################################
#Genotype data QCing
library(ggrepel)
pd.model$Race <- as.factor(pd.model$Race)

theme_set(theme_bw(base_size=18) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5)) )

snp1 <- ggplot(data = pd.model, aes(x = snpPC1, y = snpPC2, colour = Race))
snp1 <- snp1 + geom_point() +
			    labs(x = "MDS Dim 1",
					 y = "MDS Dim 2") +  
				scale_colour_grey() +
				theme(legend.position="none")+
				 geom_text_repel(data=dplyr::filter(pd.model, snpPC2< -0.05), aes(label=BrNum),colour="black")

snp2 <- ggplot(data = pd.model, aes(x = snpPC2, y = snpPC3, colour = Race))
snp2 <- snp2 + geom_point() +
			   labs(x = "MDS Dim 2",
					 y = "MDS Dim 3") +
				scale_colour_grey() +
				theme(legend.position="none")+
			   geom_text_repel(data=dplyr::filter(pd.model, snpPC2< -0.05), aes(label=BrNum),colour="black")

snp3 <- ggplot(data = pd.model, aes(x = snpPC3, y = snpPC4, colour = Race))
snp3 <- snp3 + geom_point() +
			   labs(x = "MDS Dim 3",
					 y = "MDS Dim 4") +
				scale_colour_grey() +
				theme(legend.position="none")+
			   geom_text_repel(data=dplyr::filter(pd.model, snpPC3 < -0.05 | snpPC4 < -0.03), aes(label=BrNum),colour="black")
		   
snp4 <- ggplot(data = pd.model, aes(x = snpPC4, y = snpPC5, colour = Race))
snp4 <- snp4 + geom_point() +
			   labs(x = "MDS Dim 4",
					 y = "MDS Dim 5") +
				scale_colour_grey() +
				theme(legend.position="bottom")+
			   geom_text_repel(data=dplyr::filter(pd.model, snpPC4 < -0.03), aes(label=BrNum),colour="black")
#			   geom_label(data = pd.model[pd.model$snpPC4 < -0.03,],
#							aes(label = BrNum, fill = factor(Race)), colour = "white", fontface = "bold")



################################### Supplemental Figure S5B
snp_misidentifier <- ggplot(data = pd.model, aes(x=Race, y = snpPC1, colour=Race) )
snp_misidentifier <- snp_misidentifier + 
		geom_boxplot(outlier.colour = NA, alpha = 0.1, aes(fill = Race)) +
		geom_point(alpha = 1, aes(colour = Race), position = position_jitter(width = 0.2)) +
		theme_bw(base_size = 18) + 
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
			  legend.position='none') +
			  scale_colour_grey() +
			  scale_fill_grey() +
		labs( x = "Race", y="MDS Dim 1" ) +
		geom_text_repel(data=dplyr::filter(pd.model, snpPC1 > 0 & Race=="CAUC"), aes(label=BrNum),colour="black")
pdf('MDS Genotype Exploration.pdf',width=11,height=8.5)
snp1
snp2
snp3
snp4
snp_misidentifier
dev.off()

### 
library(gridExtra)
library(grid)			   
snp_pca_boxplot <- arrangeGrob(snp1, snp2, snp3,snp4, ncol=2,nrow=2) #generates g
ggsave(snp_pca_boxplot, file="/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/Paper/pcaSnp_boxplot_panel_fetal.pdf", height=11,width=8.5) 	 
ggsave(snp_misidentifier, file="/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/Paper/misidentifed_person_fetal.pdf", height=11,width=8.5) 	 


####################################
setwd('../')

### Sex check
#see: http://www.nature.com/articles/ncomms3771
#Widespread sex differences in gene expression and splicing in the adult human brain
#c(,'DDX3Y')
sex_DEG <- c('NLGN4Y',
'USP9Y',
'DDX3Y', ##
'XIST', ##
'KDM5D',
'UTY',
'RPS4Y1', ##
'CYorf15B',
'EIF1AY', ##
'ZFY', ##
'TTTY14',
'PRKY',
'TMSB4Y',
'ZFX',
'DDX3X',
'KDM5C',
'HDHD1A',
'TSIX',
'DYZ1L25',
'KDM6A',
'ABCA6',
'HAND2')

sex_DEG <- c('DDX3Y', ##
'XIST', ##
'RPS4Y1', ##
'EIF1AY', ##
'ZFY')


library(biomaRt)
Human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
sex_check <- getBM(c('hgnc_symbol', 'ensembl_gene_id', 'start_position','end_position'), filters = 'hgnc_symbol', values = sex_DEG, mart = Human)
dat2 <- 
geneRpkm.model

geneCounts.model[sex_check[,'ensembl_gene_id'] == rownames(geneCounts.model), pd.model$Sex=="M"]


colnames(geneRpkm.model)==pd.model$RNum

dat2 = cbind(pd.model,t(geneRpkm.model[sex_check$ensembl_gene_id ,]) )
tmp <-dat2[,c('sex',sex_check$ensembl_gene_id)]
#tmp <- reshape2::melt(tmp, id.vars="sex")


sex_plot <- ggplot(data = tmp, aes(x= sex, y=ENSG00000067048)) 
sex_plot <- sex_plot + geom_boxplot(outlier.colour = NA, alpha = 0.1, aes(fill = sex)) + 
		geom_point(alpha = 1, aes(colour = sex), position = position_jitter(width = 0.2)) +
		theme_bw(base_size = 18) + 
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
			  legend.position="none") +
		scale_colour_brewer(palette="Set1") +
		scale_fill_brewer(palette="Set1") +
		labs( x = "Sex", y="RPKM", title = sex_check[sex_check$ensembl_gene_id=="ENSG00000067048",'hgnc_symbol'] )
ggsave(sex_plot, file="ENSG00000067048.png",
		   type = "cairo-png",
           height=5.5,width=9,dpi=600,units="in")

sex_plot <- ggplot(data = tmp, aes(x= sex, y=ENSG00000198692)) 
sex_plot <- sex_plot + geom_boxplot(outlier.colour = NA, alpha = 0.1, aes(fill = sex)) + 
		geom_point(alpha = 1, aes(colour = sex), position = position_jitter(width = 0.2)) +
		theme_bw(base_size = 18) + 
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
			  legend.position="none") +
		scale_colour_brewer(palette="Set1") +
		scale_fill_brewer(palette="Set1") +
		labs( x = "Sex", y="RPKM", title = sex_check[sex_check$ensembl_gene_id=="ENSG00000198692",'hgnc_symbol'] )
ggsave(sex_plot, file="ENSG00000198692.png",
		   type = "cairo-png",
           height=5.5,width=9,dpi=600,units="in")

sex_plot <- ggplot(data = tmp, aes(x= sex, y=ENSG00000129824)) 
sex_plot <- sex_plot + geom_boxplot(outlier.colour = NA, alpha = 0.1, aes(fill = sex)) + 
		geom_point(alpha = 1, aes(colour = sex), position = position_jitter(width = 0.2)) +
		theme_bw(base_size = 18) + 
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
			  legend.position="none") +
		scale_colour_brewer(palette="Set1") +
		scale_fill_brewer(palette="Set1") +
		labs( x = "Sex", y="RPKM", title = sex_check[sex_check$ensembl_gene_id=="ENSG00000129824",'hgnc_symbol'] )
ggsave(sex_plot, file="ENSG00000129824.png",
		   type = "cairo-png",
           height=5.5,width=9,dpi=600,units="in")

sex_plot <- ggplot(data = tmp, aes(x= sex, y=ENSG00000229807)) 
sex_plot <- sex_plot + geom_boxplot(outlier.colour = NA, alpha = 0.1, aes(fill = sex)) + 
		geom_point(alpha = 1, aes(colour = sex), position = position_jitter(width = 0.2)) +
		theme_bw(base_size = 18) + 
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
			  legend.position="none") +
		scale_colour_brewer(palette="Set1") +
		scale_fill_brewer(palette="Set1") +
		labs( x = "Sex", y="RPKM", title = sex_check[sex_check$ensembl_gene_id=="ENSG00000229807",'hgnc_symbol'] )
ggsave(sex_plot, file="ENSG00000229807.png",
		   type = "cairo-png",
           height=5.5,width=9,dpi=600,units="in")
		   
sex_plot <- ggplot(data = tmp, aes(x= sex, y=ENSG00000067646)) 
sex_plot <- sex_plot + geom_boxplot(outlier.colour = NA, alpha = 0.1, aes(fill = sex)) + 
		geom_point(alpha = 1, aes(colour = sex), position = position_jitter(width = 0.2)) +
		theme_bw(base_size = 18) + 
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
			  legend.position="none") +
		scale_colour_brewer(palette="Set1") +
		scale_fill_brewer(palette="Set1") +
		labs( x = "Sex", y="RPKM", title = sex_check[sex_check$ensembl_gene_id=="ENSG00000067646",'hgnc_symbol'] )
ggsave(sex_plot, file="ENSG00000067646.png",
		   type = "cairo-png",
           height=5.5,width=9,dpi=600,units="in")
		   
################### Tableone		   