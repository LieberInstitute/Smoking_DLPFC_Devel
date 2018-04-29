setwd('/users/ssemick/Tox_Expr/Adult_Smoking/RNA_Seq')
load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rawCounts_szControlDLPFC.rda') #load this file for the mapObjects 
load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda') #load rpkm data to filter junctions before DESeq2 pipeline--speed up computation!
######
library(LIBDpheno)
library(ggplot2)
library(reshape2)
#####################
#Add in LIMS information to data on cluster which isn't otherwise available
LIMS = toxicant[[1]] 
pd.model = pd[pd$Dx =="Control",]
pd.model = pd.model[pd.model$Age >16,]
id <- match( brnum(pd.model$BrNum), LIMS$brnumerical  )
id <- id[!is.na(id)]
pd.model <- cbind( LIMS[id,], pd.model)
pd.model$source <- droplevels(pd.model$source)
pd.model <- pd.model[,colSums(sapply(pd.model, is.na))!=nrow(pd.model)]
pd.model <- pd.model[ ,!colnames(pd.model) %in% 
						colnames(pd.model[,sapply(pd.model, is.logical)])[!sapply( pd.model[,sapply(pd.model, is.logical)], any, na.rm=TRUE)] ] #drop logicals without any trues
#seeing smoking and toxicant status
table(tox=pd.model$cotinine | pd.model$nicotine, smoking=pd.model$smoking) 

pd.model$modelGroup <- pd.model$cotinine| pd.model$nicotine 
pd.model <- pd.model[!is.na(pd.model$modelGroup),]
pd.model$modelGroup <- as.factor(pd.model$modelGroup)
pd.model$modelGroup <- plyr::revalue(pd.model$modelGroup, c("TRUE" = "Smoker", "FALSE" = "Non-Smoker"))
pd.model = pd.model[-which(pd.model$smoking & pd.model$modelGroup =="Non-Smoker"),]

pd.model <- pd.model[!pd.model$BrNum %in% c("Br1179","Br1105","Br2267"),] #dropping outliers on SNP PCs
#dropping those with history of smoking but no toxicant history from analysis
######### 
#Boxplots of nicotine and cotinine levels for non-psychiatric adults
setwd('plots')
setwd('Adult')
nicotine_comments = pd.model$nicotine_comments[which(pd.model$modelGroup=="Smoker")] #what are the comments concerning nicotine in the blood of positive controls
nicotine_brain_comments <- nicotine_comments[ grepl('ng/g', nicotine_comments ) ] #comments that have nicotine in brain information
nicotine_blood_comments <- nicotine_comments[ grepl('ng/mL', nicotine_comments ) ] #comments that have nicotine in blood information

###### Nicotine Brain Numbers
pd.model$nicotine_brain_number <- NA
pd.model$nicotine_brain_number <- ifelse( grepl('ng/g', pd.model$nicotine_comments), pd.model$nicotine_comments, NA)
pd.model$nicotine_brain_number <- gsub(';.*', '', pd.model$nicotine_brain_number, ignore.case = TRUE) #removing all text after ;
pd.model$nicotine_brain_number <- gsub("[^0-9|.]", "", pd.model$nicotine_brain_number)
pd.model$nicotine_brain_number <- as.numeric(pd.model$nicotine_brain_number)
#pd.model$nicotine_brain_number[which(!pd.model$nicotine & !df$cotinine & df$primarydx=='Control'& grepl('CB', df$nicotine_comments))] <- 0

######Cotinine Brain Number
pd.model$cotinine_brain_number <- NA
pd.model$cotinine_brain_number <- ifelse( grepl('ng/g', pd.model$nicotine_comments), pd.model$nicotine_comments, NA)
pd.model$cotinine_brain_number <- gsub('.*cotinine', '', pd.model$cotinine_brain_number, ignore.case = TRUE) #removing all text after ;
pd.model$cotinine_brain_number <- gsub('nicotine.*', '', pd.model$cotinine_brain_number, ignore.case = TRUE) #removing all text after ;
pd.model$cotinine_brain_number <- gsub("[^0-9|.]", "", pd.model$cotinine_brain_number)
pd.model$cotinine_brain_number <- as.numeric(pd.model$cotinine_brain_number)
#pd.model$nicotine_brain_number[which(!pd.model$nicotine & !df$cotinine & df$primarydx=='Control'& grepl('CB', df$nicotine_comments))] <- 0

###### Nicotine Blood Numbers
pd.model$nicotine_blood_number <- NA
pd.model$nicotine_blood_number <- ifelse(grepl('ng/mL', pd.model$nicotine_comments), pd.model$nicotine_comments, NA)
pd.model$nicotine_blood_number <- gsub('cotinine.*', '', pd.model$nicotine_blood_number, ignore.case = TRUE) #removing all text after ng/mL
pd.model$nicotine_blood_number <- gsub('.*Nicotine', '', pd.model$nicotine_blood_number, ignore.case = TRUE) #removing all text before Nicotine
pd.model$nicotine_blood_number <- gsub("[^0-9|.]", "", pd.model$nicotine_blood_number)
pd.model$nicotine_blood_number <- as.numeric(pd.model$nicotine_blood_number)

###### Cotinine Blood Concentrations
pd.model$cotinine_blood_number <- NA
pd.model$cotinine_blood_number <- ifelse(grepl('ng/mL', pd.model$nicotine_comments), pd.model$nicotine_comments, NA)
pd.model$cotinine_blood_number <- gsub('.*cotinine', '', pd.model$cotinine_blood_number, ignore.case = TRUE) #removing all text before ; (usually nicotine said first)
pd.model$cotinine_blood_number <- gsub("[^0-9|.]", "", pd.model$cotinine_blood_number)
pd.model$cotinine_blood_number <- as.numeric(pd.model$cotinine_blood_number)

####plotting
smoke_conc_melted <- melt(pd.model[,c("BrNum","smoking","nicotine_brain_number","cotinine_brain_number","nicotine_blood_number","cotinine_blood_number")],id.vars=c("smoking","BrNum"))
smoke_conc_melted$Tox_Source <- ifelse(grepl("brain",smoke_conc_melted$variable),"Brain", "Blood")
smoke_conc_melted$Toxicant <- ifelse(grepl("nicotine",smoke_conc_melted$variable),"Nicotine", "Cotinine")

fun_length <- function(x){
  return(data.frame(y=quantile(x,na.rm=TRUE)[1],label= paste0("n=", length(x))))
}

smoke_conc_brain <- ggplot (data = subset(smoke_conc_melted,Tox_Source=="Brain"), aes(x = Toxicant, y = value))
smoke_conc_brain <- smoke_conc_brain + geom_boxplot(outlier.colour = NA, alpha = 0.1) + 
		geom_point(alpha = 1, position = position_jitter(width = 0.2)) +
		theme_bw(base_size = 18) +
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
			  legend.position="none") + 
			  scale_fill_grey() +scale_colour_grey()+
		labs(x = 'Biomarker', y="Concentration in Brain (ng/g)")

smoke_conc_blood <- ggplot (data = subset(smoke_conc_melted,Tox_Source=="Blood"), aes(x = Toxicant, y = value))
smoke_conc_blood <- smoke_conc_blood + geom_boxplot(outlier.colour = NA, alpha = 0.1) + 
		geom_point(alpha = 1, position = position_jitter(width = 0.2)) +
		theme_bw(base_size = 18) +
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
			  legend.position="none") +
		scale_fill_grey() + scale_colour_grey()+
		labs(x = 'Biomarker', y= "Concentration") +
		labs(y="Concentration in Blood (ng/mL)")  +
	geom_hline(yintercept = 200, linetype=3) + 
    annotate("text", 2.3, 216, label = "200ng/mL", size=6)
				   			   
###### stats
cor.test(pd.model$pmi, pd.model$cotinine_blood_number, use="complete.obs")
cor.test(pd.model$RIN, pd.model$cotinine_blood_number, use="complete.obs")
cor.test(pd.model$nicotine_blood_number, pd.model$cotinine_blood_number, use="complete.obs")

summary(lm(cotinine_blood_number~nicotine_blood_number+ph,data = pd.model))

cor.test(pd.model$ph, pd.model$nicotine_blood_number / pd.model$cotinine_blood_number, use="complete.obs")
cor.test(pd.model$pmi, pd.model$nicotine_blood_number / pd.model$cotinine_blood_number, use="complete.obs")
cor.test(pd.model$totalAssignedGene, pd.model$nicotine_blood_number / pd.model$cotinine_blood_number, use="complete.obs")
cor.test(pd.model$mitoRate,pd.model$nicotine_blood_number / pd.model$cotinine_blood_number, use="complete.obs")

smoke_conc_casted <- dcast(smoke_conc_melted[,!colnames(smoke_conc_melted)%in%"variable"],  BrNum+smoking +Tox_Source~Toxicant,value.var="value")
smoke_conc_casted$smoking[is.na(smoke_conc_casted$smoking)]<- "Unknown"
smoke_conc_casted$smoking <- plyr::revalue(as.factor(smoke_conc_casted$smoking), c("TRUE"="History", "FALSE"="No History") )
nic_cot_blood = ggplot(data=subset(smoke_conc_casted,!is.na(Cotinine) & !is.na(Nicotine) ), aes(x=Nicotine, y=Cotinine))
nic_cot_blood = nic_cot_blood + geom_point() +
							   theme_bw(base_size = 18) +
							   facet_wrap(~Tox_Source,scales="free")+
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
		legend.justification = c(1, 0), 
		legend.position = c(1, 0),
		legend.background = element_rect(fill=NA, size=.5, linetype="solid")) + 
		scale_colour_brewer(palette = "Set1") +		
		scale_fill_brewer(palette = "Set1")  +
			  labs(x="Nicotine Concentration",
				   y = "Cotinine Concentration")	
pdf(file="Nicotine_Cotinine_Distributions_QC.pdf",height=8.5,width=11)
smoke_conc_blood
smoke_conc_brain
nic_cot_blood
dev.off()
library(gridExtra)		
smoke_conc_blood <- smoke_conc_blood + theme_bw(base_size = 8) + theme(plot.title = element_text(hjust = 0.5),
														panel.grid.major = element_blank(),
														panel.grid.minor = element_blank(),
														legend.position="none") +
														ggtitle("Blood")
smoke_conc_brain <- smoke_conc_brain + theme_bw(base_size = 8) + theme(plot.title = element_text(hjust = 0.5),
														panel.grid.major = element_blank(),
														panel.grid.minor = element_blank(),
														legend.position="none") +
														ggtitle("Brain")

adult_smoking_conc_plot <- arrangeGrob(smoke_conc_blood, smoke_conc_brain, ncol=2,nrow=1)
		ggsave(adult_smoking_conc_plot, file = "Adult Smoking Concentrations.png",
		   type = "cairo-png",
           height=5.5,width=9,dpi=600,units="in")				   				
ggsave(nic_cot_blood, file = "Nicotine Cotinine Scatter.png",
		   type = "cairo-png",
           height=5.5,width=9,dpi=600,units="in")				   
######### Tableone
setwd('/users/ssemick/Tox_Expr/Adult_Smoking/RNA_Seq/csvs/Adult')
library(tableone)
pd.model$source <- droplevels(pd.model$source)
adult_smoking <- CreateTableOne(vars = c("Age","Race","Sex","source", "RIN","ph","PMI","mitoRate", "smoking", "cotinine_blood_number", "nicotine_blood_number"), strata = "modelGroup" , data = pd.model)
table1 <- print(adult_smoking, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
## Save to a CSV file
write.csv(table1, file = "Adult_Smoking_TableI_NoOutliersDropped.csv")

###############
setwd('/users/ssemick/Tox_Expr/Adult_Smoking/RNA_Seq/plots/Adult')
pd.model <- pd.model[!pd.model$BrNum %in% c('Br1179','Br1105','Br2267') ,]
######### Age
t.test(Age~modelGroup, data = pd.model) #NS
a		   
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
		scale_fill_brewer(palette = "Set1") + labs(x="Smoking Status", y = "RIN")
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
		scale_fill_brewer(palette = "Set1") + labs(x="Smoking Status", y="PMI")
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
		scale_fill_brewer(palette = "Set1") + labs(x="Smoking Status", y= "Mitochondrial Mapping Rate")
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
		scale_fill_brewer(palette = "Set1") + labs(x="Smoking Status", y="Gene Assignment Rate")
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
		scale_fill_brewer(palette = "Set1") + labs(x="Smoking Status", y = "Mapping Rate")		   
########################## Ancestry PCs
t.test(snpPC1~modelGroup, data = pd.model) #Ns .09
t.test(snpPC2~modelGroup, data = pd.model) #Ns
t.test(snpPC3~modelGroup, data = pd.model) # .06 barely ns
t.test(snpPC4~modelGroup, data = pd.model) #Ns
t.test(snpPC5~modelGroup, data = pd.model) #Ns
###############################
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
ggsave( Demo_EDA, file="Exploratory_Data_Analysis_Confounds_v_Smoking.pdf",
           height=8.5,width=11)									
##############################
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
ggsave(rpkm_density, file = "Adult RPKM Density.pdf",
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
#dir.create('Gene_PCA')
#setwd('Gene_PCA')

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
PC1_2 <- ggplot(data = dat1, aes(x = PC1, y =PC2) )
PC1_2 <- PC1_2 + geom_point(aes(colour=modelGroup), alpha = 1) +
				 labs(x = paste0("PC1: ", getPcaVars(pca)[1], "% Var Expl"),
					  y = paste0("PC2: ", getPcaVars(pca)[2], "% Var Expl"),
					  colour = "Smoking") +
				scale_colour_brewer(palette="Set1") + theme(legend.position='none',legend.background = element_rect(colour='black',fill=NA, size=.5, linetype="solid")) 
#ggsave(PC1_2,file = "Adult_PC1_2 Smoking Highlighted PFC.pdf",
#           height=8.5,width=11)									
						

PC2_3 <- ggplot(data = dat1, aes(x = PC2, y =PC3) )
PC2_3 <- PC2_3 + geom_point(aes(colour=modelGroup), alpha = 1.0) +
				 labs(x = paste0("PC2: ", getPcaVars(pca)[2], "% Var Expl"),
					  y = paste0("PC3: ", getPcaVars(pca)[3], "% Var Expl"),
					  colour = "Smoking") +
				scale_colour_brewer(palette="Set1") + theme(legend.position='none',
		legend.background = element_rect(colour='black',fill=NA, size=.5, linetype="solid")) 
#ggsave(PC2_3,file = "Adult_PC2_3 Smoking Highlighted PFC.pdf",
#           height=8.5,width=11)									
		   
PC3_4 <- ggplot(data = dat1, aes(x = PC3, y =PC4) )
PC3_4 <- PC3_4 + geom_point(aes(colour=modelGroup), alpha = 1.0) +
				 labs(x = paste0("PC3: ", getPcaVars(pca)[3], "% Var Expl"),
					  y = paste0("PC4: ", getPcaVars(pca)[4], "% Var Expl"),
					  colour = "Smoking") +
					  scale_colour_brewer(palette="Set1")  + theme(legend.position='none',
		legend.background = element_rect(colour='black',fill=NA, size=.5, linetype="solid")) 
#ggsave(PC3_4, file = "Adult_PC3_4 Smoking Highlighted PFC.pdf",
#           height=8.5,width=11)									
								

PC4_5 <- ggplot(data = dat1, aes(x = PC4, y =PC5) )
PC4_5 <- PC4_5 + geom_point(aes(colour=modelGroup), alpha = 1.0) +
				 labs(x = paste0("PC4: ", getPcaVars(pca)[4], "% Var Expl"),
					  y = paste0("PC5: ", getPcaVars(pca)[5], "% Var Expl"),
					  colour = "Smoking") +
					  scale_colour_brewer(palette="Set1") + theme(legend.position='bottom') 
#ggsave(PC4_5, file = "Adult_PC4_5 Smoking Highlighted PFC.pdf",
#           height=8.5,width=11)									
								
### 
library(gridExtra)
library(grid)			   
exprs_pca_boxplot <- arrangeGrob(PC1_2, PC2_3, PC3_4,PC4_5, ncol=2,nrow=2) #generates g
ggsave(exprs_pca_boxplot, file="/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/Paper/pcaExprs_boxplot_panel_adult.pdf", height=11,width=8.5) 	 
								
								
########################## Supplemental Figure 3
library(gridExtra)
PC1_2 <- PC1_2 + theme_bw(base_size = 8)
PC2_3 <- PC2_3 + theme_bw(base_size = 8)
PC3_4 <- PC3_4 + theme_bw(base_size = 8)
PC3_4 <- PC3_4 + theme_bw(base_size = 8)
PC4_5 <- PC4_5 + theme_bw(base_size = 8)
top5_PCs <- arrangeGrob(PC1_2, PC2_3, PC3_4, PC4_5, ncol=2,nrow=2, top="No Apparent Gene Expression Outliers (Adults)")
ggsave( top5_PCs, file="top5PCs_Adult.png",
		   type = "cairo-png",
           height=5.5,width=9,dpi=600,units="in")									
		   
##########################
#Check PC1
PC = "PC1"
t.test(dat1[,PC] ~ dat1$modelGroup)
t.test(dat1[,PC] ~ as.factor(dat1$Sex))
t.test(dat1[,PC] ~ as.factor(dat1$Race))
cor.test(dat1[,PC], dat1$Age) 
cor.test(dat1[,PC], dat1$RIN) #sig
cor.test(dat1[,PC], dat1$mitoRate) #sig
cor.test(dat1[,PC], dat1$mitoMapped) #sig
cor.test(dat1[,PC], dat1$mappingRate)
cor.test(dat1[,PC], dat1$totalAssignedGene) #sig
cor.test(dat1[,PC], dat1$totalMapped)
cor.test(dat1[,PC], dat1$PMI)

PC1_RIN <- ggplot(data = dat1, aes(x=PC1, y = RIN) )
PC1_RIN <- PC1_RIN + geom_point() + 
				 geom_smooth(method="lm", se =FALSE) + 
				 scale_colour_brewer(palette="Set1") +
				 theme_bw(base_size = 18) + 
				 labs(x = paste0("PC1: ", getPcaVars(pca)[1], "% Var Expl"),
					  y = "RIN",
					  title = paste0("PC1: ", "; P = ",
					  signif(cor.test(dat1[,PC], dat1$RIN)$p.value,3) ))
ggsave(PC1_RIN, file = "PC1 Associates with RIN.png",
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
cor.test(dat1[,PC], dat1$mappingRate) ##
cor.test(dat1[,PC], dat1$totalAssignedGene)
cor.test(dat1[,PC], dat1$totalMapped)
cor.test(dat1[,PC], dat1$PMI)

PC2_mitoRate <- ggplot(data = dat1, aes(x=PC2, y = mitoRate) )
PC2_mitoRate <- PC2_mitoRate + geom_point() + 
				 scale_colour_brewer(palette="Set1") +
				 theme_bw(base_size = 18) + 
				 geom_smooth(method="lm", se =FALSE) + 
				 labs(x = paste0("PC2: ", getPcaVars(pca)[2], "% Var Expl"),
					  y = "Mito Mapping Rate",
					  title = paste0("PC2: mitoRate", "; P = ",
					  signif(cor.test(dat1[,PC], dat1$mitoRate)$p.value,3) ))
ggsave(PC2_mitoRate, file = "PC2 Associates with mitoRate.png",
		   type = "cairo-png",
           height=5.5,width=9,dpi=600,units="in")									
####################################
#Check PC3
PC = "PC3"
t.test(dat1[,PC] ~ dat1$modelGroup) #marginally signficant
t.test(dat1[,PC] ~ as.factor(dat1$Sex))
t.test(dat1[,PC] ~ as.factor(dat1$Race))#marginal
cor.test(dat1[,PC], dat1$Age)
cor.test(dat1[,PC], dat1$RIN)
cor.test(dat1[,PC], dat1$mitoRate)
cor.test(dat1[,PC], dat1$mitoMapped)
cor.test(dat1[,PC], dat1$mappingRate)
cor.test(dat1[,PC], dat1$totalAssignedGene)#
cor.test(dat1[,PC], dat1$totalMapped)
cor.test(dat1[,PC], dat1$PMI)

PC3_mitoRate <- ggplot(data = dat1, aes(x=PC3, y = mitoRate) )
PC3_mitoRate <- PC3_mitoRate + geom_point() + 
				 geom_smooth(method="lm", se =FALSE) + 
				 scale_colour_brewer(palette="Set1") +
				 theme_bw(base_size = 18) + 
				 labs(x = paste0("PC3: ", getPcaVars(pca)[3], "% Var Expl"),
					  y = "Mito Mapping Rate",
					  title = paste0("PC3: mitoRate", "; P = ",
					  signif(cor.test(dat1[,PC], dat1$mitoRate)$p.value,3) ))
ggsave(PC3_mitoRate, file = "PC3 Associates with Mito Mapping Rate.png",
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

PC4_TAG <- ggplot(data = dat1, aes(x=PC4, y = totalAssignedGene) )
PC4_TAG <- PC4_TAG + geom_point() + 
				 geom_smooth(method="lm", se =FALSE) +
				 scale_colour_brewer(palette="Set1") +
				 theme_bw(base_size = 18) + 				 
				 labs(x = paste0("PC4: ", getPcaVars(pca)[4], "% Var Expl"),
					  y = "Gene Assignment",
					  title = paste0("PC4: Gene Assignment", "; P = ",
					  signif(cor.test(dat1[,PC], dat1$totalAssignedGene)$p.value,3) ))
ggsave(PC4_TAG, file = "PC4 Associates with Total Assigned Gene.png",
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

PC5_RIN <- ggplot(data = dat1, aes(x=PC5, y = RIN) )
PC5_RIN <- PC5_RIN + geom_point() + 
				 geom_smooth(method="lm", se =FALSE) +
				 scale_colour_brewer(palette="Set1") +
				 theme_bw(base_size = 18) + 				 
				 labs(x = paste0("PC5: ", getPcaVars(pca)[5], "% Var Expl"),
					  y = "RIN",
					  title = paste0("PC5: RIN", "; P = ",
					  signif(cor.test(dat1[,PC], dat1$RIN)$p.value,3) ))
ggsave(PC5_RIN, file = "PC5 Associates with RIN.png",
		   type = "cairo-png",
           height=5.5,width=9,dpi=600,units="in")					
#doesn't strongly associate with anything
library(gridExtra)
PC1_RIN <- PC1_RIN + theme_bw(base_size = 8) + ggtitle("PC1: RIN")
PC2_mitoRate <- PC2_mitoRate + theme_bw(base_size = 8) + ggtitle("PC2: mitoRate")
PC3_mitoRate <- PC3_mitoRate + theme_bw(base_size = 8) + ggtitle("PC3: mitoRate")
PC4_TAG <- PC4_TAG + theme_bw(base_size = 8) + ggtitle("PC4: Gene Assignment")
PC5_RIN <- PC5_RIN + theme_bw(base_size = 8) + ggtitle("PC5: RIN")
top5_PCs <- arrangeGrob(PC1_RIN, PC2_mitoRate, PC3_mitoRate, PC4_TAG, PC5_RIN, ncol=5,nrow=1, top="First Five Expression PCs Associate with RNA Quality (Adults)")
ggsave( top5_PCs, file="top5PCs_Adult_pd_associations.png",
		   type = "cairo-png",
           height=5.5,width=9,dpi=600,units="in")									

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
				theme(legend.position="bottom")+
				 geom_text_repel(data=dplyr::filter(pd.model, snpPC2> 0.05), aes(label=BrNum),colour="black")

snp2 <- ggplot(data = pd.model, aes(x = snpPC2, y = snpPC3, colour = Race))
snp2 <- snp2 + geom_point() +
			   labs(x = "MDS Dim 2",
					 y = "MDS Dim 3") +
				scale_colour_grey() +
				theme(legend.position="bottom")+
			   geom_text_repel(data=dplyr::filter(pd.model, snpPC2> 0.05), aes(label=BrNum),colour="black")

snp3 <- ggplot(data = pd.model, aes(x = snpPC3, y = snpPC4, colour = Race))
snp3 <- snp3 + geom_point() +
			   labs(x = "MDS Dim 3",
					 y = "MDS Dim 4") +
				scale_colour_grey() +
				theme(legend.position="bottom")+
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
pdf('MDS Genotype Exploration.pdf')
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
ggsave(snp_pca_boxplot, file="/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/Paper/pcaSnp_boxplot_panel_adult.pdf", height=11,width=8.5) 	 
		   
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
dat2 <- geneRpkm.model

geneCounts.model[sex_check[,'ensembl_gene_id'] == rownames(geneCounts.model), pd.model$Sex=="M"]


colnames(geneRpkm.model)==pd.model$RNum

dat2 = cbind(pd.model,t(geneRpkm.model[sex_check$ensembl_gene_id ,]) )
tmp <-dat2[,c('sex',sex_check$ensembl_gene_id)]
#tmp <- reshape2::melt(tmp, id.vars="sex")


sex_plot1 <- ggplot(data = tmp, aes(x= sex, y=ENSG00000067048)) 
sex_plot1 <- sex_plot1 + geom_boxplot(outlier.colour = NA, alpha = 0.1, aes(fill = sex)) + 
		geom_point(alpha = 1, aes(colour = sex), position = position_jitter(width = 0.2)) +
		theme_bw(base_size = 18) + 
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
			  legend.position="none") +
		scale_colour_brewer(palette="Set1") +
		scale_fill_brewer(palette="Set1") +
		labs( x = "Sex", y="RPKM", title = sex_check[sex_check$ensembl_gene_id=="ENSG00000067048",'hgnc_symbol'] )
ggsave(sex_plot1, file="ENSG00000067048.png",
		   type = "cairo-png",
           height=5.5,width=9,dpi=600,units="in")

sex_plot2 <- ggplot(data = tmp, aes(x= sex, y=ENSG00000198692)) 
sex_plot2 <- sex_plot2 + geom_boxplot(outlier.colour = NA, alpha = 0.1, aes(fill = sex)) + 
		geom_point(alpha = 1, aes(colour = sex), position = position_jitter(width = 0.2)) +
		theme_bw(base_size = 18) + 
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
			  legend.position="none") +
		scale_colour_brewer(palette="Set1") +
		scale_fill_brewer(palette="Set1") +
		labs( x = "Sex", y="RPKM", title = sex_check[sex_check$ensembl_gene_id=="ENSG00000198692",'hgnc_symbol'] )
ggsave(sex_plot2, file="ENSG00000198692.png",
		   type = "cairo-png",
           height=5.5,width=9,dpi=600,units="in")

sex_plot3 <- ggplot(data = tmp, aes(x= sex, y=ENSG00000129824)) 
sex_plot3 <- sex_plot3 + geom_boxplot(outlier.colour = NA, alpha = 0.1, aes(fill = sex)) + 
		geom_point(alpha = 1, aes(colour = sex), position = position_jitter(width = 0.2)) +
		theme_bw(base_size = 18) + 
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
			  legend.position="none") +
		scale_colour_brewer(palette="Set1") +
		scale_fill_brewer(palette="Set1") +
		labs( x = "Sex", y="RPKM", title = sex_check[sex_check$ensembl_gene_id=="ENSG00000129824",'hgnc_symbol'] )
ggsave(sex_plot3, file="ENSG00000129824.png",
		   type = "cairo-png",
           height=5.5,width=9,dpi=600,units="in")

sex_plot4 <- ggplot(data = tmp, aes(x= sex, y=ENSG00000229807)) 
sex_plot4 <- sex_plot4 + geom_boxplot(outlier.colour = NA, alpha = 0.1, aes(fill = sex)) + 
		geom_point(alpha = 1, aes(colour = sex), position = position_jitter(width = 0.2)) +
		theme_bw(base_size = 18) + 
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
			  legend.position="none") +
		scale_colour_brewer(palette="Set1") +
		scale_fill_brewer(palette="Set1") +
		labs( x = "Sex", y="RPKM", title = sex_check[sex_check$ensembl_gene_id=="ENSG00000229807",'hgnc_symbol'] )
ggsave(sex_plot4, file="ENSG00000229807.png",
		   type = "cairo-png",
           height=5.5,width=9,dpi=600,units="in")
		   
sex_plot5 <- ggplot(data = tmp, aes(x= sex, y=ENSG00000067646)) 
sex_plot5 <- sex_plot5 + geom_boxplot(outlier.colour = NA, alpha = 0.1, aes(fill = sex)) + 
		geom_point(alpha = 1, aes(colour = sex), position = position_jitter(width = 0.2)) +
		theme_bw(base_size = 18) + 
		theme(plot.title = element_text(hjust = 0.5),
			  panel.grid.major = element_blank(),
			  panel.grid.minor = element_blank(),
			  legend.position="none") +
		scale_colour_brewer(palette="Set1") +
		scale_fill_brewer(palette="Set1") +
		labs( x = "Sex", y="RPKM", title = sex_check[sex_check$ensembl_gene_id=="ENSG00000067646",'hgnc_symbol'] )
ggsave(sex_plot5, file="ENSG00000067646.png",
		   type = "cairo-png",
           height=5.5,width=9,dpi=600,units="in")
########################## Supplemental Figure 4
library(gridExtra)
sex_plot1 <- sex_plot1 + theme_bw(base_size = 8) + theme(plot.title = element_text(hjust = 0.5),
														panel.grid.major = element_blank(),
														panel.grid.minor = element_blank(),
														legend.position="none")  
sex_plot2 <- sex_plot2 + theme_bw(base_size = 8) + theme(plot.title = element_text(hjust = 0.5),
														panel.grid.major = element_blank(),
														panel.grid.minor = element_blank(),
														legend.position="none")
sex_plot3 <- sex_plot3 + theme_bw(base_size = 8) + theme(plot.title = element_text(hjust = 0.5),
														panel.grid.major = element_blank(),
														panel.grid.minor = element_blank(),
														legend.position="none")
sex_plot4 <- sex_plot4 + theme_bw(base_size = 8) + theme(plot.title = element_text(hjust = 0.5),
														panel.grid.major = element_blank(),
														panel.grid.minor = element_blank(),
														legend.position="none")
sex_plot5 <- sex_plot5 + theme_bw(base_size = 8) + theme(plot.title = element_text(hjust = 0.5),
														panel.grid.major = element_blank(),
														panel.grid.minor = element_blank(),
														legend.position="none")

top5_sexplots <- arrangeGrob(sex_plot4, sex_plot1, sex_plot2, sex_plot3, sex_plot5, ncol=5,nrow=1, top="No Apparent Sex Swaps (Adults)")
ggsave( top5_sexplots, file="top5_sexplots_Adult.png",
		   type = "cairo-png",
           height=5.5,width=9,dpi=600,units="in")									
		   
		   
		   