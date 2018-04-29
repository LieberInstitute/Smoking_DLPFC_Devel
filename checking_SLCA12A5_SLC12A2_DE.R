##
library(LIBDpheno)
library(ggplot2)
theme_set(theme_bw(base_size=14) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"				 ))
## Checking differential expression in prenatal PFC in smokers vs non smokers in SLC12A5, SLC12A2 
load('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/allFeatures_bound.rda')
colnames(all_bound)[8:13] <- paste0(colnames(all_bound)[8:13],"_Interaction")
all_bound=all_bound[,c("Symbol","Type","P.Value_Interaction","P.Value_Prenatal","P.Value_Adult","adj.P.Val_Interaction","adj.P.Val_Prenatal","adj.P.Val_Adult", "logFC_Prenatal","logFC_Adult")]
foi <- all_bound[all_bound$Symbol %in% c("CHRNA7","SLC12A5","SLC12A2"),]
#write.csv(foi,file='/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/csvs/Berg_gene_transcriptFeature_stats.csv', row.names=FALSE)


##
load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda') #rpkm data
load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/processed_covMat.rda')
regionMat <- round(regionMat)

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

#Get plotting data
plot_dat <- cbind(pd.model, 
				  t(geneRpkm[rownames(geneRpkm)%in%rownames(foi),pd.model$RNum]), 
				  t(exonRpkm[rownames(exonRpkm)%in%rownames(foi),pd.model$RNum]),
				  t(jRpkm[rownames(jRpkm)%in%rownames(foi),pd.model$RNum]),
				  t(regionMat[rownames(regionMat)%in%rownames(foi),pd.model$RNum]) )
#Get map data				  


##  Create interaction boxplots
geneMap$Type = "Gene"
exonMap$Type = "Exon"
jMap$Type = "Jxn"
regions$Type="ER"
map = data.frame(Rownames=c(rownames(geneMap),rownames(exonMap),names(jMap),names(regions) ), 
Symbol = c(geneMap$Symbol,exonMap$Symbol,jMap$newGeneSymbol,regions$nearestSymbol ), 
Type = c(geneMap$Type, exonMap$Type, jMap$Type, regions$Type) )

pdf('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/plots/Berg_gene_transcriptFeature_boxplots.pdf',height=11,width=8.5)
for (i in  rownames(foi[order(foi$P.Value_Interaction, decreasing=F), ]) ) {
orig_names = colnames(plot_dat)
dat=plot_dat[,c('modelGroup','AgeGroup',i)]
colnames(dat)[3] <- 'transcriptFeature'


### Change column name for ggplot2 to work

### custom_title = paste0(i, "\n p=",as.character(signif(meqtlBestCpG$pvalue[match(i, meqtlBestCpG$UniqueID)],3)) ) #custom title

a = ggplot(dat, aes(x = AgeGroup, y = log2(transcriptFeature+1), fill=modelGroup)) +
        geom_boxplot(outlier.colour = NA, alpha = 0, col='black')  + 
		geom_point(aes(col=`modelGroup`),position = position_jitterdodge(jitter.width=0.2,dodge.width=.85)) + 
		labs(y=paste0('log2(',i,'+1)'), title = paste0(map[match(i, map$Rownames),'Symbol'], ": " ,map[match(i, map$Rownames),'Type'], "\n Interaction p=", as.character(signif(foi[i,'P.Value_Interaction'],3))) )  + 
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + 
			theme(legend.position='bottom')
print(a)		
	
}
dev.off()

### Check ratio of SLC12A2/SLC12A5
fetal_plot_dat=plot_dat[plot_dat$Age <0, ]
fetal_plot_dat[,'A2/A5'] = log2(fetal_plot_dat$ENSG00000064651+1) - log2(fetal_plot_dat$ENSG00000124140+1)

a = ggplot(fetal_plot_dat, aes(x = modelGroup, y = `A2/A5`, fill=modelGroup)) +
        geom_boxplot(outlier.colour = NA, alpha = 0, col='black')  + 
		geom_point(aes(col=`modelGroup`),position = position_jitter(w = 0.1, h = 0)) + 
		labs(y="log2( SLC12A2 RPKM+1) - log2( SLC12A5 RPKM+1)")   + 
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1")

ggsave(a, filename='/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/plots/SLC12A2_by_SLC12A5_ratio.pdf')		
t.test(fetal_plot_dat$`A2/A5` ~fetal_plot_dat$modelGroup)
t.test(`A2/A5` ~modelGroup, data = fetal_plot_dat[-which.min(fetal_plot_dat$`A2/A5`),])
cor.test(fetal_plot_dat$`A2/A5`,fetal_plot_dat$Age)

summary(lm(`A2/A5` ~modelGroup +snpPC1+Age+Sex+RIN+mitoMapped+totalAssignedGene, data = fetal_plot_dat) )

summary(lm(`A2/A5` ~modelGroup +snpPC1+Age+Sex+RIN+mitoMapped+totalAssignedGene, data = fetal_plot_dat[-which.min(fetal_plot_dat$`A2/A5`),]) )

# plot of age effect for ratio
b = ggplot(fetal_plot_dat, aes(x = Age, y = `A2/A5`, fill=modelGroup)) +
		geom_point(aes(col=`modelGroup`)) + 
		labs(x='Age', y="log2( SLC12A2 RPKM+1) - log2( SLC12A5 RPKM+1)")   + 
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + theme(legend.position='bottom')
ggsave(b, filename='/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/plots/SLC12A2_by_SLC12A5_ratio_over_prenatal_age.pdf')		
		
# plot of age effect for ratio
c = ggplot(fetal_plot_dat[-which.max(fetal_plot_dat$`Age`),], aes(x = Age, y = `A2/A5`, fill=modelGroup, col=`modelGroup`)) +
		geom_point() + geom_smooth(method='lm',se=FALSE) + 
		labs(x='Age', y="log2( SLC12A2 RPKM+1) - log2( SLC12A5 RPKM+1)")   + 
		scale_colour_brewer(palette = "Set1") +
		scale_fill_brewer(palette = "Set1") + theme(legend.position='bottom')
ggsave(c, filename='/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/plots/SLC12A2_by_SLC12A5_ratio_over_prenatal_ageXsmoking.pdf')		
		
summary(lm(`A2/A5` ~modelGroup +snpPC1+Age+Sex+RIN+mitoMapped+totalAssignedGene+modelGroup*Age, data = fetal_plot_dat) )
summary(lm(`A2/A5` ~modelGroup +snpPC1+Age+Sex+RIN+mitoMapped+totalAssignedGene+modelGroup*Age, data = fetal_plot_dat[-which.min(fetal_plot_dat$`A2/A5`),]) )
