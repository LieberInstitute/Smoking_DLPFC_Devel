#Interaction Effect Results
setwd('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq')
load('rdas/Genes_DE_sva.rda')
load('rdas/Exons_DE_sva.rda')
load('rdas/Jxns_DE_sva.rda')
load('rdas/ERs_DE_sva.rda')
load('rdas/interactionModel_allFeatures_DE_sva.rda')
theme_set(theme_bw(base_size=14) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"				 )) 

#### Volcano plot
############################### Facet Plot
g_results = dplyr::mutate(as.data.frame(intGene), sig=ifelse(intGene$adj.P.Val<0.10, "FDR<0.10", "FDR>0.10"))
g_results$HGNC_Symbol <- g_results$Symbol
g_results <- g_results[,c("HGNC_Symbol","P.Value","logFC","sig","adj.P.Val","t")]
g_results$Type = "Gene"

e_results = dplyr::mutate(as.data.frame(intExon), sig=ifelse(intExon$adj.P.Val<0.10, "FDR<0.10", "FDR>0.10"))
e_results$HGNC_Symbol <- e_results$Symbol
e_results <- e_results[,c("HGNC_Symbol","P.Value","logFC","sig","adj.P.Val","t")]
e_results$Type = "Exon"

j_results = dplyr::mutate(as.data.frame(intJxn), sig=ifelse(intJxn$adj.P.Val<0.10, "FDR<0.10", "FDR>0.10"))
j_results$HGNC_Symbol <- j_results$newGeneSymbol
j_results <- j_results[,c("HGNC_Symbol","P.Value","logFC","sig","adj.P.Val","t")]
j_results$Type = "Exon-Exon Junction"

er_results = dplyr::mutate(as.data.frame(intER), sig=ifelse(intER$adj.P.Val<0.10, "FDR<0.10", "FDR>0.10"))
er_results$HGNC_Symbol <- er_results$nearestSymbol
er_results <- er_results[,c("HGNC_Symbol","P.Value","logFC","sig","adj.P.Val","t")]
er_results$Type = "Expressed Region"

all_results <- rbind(g_results, e_results, j_results, er_results)
library(dplyr)
all_results$Type <- factor(all_results$Type, levels = c("Gene","Exon", "Exon-Exon Junction", "Expressed Region") ) 

all_volc <- ggplot(all_results, aes(logFC, -log10(P.Value))) +
  geom_point(aes(col=sig)) +
  facet_wrap(~ Type, ncol= 4, nrow=1) +
  theme(legend.justification = c(1, 1), 
		legend.position = c(1, 1),
		legend.background = element_rect(fill=NA, size=.5, linetype="solid")) +
  labs(x = expression(log[2]~"fold change"),
	   y = expression(-log[10]~"(pvalue)")) + 
  scale_color_manual(values=c("red", "grey")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, round(max(-log10(all_results$P.Value))+.5) ) ) +
  theme(legend.background = element_rect(colour = "black"),
        legend.title=element_blank()) #no legend title 
#all_volc <- all_volc + geom_text_repel(data=dplyr::filter(all_results, padj<0.05), aes(label=HGNC_Symbol))
  
#er_volc <- er_volc + geom_text_repel(data=dplyr::filter(results, padj<0.05), aes(label=nearestSymbol))
ggsave(all_volc, file = "/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/Paper/Volcano Plots for Interaction All Features.pdf",
           height=8.5,width=11,units="in")		   
		   
############## Venn Diagram with EntrezIDs ##############
library(VennDiagram)
load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda') #rpkm data
load('/users/ssemick/Tox_Expr/Adult_Smoking/RNA_Seq/rdas/Adult/ER Annotation.rda')
ann$RowName <- rownames(ER)
intER$EntrezID <- as.character(ann[match(rownames(intER),ann$RowName) ,'Entrez'])
intER$EnsemblGeneID <- rownames(geneMap)[match(intER$EntrezID, geneMap$EntrezID)]
intJxn$EntrezID <- geneMap[intJxn[,'newGeneID'], c("EntrezID")]

### Gene lists for each feature type (entrezID)
intGene_DE = unique(intGene[intGene$adj.P.Val<0.1, 'EntrezID' ])
intGene_DE = intGene_DE[!is.na(intGene_DE)]

intExon_DE = unique(intExon[intExon$adj.P.Val<0.1, 'EntrezID' ])
intExon_DE = intExon_DE[!is.na(intExon_DE)]

intJxn_DE = unique(intJxn[intJxn$adj.P.Val<0.1, 'EntrezID' ])
intJxn_DE = intJxn_DE[!is.na(intJxn_DE)]

intER_DE = unique(intER[intER$adj.P.Val<0.1, 'EntrezID' ])
intER_DE = intER_DE[!is.na(intER_DE)]

univ = unique(c(intER$EntrezID, intJxn$EntrezID, intExon$EntrezID, intGene$EntrezID ) )
univ = univ[!is.na(univ)]
#The goal of the Venn Diagram is to count how many words are common between SNP_pop_1 and SNP_pop_2, between SNP_pop_1 and SNP_pop_3 and so on...
#The venn.diagram function do it automatically and draw it! (you will get a png file in your current working directory)
library(VennDiagram)
DE_List = list(intGene_DE, intExon_DE, intJxn_DE, intER_DE )
##### Drawing the venn diagram   
entrez_venn = venn.diagram(x = DE_List,
							category.names = c("Genes" , "Exons" , "Junctions", "Expressed Regions"),
							filename = NULL,
							fill = c('red', 'blue', 'green', 'yellow') )
pdf(file='/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/Paper/interaction_feature_overlap_EntrezID.pdf')
    grid.draw(entrez_venn)
dev.off()

		
############## Venn Diagram with EnsemblGeneIDs ##############

### Gene lists for each feature type (ensemblGeneId)
load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda') #rpkm data
intJxn$EntrezID <- geneMap[intJxn[,'newGeneID'], c("EntrezID")]
load('/users/ssemick/Tox_Expr/Adult_Smoking/RNA_Seq/rdas/Adult/ER Annotation.rda')
ann$RowName <- rownames(ER)
intER$EntrezID <- as.character(ann[match(rownames(intER),ann$RowName) ,'Entrez'])
intER$EnsemblGeneID <- rownames(geneMap)[match(intER$EntrezID, geneMap$EntrezID)]
#> table(is.na(intER[intER$adj.P.Val<0.1, 'EnsemblGeneID' ]))
#
#FALSE  TRUE
#  553    25

intGene_DE = unique( rownames(intGene[intGene$adj.P.Val<0.1, ]) )
intGene_DE = intGene_DE[!is.na(intGene_DE)]

intExon_DE = unique(intExon[intExon$adj.P.Val<0.1, 'Geneid' ])
intExon_DE = intExon_DE[!is.na(intExon_DE)]

intJxn_DE = unique(intJxn[intJxn$adj.P.Val<0.1, 'newGeneID' ])
intJxn_DE = intJxn_DE[!is.na(intJxn_DE)]

intER_DE = unique(intER[intER$adj.P.Val<0.1, 'EnsemblGeneID' ])
intER_DE = intER_DE[!is.na(intER_DE)]
		   
DE_List = list(intGene_DE, intExon_DE, intJxn_DE, intER_DE )
##### Drawing the venn diagram   
ensembl_venn = venn.diagram(x = DE_List,
							category.names = c("Genes" , "Exons" , "Junctions", "Expressed Regions"),
							filename = NULL,
							fill = c('red', 'blue', 'green', 'yellow') )
pdf(file='/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/Paper/interaction_feature_overlap_EnsemblGeneIDs.pdf')
    grid.draw(ensembl_venn)
dev.off()