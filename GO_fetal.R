#qsub -l bluejay -l mf=30G,h_vmem=50G,h_stack=256M -o logs -e logs -cwd -b y R CMD BATCH --no-save GO_fetal.R
#allFeatures fetal smoking Gene Ontology
library('clusterProfiler')
library('org.Hs.eg.db')
library('WriteXLS')
library('ReactomePA')
library('DOSE')
library('bumphunter')
library('TxDb.Hsapiens.UCSC.hg19.knownGene')

setwd('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq')
load('rdas/Genes_DE_sva.rda')
load('rdas/Exons_DE_sva.rda')
load('rdas/Jxns_DE_sva.rda')
load('rdas/ERs_DE_sva.rda')

load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rawCounts_szControlDLPFC.rda') #load this file for the mapObjects 

#add annotation for ERs
load('/users/ssemick/Tox_Expr/Adult_Smoking/RNA_Seq/rdas/Adult/ER Annotation.rda')
ann$RowName <- rownames(ER)
fetalER$EntrezID <- as.character(ann[match(rownames(fetalER),ann$RowName) ,'Entrez'])
#genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
#ann <- matchGenes(fetalER, subject = genes)
#fetalER$EntrezID <- as.character(ann$Entrez)
#save(ER,ann, file = "/users/ssemick/Tox_Expr/Adult_Smoking/RNA_Seq/rdas/limma/ER_Annotation.rda")
# Annotation for Jxn
fetalJxn$EntrezID <- geneMap[fetalJxn[,'newGeneID'], c("EntrezID")]
##############
univ = c(fetalGene$EntrezID, fetalExon$EntrezID, fetalJxn$EntrezID, fetalER$EntrezID) #, fetalER$EntrezID
univ = as.character(unique(univ[!is.na(univ)]))
#################
##Fixing type
fetalGene$type <- "Gene"
fetalExon$type <- "Exon"
fetalJxn$type <- "Jxn"
fetalER$type <- "ER"

## up and down stats
statList = list(Gene = fetalGene, Exon = fetalExon, Junction = fetalJxn, ER = fetalER ) #
sigStatList = lapply(statList, function(x) 
  x[which(x$P.Value < 1e-3), c("logFC", "adj.P.Val", "EntrezID", "type")] )
sigStats = do.call("rbind", sigStatList)
sigStats$Sign = sign(sigStats$logFC)
sigStats = sigStats[!is.na(sigStats$EntrezID),]
sigStats = sigStats[!duplicated(sigStats[,c("EntrezID", "type")]),]
sigStats$EntrezID = as.character(sigStats$EntrezID)
## by sign
sigStats$lab =paste0(sigStats$Sign, "_", sigStats$type)
tapply(sigStats$EntrezID,sigStats$lab,length)
gList = split(sigStats$EntrezID, sigStats$lab)

gList$`-1_All` <- unique(c(gList[['-1_Exon']], gList[['-1_Gene']], gList[['-1_Jxn']], gList[['-1_ER']] ))

gList$`1_All` <- unique(c(gList[['1_Exon']], gList[['1_Gene']], gList[['1_Jxn']], gList[['1_ER']] ))
##
setwd('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq')
save(statList,sigStats,gList,univ, file = '/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/Prenatal_e3_GO_Split.rda')
##
compareKegg = compareCluster(gList, fun = "enrichKEGG",organism = "human", 
                             universe = univ, qvalueCutoff = 0.2, pvalueCutoff = 0.05)
compareGoMf = compareCluster(gList, fun = "enrichGO", universe = univ, ont = "MF",OrgDb=org.Hs.eg.db, 
                             qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
compareGoBp = compareCluster(gList, fun = "enrichGO", universe = univ, ont = "BP",OrgDb=org.Hs.eg.db, 
                             qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
compareGoCc = compareCluster(gList, fun = "enrichGO", universe = univ, ont = "CC",OrgDb=org.Hs.eg.db, 
                             qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
#compareDO   = compareCluster(gList, fun = "enrichDO", universe = univ,
#							 qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
compareReact= compareCluster(gList, fun="enrichPathway", universe = univ, 
							 qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
compareGo = lapply(list(compareGoMf=compareGoMf,
compareGoBp=compareGoBp,
compareGoCc=compareGoCc),simplify)
#fun One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway"

pdf(file = '/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/plots/Prenatal_e3_GO_Split.pdf',width = 8.5, height=11)
plot(compareKegg,colorBy="qvalue", font.size =10,title = "KEGG Terms DE Group Comparisons")
plot(compareGo[[1]],colorBy="qvalue", font.size =10, title = "GO-MF Terms DE Group Comparisons")
plot(compareGo[[2]],colorBy="qvalue", font.size =10, title = "GO-BP Terms DE Group Comparisons")
plot(compareGo[[3]],colorBy="qvalue", font.size =10, title = "GO-CC Terms DE Group Comparisons")
#plot(compareDO,colorBy="qvalue", font.size =10,title = "DO Terms DE Group Comparisons")
plot(compareReact,colorBy="qvalue", font.size =10,title = "Reactome Terms DE Group Comparisons")
dev.off()

#compareDO, not included
#, compareDO = compareDO
save(compareKegg, compareGo, compareReact, file = '/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/Prenatal_e3_Enrichment_Split.rda' )
openxlsx::write.xlsx(lapply(c(compareGo,compareKegg=compareKegg, compareReact = compareReact, compareDO=compareDO ), summary),
file = '/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/csvs/Prenatal_e3_Enrichment_Split.xlsx')

###########
#no up down
sigStats$lab = sigStats$type
gList = split(sigStats$EntrezID, sigStats$lab)
gList$`All` <- unique(c(gList[['Exon']], gList[['Gene']], gList[['Jxn']], gList[['ER']]))

## compare clusters 
compareKegg = compareCluster(gList, fun = "enrichKEGG",organism = "human", 
                             universe = univ, qvalueCutoff = 0.2, pvalueCutoff = 0.05)
compareGoMf = compareCluster(gList, fun = "enrichGO", universe = univ, ont = "MF",OrgDb=org.Hs.eg.db, 
                             qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
compareGoBp = compareCluster(gList, fun = "enrichGO", universe = univ, ont = "BP",OrgDb=org.Hs.eg.db, 
                             qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
compareGoCc = compareCluster(gList, fun = "enrichGO", universe = univ, ont = "CC",OrgDb=org.Hs.eg.db, 
                             qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
#compareDO   = compareCluster(gList, fun = "enrichDO", universe = univ,
#							 qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
#compareReact= compareCluster(gList, fun="enrichPathway", universe = univ,  qvalueCutoff = 0.2, pvalueCutoff = 0.05, readable = TRUE)
compareGo = lapply(list(compareGoMf=compareGoMf,compareGoBp=compareGoBp,compareGoCc=compareGoCc),simplify)

pdf(file = '/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/plots/Prenatal_e3_GO_Unsplit.pdf',width = 8.5, height=11)
plot(compareKegg,colorBy="qvalue", font.size =10,title = "KEGG Terms DE Group Comparisons")
plot(compareGo[[1]],colorBy="qvalue", font.size =10, title = "GO-MF Terms DE Group Comparisons")
plot(compareGo[[2]],colorBy="qvalue", font.size =10, title = "GO-BP Terms DE Group Comparisons")
plot(compareGo[[3]],colorBy="qvalue", font.size =10, title = "GO-CC Terms DE Group Comparisons")
#plot(compareDO,colorBy="qvalue", font.size =10,title = "DO Terms DE Group Comparisons")
#plot(compareReact,colorBy="qvalue", font.size =10,title = "Reactome Terms DE Group Comparisons")
dev.off()
#compareDO,
save(compareKegg, compareGo, file = '/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/rdas/Prenatal_e3_Enrichment_Unsplit.rda')
openxlsx::write.xlsx(lapply(c(compareGo,compareKegg=compareKegg ), summary),file = '/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/csvs/Prenatal_e3_Enrichment_Unsplit.xlsx')