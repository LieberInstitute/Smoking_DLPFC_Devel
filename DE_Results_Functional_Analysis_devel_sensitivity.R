#Results Analysis
setwd('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq')
load('rdas/developmentModel_allFeatures_DE_sva_nonSmokeExposed.rda')

TESTED_UNIV <- unique(c(mainGene$Symbol,mainExon$Symbol, mainJxn$newGeneSymbol, mainER$nearestSymbol ) )
DE <- unique(c(mainGene[mainGene$adj.P.Val<0.1,'Symbol'],mainExon[mainExon$adj.P.Val<0.1,'Symbol'], mainJxn[mainJxn$adj.P.Val<0.1,'newGeneSymbol'], mainER[mainER$adj.P.Val<0.1,'nearestSymbol'] ))
NOT_DE <- TESTED_UNIV[!TESTED_UNIV %in% DE]


autism_SFARI = read.csv('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/csvs/sets/SFARI-Gene_genes_export30-08-2017.csv')
autism_SFARI_expressed = autism_SFARI[autism_SFARI$gene.symbol %in%  mainGene$Symbol,]

autism_autdb = read.csv('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/csvs/sets/AutDb_gene-summary.csv')
autism_autdb_expressed = autism_autdb[autism_autdb$Gene.Symbol %in% mainGene$Symbol,]

############### Testing additional enrichment ###############
aej_sets = openxlsx::read.xlsx('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/csvs/sets/Supplementary Tables for paper.Birnbaum November 2013.AJP.xlsx')
aej_sets_expressed = aej_sets[aej_sets$Gene.Symbol %in% TESTED_UNIV, ]
splitSets = split(aej_sets_expressed, aej_sets_expressed$Gene.Set)
splitSets$SFARI = data.frame(Gene.Symbol = autism_SFARI_expressed$gene.symbol)
splitSets$AutDb = autism_autdb_expressed

## Devel Main
develEnrich = sapply(splitSets, function(x) {
DE_OVERLAP = c( sum(DE%in% x$Gene.Symbol),sum(!DE%in% x$Gene.Symbol))
NOT_DE_OVERLAP= c(sum(NOT_DE%in% x$Gene.Symbol), sum(!NOT_DE%in% x$Gene.Symbol) )
enrich_table = cbind(DE_OVERLAP, NOT_DE_OVERLAP)
res = fisher.test(enrich_table)
dat=c(res$p.value, 
	  res$estimate, 
	  nGeneEnriched =enrich_table[1] , 
	  nGeneSet = enrich_table[1] +enrich_table[3],  
	  percentEnrich = (enrich_table[1] / enrich_table[2])*100, 
	  backgroundEnrich = (enrich_table[3] / enrich_table[4] ) *100,
	  geneSymbolsEnriched = paste( DE[DE%in% x$Gene.Symbol] , collapse=";") )

names(dat) <- c("P.Value","Odds Ratio", "nGeneEnriched", "nGeneSet", "percentEnrich", "percentBackground","geneSymbolsEnriched")
return(dat)
})
write.csv(t(develEnrich), file='/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/csvs/developmentalTargetedEnrichmentAnalysis.csv' ,row.names=TRUE)
