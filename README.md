Code used to generate results for the paper: *Developmental effects of maternal smoking during pregnancy on the human frontal cortex transcriptome*

## Citation
Coming soon

## Brief analysis summary
We used a large cohort of prenatal (N=33) and adult (N=207) postmortem dorsolateral prefrontal cortex (DLPFC) samples to identify genes associated with exposure of smoking. These genes differed between prenatal and adult life, suggesting smoking exerts distinct transcriptomic effects dependent upon the developmental stage of the human cortex. Furthermore, these genes affected by smoking during prenatal life and in a developmentally-dependent fashion were enriched for the gene signatures of neuropsychiatric disease, especially autism spectrum disorder. Overall, these results support the possibility of maternal smoking during pregnancy as conferring molecular changes that may increase risk for neuropsychiatric disease later in the offspring's life.

## Script summary

#### Exploratory data analysis
* *adult_smoking_EDA_v2.R*: Exploring RNAseq data for adult cohort.
* *fetal_smoking_EDA_v2.R*: Exploring RNAseq data for prenatal cohort.
* *summary_of_sample_characteristics_Table1.R*: Summarising demographic information for prenatal and adult cohorts (Table 1).

#### Differential expression analysis runs
* *limma_differential_expression_analyses_sva_cluster.R*: Testing for differential expression between smoking exposed and unexposed prenatal samples as well between adult non-smokers and smokers.
* *limma_differential_expression_analyses_sva_cluster_interaction.R*: Testing for development-dependent effects of smoke exposure.
* *limma_differential_expression_adult_schizo_replication.R*: Testing for replication of adult smoking effects using a seperate cohort of patients diagnosed with schizophrenia.
* *limma_differential_expression_adult_sensitivity_analysis.R*: Ordinal sensitivity model for smoking effects in adults.
* *limma_differential_expression_analyses_sva_cluster_AA_only_prenatal_sensitivity.R*: Sensitivity analysis of prenatal differential expression analysis using only African-american samples (no Caucasian samples, N=3).
* *limma_differential_expression_analyses_sva_cluster_developmental_nonSmokers_reviewer01_suggestion.R*: Running differential expression analysis to delineate developmental effects (no smoke-exposure samples only).
* *limma_differential_expression_analyses_sva_cluster_secondaryAnalysis_review03Suggestion_noSVA.R*: Running differential expression analyses without surrogate variable analysis (SVA) as a sensitivity model.

#### Results summary

###### General output and tables
* *mainEffect_smoking_results.R*: Summarising effects of smoking exposure on gene in expression in prenatal and adult -only cohorts.
* *interaction_results_table.R*: Summarising development-dependent smoking exposure effects into a table. 
* *interactionEffect_smoking_results.R*: Summarising development-dependent smoking exposure effects.
* *ordinal_sensitivity_results.R*: Asessing the effect of smoking in adults using a smoking-intensity ordinal sensitivity model.
* *checking_SLCA12A5_SLC12A2_DE.R*: Post-hoc analyses of genes SLC12A5 and SLC12A2.

###### Visualizations
* *CEP85_junc_plots.R*: Visualizing significantly differentially expressed *CEP85* exon-exon junctions.
* *MARCO_gene_plots.R*: Visualizing differential expression of *MARCO*.
* *Prenatal_gene_plots.R*: Boxplots of significant differentially expressed genes in the prenatal cohort.
* *interaction_genes_boxplots.R*: Boxplots of significant development-dependent smoking-exposure genes.

#### Gene ontology analyses
* *GO_adult.R*: Gene ontology analysis of differentially expressed genes in adults.
* *GO_fetal.R*: Gene ontology od differentially expressed genes in prenatal samples.
* *gene_length_bias_enrichment_analysis.R*: Assessing whether gene length associates with differential gene expression and hence influences gene ontology results.

#### Functional enrichment analysis
* *autism_geschwind_enrichmentTest.R*: Testing for enrichment of differentially expressed genes for ASD-associated gene expression signatures. 
* *autism_geschwind_tstat_cutTable.R*: Pulling ASD statistics for differentially expressed genes.
* *DE_Results_Functional_Analysis_devel_sensitivity.R*: Assessing developmentally regulated genes.
* *schizophrenia_LIBD_enrichmentTest.R*: Testing for enrichment of differentially expressed genes for schizophrenia-associated gene expression signatures.
* *visualizing_prenatal_DE_and_autism_geschwind_DE.R*: Visualizations of genes associated with both smoking exposure and ASD.