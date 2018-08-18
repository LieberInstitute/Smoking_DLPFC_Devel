Code used to generate results for the paper: *Developmental effects of maternal smoking during pregnancy on the human frontal cortex transcriptome*

## Citation
Coming soon

## Script summary

#### Exploratory data analysis
* *adult_smoking_EDA_v2.R*: Exploring RNAseq data for adult cohort.
* *fetal_smoking_EDA_v2.R*: Exploring RNAseq data for prenatal cohort.
* *summary_of_sample_characteristics_Table1.R*: Summarising demographic information for prenatal and adult cohorts (Table 1).

#### Differential expression analysis runs
* *limma_differential_expression_analyses_sva_cluster.R*:
* *limma_differential_expression_analyses_sva_cluster_interaction.R*:
* *limma_differential_expression_adult_schizo_replication.R*:
* *limma_differential_expression_adult_sensitivity_analysis.R*:
* *limma_differential_expression_analyses_sva_cluster_AA_only_prenatal_sensitivity.R*:
* *limma_differential_expression_analyses_sva_cluster_developmental_nonSmokers_reviewer01_suggestion.R*
* *limma_differential_expression_analyses_sva_cluster_secondaryAnalysis_review03Suggestion_noSVA.R*

#### Results summary

###### General output and tables
* *mainEffect_smoking_results.R*
* *interaction_results_table.R*
* *interactionEffect_smoking_results.R*
* *ordinal_sensitivity_results.R*
* *checking_SLCA12A5_SLC12A2_DE.R*

###### Visualizations
* *CEP85_junc_plots.R*
* *MARCO_gene_plots.R*
* *Prenatal_gene_plots.R*
* *interaction_genes_boxplots.R*

#### Gene ontology analyses
* *GO_adult.R*
* *GO_fetal.R*
* *gene_length_bias_enrichment_analysis.R*

#### Functional enrichment analysis
* *autism_geschwind_enrichmentTest.R*
* *autism_geschwind_tstat_cutTable.R*
* *DE_Results_Functional_Analysis_devel_sensitivity.R*
* *schizophrenia_LIBD_enrichmentTest.R*
* *visualizing_prenatal_DE_and_autism_geschwind_DE.R*