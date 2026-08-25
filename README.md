# Single Nucleotide Variant and Gene-Based Analysis Pipeline

This repository contains the scripts used for genome-wide single nucleotide variant (SNV) and gene-based analyses of metabolic and sensory age-related hearing loss (ARHL) phenotypes in the Canadian Longitudinal Study on Aging (CLSA).

## Pipeline

The main analysis workflow includes:

1. Genotype quality control and data preparation
2. Heritability estimation
3. Genome-wide association analysis
4. Genome build conversion
5. X-chromosome variant- and gene-based analyses
6. HLA analysis
7. Fine-mapping
8. Conditional and interaction analyses
9. Genetic correlation analysis
10. Candidate variant annotation
11. Statistical analyses and data visualization

# SNV and Gene-Based Analysis Pipeline

This repository contains the scripts used for genome-wide SNV and gene-based analyses of metabolic and sensory age-related hearing loss (ARHL) phenotypes in the Canadian Longitudinal Study on Aging (CLSA).

## Scripts

### 00_prepare_imp_qc.sh
Prepares the imputed genotype data for quality control and downstream analyses.

### 01a_clsaARHL_prep_plink.sh
Prepares the primary PLINK files for the metabolic and sensory phenotypes.

### 01b_prep_plink_files_cluster4.sh
Prepares ancestry-specific PLINK files for participants in CLSA PCA cluster 4.

### 02_heritability_test_greml.sh
Estimates SNV-based heritability for the metabolic and sensory phenotypes using GCTA-GREML.

### 03a_gwas.sh
Performs genome-wide SNV association analyses for the metabolic and sensory phenotypes.

### 03b_gwas_interaction.sh
Performs targeted interaction analyses for selected SNVs.

### 04_liftdown.sh
Converts SNV coordinates from GRCh38 to GRCh37/hg19 for downstream analyses.

### 05_format_xchr.bgen_to_plink.sh
Converts the imputed X-chromosome BGEN data to PLINK format.

### 06_prep_xchr_plink_files.sh
Prepares the X-chromosome PLINK files for association analysis.

### 07_liftover_xchr_bed38_to_bed19.sh
Converts X-chromosome variant coordinates from GRCh38 to GRCh37/hg19.

### 07a_format_xchr_bim_to_bed.R
Formats X-chromosome variant positions for genome build conversion.

### 08a_xwas.sh
Performs X-chromosome association analysis using XWAS.

### 08b_extract_xchr_snp_pval.sh
Extracts SNV association P values from the XWAS results.

### 08c_extract_xchr_genes.sh
Extracts genes for X-chromosome gene-based association analysis.

### 08d_xchr_gene_based_test.sh
Performs gene-based association testing for the X chromosome.

### 08e_clean_xchr_gene_results.sh
Cleans and formats the X-chromosome gene-based association results.

### 08f_run_QC_noPCA_noRelated.sh
Runs the modified XWAS quality-control procedure without PCA or relatedness filtering.

### 09a_prepare_hla_hped.R
Prepares the HLA genotype data in HPED format for HLA analysis.

### 09b_HATK_hla.sh
Performs HLA association analysis using HATK.

### 10a_finemapping.sh
Performs fine-mapping of the genome-wide association results using FINEMAP.

### 10b_polyfun.sh
Performs functionally informed fine-mapping using PolyFun.

### 10c_finemap_PIP_plot.R
Processes the fine-mapping results and generates plots of posterior inclusion probabilities.

### 10d_polyfun_effect_size.R
Compares effect sizes of fine-mapped variants between the metabolic and sensory phenotypes.

### 11_extract_snp_covar_plink.sh
Extracts selected SNVs for use as covariates in conditional association analyses.

### 12a_download_refLD.sh
Downloads the reference LD scores required for LD score regression.

### 12b_ldsc_genetic_correlation.sh
Estimates the genetic correlation between the metabolic and sensory phenotypes using LD score regression.

### 13_calc_maf.sh
Calculates minor allele frequencies for selected variants.

### 14a_cadd_select_candidate_variants.sh
Selects candidate variants for functional annotation using CADD scores.

### 14b_snps_within_deafness_genes_caddScores.R
Identifies SNVs within deafness-associated genes and processes their CADD scores.

### clinical_descriptive_stats_table.R
Generates descriptive statistics for the study cohort.

### phenotype_age_sex_analysis.R
Examines the relationships between the metabolic and sensory phenotypes, age, and sex.

### select_clinical_covars.R
Evaluates clinical variables for inclusion as covariates in the association analyses.

### miami_plots.R
Generates Miami plots for the SNV and gene-based association results.
