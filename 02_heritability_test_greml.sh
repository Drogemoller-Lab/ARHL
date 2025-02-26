#!/bin/bash

DATA=/home/projects/hearing_loss/clsa_v2/plink_output
DATA2=/home/projects/hearing_loss/clsa_v2/plink_files
GREML=/home/projects/hearing_loss/clsa_v2/greml
GCTA=/opt/gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static

# This script: 
# GCTA-GRM: estimating genetic relatedness from SNPs
# GCTA-GREML: Estimates variance explained by all the SNPs

for sample in "met_better" "met_worse" "sen_better" "sen_worse"
do
	$GCTA --bfile $DATA/clsa_imp_v3_${sample} --autosome\
		--make-grm --out $GREML/clsa_imp_v3_${sample} --thread-num 16
	$GCTA --grm $GREML/clsa_imp_v3_${sample} --pheno $DATA2/updated_phen_${sample}.txt \
		--reml --qcovar $DATA2/clsa_pca10.txt\
		--out $GREML/clsa_imp_v3_${sample} --thread-num 16
done


