#/bin/bash

DATA=/home/projects/hearing_loss/clsaARHL_SA/plinkFiles
COVAR=/home/projects/hearing_loss/clsaARHL_SA/clsaFiles
GWAS=/home/projects/hearing_loss/clsaARHL_SA/gwas
PLINK=/opt/plink-1.09/plink


#extract chr5:73780686:C:A for the metabolic phenotype
$PLINK --bfile $DATA/met_better --extract $COVAR/chr5_73746407.txt \
	--recode --out $DATA/met_chr5_73746407


# Perform association testing
#covar file header FID IID sen.better.rank met.better.rank sex age diabetes hypertension PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 chr5_73746407

 $PLINK --bfile $DATA/met_better \
        --linear --ci 0.95 --covar $COVAR/covar_chr5_73746407.txt \
        --covar-name sen.better.rank, age, diabetes, hypertension, PC1-PC10, chr5_73746407 \
        --out $GWAS/met_better_gwas_chr5Miss_v2

