#/bin/bash

DATA=/home/projects/archive/CLSA/VADA/plink_files
INPUT=/home/projects/hearing_loss/clsaARHL_SA/clsaFiles
OUTPUT=/home/projects/hearing_loss/clsaARHL_SA/plinkFiles
PLINK=/opt/plink-1.09/plink



#Filtering the formatted bfiles so that the sex and phenotype column is correct
#Then filtering by MAF 0.01 so we only have biallelic SNPS

$PLINK --bfile $DATA/clsa_imp_v3 --remove $INPUT/comb_axiome_kinship_sex_mismatch_to_remove.txt\
	--update-sex $INPUT/updated_sex_v2.txt --exclude $INPUT/imp_markers_to_exclude.txt\
	--maf 0.01 --geno 0.05 --hwe 0.000001 --mind 0.05 \
     	--make-bed --out $OUTPUT/clsa_imp_v3_clean


#create bed files for the different phenotypes
$PLINK --bfile $OUTPUT/clsa_imp_v3_clean\
        --pheno $INPUT/phen_met_better.txt \
        --make-bed --out $OUTPUT/met_better

$PLINK --bfile $OUTPUT/clsa_imp_v3_clean\
        --pheno $INPUT/phen_met_worse.txt \
        --make-bed --out $OUTPUT/met_worse


$PLINK --bfile $OUTPUT/clsa_imp_v3_clean\
        --pheno $INPUT/phen_sen_better.txt \
        --make-bed --out $OUTPUT/sen_better

$PLINK --bfile $OUTPUT/clsa_imp_v3_clean\
        --pheno $INPUT/phen_sen_worse.txt \
        --make-bed --out $OUTPUT/sen_worse

