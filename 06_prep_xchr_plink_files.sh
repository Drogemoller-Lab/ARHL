#/bin/bash

DATA=/home/projects/archive/CLSA/raw_data/genetic_data
XCHR=/home/projects/hearing_loss/clsaARHL_SA/xchr
FILES=/home/projects/hearing_loss/clsaARHL_SA/clsaFiles
PLINK=/opt/plink-1.09/plink

#imputed markers to be removed will be saved in a separate file
zcat $DATA/clsa_mfi_23_v3.txt.gz | awk ' $5 <= 0.8 || $7 <= 0.8 {print $1}' | awk 'NR>1'  > $FILES/xchr_imp_markers_to_exclude.txt

#Filtering the formatted bfiles so that the sex is correct 
#remove axiome lines and related samples

$PLINK --bfile $XCHR/clsa_imp_23_v3 --remove $FILES/comb_axiome_kinship_sex_mismatch_to_remove.txt\
        --update-sex $FILES/updated_sex_v2.txt --exclude $FILES/xchr_imp_markers_to_exclude.txt\
        --make-bed --out $XCHR/clsa_imp_x_v3_clean

#create bed files for the metabolic and sensory continuous phenotypes (20,332 samples)

for phenotype in "met_better" "sen_better"
do
$PLINK --bfile $XCHR/clsa_imp_x_v3_clean\
        --pheno $FILES/phen_${phenotype}.txt \
        --make-bed --out $XCHR/xchr_${phenotype}
done


