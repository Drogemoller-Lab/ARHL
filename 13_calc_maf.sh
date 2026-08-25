OUTPUT=/home/projects/hearing_loss/clsaARHL_SA/clsaFiles
INPUT=/home/projects/archive/previous_projects/
PLINK=/opt/plink-1.09/plink


#calcualte maf
#$PLINK --bfile $INPUT/plinkFiles/clsa_imp_v3_clean --freq --out $OUTPUT/clsa_imp_v3_clean_maf

#for the X CHR
$PLINK --bfile $INPUT/xchr/clsa_imp_x_v3_clean --freq --out $OUTPUT/clsa_imp_x_v3_clean_maf


