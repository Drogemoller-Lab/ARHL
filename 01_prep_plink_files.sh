#/bin/bash

DATA=/home/projects/archive/CLSA/VADA/plink_files
INPUT=/home/projects/hearing_loss/clsaARHL_SA/clsaFiles
OUTPUT=/home/projects/archive/previous_projects/plinkFiles
PLINK=/opt/plink-1.09/plink



#Filtering the formatted bfiles so that the sex and phenotype column is correct
#Then filtering by MAF 0.01 so we only have biallelic SNPS
#$PLINK --bfile $DATA/clsa_imp_v3 --remove $INPUT/comb_axiome_kinship_sex_mismatch_to_remove.txt\
#	--update-sex $INPUT/updated_sex_v2.txt --exclude $INPUT/imp_markers_to_exclude.txt\
#	--maf 0.01 --geno 0.05 --hwe 0.000001 --mind 0.05 \
#     	--make-bed --out $OUTPUT/clsa_imp_v3_clean
#
#
##create bed files for the different phenotypes
#for sample in "met_better" "met_worse" "sen_better" "sen_worse"
#do
#	$PLINK --bfile $OUTPUT/clsa_imp_v3_clean\
#		--pheno $INPUT/phen_${sample}.txt \
#		--make-bed --out $OUTPUT/${sample}
#done
#
##creat bed files for the European ancestry clusters then add the 2 phenotypes
#$PLINK --bfile $OUTPUT/clsa_imp_v3_clean --keep $INPUT/cluster4_id.txt \
#        --make-bed --out $OUTPUT/clean_cluster4
#
#
##create bed files for the different phenotypes
#for sample in "met_better" "sen_better"
#do
#        $PLINK --bfile $OUTPUT/clean_cluster4\
#                --pheno $INPUT/phen_${sample}.txt \
#                --make-bed --out $OUTPUT/${sample}_cluster4
#done


#creat bed files for the European ancestry clusters then add the 2 phenotypes
$PLINK --bfile $OUTPUT/clsa_imp_v3_clean --keep $INPUT/cluster5_id.txt \
        --make-bed --out $OUTPUT/clean_cluster5


#create bed files for the different phenotypes
for sample in "met_better" "sen_better"
do
        $PLINK --bfile $OUTPUT/clean_cluster5\
                --pheno $INPUT/phen_${sample}.txt \
                --make-bed --out $OUTPUT/${sample}_cluster5
done



