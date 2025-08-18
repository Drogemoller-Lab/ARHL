#/bin/bash

DATA=/home/projects/archive/previous_projects/plinkFiles/
COVAR=/home/projects/hearing_loss/clsaARHL_SA/clsaFiles
PLINK=/opt/plink-1.09/plink
GWAS=/home/projects/hearing_loss/clsaARHL_SA/gwas
# Perform association testing
#covar files header 
#FID IID sen.better.rank met.better.rank sex age diabetes hypertension PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
#FID IID sen.better.rank met.better.rank sex age diabetes hypertension ePC1 ePC2 ePC3 ePC4 ePC5 ePC6 ePC7 ePC8 ePC9 ePC10

#$PLINK --bfile $DATA/met_better \
#       --linear --ci 0.95 --covar $COVAR/covar_all_v2.txt \
#       --covar-name sen.better.rank, age, diabetes-PC10\
#       --out $GWAS/met_better_gwas
#
#$PLINK --bfile $DATA/sen_better \
#	--linear  --ci 0.95 --covar $COVAR/covar_all_v2.txt \
#	--covar-name met.better.rank-PC10\
#	--out $GWAS/sen_better_gwas
#
#for cluster 4 (european)
$PLINK --bfile $DATA/met_better_cluster4 \
       --linear --ci 0.95 --covar $COVAR/covar_euro.txt \
       --covar-name age-ePC10 \
       --out $GWAS/met_better_cluster4_gwas_noSen

$PLINK --bfile $DATA/sen_better_cluster4 \
       --linear  --ci 0.95 --covar $COVAR/covar_euro.txt \
       --covar-name sex-ePC10\
       --out $GWAS/sen_better_cluster4_gwas_NoMet


#for cluster 5
#for cluster 4 (european)
#$PLINK --bfile $DATA/met_better_cluster5 \
#       --linear --ci 0.95 --covar $COVAR/covar_all_v2.txt \
#       --covar-name sen.better.rank, age, diabetes-PC10 \
#       --out $GWAS/met_better_cluster5_gwas
#
#$PLINK --bfile $DATA/sen_better_cluster5 \
#       --linear  --ci 0.95 --covar $COVAR/covar_all_v2.txt \
#       --covar-name met.better.rank-PC10\
#       --out $GWAS/sen_better_cluster5_gwas
#
#sex-stratified
#metabolic
#
#for sex in "males" "females"
#do
#	$PLINK --bfile $DATA/met_better --filter-${sex} \
#		--linear --ci 0.95 --covar $COVAR/covar_all_v2.txt \
#       		--covar-name sen.better.rank, age, diabetes, hypertension, PC1-PC10\
#       		--out $GWAS/met_better_${sex}_gwas
#done
#
##sensory
#for sex in "males" "females"
#do
#	$PLINK --bfile $DATA/sen_better --filter-${sex} \
#		--linear  --ci 0.95 --covar $COVAR/covar_all_v2.txt \
#       		--covar-name met.better.rank, age, diabetes, hypertension, PC1-PC10\
#       		--out $GWAS/sen_better_${sex}_gwas
#done
#
#
##sort the data 
for phenotype in "sen_better_cluster4_gwas_NoMet" "met_better_cluster4_gwas_noSen" 
do
	sort -g -k12 $GWAS/${phenotype}.assoc.linear | awk '$12 != "NA"'|awk '$12 != 0'|awk 'NR==1 || $5 == "ADD"' > $GWAS/${phenotype}_sorted_noNA_add
done


