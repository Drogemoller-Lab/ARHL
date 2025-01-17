#/bin/bash

DATA=/home/projects/CLSA/raw_data/genetic_data
GWAS=/home/projects/hearing_loss/clsa_v2/plink_files

# zcat to unzip the imputation files
#filter for markers with quality (Rsq) > 0.8 for both batches (column 5 and 7)
#delete the header (NR)
# combine the files
for chr in {1..22}
do
zcat $DATA/clsa_mfi_${chr}_v3.txt.gz | awk ' $5 <= 0.8 || $7 <= 0.8 {print $1}' | awk 'NR>1'  > $GWAS/clsa_imp_${chr}_v3.txt
cat $GWAS/clsa_imp_${chr}_v3.txt >> $GWAS/imp_markers_to_exclude.txt
done


