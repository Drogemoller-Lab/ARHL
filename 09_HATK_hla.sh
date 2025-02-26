#/bin/bash
DATA=/home/projects/hearing_loss/clsaARHL_SA/hla
PLINK=/opt/plink-1.09/plink
FILES=/home/projects/hearing_loss/clsaARHL_SA/plinkFiles
COVAR=/home/projects/hearing_loss/clsaARHL_SA/clsaFiles

#extract chromosome # 6
#metabolic
$PLINK --bfile $FILES/met_better --chr 6 \
	--make-bed --out $DATA/met_better_chr6
#sensory
$PLINK --bfile $FILES/sen_better --chr 6 \
        --make-bed --out $DATA/sen_better_chr6

# run the hatk tool in the originat directory where the HATK.py is present (use cd in the code)
#IMGT2SNP: this is to generate dictionary and file formats that can be used in the next steps. this tool uses IPD-IMGT/HLA database to generate the latest nomenclature
#this tool needs conda environment
cd /home/ahmeds26@med.umanitoba.ca/HATK
conda activate HATK

python3.8 HATK.py \
    	--imgt2seq \
    	--hg 38 \
    	--imgt 3490 \
  	--out $DATA/met.hg38.imgt3490 \
    	--imgt-dir /home/ahmeds26@med.umanitoba.ca/HATK/IMGTHLA \
    	--multiprocess 2

#this step is to update the hla nomenclature of the .hped file
#input the .hat from the first step
#output the .chped 
python3.8 HATK.py \
    --nomencleaner \
    --hat $DATA/HLA_ALLELE_TABLE.imgt3490.hat \
    --hped $DATA/met_better_hla_noDRB345.hped \
    --out $DATA/met_better_hla_noDRB345_named \
    --imgt 3490  

#convert to plink bed files
#use the .chped 

python3.8 HATK.py \
   --bmarkergenerator \
    --chped $DATA/met_better_hla_noDRB345_named.chped \
    --out $DATA/met_better_hla_noDRB345_named_hatk \
    --hg 38 \
    --dict-AA $DATA/HLA_DICTIONARY_AA.hg38.imgt3490 \
    --dict-SNPS $DATA/HLA_DICTIONARY_SNPS.hg38.imgt3490


#perform GWAS
$PLINK --bfile $DATA/met_better_hla_noDRB345_named_hatk \
	--linear --ci 0.95 --covar $COVAR/covar_all_v2.txt \
	--covar-name sen.better.rank, age, diabetes, hypertension, PC1-PC10 \
	--out $DATA/met_better_hla_gwas
#sex-stratified
$PLINK --bfile $DATA/met_better_hla_noDRB345_named_hatk --filter-males\
       --linear --ci 0.95 --covar $COVAR/covar_all_v2.txt \
        --covar-name sen.better.rank, age, diabetes, hypertension, PC1-PC10 \
        --out $DATA/met_better_hla_males_gwas

$PLINK --bfile $DATA/met_better_hla_noDRB345_named_hatk --filter-females\
        --linear --ci 0.95 --covar $COVAR/covar_all_v2.txt \
        --covar-name sen.better.rank, age, diabetes, hypertension, PC1-PC10 \
        --out $DATA/met_better_hla_females_gwas

sort -g -k12 $DATA/met_better_hla_gwas.assoc.linear | awk '$12 != "NA"'|awk '$12 != 0'|awk 'NR==1 || $5 == "ADD"' > $DATA/met_better_hla_gwas_sorted_noNA_add


sort -g -k12 $DATA/met_better_hla_males_gwas.assoc.linear | awk '$12 != "NA"'|awk '$12 != 0'|awk 'NR==1 || $5 == "ADD"' > $DATA/met_better_hla_males_gwas_sorted_noNA_add

sort -g -k12 $DATA/met_better_hla_females_gwas.assoc.linear | awk '$12 != "NA"'|awk '$12 != 0'|awk 'NR==1 || $5 == "ADD"' > $DATA/met_better_hla_females_gwas_sorted_noNA_add


######################
#for the sensory phenotype
#python3.8 HATK.py \
#       --imgt2seq \
#       --hg 38 \
#       --imgt 3490 \
#       --out $DATA/sen.hg38.imgt3490 \
#       --imgt-dir /home/ahmeds26@med.umanitoba.ca/HATK/IMGTHLA \
#       --multiprocess 2

#this step is to update the hla nomenclature of the .hped file
#input the .hat from the first step
#output the .chped
#python3.8 HATK.py \
#    --nomencleaner \
#    --hat $DATA/HLA_ALLELE_TABLE.imgt3490.hat \
#    --hped $DATA/sen_better_hla_noDRB345.hped\
#    --out $DATA/sen_better_hla_noDRB345_named \
#    --imgt 3490


#convert to plink bed files
#use the .chped
# 
# python3.8 HATK.py \
#     --bmarkergenerator \
#     --chped $DATA/sen_better_hla_noDRB345_named.chped \
#     --out $DATA/sen_better_hla_noDRB345_named_hatk \
#     --hg 38 \
#     --dict-AA $DATA/HLA_DICTIONARY_AA.hg38.imgt3490 \
#     --dict-SNPS $DATA/HLA_DICTIONARY_SNPS.hg38.imgt3490
# 
# 
# #perform GWAS
# $PLINK --bfile $DATA/sen_better_hla_noDRB345_named_hatk \
# 	--linear --ci 0.95 --covar $COVAR/covar_all_v2.txt \
# 	--covar-name met.better.rank, age, sex, diabetes, hypertension, PC1-PC10 \
# 	--out $DATA/sen_better_hla_gwas
# 
# #sex-stratified sensory
# $PLINK --bfile $DATA/sen_better_hla_noDRB345_named_hatk --filter-males \
#         --linear --ci 0.95 --covar $COVAR/covar_all_v2.txt \
#         --covar-name met.better.rank, age, diabetes, hypertension, PC1-PC10 \
#         --out $DATA/sen_better_hla_males_gwas
# $PLINK --bfile $DATA/sen_better_hla_noDRB345_named_hatk --filter-females \
#         --linear --ci 0.95 --covar $COVAR/covar_all_v2.txt \
#         --covar-name met.better.rank, age, diabetes, hypertension, PC1-PC10 \
#         --out $DATA/sen_better_hla_females_gwas
# 
# sort -g -k12 $DATA/sen_better_hla_gwas.assoc.linear | awk '$12 != "NA"'|awk '$12 != 0'|awk 'NR==1 || $5 == "ADD"' > $DATA/sen_better_hla_gwas_sorted_noNA_add
# 
# sort -g -k12 $DATA/sen_better_hla_males_gwas.assoc.linear | awk '$12 != "NA"'|awk '$12 != 0'|awk 'NR==1 || $5 == "ADD"' > $DATA/sen_better_hla_males_gwas_sorted_noNA_add
# 
# sort -g -k12 $DATA/sen_better_hla_females_gwas.assoc.linear | awk '$12 != "NA"'|awk '$12 != 0'|awk 'NR==1 || $5 == "ADD"' > $DATA/sen_better_hla_females_gwas_sorted_noNA_add
# 


