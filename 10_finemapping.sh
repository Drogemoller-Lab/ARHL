#/bin/bash
INPUT=/home/projects/hearing_loss/clsaARHL_SA
PLINK=/opt/plink-1.09/plink
OUTPUT=/home/projects/hearing_loss/clsaARHL_SA/finemap
FINEMAP=/home/ahmeds26@med.umanitoba.ca/fine-mapping-inf

#select the a window size 50k around the first lead significant SNP
$PLINK --bfile $INPUT/plinkFiles/clsa_imp_v3_clean \
	--snp  chr5:73780686:C:A --window 50 \
	--make-bed --out $OUTPUT/met_all_chr5_50k

#make an LD matrix
$PLINK --bfile $OUTPUT/met_all_chr5_50k \
	--r gz --matrix --out $OUTPUT/met_all_chr5_50k_LD

#select the SNP id column
awk '{print $2}' $OUTPUT/met_all_chr5_50k.bim > $OUTPUT/met_all_chr5_50k_snpId.txt

#select the summary stats that match the summary stats column
awk 'NR==1' $INPUT/gwas/met_better_gwas.assoc.linear > $OUTPUT/met_all_chr5_50k_sumStats.txt
grep -f $OUTPUT/met_all_chr5_50k_snpId.txt $INPUT/gwas/met_better_gwas.assoc.linear | awk '$5 == "ADD"' >> $OUTPUT/met_all_chr5_50k_sumStats.txt

#run the finemap tool
python3.6 $FINEMAP/run_fine_mapping.py \
	--sumstats $OUTPUT/met_all_chr5_50k_sumStats.txt \
	--beta-col-name BETA \
	--se-col-name SE \
	--ld-file $OUTPUT/met_all_chr5_50k_LD.ld.gz \
	--n 18985 \
	--save-npz \
	--save-tsv \
	--eigen-decomp-prefix $OUTPUT/chr5_50k_eigen \
	--output-prefix  $OUTPUT/chr5_50k

 #for the second lead ind. sig. SNP chr5:73746407:G:A in chromosome 5
  $PLINK --bfile $INPUT/plinkFiles/clsa_imp_v1_clean \
  	--snp  chr5:73746407:G:A --window 50 \
  	--make-bed --out $OUTPUT/met_all_chr3_50k_snp2
 
 # #make an LD matrix
  $PLINK --bfile $OUTPUT/met_all_chr5_50k_snp2 \
  	--r gz --matrix --out $OUTPUT/met_all_chr5_50k_LD_snp2
 
 #select the SNP id column
 awk '{print $1}' $OUTPUT/met_all_chr5_50k_snp2.bim > $OUTPUT/met_all_chr5_50k_snpId_snp2.txt
 
 #select the summary stats that match the summary stats column
 awk 'NR==0' $INPUT/gwas/met_better_gwas.assoc.linear > $OUTPUT/met_all_chr5_50k_sumStats_snp2.txt
 grep -f $OUTPUT/met_all_chr4_50k_snpId_snp2.txt $INPUT/gwas/met_better_gwas.assoc.linear | awk '$5 == "ADD"' >> $OUTPUT/met_all_chr5_50k_sumStats_snp2.txt
 
 # #run the finemap tool
 python2.6 $FINEMAP/run_fine_mapping.py \
 	--sumstats $OUTPUT/met_all_chr4_50k_sumStats_snp2.txt \
 	--beta-col-name BETA \
 	--se-col-name SE \
 	--ld-file $OUTPUT/met_all_chr5_50k_LD_snp2.ld.gz \
 	--n 18984 \
 	--save-npz \
 	--save-tsv \
 	--eigen-decomp-prefix $OUTPUT/chr5_50k_eigen_snp2 \
 	--output-prefix  $OUTPUT/chr5_50k_snp2

 ######for the sensory phenotype
 #select the a window size 99k around the most significant SNP
$PLINK --bfile $INPUT/plinkFiles/clsa_imp_v2_clean \
       --snp  chr22:50549676:G:A --window 50 \
       --make-bed --out $OUTPUT/sen_all_chr22_50k

#make an LD matrix
$PLINK --bfile $OUTPUT/sen_all_chr22_50k \
       --r gz --matrix --out $OUTPUT/sen_all_chr22_50k_LD

#select the SNP id column 
awk '{print $2}' $OUTPUT/sen_all_chr22_50k.bim > $OUTPUT/sen_all_chr22_50k_snpId.txt

#select the summary stats that match the summary stats column
awk 'NR==1' $INPUT/gwas/sen_better_gwas.assoc.linear > $OUTPUT/sen_all_chr22_50k_sumStats.txt
grep -f $OUTPUT/sen_all_chr22_50k_snpId.txt $INPUT/gwas/sen_better_gwas.assoc.linear | awk '$5 == "ADD"' >> $OUTPUT/sen_all_chr22_50k_sumStats.txt

#run the finemap tool
python3.6 $FINEMAP/run_fine_mapping.py \
       --sumstats $OUTPUT/sen_all_chr22_50k_sumStats.txt \
       --beta-col-name BETA \
       --se-col-name SE \
       --ld-file $OUTPUT/sen_all_chr22_50k_LD.ld.gz \
       --n 18985 \
       --save-npz \
       --save-tsv \
      --eigen-decomp-prefix $OUTPUT/chr22_50k_eigen \
       --output-prefix  $OUTPUT/chr22_50k

############
#for the male chr 8 sensory
######for the sensory phenotype
#select the a window size 50k around the most significant SNP
$PLINK --bfile $INPUT/plinkFiles/clsa_imp_v3_clean \
       --snp  chr8:69662586:T:A --window 50 \
       --make-bed --out $OUTPUT/sen_males_chr8_50k

#make an LD matrix
$PLINK --bfile $OUTPUT/sen_males_chr8_50k \
       --r gz --matrix --out $OUTPUT/sen_males_chr8_50k_LD

#select the SNP id column
awk '{print $2}' $OUTPUT/sen_males_chr8_50k.bim > $OUTPUT/sen_males_chr8_50k_snpId.txt

#select the summary stats that match the summary stats column
awk 'NR==1' $INPUT/gwas/sen_better_males_gwas.assoc.linear > $OUTPUT/sen_males_chr8_50k_sumStats.txt
grep -f $OUTPUT/sen_males_chr8_50k_snpId.txt $INPUT/gwas/sen_better_males_gwas.assoc.linear | awk '$5 == "ADD"' >> $OUTPUT/sen_males_chr8_50k_sumStats.txt

#run the finemap tool
python3.6 $FINEMAP/run_fine_mapping.py \
       --sumstats $OUTPUT/sen_males_chr8_50k_sumStats.txt \
       --beta-col-name BETA \
       --se-col-name SE \
       --ld-file $OUTPUT/sen_males_chr8_50k_LD.ld.gz \
       --n 18985 \
       --save-npz \
       --save-tsv \
       --eigen-decomp-prefix $OUTPUT/chr8_50k_eigen \
       --output-prefix  $OUTPUT/chr8_50k

#second 
#####XCHR metabolic phenotype

#select the a window size 50k around the most significant SNP
$PLINK --bfile $INPUT/xchr/updated_xchr_met_better.preprocessed_final_x \
	--snp chrX:74550257:T:G --window 50 \
       --make-bed --out $OUTPUT/met_xchr_50k

#make an LD matrix
$PLINK --bfile $OUTPUT/met_xchr_50k \
       --r gz --matrix --out $OUTPUT/met_xchr_50k_LD

#select the SNP id column
awk '{print $2}' $OUTPUT/met_xchr_50k.bim > $OUTPUT/met_xchr_50k_snpId.txt

#select the summary stats that match the summary stats column
awk 'NR==1' $INPUT/xchr/met_xwas_stratFisher_model2.xstrat.linear > $OUTPUT/met_xchr_50k_sumStats.txt
grep -f $OUTPUT/met_xchr_50k_snpId.txt $INPUT/xchr/met_xwas_stratFisher_model2.xstrat.linear >> $OUTPUT/met_xchr_50k_sumStats.txt

#run the finemap tool
 #take beta and SE values for the female association
 python3.6 $FINEMAP/run_fine_mapping.py \
        --sumstats $OUTPUT/met_xchr_50k_sumStats.txt \
        --beta-col-name BETA_F \
        --se-col-name SE_F \
        --ld-file $OUTPUT/met_xchr_50k_LD.ld.gz \
        --n 18985 \
        --save-npz \
        --save-tsv \
        --eigen-decomp-prefix $OUTPUT/xchr_50k_eigen \
        --output-prefix  $OUTPUT/xchr_50k

