#/bin/bash
FILES=/home/projects/hearing_loss/clsaARHL_SA/clsaFiles
DATA=/home/projects/hearing_loss/clsaARHL_SA/xchr
PLINK=/opt/plink-1.09/plink
XWAS=/opt/github_repository/xwas-3.0-master/

#extract the SNPs that has been successfully lifted from hg38 to hg19
$PLINK --bfile $DATA/xchr_met_better --extract \
	$DATA/old_snp_ids.txt --make-bed --out $DATA/updated_xchr_met_better

#replace the old with the new bim 

mv $DATA/updated_xchr_met_better.bim $DATA/old_xchr_met_better.bim
cp $DATA/updated_xchr_hg19.bim $DATA/updated_xchr_met_better.bim


#conduct XWAS qc 
#use -s flag to skip IBD relatedness filtering

cd /home/projects/hearing_loss/clsaARHL_SA/xchr
#FID IID sen.better.rank met.better.rank sex age diabetes hypertension PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
#copy covariate file and generate two files for metabolic and sensory
cp /home/projects/hearing_loss/clsaARHL_SA/clsaFiles/covar_all_v2.txt .

##################for the metabolic phenotype
#Remove metabolic phenotype and sex from the covariate file
awk '{$4=$5=""; print $0}' covar_all_v2.txt > covar_met.txt
#run X-chr QC
./run_QC_noPCA_noRelated.sh -s -g params_file.txt

#run sex stratified XWAS with multi male models
	--xwas --strat-sex --fishers --multi-xchr-model --covar $DATA/covar_met.txt \
	--out $DATA/met_xwas_stratFisher



##################for the sensory phenotype
#extract the SNPs that has been successfully lifted from hg38 to hg19
$PLINK --bfile $DATA/xchr_sen_better --extract \
       $DATA/old_snp_ids.txt --make-bed --out $DATA/updated_xchr_sen_better

#replace the old with the new bim 

mv $DATA/updated_xchr_sen_better.bim $DATA/old_xchr_sen_better.bim
cp $DATA/updated_xchr_hg19.bim $DATA/updated_xchr_sen_better.bim


#conduct XWAS qc 
#use -s flag to skip IBD relatedness filtering

cd /home/projects/hearing_loss/clsaARHL_SA/xchr


#FID IID sen.better.rank met.better.rank sex age diabetes hypertension PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
#for the sensory phenotype remove sensory phenotype column and sex from the covaiates file
awk '{$3=$5=""; print $0}' covar_all_v2.txt > covar_sen.txt
#run X-chr QC
./run_QC_noPCA_noRelated.sh -s -g params_file.txt

# #run sex stratified XWAS with multi male models
 $XWAS/bin/xwas --noweb --bfile $DATA/updated_xchr_sen_better.preprocessed_final_x\
        --xwas --strat-sex --fishers --multi-xchr-model --covar $DATA/covar_sen.txt \
        --out $DATA/sen_xwas_stratFisher
#########
