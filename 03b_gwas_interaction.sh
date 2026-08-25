DATA=/home/projects/archive/previous_projects/plinkFiles/
COVAR=/home/projects/hearing_loss/clsaARHL_SA/clsaFiles
PLINK=/opt/plink-1.09/plink
GWAS=/home/projects/hearing_loss/clsaARHL_SA/gwas

# Targeted SNP interaction analysis for selected GWAS variants.
# Uses PLINK linear interaction models with the GWAS covariates.

# Perform association testing
#covar files header 
#FID IID sen.better.rank met.better.rank sex age diabetes hypertension PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
#FID IID sen.better.rank met.better.rank sex age diabetes hypertension ePC1 ePC2 ePC3 ePC4 ePC5 ePC6 ePC7 ePC8 ePC9 ePC10

#get the variant that i want to investigate chr8:69662586:T:A
$PLINK --bfile $DATA/sen_better \
       --snp chr8:69662586:T:A \
--make-bed \
--out $GWAS/sen_better_chr8_69662586

$PLINK --bfile $GWAS/sen_better_chr8_69662586 \
       --linear interaction --ci 0.95 --covar $COVAR/covar_all_v2.txt \
       --parameters 1-16,18 \
       --covar-name met.better.rank-PC10 \
       --out $GWAS/sen_better_gwas_chr8_69662586
$PLINK --bfile $DATA/sen_better \
        --snp chr8:69654415:G:T \
        --make-bed \
        --out $GWAS/sen_better_chr8_69654415

$PLINK --bfile $GWAS/sen_better_chr8_69654415 \
        --linear interaction --ci 0.95 --covar $COVAR/covar_all_v2.txt \
        --parameters 1-16,18 \
        --covar-name met.better.rank-PC10 \
        --out $GWAS/sen_better_gwas_chr8_69654415

$PLINK --bfile $DATA/sen_better \
       --snp chr22:50549676:G:A \
      --make-bed \
      --out $GWAS/sen_better_chr22_50549676

$PLINK --bfile $GWAS/sen_better_chr22_50549676 \
       --linear interaction --ci 0.95 --covar $COVAR/covar_all_v2.txt \
       --parameters 1-16,18 \
       --covar-name met.better.rank-PC10 \
       --out $GWAS/sen_better_gwas_chr22_50549676

