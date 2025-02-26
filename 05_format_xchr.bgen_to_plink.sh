#/bin/bash

DATA=/home/projects/archive/CLSA/raw_data/genetic_data
XCHR=/home/projects/hearing_loss/clsaARHL_SA/xchr
PLINK=/opt/plink2/plink2

# Convert to plink format

$PLINK --bgen $DATA/clsa_imp_23_v3.bgen ref-first \
        --sample $XCHR/clsa_imp_v3_updated.sample \
        --allow-extra-chr --max-alleles 2  --make-bed \
        -out $XCHR/clsa_imp_23_v3

