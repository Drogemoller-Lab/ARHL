#/bin/bash

CHAIN=/home/projects/lab/reference_data/metaxcan/data/liftover
LIFTOVER=/home/projects/hearing_loss/scripts
DATA=/home/projects/hearing_loss/clsaARHL_SA/xchr

$LIFTOVER/liftOver  $DATA/xchr_hg38.bed $CHAIN/hg38ToHg19.over.chain.gz $DATA/xchr_hg19_converted.bed $DATA/unmapped
