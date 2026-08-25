GWAS_TOOLS=/opt/github_repository/summary-gwas-imputation/src
DATA=/home/projects/hearing_loss/clsaARHL_SA/polyfun_noopp
OUTPUT=/home/projects/hearing_loss/clsaARHL_SA/polyfun_noopp
REF=/home/projects/archive/lab/reference_data/metaxcan/data/liftover

#liftover the GWAS summary stats from genome build hg38 to hg19

mkdir -p $OUTPUT

for stats in "sen_cluster4" "met_cluster4"
do
python3.8 $GWAS_TOOLS/gwas_parsing.py \
	-gwas_file $DATA/${stats}_to_liftover.txt \
	-liftover $REF/hg38ToHg19.over.chain.gz \
        -output_column_map CHR chromosome \
        -output_column_map BP position \
        -output_column_map A1 effect_allele \
        -output_column_map A2 non_effect_allele \
        -output_column_map BETA effect_size \
        -output_column_map P pvalue \
        -output_column_map NMISS sample_size \
        -output_column_map SE standard_error \
        --enforce_numeric_columns \
        --force_special_handling \
        -input_pvalue_fix 0 \
	--chromosome_format \
	-output_order panel_variant_id chromosome position effect_allele non_effect_allele \
	pvalue effect_size standard_error zscore sample_size\
	-output $OUTPUT/${stats}_lifted.txt
done


