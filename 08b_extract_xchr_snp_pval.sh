#/bin/bash
XCHR=/home/projects/archive/previous_projects/xchr
XCHR_GENE=/home/projects/hearing_loss/clsaARHL_SA/xchGeneBased/

awk 'NR > 1 {print $2, $8}' ${XCHR}/met_xwas_stratFisher_model2.xstrat.linear > ${XCHR_GENE}/met_male_xwas_snp_pval.txt

awk 'NR > 1 {print $2, $11}' ${XCHR}/met_xwas_stratFisher_model2.xstrat.linear > ${XCHR_GENE}/met_female_xwas_snp_pval.txt

awk 'NR > 1 {print $2, $8}' ${XCHR}/sen_xwas_stratFisher_model2.xstrat.linear > ${XCHR_GENE}/sen_male_xwas_snp_pval.txt

awk 'NR > 1 {print $2, $11}' ${XCHR}/sen_xwas_stratFisher_model2.xstrat.linear > ${XCHR_GENE}/sen_female_xwas_snp_pval.txt

