#/bin/bash

XCHR=/home/projects/archive/previous_projects/xchr
CADD=/home/projects/hearing_loss/clsaARHL_SA/cadd
HARM=/home/projects/hearing_loss/clsaARHL_SA/gwasFumaInput

#all these coordinates are in hg19
#metabolic
awk  '$1 == 5 && $2 >= 72921983 && $2 <= 73237818 {print $1, $2, $3, $4, $5}' $HARM/met_all_fuma_input.txt  > $CADD/met_arhgef28_var2.txt
awk  '$1 == 11 && $2 >= 9406169 && $2 <= 9469673 {print $1, $2, $3, $4, $5}' $HARM/met_all_fuma_input.txt > $CADD/met_ipo7_var2.txt
awk  '$1 == 16 && $2 >= 31191431 && $2 <= 31203127 {print $1, $2, $3, $4, $5}' $HARM/met_all_fuma_input.txt  > $CADD/met_fus_var2.txt

#sensory
awk  '$1 == 6 && $2 >= 45868045 && $2 <= 46048132 {print $1, $2, $3, $4, $5}' $HARM/sen_all_fuma_input.txt  > $CADD/sen_clic5_var2.txt
awk  '$1 == 17 && $2 >= 18012020 && $2 <= 18083116 {print $1, $2, $3, $4, $5}' $HARM/sen_all_fuma_input.txt  > $CADD/sen_myo15a_var2.txt
awk  '$1 == 17 && $2 >= 17991200 && $2 <= 18011285 {print $1, $2, $3, $4, $5}' $HARM/sen_all_fuma_input.txt  > $CADD/sen_drg2_var2.txt
awk  '$1 == 22 && $2 >= 50986462 && $2 <= 50989451 {print $1, $2, $3, $4, $5}' $HARM/sen_all_fuma_input.txt  > $CADD/sen_klhdc7b_var2.txt

#X chromosome
#metabolic
awk  '$1 == 23 && $3 >= 129115083 && $3 <= 129192058 {print $1, $3, $4, $5, $14}' $XCHR/met_xwas_stratFisher_model2.xstrat.linear  > $CADD/met_bcorl1_var2.txt

awk  '$1 == 23 && $3 >=  129198849 && $3 <= 129244691 {print $1, $3, $4, $5, $14}' $XCHR/met_xwas_stratFisher_model2.xstrat.linear  > $CADD/met_elf4_var2.txt

awk  '$1 == 23 && $3 >= 28605516 && $3 <= 29974840 {print $1, $3, $4, $5, $14}' $XCHR/met_xwas_stratFisher_model2.xstrat.linear  > $CADD/met_il1rapl1_var2.txt

awk  '$1 == 23 && $3 >= 28605516 && $3 <= 29974840 {print $1, $3, $4, $5, $14}' $XCHR/met_xwas_stratFisher_model2.xstrat.linear  > $CADD/met_il1rapl1_var2.txt

#sensory
awk  '$1 == 23 && $3 >= 129115083 && $3 <= 129192058 {print $1, $3, $4, $5, $14}' $XCHR/sen_xwas_stratFisher_model2.xstrat.linear  > $CADD/sen_bcorl1_var2.txt

awk  '$1 == 23 && $3 >= 122734412 && $3 <= 122866906 {print $1, $3, $4, $5, $14}' $XCHR/sen_xwas_stratFisher_model2.xstrat.linear  > $CADD/sen_thoc2_var2.txt

awk  '$1 == 23 && $3 >= 45007619 && $3 <= 45060146 {print $1, $3, $4, $5, $14}' $XCHR/sen_xwas_stratFisher_model2.xstrat.linear  > $CADD/sen_cxorf36_var2.txt

awk  '$1 == 23 && $3 >= 38128416 && $3 <= 38186817 {print $1, $3, $4, $5, $14}' $XCHR/sen_xwas_stratFisher_model2.xstrat.linear  > $CADD/sen_rpgr_var2.txt










