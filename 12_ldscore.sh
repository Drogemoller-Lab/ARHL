DATA=/home/projects/hearing_loss/clsaARHL_SA/ldsc
LDSC=/home/ahmeds26@med.umanitoba.ca/ldsc

#reformat summary stats

#remove the original CHR, BP, A1, A2, TEST columns
cut --complement -f1-6 $DATA/sumstats/met_eur_gwas_tidy_rsID_noSen > $DATA/sumstats/met_eur_gwas_tidy_rsID_noSen_updated
cut --complement -f1-6 $DATA/sumstats/sen_eur_gwas_tidy_rsID_noMet > $DATA/sumstats/sen_eur_gwas_tidy_rsID_noMet_updated

$LDSC/munge_sumstats.py \
	--sumstats $DATA/sumstats/met_eur_gwas_tidy_rsID_noSen_updated \
	--N 17531 \
	--snp rsid \
	--a1 allele_1 \
	--a2 allele_0 \
	--p P \
	--signed-sumstats STAT,0 \
	--out $DATA/sumstats/met_munge_noSen

$LDSC/munge_sumstats.py \
        --sumstats $DATA/sumstats/sen_eur_gwas_tidy_rsID_noMet_updated \
        --N 17531 \
        --snp rsid \
        --a1 allele_1 \
        --a2 allele_0 \
        --p P \
        --signed-sumstats STAT,0 \
        --out $DATA/sumstats/sen_munge_noMet


#calculate genetic correlation
$LDSC/ldsc.py \
	--rg $DATA/sumstats/met_munge_noSen.sumstats.gz,$DATA/sumstats/sen_munge_noMet.sumstats.gz \
	--ref-ld-chr $DATA/LDscores/ \
	--w-ld-chr $DATA/LDweights/ \
	--out $DATA/corr_met_sen_noOpp



