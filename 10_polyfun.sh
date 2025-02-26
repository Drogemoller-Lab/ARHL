POLYFUN=/home/software/polyfun/
FINEMAP=/home/ahmeds26@med.umanitoba.ca/finemap_v1.4.1_x86_64/
DATA=/home/projects/hearing_loss/clsaARHL_SA/polyfun
PLINK=/opt/plink-1.09/plink
FILES=/home/projects/hearing_loss/clsaARHL_SA/plinkFiles

#run polyfun conda environment
#run the munging script

for phen in "met" "sen"
do
	for sex in "male" "female"
	do
		python3.8 $POLYFUN/munge_polyfun_sumstats.py \
			--sumstats $DATA/${phen}_cluster4_${sex}_lifted_hrm.txt \
			--out $DATA/${phen}_cluster4_${sex}_munged.parquet \
			--min-info 0 \
			--min-maf 0
		python3.8 $POLYFUN/extract_snpvar.py \
                      --sumstats $DATA/${phen}_cluster4_${sex}_munged.parquet \
                      --out $DATA/${phen}_cluster4_${sex}_snps_with_var.gz \
                      --allow-missing
	done
done

############
#run genome-wide polyfun at 1e-5

#run the finemapping
for phen in  "sen" "met"
do
	for tool in "finemap" "susie"
	do
		for sex in "male" "female"
		do
			#Genome-wide fine-mapping step 1: Creating region-specific jobs
			python3.8 $POLYFUN/create_finemapper_jobs.py\
				--sumstats $DATA/${phen}_cluster4_${sex}_snps_with_var.gz \
				--n 17531 \
				--pvalue-cutoff 1e-5 \
				--finemap-exe $FINEMAP/finemap_v1.4.1_x86_64 \
				--method ${tool} \
				--max-num-causal 5 \
				--out-prefix ${phen}_cluster4_${sex}_polyfun_all_${tool} \
				--jobs-file $DATA/${phen}_cluster4_${sex}_polyfun_all_${tool}_jobs.txt
			#Genome-wide fine-mapping step 2: Invoking the region-specific jobs
			chmod +x $DATA/${phen}_cluster4_${sex}_polyfun_all_${tool}_jobs.txt
			$DATA/${phen}_cluster4_${sex}_polyfun_all_${tool}_jobs.txt
			
			#Genome-wide fine-mapping step 3: Aggregating the results
			python3.8 $POLYFUN/aggregate_finemapper_results.py \
				--out-prefix $DATA/${phen}_cluster4_${sex}_polyfun_all_${tool} \
				--sumstats $DATA/${phen}_cluster4_${sex}_snps_with_var.gz \
				--out $DATA/${phen}_cluster4_${sex}_polyfun_${tool}_agg.txt.gz \
				--pvalue-cutoff 1e-5
		done
	done
done

################################3

#run polyfun on sensory male 
#run the munging script

python3.8 $POLYFUN/munge_polyfun_sumstats.py \
	--sumstats $DATA/sen_cluster4_male_lifted_hrm.txt \
        --out $DATA/male_cluster4_male_munged.parquet \
        --min-info 0 \
        --min-maf 0
python3.8 $POLYFUN/extract_snpvar.py \
	--sumstats $DATA/sen_cluster4_male_munged.parquet \
        --out $DATA/sen_cluster4_male_snps_with_var.gz \
        --allow-missing

#Genome-wide fine-mapping step 1: Creating region-specific jobs
python3.8 $POLYFUN/create_finemapper_jobs.py\
	--sumstats $DATA/sen_cluster4_male_snps_with_var.gz \
	--n 9191 \
	--pvalue-cutoff 1e-5 \
        --method susie \
        --max-num-causal 5 \
        --out-prefix sen_cluster4_male_polyfun_all_susie \
        --jobs-file $DATA/sen_cluster4_male_polyfun_all_susie_jobs.txt
#Genome-wide fine-mapping step 2: Invoking the region-specific jobs
chmod +x $DATA/sen_cluster4_male_polyfun_all_susie_jobs.txt
$DATA/sen_cluster4_male_polyfun_all_susie_jobs.txt

#Genome-wide fine-mapping step 3: Aggregating the results
python3.8 $POLYFUN/aggregate_finemapper_results.py \
	--out-prefix $DATA/sen_cluster4_male_polyfun_all_susie \
        --sumstats $DATA/sen_cluster4_male_snps_with_var.gz \
        --out $DATA/sen_cluster4_male_polyfun_susie_agg.txt.gz \
        --pvalue-cutoff 1e-5




