#!/bin/sh

#/home/projects/hearing_loss/clsaARHL_SA/xchGeneBased/mod_gene_based_test.sh params_file_genes_met.txt

#/home/projects/hearing_loss/clsaARHL_SA/xchGeneBased/mod_gene_based_test.sh params_file_genes_met_male.txt

#/home/projects/hearing_loss/clsaARHL_SA/xchGeneBased/mod_gene_based_test.sh params_file_genes_met_female.txt


#/home/projects/hearing_loss/clsaARHL_SA/xchGeneBased/mod_gene_based_test.sh params_file_genes_sen.txt


#/home/projects/hearing_loss/clsaARHL_SA/xchGeneBased/mod_gene_based_test.sh params_file_genes_sen_male.txt


#/home/projects/hearing_loss/clsaARHL_SA/xchGeneBased/mod_gene_based_test.sh params_file_genes_sen_female.txt

for phen in "met" "sen" "met_male" "met_female" "sen_male" "sen_female"
do
sed -e "s/\r//g" ${phen}_gene_test_result_v2.txt.sort > ${phen}_gene_test_result_v2.txt.sort.clean

done



