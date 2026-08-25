#/bin/bash
GENES=/home/projects/hearing_loss/clsaARHL_SA/xchGeneBased

awk '$2 == "X"' ${GENES}/magma_genes.txt | awk 'NR > 1 {print 23, $3, $4, $10}' > ${GENES}/xchr_magma_genes_v2.txt

