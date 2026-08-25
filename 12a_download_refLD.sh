#/bin/bash

#download LD files
#wget https://zenodo.org/records/10515792/files/1000G_Phase3_baselineLD_v2.2_ldscores.tgz
#tar -zxvf ../1000G_Phase3_baselineLD_v2.2_ldscores.tgz

#download weight files

wget https://zenodo.org/records/10515792/files/1000G_Phase3_weights_hm3_no_MHC.tgz

mkdir -p LDweights
tar --strip-components=1 -zxvf 1000G_Phase3_weights_hm3_no_MHC.tgz -C LDweights


