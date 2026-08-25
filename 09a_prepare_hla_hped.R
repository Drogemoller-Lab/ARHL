library(tidyr)

#needs tidyr v1.0.0+
#Remove comment lines above header before running this script
#assume header as below separated by comma
#ADM_GWAS_COM,Locus,Allele1,Allele1_Probability,Allele2,Allele2_Probability,Combined_Probability

hla_CLSA <- read.table(file="/home/projects/archive/CLSA/raw_data/genetic_data/clsa_hla_v3.csv", 
                       sep =",", header = T,
                       colClasses=c('numeric','character','character','numeric','character','numeric','numeric'))

#only keep the alleles
hla_CLSA_long_2alleles <- hla_CLSA[,c(1,2,3,5)]

#fill missing values using 0 as compatible with PLINK ped format
hla_CLSA_wide_2alleles<-pivot_wider(data = hla_CLSA_long_2alleles, 
                                    id_cols = ADM_GWAS_COM, 
                                    names_from = Locus, 
                                    values_from = c("Allele1", "Allele2"),
                                    names_vary = "slowest",
                                    values_fill = "0")

#import fam file of chr6
met_chr6_fam <- read.table("/home/projects/hearing_loss/clsaARHL_SA/hla/met_better_chr6.fam")
colnames(met_chr6_fam)[c(1,2,5,6)] <- c("FID", "IID", "sex", "phenotype")


#merge with the .fam file to generate .hped format
met_hla_CLSA_wide_2alleles_hped <- merge(met_chr6_fam, hla_CLSA_wide_2alleles, by.x = "FID", by.y = "ADM_GWAS_COM")

#remove DRG3,4,and 5
met_hla_CLSA_wide_2alleles_hped_no_drb345 <- met_hla_CLSA_wide_2alleles_hped[, -c(19:24)]
  
#save file
write.table(met_hla_CLSA_wide_2alleles_hped_no_drb345, "/home/projects/hearing_loss/clsaARHL_SA/hla/met_better_hla_noDRB345.hped",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

############

#for the sensory phenotype
sen_chr6_fam <- read.table("/home/projects/hearing_loss/clsaARHL_SA/hla/sen_better_chr6.fam")
colnames(sen_chr6_fam)[c(1,2,5,6)] <- c("FID", "IID", "sex", "phenotype")

#merge with the .fam file to generate .hped format
sen_hla_CLSA_wide_2alleles_hped <- merge(sen_chr6_fam, hla_CLSA_wide_2alleles, by.x = "FID", by.y = "ADM_GWAS_COM")

#remove DRG3,4,and 5
sen_hla_CLSA_wide_2alleles_hped_no_drb345 <- sen_hla_CLSA_wide_2alleles_hped[, -c(19:24)]

#save file
write.table(sen_hla_CLSA_wide_2alleles_hped_no_drb345, "/home/projects/hearing_loss/clsaARHL_SA/hla/sen_better_hla_noDRB345.hped",
            quote = FALSE, row.names = FALSE, col.names = FALSE)



