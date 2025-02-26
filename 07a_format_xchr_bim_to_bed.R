library(data.table)
library(dplyr)
# format .bim to ucsc .bed file
#read chromosome X .bim file 
xchr_bim <- fread("/home/projects/hearing_loss/clsa_v2/xchr/xchr_v2/clsa_imp_x_v3_clean_v2.bim.gz", header =FALSE)

#select columns of chromosome and position
xchr_bed <- select(xchr_bim, V1, V4, V2)
#rename columns
colnames(xchr_bed) <- c("CHR", "start_hg38", "snp_id")
# add +1 to the position in the third position
#here we need to prevent scientific notation 
xchr_bed$end_hg38 <- as.integer(xchr_bed$start_hg38 + 1)
#replace 23 with X
xchr_bed$CHR <- gsub(23, "chrX", xchr_bed$CHR)
#rearrange the columns
xchr_bed <- select(xchr_bed, CHR, start_hg38, end_hg38, snp_id)

#save as .bed file

write.table(xchr_bed, "/home/projects/hearing_loss/clsa_v2/xchr/xchr_v2/xchr_hg38.bed",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
######################
#run the LiftOver script on bash
#######################
#convert the lifted .bed back to .bim
xchr_bed_hg19 <- fread("/home/projects/hearing_loss/clsa_v2/xchr/xchr_v2/xchr_hg19_converted.bed")
#rename columns 
colnames(xchr_bed_hg19) <- c("CHR", "start_hg19", "end_hg19", "snp_id")
#merge with older bed file by SNP identifier and chromosome to keep only snps on X chromosome
xchr_bed_hg19.merge <- merge(xchr_bed, xchr_bed_hg19, by = c("CHR", "snp_id"))
#save the old SNP ids to keep
old_snp_ids <- as.data.frame(xchr_bed_hg19.merge$snp_id)
write.table(old_snp_ids, "/home/projects/hearing_loss/clsa_v2/xchr/xchr_v2/old_snp_ids.txt", quote = F, col.names = F, row.names = F)

#generate new SNP ids using the start position from the positions on the lifted SNPs
#first separate the identifier 
xchr_bed_hg19.separated <- separate(xchr_bed_hg19.merge, col = snp_id, into = c("CHR.x", "BP.x", "major", "minor"), sep = ":")
#paste with the new values 
xchr_bed_hg19.separated$snp_id_new <- with(xchr_bed_hg19.separated, paste(CHR,start_hg19, major, minor, sep=":"))
#add 0 value for position in morgans
xchr_bed_hg19.separated$V3 <- 0
#select columns
xchr_hg19.bim <- select(xchr_bed_hg19.separated, CHR, snp_id_new, V3, start_hg19, minor, major)
#rename chromosome to 23
xchr_hg19.bim$CHR <- gsub("chrX", 23, xchr_hg19.bim$CHR)
#mark duplicates in new column
#xchr_hg19.bim$dup <- duplicated(xchr_hg19.bim$snp_id_new)
#tag the duplicated markers to be excluded
#xchr_hg19.bim$snp_id_nodups <- ifelse(xchr_hg19.bim$dup=="TRUE", gsub("^", "dup", xchr_hg19.bim$snp_id_new), xchr_hg19.bim$snp_id_new)
#select relevant columns
#xchr_hg19_nodups.bim <- select(xchr_hg19.bim, CHR, snp_id_nodups, V3, start_hg19, minor, major)
#save the updated .bim file
write.table(xchr_hg19.bim, "/home/projects/hearing_loss/clsa_v2/xchr/xchr_v2/updated_xchr_hg19.bim",
            quote = FALSE, col.names = FALSE, row.names = FALSE)

#save duplicated snps
#duplicates <- as.data.frame(xchr_hg19_nodups.bim[grep("^dup", xchr_hg19_nodups.bim$snp_id_nodups), 2])

#write.table(duplicates, "/home/projects/hearing_loss/clsa_v2/xchr/duplicate_snps_to_exclude.txt",
#            quote = FALSE, col.names = FALSE, row.names = FALSE)


#mark SNP names > 35 to be removed for converf and smartpca
#count the number of characters in snp column
#xchr_hg19_nodups.bim <- fread("/home/projects/hearing_loss/clsa_v2/xchr/updated_xchr_hg19.bim")
xchr_hg19_nodups.bim$count <- nchar(as.character(xchr_hg19_nodups.bim$V2))
long_snps_toremove <- subset(xchr_hg19_nodups.bim, count > 35)
long_snps_toremove <- select(long_snps_toremove, V2)

write.table(long_snps_toremove, "/home/projects/hearing_loss/clsa_v2/xchr/long_snps_toremove.txt",
            quote = F, col.names = F, row.names = F)



#then replaced the old bim with the new bim



