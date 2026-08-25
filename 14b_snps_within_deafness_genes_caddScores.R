#load libraries
#library(usethis)
#library(devtools)
library(data.table)
#install.packages("remotes")
#remotes::install_github("juliedwhite/miamiplot")
library(miamiplot)
library(tidyverse)
library(hrbrthemes)
library(broom)
library(gridExtra)
library(ggpubr)

if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")


#########
# QQplot for SNPs within these genes
met_gwas_hg19 <- fread ("/home/projects/hearing_loss/clsaARHL_SA/gwasFumaInput/met_all_fuma_input.txt")
sen_gwas_hg19 <- fread ("/home/projects/hearing_loss/clsaARHL_SA/gwasFumaInput/sen_all_fuma_input.txt")

gene_list_magma <- fread ("/home/projects/hearing_loss/clsaARHL_SA/MAGMA_output/genes_position_magma.txt")
gene_list_magma_auto <- subset(gene_list_magma, CHR != "X")
gene_list_magma_auto$CHR <- as.integer(gene_list_magma_auto$CHR)

# Function to check if a SNP is within 50,000 bp of a gene's boundaries
is_within_gene <- function(chr, pos, gene_df, window = 50000) {
  gene_names <- gene_df %>%
    filter(CHR == chr, (START - window) <= pos, (STOP + window) >= pos) %>%
    select(SYMBOL) %>%
    pull()
  
  if (length(gene_names) == 0) {
    return(NA)
  } else {
    return(paste(gene_names, collapse = ","))
  }
}


# Apply the function to each SNP
met_snps_within_genes <- met_gwas_hg19 %>%
  rowwise() %>%
  mutate(gene = is_within_gene(chromosome, position, gene_list_magma_auto)) %>%
  filter(!is.na(gene)) %>%
  ungroup()

sen_snps_within_genes <- sen_gwas_hg19 %>%
  rowwise() %>%
  mutate(gene = is_within_gene(chromosome, position, gene_list_magma_auto)) %>%
  filter(!is.na(gene)) %>%
  ungroup()

# # Write the result
# write.table(met_snps_within_genes, "/home/projects/hearing_loss/clsaARHL_SA/MAGMA_output/met_snps_within_genes_window.txt", quote = FALSE, row.names = FALSE)
# write.table(sen_snps_within_genes, "/home/projects/hearing_loss/clsaARHL_SA/MAGMA_output/sen_snps_within_genes_window.txt", quote = FALSE, row.names = FALSE)

###############
# snps within genes for the X chromosome
metabolic_xchr <- fread(file = "/home/projects/archive/previous_projects/xchr/met_xwas_stratFisher_model2.xstrat.linear")
sensory_xchr <- fread(file = "/home/projects/archive/previous_projects/xchr/sen_xwas_stratFisher_model2.xstrat.linear")
gene_list_magma_x <- subset(gene_list_magma, CHR == "X")
gene_list_magma_x$CHR <- as.integer(gsub("X", 23, gene_list_magma_x$CHR))



# Apply the function to each SNP
met_snps_within_genes_x <- metabolic_xchr %>%
  rowwise() %>%
  mutate(gene = is_within_gene(CHR, BP, gene_list_magma_x)) %>%
  filter(!is.na(gene)) %>%
  ungroup()

sen_snps_within_genes_x <- sensory_xchr %>%
  rowwise() %>%
  mutate(gene = is_within_gene(CHR, BP, gene_list_magma_x)) %>%
  filter(!is.na(gene)) %>%
  ungroup()

#write the files
#write.table(met_snps_within_genes_x, "/home/projects/hearing_loss/clsaARHL_SA/MAGMA_output/met_snps_within_genes_window_x.txt", quote = FALSE, row.names = FALSE)
#write.table(sen_snps_within_genes_x, "/home/projects/hearing_loss/clsaARHL_SA/MAGMA_output/sen_snps_within_genes_window_x.txt", quote = FALSE, row.names = FALSE)

#format data frames to prepare for row bind
met_snps_within_genes <- fread("/home/projects/hearing_loss/clsaARHL_SA/MAGMA_output/met_snps_within_genes_window.txt")
sen_snps_within_genes <- fread("/home/projects/hearing_loss/clsaARHL_SA/MAGMA_output/sen_snps_within_genes_window.txt")
met_snps_within_genes_x <- fread("/home/projects/hearing_loss/clsaARHL_SA/MAGMA_output/met_snps_within_genes_window_x.txt")
sen_snps_within_genes_x <- fread("/home/projects/hearing_loss/clsaARHL_SA/MAGMA_output/sen_snps_within_genes_window_x.txt")

met_snps_within_genes.format <- met_snps_within_genes |> 
  dplyr::select(chromosome, position, effect_allele, non_effect_allele, pvalue, gene) |> 
  dplyr::rename(CHR = chromosome,
                BP = position,
                A1 = effect_allele,
                A2 = non_effect_allele,
                P = pvalue)

sen_snps_within_genes.format <- sen_snps_within_genes |> 
  dplyr::select(chromosome, position, effect_allele, non_effect_allele, pvalue, gene) |> 
  dplyr::rename(CHR = chromosome,
                BP = position,
                A1 = effect_allele,
                A2 = non_effect_allele,
                P = pvalue)

met_snps_within_genes_x.format <- met_snps_within_genes_x |> 
  dplyr::select(CHR, BP, A1, A2, P_comb_Fisher, gene) |> 
  dplyr::rename(P = P_comb_Fisher)

sen_snps_within_genes_x.format <- sen_snps_within_genes_x |> 
  dplyr::select(CHR, BP, A1, A2, P_comb_Fisher, gene) |> 
  dplyr::rename(P = P_comb_Fisher)

#####
#combine the x chromosome data 
met_snps_within_genes.comb <- rbind (met_snps_within_genes.format, met_snps_within_genes_x.format)
sen_snps_within_genes.comb <- rbind(sen_snps_within_genes.format, sen_snps_within_genes_x.format)

###############

dfn_genes <- fread("/home/projects/hearing_loss/clsaARHL_SA/MAGMA_output/deafness_genes.txt")

#separate multiple genes
met_snps_within_genes.comb <- separate(met_snps_within_genes.comb, col = gene, into = c("gene1", "gene2", "gene3", "gene4"), sep = ",")
sen_snps_within_genes.comb <- separate(sen_snps_within_genes.comb, col = gene, into = c("gene1", "gene2", "gene3", "gene4"), sep = ",")

met_snps_within_genes.comb$dfn <- ifelse(met_snps_within_genes.comb$gene1 %in% dfn_genes$gene, met_snps_within_genes.comb$gene1,
                                         ifelse (met_snps_within_genes.comb$gene2 %in% dfn_genes$gene, met_snps_within_genes.comb$gene2, 
                                                 ifelse(met_snps_within_genes.comb$gene3 %in% dfn_genes$gene, met_snps_within_genes.comb$gene3,
                                                        ifelse(met_snps_within_genes.comb$gene4 %in% dfn_genes$gene, met_snps_within_genes.comb$gene4, NA))))

met_snps_within_genes.comb$deafness <- ifelse(is.na(met_snps_within_genes.comb$dfn), "non_deafness_gene", "deafness_gene")

sen_snps_within_genes.comb$dfn <- ifelse(sen_snps_within_genes.comb$gene1 %in% dfn_genes$gene, sen_snps_within_genes.comb$gene1,
                                         ifelse (sen_snps_within_genes.comb$gene2 %in% dfn_genes$gene, sen_snps_within_genes.comb$gene2, 
                                                 ifelse(sen_snps_within_genes.comb$gene3 %in% dfn_genes$gene, sen_snps_within_genes.comb$gene3,
                                                        ifelse(sen_snps_within_genes.comb$gene4 %in% dfn_genes$gene, sen_snps_within_genes.comb$gene4, NA))))

sen_snps_within_genes.comb$deafness <- ifelse(is.na(sen_snps_within_genes.comb$dfn), "non_deafness_gene", "deafness_gene")
#metabolic phenotype
#if Gene 2 is TRIOBP, replace NOL12 with TRIOBP
met_snps_within_genes.comb$final_gene <- ifelse(!is.na(met_snps_within_genes.comb$gene2) & met_snps_within_genes.comb$gene2 == met_snps_within_genes.comb$dfn, 
                                                met_snps_within_genes.comb$gene2, 
                                                met_snps_within_genes.comb$gene1)

sen_snps_within_genes.comb$final_gene <- ifelse(!is.na(sen_snps_within_genes.comb$gene2) & sen_snps_within_genes.comb$gene2 == sen_snps_within_genes.comb$dfn, 
                                                sen_snps_within_genes.comb$gene2, 
                                                sen_snps_within_genes.comb$gene1)
# Calculate expected -log10(p-values) under null hypothesis
met_snps_within_genes.comb <- met_snps_within_genes.comb %>%
  arrange(P) %>%
  mutate(expected = -log10(ppoints(n())),
         observed = -log10(P))

sen_snps_within_genes.comb <- sen_snps_within_genes.comb %>%
  arrange(P) %>%
  mutate(expected = -log10(ppoints(n())),
         observed = -log10(P))


#calculate observed and expected values for only deafness genes
met_snps_within_genes_deafness <- filter(met_snps_within_genes.comb, deafness == "deafness_gene")
sen_snps_within_genes_deafness <- filter(sen_snps_within_genes.comb, deafness == "deafness_gene")

met_snps_within_genes_deafness <- met_snps_within_genes_deafness %>%
  arrange(P) %>%
  mutate(expected = -log10(ppoints(n())),
         observed = -log10(P))

sen_snps_within_genes_deafness <- sen_snps_within_genes_deafness %>%
  arrange(P) %>%
  mutate(expected = -log10(ppoints(n())),
         observed = -log10(P))

#label the top snp for each gene
met_snps_within_genes.comb.top <- met_snps_within_genes.comb %>% 
  group_by(final_gene) %>% 
  slice_min(P) 
met_snps_within_genes_deafness.top <- met_snps_within_genes_deafness %>% 
  group_by(final_gene) %>% 
  slice_min(P) 
sen_snps_within_genes.comb.top <- sen_snps_within_genes.comb %>% 
  group_by(final_gene) %>% 
  slice_min(P) 
sen_snps_within_genes_deafness.top <- sen_snps_within_genes_deafness %>% 
  group_by(final_gene) %>% 
  slice_min(P) 

# Create the QQ plot using ggplot2
library(ggplot2)
library(ggrepel)

#combine the the two dataframes
# Combine the datasets
met_snps_within_genes.comb$source <- "Variants within all genes"
met_snps_within_genes_deafness$source <- "Variants within deafness genes"

# Combine both dataframes into one
met_combined_df <- bind_rows(met_snps_within_genes.comb, met_snps_within_genes_deafness)

# Combine the datasets
sen_snps_within_genes.comb$source <- "Variants within all genes"
sen_snps_within_genes_deafness$source <- "Variants within deafness genes"

# Combine both dataframes into one
sen_combined_df <- bind_rows(sen_snps_within_genes.comb, sen_snps_within_genes_deafness)



# Create the ggplot
met_qqplot <- ggplot(met_combined_df) +
  geom_point(aes(x = expected, y = observed, color = source), size = 1) +
  geom_abline(intercept = 0, slope = 1, col = "black", linetype = "dashed") +
  scale_color_manual(values = c("Variants within all genes" = "steelblue", "Variants within deafness genes" = "grey")) +
  
  scale_y_continuous(limits = c(0, 14)) +
  labs(title = "Metabolic phenotype", x = "Expected -log10(p)") +
  theme_minimal() + theme(legend.title = element_blank(), axis.title.x = element_text(size = 8),
                          axis.title.y = element_blank(), plot.title = element_text(size = 10, hjust = 0.5, vjust = 0.5),
                          legend.text = element_text(size = 7))

sen_qqplot <- ggplot(sen_combined_df) +
  geom_point(aes(x = expected, y = observed, color = source), size = 1) +
  geom_abline(intercept = 0, slope = 1, col = "black", linetype = "dashed") +
  scale_color_manual(values = c("Variants within all genes" = "darkred", "Variants within deafness genes" = "grey")) +
  scale_y_continuous(limits = c(0, 14)) +
  labs(title = "Sensory phenotype", x = "Expected -log10(p)") +
  theme_minimal() + theme(legend.title = element_blank(), axis.title.x = element_text(size = 8),
                          axis.title.y = element_blank(), plot.title = element_text(size = 10, hjust = 0.5, vjust = 0.5),
                          legend.text = element_text(size = 7))


# Combine plots
combined_qqplot <- ggarrange(met_qqplot, sen_qqplot, ncol = 2, common.legend = F,
                             label.x = "Expected -log10(p)", legend = "bottom")

annotated_qqplot <- annotate_figure(combined_qqplot, 
                                    left = text_grob("Observed -log10(p)", rot = 90, vjust = 1, size = 8)) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))

ggplot2::ggsave("/home/projects/hearing_loss/clsaARHL_SA/annotated_qqplot_window.tiff",
                height=4, width=6, units='in', dpi=300)

#calculate Z scores
# Filter invalid p-values first
met_snps_within_genes.comb <- met_snps_within_genes.comb %>%
  filter(!is.na(P), P > 0, P <= 1)

sen_snps_within_genes.comb <- sen_snps_within_genes.comb %>%
  filter(!is.na(P), P > 0, P <= 1)

met_snps_within_genes.comb$z <- qnorm(1 - met_snps_within_genes.comb$P / 2)
met_snps_within_genes.comb$z2 <- met_snps_within_genes.comb$z^2

sen_snps_within_genes.comb$z <- qnorm(1 - sen_snps_within_genes.comb$P / 2)
sen_snps_within_genes.comb$z2 <- sen_snps_within_genes.comb$z^2

# Split by deafness gene status 
met_deafness     <- met_snps_within_genes.comb %>% filter(deafness == "deafness_gene")
met_non_deafness <- met_snps_within_genes.comb %>% filter(deafness == "non_deafness_gene")

sen_deafness     <- sen_snps_within_genes.comb %>% filter(deafness == "deafness_gene")
sen_non_deafness <- sen_snps_within_genes.comb %>% filter(deafness == "non_deafness_gene")

# Wilcoxon test: are deafness genes more significant?
met_wilcox <- wilcox.test(met_deafness$z2, met_non_deafness$z2, alternative = "greater")
sen_wilcox <- wilcox.test(sen_deafness$z2, sen_non_deafness$z2, alternative = "greater")

# ── Mean Z² per group ────────────────────────────────────────────────
met_mean_Z2_deaf  <- mean(met_deafness$z2,     na.rm = TRUE)
met_mean_Z2_other <- mean(met_non_deafness$z2, na.rm = TRUE)

sen_mean_Z2_deaf  <- mean(sen_deafness$z2,     na.rm = TRUE)
sen_mean_Z2_other <- mean(sen_non_deafness$z2, na.rm = TRUE)

cat(
  "Metabolic phenotype\n",
  "Wilcox P =", signif(met_wilcox$p.value, 3), "\n",
  "Mean Z² deafness =", round(met_mean_Z2_deaf, 2), "\n",
  "Mean Z² other =", round(met_mean_Z2_other, 2), "\n\n"
)

cat(
  "Sensory phenotype\n",
  "Wilcox P =", signif(sen_wilcox$p.value, 3), "\n",
  "Mean Z² deafness =", round(sen_mean_Z2_deaf, 2), "\n",
  "Mean Z² other =", round(sen_mean_Z2_other, 2), "\n"
)
###################
#to extract the SNPs that deviated from the null line
#get SNPs at p < 1e-5 for SNPs within all genes and p < 1e-3 for SNPs within deafness genes
#for metabolic phenotype
#for SNPs in all genes
met_snps_within_genes.1e5 <- filter(met_snps_within_genes.comb, P < 1e-5)
met_snps_within_genes.names <- unique(met_snps_within_genes.1e5$final_gene)
#for SNPs in deafness genes
met_snps_within_genes_deafness.1e3 <- filter(met_snps_within_genes_deafness, P < 1e-3)
met_snps_within_genes_deafness.names <- unique(met_snps_within_genes_deafness.1e3$final_gene)

#for sensory phenotype
#for SNPs in all genes
sen_snps_within_genes.1e5 <- filter(sen_snps_within_genes.comb, P < 1e-5)
sen_snps_within_genes.names <- unique(sen_snps_within_genes.1e5$final_gene)
#for SNPs in deafness genes
sen_snps_within_genes_deafness.1e3 <- filter(sen_snps_within_genes_deafness, P < 1e-3)
sen_snps_within_genes_deafness.names <- unique(sen_snps_within_genes_deafness.1e3$final_gene)

intersect(met_snps_within_genes.names, sen_snps_within_genes.names)
intersect(met_snps_within_genes_deafness.names, sen_snps_within_genes_deafness.names)

#save the files to submit to CADD
write.table(met_snps_within_genes.1e5, "/home/projects/hearing_loss/clsaARHL_SA/cadd/input/met_snps_within_genes_window.1e5.txt", quote = F, row.names = F)
write.table(sen_snps_within_genes.1e5, "/home/projects/hearing_loss/clsaARHL_SA/cadd/input/sen_snps_within_genes_window.1e5.txt", quote = F, row.names = F)
write.table(met_snps_within_genes_deafness.1e3, "/home/projects/hearing_loss/clsaARHL_SA/cadd/input/met_snps_within_genes_deafness_window.1e3.txt", quote = F, row.names = F)
write.table(sen_snps_within_genes_deafness.1e3, "/home/projects/hearing_loss/clsaARHL_SA/cadd/input/sen_snps_within_genes_deafness_window.1e3.txt", quote = F, row.names = F)

###########################################
#add the CADD results to the tables of deafness genes
met_snps_within_genes_deafness.cadd_out <- fread("/home/projects/hearing_loss/clsaARHL_SA/cadd/output/met_snps_within_genes_deafness_window")
met_snps_within_genes_deafness.cadd_out <- select(met_snps_within_genes_deafness.cadd_out, Chrom, Pos,
                                                  ConsDetail, protPos, oAA, nAA, protPos, SIFTcat, SIFTval, PolyPhenCat, PolyPhenVal, PHRED)
#merge with GWAS data
met_snps_within_genes_deafness.cadd_out <- merge(select(met_snps_within_genes_deafness.1e3, CHR, BP, A1, A2, P, final_gene), met_snps_within_genes_deafness.cadd_out,
                                                 by.x= c("CHR", "BP"), by.y = c("Chrom", "Pos"), all.x = TRUE)
met_snps_within_genes_deafness.cadd_out <- dplyr::rename(met_snps_within_genes_deafness.cadd_out, met_P=P)

sen_snps_within_genes_deafness.cadd_out <- fread("/home/projects/hearing_loss/clsaARHL_SA/cadd/output/sen_snps_within_genes_deafness_window")
sen_snps_within_genes_deafness.cadd_out <- select(sen_snps_within_genes_deafness.cadd_out, Chrom, Pos,
                                                  ConsDetail, protPos, oAA, nAA, protPos, SIFTcat, SIFTval, PolyPhenCat, PolyPhenVal, PHRED)
#merge with GWAS data
sen_snps_within_genes_deafness.cadd_out <- merge(select(sen_snps_within_genes_deafness.1e3, CHR, BP, A1, A2, P, final_gene), sen_snps_within_genes_deafness.cadd_out,
                                                 by.x= c("CHR", "BP"), by.y = c("Chrom", "Pos"), all.x = TRUE)
sen_snps_within_genes_deafness.cadd_out <- dplyr::rename(sen_snps_within_genes_deafness.cadd_out, sen_P=P)


#######
#window tables
met_snps_within_deafness_genes <- read.table("/home/projects/hearing_loss/clsaARHL_SA/cadd/output/output_manual/met_deafness_genes.merge.txt", header = TRUE)
sen_snps_within_deafness_genes <- read.table("/home/projects/hearing_loss/clsaARHL_SA/cadd/output/output_manual/sen_deafness_genes.merge.txt", header = TRUE)

met_snps_within_deafness_genes <- dplyr::rename(met_snps_within_deafness_genes,
                                                met_P = P,
                                                variant_conseq = ConsDetail,
                                                gene = dfn)
sen_snps_within_deafness_genes <- dplyr::rename(sen_snps_within_deafness_genes,
                                                sen_P = P,
                                                variant_conseq = ConsDetail,
                                                gene = dfn)
#get the pvalue for the opposite phenotype
met_snps_within_deafness_genes.annotated <- merge(met_snps_within_deafness_genes, select(sen_all_hg19, chromosome, position, pvalue),
                                                  by.x = c("CHR", "BP"), by.y = c("chromosome", "position"))
sen_snps_within_deafness_genes.annotated <- merge(sen_snps_within_deafness_genes, select(met_all_hg19, chromosome, position, pvalue),
                                                  by.x = c("CHR", "BP"), by.y = c("chromosome", "position"))
#rename and rearrange columns
met_snps_within_deafness_genes.annotated.2 <- met_snps_within_deafness_genes.annotated |> 
  dplyr::rename(sen_P = pvalue) |> 
  select(CHR, BP, A1, A2, variant_conseq, PHRED, met_P, sen_P, gene)
sen_snps_within_deafness_genes.annotated.2 <- sen_snps_within_deafness_genes.annotated |> 
  dplyr::rename(met_P = pvalue) |> 
  select(CHR, BP, A1, A2, variant_conseq, PHRED, met_P, sen_P, gene)


#get the full list
snps_within_genes_deafness.annotated <- rbind(met_snps_within_deafness_genes.annotated.2, sen_snps_within_deafness_genes.annotated.2)
#save it
write.table(snps_within_genes_deafness.annotated, "/home/projects/hearing_loss/clsaARHL_SA/cadd/output/snps_within_genes_deafness_annotated_window_v2.txt",
            quote = F, row.names = F)
#####################################
#save the files
write.table(met_snps_within_genes_deafness.cadd_out, "/home/projects/hearing_loss/clsaARHL_SA/cadd/output/met_snps_within_genes_deafness_annotated.txt",
            quote = F, row.names = F)

write.table(sen_snps_within_genes_deafness.cadd_out, "/home/projects/hearing_loss/clsaARHL_SA/cadd/output/sen_snps_within_genes_deafness_annotated.txt",
            quote = F, row.names = F)
write.table(snps_within_genes_deafness.cadd_out, "/home/projects/hearing_loss/clsaARHL_SA/cadd/output/snps_within_genes_deafness_annotated.txt",
            quote = F, row.names = F)
###########################################



