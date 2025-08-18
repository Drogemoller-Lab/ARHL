#load libraries
#library(usethis)
#library(devtools)
library(data.table)
library(tidyverse)
library(hrbrthemes)
library(broom)
library(gridExtra)
library(ggpubr)

if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")


###################################

#get the top SNPs from finemapping results
##finemap results
#for the genome-wide metabolic
met_stats_all_finemap <- fread("/home/projects/hearing_loss/clsaARHL_SA/polyfun/met_cluster4_polyfun_finemap_agg.txt.gz")
met_stats_all_finemap_2 <- separate(met_stats_all_finemap, col = CREDIBLE_SET, into = c("CHR.x", "POSITION", "SET"), sep = ":")
met_stats_all_finemap_credible <- filter(met_stats_all_finemap_2, PIP >= 0.5 | SNP == "rs6453022"& SET == 1)
met_stats_all_finemap_credible$finemap_group <- "Variants significant in metabolic phenotype"
#merge with gwas results (hg19)
met_all_hg19 <- fread ("/home/projects/hearing_loss/clsaARHL_SA/gwasFumaInput/met_all_fuma_input.txt")
met_stats_all_finemap_credible.merge <- merge(met_stats_all_finemap_credible, met_all_hg19, 
                                              by.x = c("CHR", "BP"), by.y = c("chromosome", "position"))

#for the genome-wide sensory
sen_stats_all_finemap <- fread("/home/projects/hearing_loss/clsaARHL_SA/polyfun/sen_cluster4_polyfun_finemap_agg.txt.gz")
sen_stats_all_finemap_2 <- separate(sen_stats_all_finemap, col = CREDIBLE_SET, into = c("CHR.x", "POSITION", "SET"), sep = ":")
sen_stats_all_finemap_credible <- filter(sen_stats_all_finemap_2,PIP >= 0.5 & SET == 1)

sen_stats_all_finemap_credible$finemap_group <- "Variants significant in sensory phenotype"
#merge with gwas results (hg19)
sen_all_hg19 <- fread ("/home/projects/hearing_loss/clsaARHL_SA/gwasFumaInput/sen_all_fuma_input.txt")
sen_stats_all_finemap_credible.merge <- merge(sen_stats_all_finemap_credible, sen_all_hg19, 
                                              by.x = c("CHR", "BP"), by.y = c("chromosome", "position"))

#############################################
#get a list of significant SNPs
polyfun_all_sig_snps <- rbind(select(met_stats_all_finemap_credible, SNP, CHR, BP, PIP, finemap_group),
                              select(sen_stats_all_finemap_credible, SNP, CHR, BP, PIP, finemap_group))

#get the effect sizes of these SNPs from the GWAS results hg19
met_polyfun_sig_snps <- merge(polyfun_all_sig_snps, met_all_hg19, 
                              by.x = c("CHR", "BP"), by.y = c("chromosome", "position")) 

sen_polyfun_sig_snps <- merge(polyfun_all_sig_snps, sen_all_hg19, 
                              by.x = c("CHR", "BP"), by.y = c("chromosome", "position")) 

data_polyfun_sig_snps <- merge(met_polyfun_sig_snps, sen_polyfun_sig_snps, by = c("CHR", "BP", "SNP", "finemap_group",
                                                                                  "effect_allele", "non_effect_allele",
                                                                                  "PIP", "sample_size"))
#add gene names
data_polyfun_sig_snps$genes <- ifelse(data_polyfun_sig_snps$SNP == "rs6453022", "ARHGEF28", 
                                      ifelse(data_polyfun_sig_snps$SNP == "rs76510133", "EYA4",
                                             ifelse(data_polyfun_sig_snps$SNP == "rs4832151", "TCF7L1",
                                                    ifelse(data_polyfun_sig_snps$SNP == "rs7428702", "HEG1",
                                                           ifelse(data_polyfun_sig_snps$SNP == "rs17008250", "LINC02994",
                                                                  ifelse(data_polyfun_sig_snps$SNP == "rs77196213", "LINC02716",
                                                                         ifelse(data_polyfun_sig_snps$SNP == "rs56105727", "FGF14",
                                                                                ifelse(data_polyfun_sig_snps$SNP == "rs1604887", "OTX2-AS1",
                                                                                       ifelse(data_polyfun_sig_snps$SNP == "rs200974429", "LOC105370982",
                                                                                              ifelse(data_polyfun_sig_snps$SNP == "rs72628292", "FENDRR",
                                                                                                     ifelse(data_polyfun_sig_snps$SNP == "rs36100160", "CDYL", 
                                                                                                            ifelse(data_polyfun_sig_snps$SNP == "rs34610101", "CAMKK1",
                                                                                                                   ifelse(data_polyfun_sig_snps$SNP == "rs36062310", "KLHDC7B",
                                                                                                                          ifelse(data_polyfun_sig_snps$SNP == "rs80122938", "PDPN",
                                                                                                                                 ifelse(data_polyfun_sig_snps$SNP == "rs115671365", "GFOD1",
                                                                                                                                        ifelse(data_polyfun_sig_snps$SNP == "rs78762097", "THBS2",
                                                                                                                                               ifelse(data_polyfun_sig_snps$SNP == "rs759423419", "SLC13A1",
                                                                                                                                                      ifelse(data_polyfun_sig_snps$SNP == "rs11643815", "MT4",    
                                                                                                                                                             ifelse(data_polyfun_sig_snps$SNP == "rs111594798", "AFAP1L2", #keep the SNP in credible set 1
                                                                                                                                                                    NA)))))))))))))))))))

#plot the effect sizes
#add scaled effect sizes
data_polyfun_sig_snps$scaled_effect_size.x <- scale(data_polyfun_sig_snps$effect_size.x)
data_polyfun_sig_snps$scaled_effect_size.y <- scale(data_polyfun_sig_snps$effect_size.y)

#fit linear regression model
model_results <- data_polyfun_sig_snps |>
  group_by(finemap_group) |>
  do(tidy(lm(scaled_effect_size.y ~ scaled_effect_size.x, data = .), conf.int = TRUE)) |> 
  filter(term == "scaled_effect_size.x") |> 
  select(finemap_group, estimate, conf.low, conf.high, p.value)

annotation_text <- model_results %>%
  mutate(label = paste0(finemap_group, ": Slope = ", round(estimate, 2), 
                        " (", round(conf.low, 2), ", ", round(conf.high, 2), 
                        "), P value = ", format.pval(p.value, digits = 2))) %>%
  pull(label) %>%
  paste(collapse = "\n")

#plot
data_polyfun_sig_snps |>
  ggplot(aes(x = scaled_effect_size.x, y = scaled_effect_size.y, color = finemap_group)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, col = "grey", linetype = "dashed") +
  labs(x = "Metabolic (scaled effect size)", y = "Sensory (scaled effect size)", color = "Phenotype") +
  geom_smooth(aes(group = finemap_group, color = finemap_group), method = "lm", formula = y ~ x, linetype = "dashed") +
  scale_color_manual(values = c("Variants significant in metabolic phenotype" = "steelblue4", "Variants significant in sensory phenotype" = "darkred")) +
  # scale_y_continuous(limits = c(-3,3)) +
  ggrepel::geom_label_repel(aes(label = paste(SNP, "\n", genes)),
                            size = 2.5, segment.size = 0.2, point.padding = 0.3,
                            force = 2, box.padding = 0.5, max.overlaps = 50, fontface = "italic") +
  # annotate("text", x = Inf, y = Inf, label = annotation_text, hjust = 1.1, vjust = 1, 
  #          size = 3, color = "black", fontface = "italic") +
  theme_minimal()

ggplot2::ggsave("/home/projects/hearing_loss/clsaARHL_SA/plots/effect_size_all.tiff",
                height=5, width=10, units='in', dpi=300)

#####################


#plot a miami plot for the finemap results
met_stats_all_finemap_2$finemap_group <- "metabolic" 
sen_stats_all_finemap_2$finemap_group <- "sensory" 

#merge
data_all_finemap <- rbind(met_stats_all_finemap_2, sen_stats_all_finemap_2)
data_all_finemap1 <- filter(data_all_finemap, PIP > 0.05)
data_all_finemap2 <- filter(data_all_finemap, PIP > 1e-5)


#miami plot
##########
#create columns for variant consequence after fine mapping
data_all_finemap2$genes <- ifelse(data_all_finemap2$SNP == "rs6453022"  & data_all_finemap2$finemap_group == "metabolic", "ARHGEF28\n rs6453022\n Pro284Gln", 
                                  ifelse(data_all_finemap2$SNP == "rs76510133" & data_all_finemap2$finemap_group == "metabolic", "EYA4\n rs76510133 ",
                                         ifelse(data_all_finemap2$SNP == "rs4832151"  & data_all_finemap2$finemap_group == "metabolic", "TCF7L1\n rs4832151 ",
                                                ifelse(data_all_finemap2$SNP == "rs7428702"  & data_all_finemap2$finemap_group == "metabolic", "HEG1\n rs7428702 ",
                                                       ifelse(data_all_finemap2$SNP == "rs17008250" & data_all_finemap2$finemap_group == "metabolic", "LINC02994\n rs17008250 ",
                                                              ifelse(data_all_finemap2$SNP == "rs77196213" & data_all_finemap2$finemap_group == "metabolic", "LINC02716\n rs77196213 ",
                                                                     ifelse(data_all_finemap2$SNP == "rs56105727"  & data_all_finemap2$finemap_group == "metabolic", "FGF14\n rs56105727 ",
                                                                            ifelse(data_all_finemap2$SNP == "rs1604887"  & data_all_finemap2$finemap_group == "metabolic", "OTX2-AS1\n rs1604887 ",
                                                                                   ifelse(data_all_finemap2$SNP == "rs200974429"  & data_all_finemap2$finemap_group == "metabolic", "LOC105370982\n rs200974429 ",
                                                                                          ifelse(data_all_finemap2$SNP == "rs72628292" & data_all_finemap2$finemap_group == "metabolic", "FENDRR\n rs72628292 ",
                                                                                                 
                                                                                                 ifelse(data_all_finemap2$SNP == "rs36100160"  & data_all_finemap2$finemap_group == "sensory", "CDYL\n rs36100160 ", 
                                                                                                        ifelse(data_all_finemap2$SNP == "rs34610101"  & data_all_finemap2$finemap_group == "sensory", "CAMKK1\n rs34610101 ",
                                                                                                               ifelse(data_all_finemap2$SNP == "rs36062310"  & data_all_finemap2$finemap_group == "sensory", "KLHDC7B\n rs36062310 \n Val504Met",
                                                                                                                      ifelse(data_all_finemap2$SNP == "rs80122938"  & data_all_finemap2$finemap_group == "sensory", "PDPN\n rs80122938 ",
                                                                                                                             ifelse(data_all_finemap2$SNP == "rs115671365"  & data_all_finemap2$finemap_group == "sensory", "GFOD1\n rs115671365 ",
                                                                                                                                    ifelse(data_all_finemap2$SNP == "rs78762097"  & data_all_finemap2$finemap_group == "sensory", "THBS2\n rs78762097 ",
                                                                                                                                           ifelse(data_all_finemap2$SNP == "rs759423419"  & data_all_finemap2$finemap_group == "sensory", "SLC13A1\n rs759423419",
                                                                                                                                                  ifelse(data_all_finemap2$SNP == "rs11643815"  & data_all_finemap2$finemap_group == "sensory", "MT4\n rs11643815\n Gly48Asp",    
                                                                                                                                                         ifelse(data_all_finemap2$SNP == "rs111594798" & data_all_finemap2$finemap_group == "sensory", "AFAP1L2\n rs111594798", #keep the SNP in vredible set 1
                                                                                                                                                                NA)))))))))))))))))))


data_all_finemap2$conseq <- ifelse(data_all_finemap2$SNP == "rs6453022"  & data_all_finemap2$finemap_group == "metabolic", "Missense variant", 
                                   ifelse(data_all_finemap2$SNP == "rs76510133" & data_all_finemap2$finemap_group == "metabolic", "3 prime UTR variant",
                                          ifelse(data_all_finemap2$SNP == "rs4832151"  & data_all_finemap2$finemap_group == "metabolic", "Intronic variant",
                                                 ifelse(data_all_finemap2$SNP == "rs7428702"  & data_all_finemap2$finemap_group == "metabolic", "Intronic variant",
                                                        ifelse(data_all_finemap2$SNP == "rs17008250" & data_all_finemap2$finemap_group == "metabolic", "Intronic variant",
                                                               ifelse(data_all_finemap2$SNP == "rs77196213" & data_all_finemap2$finemap_group == "metabolic", "Upstream gene variant",
                                                                      ifelse(data_all_finemap2$SNP == "rs200974429"  & data_all_finemap2$finemap_group == "metabolic", "Intronic variant",
                                                                             ifelse(data_all_finemap2$SNP == "rs56105727"  & data_all_finemap2$finemap_group == "metabolic", "Intronic variant",
                                                                                    ifelse(data_all_finemap2$SNP == "rs1604887"  & data_all_finemap2$finemap_group == "metabolic", "Intronic variant",
                                                                                           ifelse(data_all_finemap2$SNP == "rs72628292" & data_all_finemap2$finemap_group == "metabolic", "Intronic variant",
                                                                                                  #sensory phenotype
                                                                                                  ifelse(data_all_finemap2$SNP == "rs36100160"  & data_all_finemap2$finemap_group == "sensory", "Intronic variant", 
                                                                                                         ifelse(data_all_finemap2$SNP == "rs34610101"  & data_all_finemap2$finemap_group == "sensory", "Intronic variant",
                                                                                                                ifelse(data_all_finemap2$SNP == "rs36062310"  & data_all_finemap2$finemap_group == "sensory", "Missense variant",
                                                                                                                       ifelse(data_all_finemap2$SNP == "rs80122938"  & data_all_finemap2$finemap_group == "sensory", "Intergenic variant",
                                                                                                                              ifelse(data_all_finemap2$SNP == "rs115671365"  & data_all_finemap2$finemap_group == "sensory", "Intronic variant",
                                                                                                                                     ifelse(data_all_finemap2$SNP == "rs78762097"  & data_all_finemap2$finemap_group == "sensory", "Regulatory region variant",
                                                                                                                                            ifelse(data_all_finemap2$SNP == "rs759423419"  & data_all_finemap2$finemap_group == "sensory", "Intergenic variant",
                                                                                                                                                   ifelse(data_all_finemap2$SNP == "rs11643815" & data_all_finemap2$finemap_group == "sensory", "Missense variant",
                                                                                                                                                          ifelse(data_all_finemap2$SNP == "rs111594798" & data_all_finemap2$finemap_group == "sensory", "Intronic variant",
                                                                                                                                                                 NA)))))))))))))))))))

##########

#create chromosome labels
labels <- c(1:22)

#color label 
variant_colors <- c(
  "Intronic variant" = "#ed8b00", 
  "Upstream gene variant" = "#7b2d26", 
  "3 prime UTR variant" = "#0f85a0", 
  "Intergenic variant" = "#006820", 
  "Missense variant" = "#dd4124", 
  "Regulatory region variant" = "#00496f"
)
# Upper plot
upper_plot <- data_all_finemap2 |> 
  mutate(CHR = factor(CHR, levels = as.character(1:22))) |> 
  filter(finemap_group == "metabolic") |> 
  ggplot(aes(x = CHR, y = PIP)) +
  geom_point(data = data_all_finemap2 |>  filter(PIP > 0.1 & finemap_group == "metabolic"), size = 3, color = "#C1C1C3", alpha = 0.5) + 
  geom_jitter(data = data_all_finemap2 |>  filter(PIP <= 0.3 & finemap_group == "metabolic"), size = 3, color = "#C1C1C3", alpha = 0.5, width = 0.4)  +
  geom_point(data = data_all_finemap2 |>  filter(PIP > 0.5 & SET == 1 & finemap_group == "metabolic"), size = 3, aes(color = as.factor(conseq))) + 
  geom_abline(intercept = 0.5, slope = 0, col = "red", linetype = "dashed") +
  scale_color_manual(values = variant_colors) + 
  scale_x_discrete(limits = as.character(1:22), position = "bottom") +
  scale_y_continuous(limits = c(0, max(data_all_finemap2$PIP)*1.05), expand = expansion(mult = c(0.02, 0))) +
  labs(x = "", y = "Metabolic (PIP)", color = "Variant Consequence") +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.title.y = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 15),
    plot.margin = margin(t = 10, b = 0, l = 10, r = 10)
  ) + ggrepel::geom_label_repel(data = data_all_finemap2 |> filter(!is.na(genes) & finemap_group == "metabolic"), aes(label = genes),
                                size = 3, segment.size = 0.2, point.padding = 0.3,
                                ylim = c(max(data_all_finemap2$PIP/2), NA), min.segment.length = 0, 
                                force = 2, box.padding = 0.5, fontface = "italic")

# Lower plot
lower_plot <- data_all_finemap2 |> 
  mutate(CHR = factor(CHR, levels = as.character(1:22))) |> 
  filter(finemap_group == "sensory") |> 
  ggplot(aes(x = CHR, y = PIP)) +
  geom_point(data = data_all_finemap2 |>  filter(PIP > 0.1 & finemap_group == "sensory"), size = 3, color = "#C1C1C3", alpha = 0.5) + 
  geom_jitter(data = data_all_finemap2 |>  filter(PIP <= 0.3 & finemap_group == "sensory"), size = 3, color = "#C1C1C3", alpha = 0.5, width = 0.4)  +
  geom_point(data = data_all_finemap2 |>  filter(PIP > 0.5 & SET == 1 & finemap_group == "sensory"), size = 3, aes(color = as.factor(conseq))) + 
  scale_color_manual(values = variant_colors) + 
  geom_abline(intercept = -0.5, slope = 0, col = "red", linetype = "dashed") +
  scale_x_discrete(limits = as.character(1:22), position = "top") +
  scale_y_reverse(limits = c(max(data_all_finemap2$PIP)*1.05, 0), expand = expansion(mult = c(0.02, 0))) +
  labs(x = "", y = "Sensory (PIP)", color = "Variant Consequence") +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.title.y = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 15),
    plot.margin = margin(t = 10, b = 0, l = 10, r = 10)
  ) + ggrepel::geom_label_repel(data = data_all_finemap2 |> filter(!is.na(genes) & finemap_group == "sensory"), aes(label = genes),
                                size = 3, segment.size = 0.2, point.padding = 0.3,
                                ylim = c(NA, -max(data_all_finemap2$PIP/2)), min.segment.length = 0, 
                                force = 2, box.padding = 0.5, fontface = "italic")
# Combine plots
ggarrange(upper_plot, lower_plot, nrow = 2, common.legend = F)
finemap_plot <- ggarrange(upper_plot, lower_plot, nrow = 2)
ggplot2::ggsave("/home/projects/hearing_loss/clsaARHL_SA/finemap_plot.tiff",
                height=3, width=7, units='in', dpi=300)

###########################################################