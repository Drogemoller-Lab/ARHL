#get the top SNPs from finemapping results
library(data.table)
library(tidyverse)
# Requires harmonized, hg19-lifted European ancestry (cluster 4) GWAS
# summary statistics generated without adjustment for the opposite phenotype.
# These files are used to obtain effect sizes for the phenotype comparison plots.

##finemap results
#for the genome-wide metabolic
met_stats_all_finemap <- fread("/home/projects/hearing_loss/clsaARHL_SA/polyfun/met_cluster4_polyfun_finemap_agg.txt.gz")
met_stats_all_finemap_2 <- separate(met_stats_all_finemap, col = CREDIBLE_SET, into = c("CHR.x", "POSITION", "SET"), sep = ":")
met_stats_all_finemap_credible <- filter(met_stats_all_finemap_2, PIP >= 0.5 | SNP == "rs6453022"& SET == 1)
met_stats_all_finemap_credible$finemap_group <- "Variants significant in metabolic phenotype"
#merge with gwas results (hg19)
met_cluster4_hg19 <- fread ("/home/projects/hearing_loss/clsaARHL_SA/polyfun_noopp/met_cluster4_lifted_hrm.txt")
met_stats_all_finemap_credible.merge <- merge(met_stats_all_finemap_credible, met_cluster4_hg19, 
                                              by = c("CHR", "BP"))

#for the genome-wide sensory
sen_stats_all_finemap <- fread("/home/projects/hearing_loss/clsaARHL_SA/polyfun/sen_cluster4_polyfun_finemap_agg.txt.gz")
sen_stats_all_finemap_2 <- separate(sen_stats_all_finemap, col = CREDIBLE_SET, into = c("CHR.x", "POSITION", "SET"), sep = ":")
sen_stats_all_finemap_credible <- filter(sen_stats_all_finemap_2,PIP >= 0.5 & SET == 1)

sen_stats_all_finemap_credible$finemap_group <- "Variants significant in sensory phenotype"
#merge with gwas results (hg19)
sen_cluster4_hg19 <- fread ("/home/projects/hearing_loss/clsaARHL_SA/polyfun_noopp/sen_cluster4_lifted_hrm.txt")
sen_stats_all_finemap_credible.merge <- merge(sen_stats_all_finemap_credible, sen_cluster4_hg19, 
                                              by= c("CHR", "BP"))

#############################################
#get a list of significant SNPs
polyfun_all_sig_snps <- rbind(dplyr::select(met_stats_all_finemap_credible, SNP, CHR, BP, PIP, finemap_group),
                              select(sen_stats_all_finemap_credible, SNP, CHR, BP, PIP, finemap_group))

#get the effect sizes of these SNPs from the GWAS results hg19
met_polyfun_sig_snps <- merge(polyfun_all_sig_snps, met_cluster4_hg19, 
                              by = c("CHR", "BP")) 

sen_polyfun_sig_snps <- merge(polyfun_all_sig_snps, sen_cluster4_hg19, 
                              by = c("CHR", "BP")) 

data_polyfun_sig_snps <- merge(met_polyfun_sig_snps, sen_polyfun_sig_snps, by = c("CHR", "BP", "SNP.x", "finemap_group",
                                                                                  "A1", "A2",
                                                                                  "PIP", "sample_size"))
colnames(data_polyfun_sig_snps)[3] <- "SNP"
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
data_polyfun_sig_snps$scaled_beta.x <- as.numeric(scale(data_polyfun_sig_snps$Beta.x))
data_polyfun_sig_snps$scaled_beta.y <- as.numeric(scale(data_polyfun_sig_snps$Beta.y))

#fit linear regression model
model_results <- data_polyfun_sig_snps |>
  group_by(finemap_group) |>
  do(tidy(lm(scaled_beta.y ~ scaled_beta.x, data = .), conf.int = TRUE)) |> 
  filter(term == "scaled_beta.x") |> 
  select(finemap_group, estimate, conf.low, conf.high, p.value)

annotation_text <- model_results %>%
  mutate(label = paste0(finemap_group, ": Slope = ", round(estimate, 2), 
                        " (", round(conf.low, 2), ", ", round(conf.high, 2), 
                        "), P value = ", format.pval(p.value, digits = 2))) %>%
  pull(label) %>%
  paste(collapse = "\n")

#plot
data_polyfun_sig_snps |>
  ggplot(aes(x = scaled_beta.x, y = scaled_beta.y, color = finemap_group)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, col = "grey", linetype = "dashed") +
  labs(x = "Metabolic (standardized effect size)", y = "Sensory (standardized effect size)", color = "Phenotype") +
  geom_smooth(aes(group = finemap_group, color = finemap_group), method = "lm", formula = y ~ x, linetype = "dashed") +
  scale_color_manual(values = c("Variants significant in metabolic phenotype" = "steelblue4", "Variants significant in sensory phenotype" = "darkred")) +
  # scale_y_continuous(limits = c(-3,3)) +
  ggrepel::geom_label_repel(aes(label = paste(SNP, "\n", genes)),
                            size = 2.5, segment.size = 0.2, point.padding = 0.3,
                            force = 2, box.padding = 0.5, max.overlaps = 50, fontface = "italic") +
  # annotate("text", x = Inf, y = Inf, label = annotation_text, hjust = 1.1, vjust = 1, 
  #          size = 3, color = "black", fontface = "italic") +
  theme_minimal()

ggplot2::ggsave("/home/projects/hearing_loss/clsaARHL_SA/plots/effect_size_noopp.tiff",
                height=7, width=10, units='in', dpi=300)

#test the interaction
interaction_model <- lm(scaled_beta.x ~ scaled_beta.y * finemap_group, data = data_polyfun_sig_snps)
summary(interaction_model)

###try all variants
#met_stats_all_finemap <- fread("/home/projects/hearing_loss/clsaARHL_SA/polyfun/met_cluster4_polyfun_finemap_agg.txt.gz")
#met_stats_all_finemap_2 <- separate(met_stats_all_finemap, col = CREDIBLE_SET, into = c("CHR.x", "POSITION", "SET"), sep = ":")
met_stats_all_finemap_credible <- filter(met_stats_all_finemap_2, PIP >= 0.5 | SNP == "rs6453022"& SET == 1)
met_stats_all_finemap_credible$finemap_group <- "Variants significant in metabolic phenotype"
#merge with gwas results (hg19)
#met_cluster4_hg19 <- fread ("/home/projects/hearing_loss/clsaARHL_SA/polyfun_noopp/met_cluster4_lifted_hrm.txt")
met_stats_all_finemap_credible.merge <- merge(met_stats_all_finemap_credible, met_cluster4_hg19, 
                                              by = c("CHR", "BP"))

#for the genome-wide sensory
#sen_stats_all_finemap <- fread("/home/projects/hearing_loss/clsaARHL_SA/polyfun/sen_cluster4_polyfun_finemap_agg.txt.gz")
#sen_stats_all_finemap_2 <- separate(sen_stats_all_finemap, col = CREDIBLE_SET, into = c("CHR.x", "POSITION", "SET"), sep = ":")
sen_stats_all_finemap_credible <- filter(sen_stats_all_finemap_2,PIP >= 0.5 & SET == 1)

sen_stats_all_finemap_credible$finemap_group <- "Variants significant in sensory phenotype"
#merge with gwas results (hg19)
#sen_cluster4_hg19 <- fread ("/home/projects/hearing_loss/clsaARHL_SA/polyfun_noopp/sen_cluster4_lifted_hrm.txt")
sen_stats_all_finemap_credible.merge <- merge(sen_stats_all_finemap_credible, sen_cluster4_hg19, 
                                              by= c("CHR", "BP"))

polyfun_all_sig_snps <- rbind(select(met_stats_all_finemap_credible, SNP, CHR, BP, PIP, finemap_group),
                              select(sen_stats_all_finemap_credible, SNP, CHR, BP, PIP, finemap_group))

#get the effect sizes of these SNPs from the GWAS results hg19
met_polyfun_sig_snps <- merge(polyfun_all_sig_snps, met_cluster4_hg19, 
                              by = c("CHR", "BP")) 

sen_polyfun_sig_snps <- merge(polyfun_all_sig_snps, sen_cluster4_hg19, 
                              by = c("CHR", "BP")) 

data_polyfun_sig_snps <- merge(met_polyfun_sig_snps, sen_polyfun_sig_snps, by = c("CHR", "BP", "SNP.x", "finemap_group",
                                                                                  "A1", "A2",
                                                                                  "PIP", "sample_size"))
colnames(data_polyfun_sig_snps)[3] <- "SNP"


#plot the effect sizes
#add scaled effect sizes
data_polyfun_sig_snps$scaled_beta.x <- as.numeric(scale(data_polyfun_sig_snps$Beta.x))
data_polyfun_sig_snps$scaled_beta.y <- as.numeric(scale(data_polyfun_sig_snps$Beta.y))

#scale errors
data_polyfun_sig_snps$scaled_SE.x <- as.numeric(scale(data_polyfun_sig_snps$SE.x))
data_polyfun_sig_snps$scaled_SE.y <- as.numeric(scale(data_polyfun_sig_snps$SE.y))

#fit linear regression model
model_results2 <- data_polyfun_sig_snps |>
  group_by(finemap_group) |>
  do(tidy(lm(scaled_beta.x ~ scaled_beta.y, data = .), conf.int = TRUE)) |> 
  filter(term == "scaled_beta.y") |> 
  select(finemap_group, estimate, conf.low, conf.high, p.value)

annotation_text <- model_results2 %>%
  mutate(label = paste0(finemap_group, ": Slope = ", round(estimate, 2), 
                        " (", round(conf.low, 2), ", ", round(conf.high, 2), 
                        "), P value = ", format.pval(p.value, digits = 2))) %>%
  pull(label) %>%
  paste(collapse = "\n")

#test the interaction
interaction_model <- lm(scaled_beta.y ~ scaled_beta.x * finemap_group, data = data_polyfun_sig_snps)
summary(interaction_model)
data_polyfun_sig_snps |>
  ggplot(aes(x = scaled_beta.x, y = scaled_beta.y, color = finemap_group)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 0, col = "grey", linetype = "dashed") +
  geom_vline(xintercept = 0, col = "grey", linetype = "dashed") +
  labs(x = "Metabolic (standardized effect size)", y = "Sensory (standardized effect size)", color = "Phenotype") +
  geom_smooth(aes(group = finemap_group, color = finemap_group), method = "lm", formula = y ~ x, linetype = "dashed", se = F) +
  scale_color_manual(values = c("Variants significant in metabolic phenotype" = "steelblue4", "Variants significant in sensory phenotype" = "darkred")) +
  ggrepel::geom_label_repel(aes(label = paste(SNP, "\n", genes)),
                            size = 2.5, segment.size = 0.2, point.padding = 0.3,
                            force = 2, box.padding = 0.5, max.overlaps = 50, fontface = "italic") +
  geom_errorbar(
    data = subset(data_polyfun_sig_snps, finemap_group == "Variants significant in sensory phenotype"),
    aes(ymin = scaled_beta.y - scaled_SE.y,
        ymax = scaled_beta.y + scaled_SE.y),
    width = 0.05,
    alpha = 0.3
  ) +
  geom_errorbarh(
    data = subset(data_polyfun_sig_snps, finemap_group == "Variants significant in metabolic phenotype"),
    aes(xmin = scaled_beta.x - scaled_SE.x,
        xmax = scaled_beta.x + scaled_SE.x),
    height = 0.05,
    alpha = 0.3
  ) +
  theme_minimal()

ggplot2::ggsave("/home/projects/hearing_loss/clsaARHL_SA/plots/effect_size_noopp_errBars.tiff",
                height=7, width=10, units='in', dpi=300)

data_polyfun_sig_snps %>%
  group_by(finemap_group) %>%
  summarize(
    correlation = cor(scaled_beta.x, scaled_beta.y),
  )

