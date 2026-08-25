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

###########################################
## ggmiami function
ggmiami_gwas2 <-function (data, split_by, split_at, chr = "chr", pos = "pos", 
                         p = "p", chr_colors = c("black", "grey"), upper_chr_colors = NULL, 
                         lower_chr_colors = NULL, upper_ylab = "-log10(p)", lower_ylab = "-log10(p)", 
                         genome_line = 5e-08, genome_line_color = "red", suggestive_line = 1e-05, 
                         suggestive_line_color = "blue", hits_label_col = NULL, hits_label = NULL, 
                         top_n_hits = 5, upper_labels_df = NULL, lower_labels_df = NULL, 
                         upper_highlight = NULL, upper_highlight_col = NULL, upper_highlight_color = "green", 
                         lower_highlight = NULL, lower_highlight_col = NULL, lower_highlight_color = "green") 
{
  plot_data <- prep_miami_data(data = data, split_by = split_by, 
                               split_at = split_at, chr = chr, pos = pos, p = p)
  # Create a vector for x-axis labels with chromosome 23 labeled as "X"
  custom_labels <- plot_data$axis$chr
  custom_labels[custom_labels == 23] <- "X"
  
  if (all(!is.null(chr_colors), any(!is.null(upper_chr_colors), 
                                    !is.null(lower_chr_colors)))) {
    stop("You have specified both chr_colors and upper_chr_colors and/or\n         lower_chr_colors. This package does not know how to use both\n         information simultaneously. Please only use one method for coloring:\n         either chr_colors, for making upper and lower plot have the same\n         colors, or upper_chr_colors + lower_chr_colors for specifying different\n         colors for upper and lower plot.")
  }
  if (upper_ylab == "-log10(p)") {
    upper_ylab <- expression("-log"[10] * "(p)")
  }
  else {
    upper_ylab <- bquote(atop(.(upper_ylab), "-log"[10] * 
                                "(p)"))
  }
  if (lower_ylab == "-log10(p)") {
    lower_ylab <- expression("-log"[10] * "(p)")
  }
  else {
    lower_ylab <- bquote(atop(.(lower_ylab), "-log"[10] * 
                                "(p)"))
  }
  upper_plot <- ggplot2::ggplot() + ggplot2::geom_point(data = plot_data$upper, 
                                                        aes(x = .data$rel_pos, y = .data$logged_p, color = as.factor(.data$chr)), 
                                                        size = 3) + ggplot2::scale_x_continuous(labels = custom_labels, 
                                                                                                breaks = plot_data$axis$chr_center, 
                                                                                                expand = ggplot2::expansion(mult = 0.01), 
                                                                                                guide = ggplot2::guide_axis(check.overlap = TRUE)) + 
    ggplot2::scale_y_continuous(limits = c(0, 12), 
                                expand = ggplot2::expansion(mult = c(0.02, 0))) + 
    ggplot2::labs(x = "", y = upper_ylab) + ggplot2::theme_classic() + 
    ggplot2::theme(legend.position = "none", axis.title.y = ggplot2::element_text(size = 20), 
                   axis.title.x = ggplot2::element_blank(), axis.text = element_text(size = 15),
                   plot.margin = ggplot2::margin(t = 10, b = 0, l = 10, r = 10))
  lower_plot <- ggplot2::ggplot() + ggplot2::geom_point(data = plot_data$lower, 
                                                        aes(x = .data$rel_pos, y = .data$logged_p, color = as.factor(.data$chr)), 
                                                        size = 3) + ggplot2::scale_x_continuous(breaks = plot_data$axis$chr_center,
                                                                                                position = "top", expand = ggplot2::expansion(mult = 0.01)) + 
    ggplot2::scale_y_reverse(limits = c(12, 0), 
                             expand = ggplot2::expansion(mult = c(0, 0.02))) + 
    ggplot2::labs(x = "", y = lower_ylab) + ggplot2::theme_classic() + 
    ggplot2::theme(legend.position = "none", axis.title.y = ggplot2::element_text(size = 20), 
                   axis.text.x = ggplot2::element_blank(), axis.text = element_text(size = 15),
                   axis.title.x = ggplot2::element_blank(), plot.margin = ggplot2::margin(t = 0, 
                                                                                          b = 10, l = 10, r = 10))
  if (all(!is.null(chr_colors), is.null(upper_chr_colors), 
          is.null(lower_chr_colors))) {
    if (length(chr_colors) == 2) {
      chr_colors <- rep(chr_colors, length.out = nrow(plot_data$axis))
    }
    else if (length(chr_colors) == nrow(plot_data$axis)) {
      chr_colors <- chr_colors
    }
    else {
      stop("The number of colors specified in {chr_colors} does not match the\n         number of chromosomes to be displayed.")
    }
    upper_plot <- upper_plot + ggplot2::scale_color_manual(values = chr_colors)
    lower_plot <- lower_plot + ggplot2::scale_color_manual(values = chr_colors)
  }
  else if (all(is.null(chr_colors), !is.null(upper_chr_colors), 
               !is.null(lower_chr_colors))) {
    if (length(upper_chr_colors) == 2) {
      upper_chr_colors <- rep(upper_chr_colors, length.out = nrow(plot_data$axis))
    }
    else if (length(upper_chr_colors) == nrow(plot_data$axis)) {
      upper_chr_colors <- upper_chr_colors
    }
    else {
      stop("The number of colors specified in {upper_chr_colors} does not match\n           the number of chromosomes to be displayed.")
    }
    if (length(lower_chr_colors) == 2) {
      lower_chr_colors <- rep(lower_chr_colors, length.out = nrow(plot_data$axis))
    }
    else if (length(lower_chr_colors) == nrow(plot_data$axis)) {
      lower_chr_colors <- lower_chr_colors
    }
    else {
      stop("The number of colors specified in {lower_chr_colors} does not match\n           the number of chromosomes to be displayed.")
    }
    upper_plot <- upper_plot + ggplot2::scale_color_manual(values = upper_chr_colors)
    lower_plot <- lower_plot + ggplot2::scale_color_manual(values = lower_chr_colors)
  }
  else if (all(is.null(chr_colors), any(is.null(upper_chr_colors), 
                                        is.null(lower_chr_colors)))) {
    stop("It looks like you've specified one of upper or lower chr colors\n         without specifying the other. This package needs both colors, unless\n         you want the upper and lower plot to have the same colors, which is\n         done using {chr_colors}.")
  }
  if (!is.null(suggestive_line)) {
    upper_plot <- upper_plot + ggplot2::geom_hline(yintercept = -log10(suggestive_line), 
                                                   color = suggestive_line_color, linetype = "dotted", 
                                                   linewidth = 0.5)
    lower_plot <- lower_plot + ggplot2::geom_hline(yintercept = -log10(suggestive_line), 
                                                   color = suggestive_line_color, linetype = "dotted", 
                                                   linewidth = 0.5)
  }
  if (!is.null(genome_line)) {
    upper_plot <- upper_plot + ggplot2::geom_hline(yintercept = -log10(genome_line), 
                                                   color = genome_line_color, linetype = "dashed", linewidth = 0.7)
    lower_plot <- lower_plot + ggplot2::geom_hline(yintercept = -log10(genome_line), 
                                                   color = genome_line_color, linetype = "dashed", linewidth = 0.7)
  }
  if (all(!is.null(hits_label_col), any(!is.null(upper_labels_df), 
                                        !is.null(lower_labels_df)))) {
    stop("You have specified both hits_label_col and a *_labels_df. This\n         package does not know how to use both information simultaneously.\n         Please only use one method for labelling: either hits_label_col (with\n         or without hits_label), or *_labels_df.")
  }
  if (all(!is.null(hits_label_col), is.null(upper_labels_df), 
          is.null(lower_labels_df))) {
    upper_labels_df <- make_miami_labels(data = plot_data$upper, 
                                         hits_label_col = hits_label_col, hits_label = hits_label, 
                                         top_n_hits = top_n_hits)
    lower_labels_df <- make_miami_labels(data = plot_data$lower, 
                                         hits_label_col = hits_label_col, hits_label = hits_label, 
                                         top_n_hits = top_n_hits)
    upper_plot <- upper_plot + ggrepel::geom_label_repel(data = upper_labels_df, 
                                                         aes(x = .data$rel_pos, y = .data$logged_p, label = .data$label), 
                                                         size = 5, segment.size = 0.2, point.padding = 0.3, 
                                                         ylim = c(plot_data$maxp/2, NA), min.segment.length = 0, 
                                                         force = 2, box.padding = 0.5, fontface = "italic")
    lower_plot <- lower_plot + ggrepel::geom_label_repel(data = lower_labels_df, 
                                                         aes(x = .data$rel_pos, y = .data$logged_p, label = .data$label), 
                                                         size = 5, segment.size = 0.2, point.padding = 0.3, 
                                                         ylim = c(NA, -(plot_data$maxp/2)), min.segment.length = 0, 
                                                         force = 2, box.padding = 0.5, fontface = "italic")
  }
  if (all(is.null(hits_label_col), !is.null(upper_labels_df))) {
    checkmate::assertNames(colnames(upper_labels_df), identical.to = c("rel_pos", 
                                                                       "logged_p", "label"))
    upper_plot <- upper_plot + ggrepel::geom_label_repel(data = upper_labels_df, 
                                                         aes(x = .data$rel_pos, y = .data$logged_p, label = .data$label), 
                                                         size = 2, segment.size = 0.2, point.padding = 0.3, 
                                                         ylim = c(plot_data$maxp/2, NA), min.segment.length = 0, 
                                                         force = 2, box.padding = 0.5, fontface = "italic",
                                                         max.overlaps = 50)
  }
  if (all(is.null(hits_label_col), !is.null(lower_labels_df))) {
    checkmate::assertNames(colnames(lower_labels_df), identical.to = c("rel_pos", 
                                                                       "logged_p", "label"))
    lower_plot <- lower_plot + ggrepel::geom_label_repel(data = lower_labels_df, 
                                                         aes(x = .data$rel_pos, y = .data$logged_p, label = .data$label), 
                                                         size = 2, segment.size = 0.2, point.padding = 0.3, 
                                                         ylim = c(NA, -(plot_data$maxp/2)), min.segment.length = 0, 
                                                         force = 2, box.padding = 0.5, fontface = "italic",
                                                         max.overlaps = 50)
  }
  if (all(!is.null(upper_highlight), !is.null(upper_highlight_col))) {
    upper_highlight_df <- highlight_miami(data = plot_data$upper, 
                                          highlight = upper_highlight, highlight_col = upper_highlight_col, 
                                          highlight_color = upper_highlight_color)
    upper_plot <- upper_plot + ggplot2::geom_point(data = upper_highlight_df, 
                                                   aes(x = .data$rel_pos, y = .data$logged_p), color = upper_highlight_df$color, 
                                                   size = 3)
  }
  if (all(!is.null(lower_highlight), !is.null(lower_highlight_col))) {
    lower_highlight_df <- highlight_miami(data = plot_data$lower, 
                                          highlight = lower_highlight, highlight_col = lower_highlight_col, 
                                          highlight_color = lower_highlight_color)
    lower_plot <- lower_plot + ggplot2::geom_point(data = lower_highlight_df, 
                                                   aes(x = .data$rel_pos, y = .data$logged_p), color = lower_highlight_df$color, 
                                                   size = 3)
  }
  gridExtra::grid.arrange(upper_plot, lower_plot, nrow = 2)
}


#####
