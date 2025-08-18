#figure 2
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


setwd("/home/projects/hearing_loss/clsaARHL_SA/gwas")
metabolic_all <- fread(file = "met_better_gwas_sorted_noNA_add.gz", header = TRUE)
sensory_all <- fread(file = "sen_better_gwas_sorted_noNA_add.gz", header = TRUE)

#add phenotype column
metabolic_all <- mutate(metabolic_all, group = "metabolic")
sensory_all <- mutate(sensory_all, group = "sensory")


#select relevant columns and add group column
metabolic_all.1 <- metabolic_all |> 
  select(CHR, BP, SNP, P, group) 

sensory_all.1 <- sensory_all|>
  select(CHR, BP, SNP, P, group) 

#add x chromosome (in hg19)
metabolic_xchr <- fread(file = "/home/projects/archive/previous_projects/xchr/met_xwas_stratFisher_model2.xstrat.linear")
sensory_xchr <- fread(file = "/home/projects/archive/previous_projects/xchr/sen_xwas_stratFisher_model2.xstrat.linear")
#select relevant columns
metabolic_xchr.1 <- metabolic_xchr |>
  select(CHR, BP, SNP, P_comb_Fisher) |>
  dplyr::rename(P = P_comb_Fisher) |>
  mutate(group = "metabolic") 

sensory_xchr.1 <- sensory_xchr |>
  select(CHR, BP, SNP, P_comb_Fisher) |>
  dplyr::rename(P = P_comb_Fisher) |>
  mutate(group = "sensory") 


data_all <- rbind(metabolic_all.1, sensory_all.1, metabolic_xchr.1, sensory_xchr.1)


#get significant snps 
sig_snps_met <- data_all |> 
  filter(group == "metabolic" & P < 5e-8)
sig_snps_met <- sig_snps_met$SNP

sig_snps_sen <- data_all |> 
  filter(group == "sensory" & P < 5e-8)
sig_snps_sen <- sig_snps_sen$SNP


data_all$rs <- ifelse(data_all$SNP == "chr5:73780686:C:A" & data_all$group == "metabolic", "rs6453022\n ARHGEF28\n Pro284Gln",
                      ifelse(data_all$SNP == "chr5:73776529:T:C" & data_all$group == "metabolic", "rs7714670\n ARHGEF28\n Trp225Arg",   
                             ifelse(data_all$SNP == "chrX:74550257:T:G" & data_all$group == "metabolic", "rs895513076\n TERF1P7", 
                                    ifelse(data_all$SNP == "chr22:50549676:G:A" & data_all$group == "sensory", "rs36062310\n KLHDC7B\n Val504Met", NA))))



rsID <- c("rs6453022\n ARHGEF28\n Pro284Gln", "rrs7714670\n ARHGEF28\n Trp225Arg", "rs895513076\n TERF1P7", "rs36062310\n KLHDC7B\n Val504Met")
ggmiami_gwas2(data = data_all, split_by = "group", split_at = "metabolic", chr = "CHR", pos = "BP", p="P",
              chr_colors = NULL, upper_ylab = "Metabolic",
              lower_ylab = "Sensory", upper_chr_colors = c("#9BB3C5", "#6B6B6B"),
              lower_chr_colors = c("#D4A4A4", "#6B6B6B"),suggestive_line = 1e-05, 
              suggestive_line_color = "black", genome_line_color = "red",
              upper_highlight_col = "SNP", upper_highlight = sig_snps_met, upper_highlight_color = "#5ced73",
              lower_highlight_col = "SNP", lower_highlight = sig_snps_sen, lower_highlight_color = "#5ced73",
              top_n_hits = 10, hits_label_col = "rs", hits_label = rsID)




#miami plot for the MAGMA results without x chromosome
met_all_magma <- fread("/home/projects/hearing_loss/clsaARHL_SA/MAGMA_output/met_all/magma.genes.out")
sen_all_magma <- fread("/home/projects/hearing_loss/clsaARHL_SA/MAGMA_output/sen_all/magma.genes.out")

#select relevant columns
met_all_magma <- met_all_magma |> 
  select(SYMBOL, CHR, START, STOP, P)
sen_all_magma <- sen_all_magma |> 
  select(SYMBOL, CHR, START, STOP, P)

#gene-based analysis from XWAS
met_all_xchr_gene <- fread ("/home/projects/hearing_loss/clsaARHL_SA/xchGeneBased/met_gene_test_result_v2.txt.sort.clean")
sen_all_xchr_gene <- fread ("/home/projects/hearing_loss/clsaARHL_SA/xchGeneBased/sen_gene_test_result_v2.txt.sort.clean")
gene_list_magma <- fread ("/home/projects/hearing_loss/clsaARHL_SA/MAGMA_output/genes_position_magma.txt")
gene_list_magma_x <- subset(gene_list_magma, CHR == "X")
gene_list_magma_x$CHR <- as.integer(gsub("X", 23 , gene_list_magma_x$CHR))

#merge the list with gene-based results
met_all_xchr_gene_merge <- merge (gene_list_magma_x, met_all_xchr_gene, by.x = "SYMBOL", by.y = "gene")
met_all_xchr_gene_merge <- met_all_xchr_gene_merge |> 
  select(SYMBOL, CHR, START, STOP, tail_p_value) |> 
  dplyr::rename(P = tail_p_value)
sen_all_xchr_gene_merge <- merge (gene_list_magma_x, sen_all_xchr_gene, by.x = "SYMBOL", by.y = "gene")
sen_all_xchr_gene_merge <- sen_all_xchr_gene_merge |> 
  select(SYMBOL, CHR, START, STOP, tail_p_value) |> 
  dplyr::rename(P = tail_p_value)

#merge x chromosome with data
met_all_magma_to_plot <- rbind(met_all_magma, met_all_xchr_gene_merge)
sen_all_magma_to_plot <- rbind(sen_all_magma, sen_all_xchr_gene_merge)

#add representative position column

met_all_magma_to_plot$BP <- round((met_all_magma_to_plot$START + met_all_magma_to_plot$STOP)/2)
sen_all_magma_to_plot$BP <- round((sen_all_magma_to_plot$START + sen_all_magma_to_plot$STOP)/2)

#add group column
met_all_magma_to_plot$group <- "metabolic"
sen_all_magma_to_plot$group <- "sensory"
#combine the data
data_all_magma_to_plot <- rbind(met_all_magma_to_plot, sen_all_magma_to_plot)

#select the significant genes 
sig_genes_met <- met_all_magma_to_plot$SYMBOL[met_all_magma_to_plot$P < 2.6e-6]
sig_genes_sen <- sen_all_magma_to_plot$SYMBOL[sen_all_magma_to_plot$P < 2.6e-6]

# #remove non-relevant genes
data_all_magma_to_plot$SYMBOL <- gsub("ZNF143", NA , data_all_magma_to_plot$SYMBOL)


ggmiami_gwas2(data = data_all_magma_to_plot, split_by = "group", split_at = "metabolic", chr = "CHR", pos = "BP", p="P",
              chr_colors = NULL, upper_ylab = "Metabolic",
              lower_ylab = "Sensory", upper_chr_colors = c("#9BB3C5", "#6B6B6B"),
              lower_chr_colors = c("#D4A4A4", "#6B6B6B"),suggestive_line = NULL, genome_line = 2.6e-6,
              suggestive_line_color = "black", genome_line_color = "red",
              upper_highlight_col = "SYMBOL", upper_highlight = sig_genes_met, upper_highlight_color = "#5ced73",
              lower_highlight_col = "SYMBOL", lower_highlight = sig_genes_sen, lower_highlight_color = "#5ced73",
              top_n_hits = 4, hits_label_col = "SYMBOL")


ggplot2::ggsave("/home/projects/hearing_loss/clsaARHL_SA/plots/magma_miami.tiff",
                height=20, width=40, units='in', dpi=300)

#####################

####updated miami function
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


