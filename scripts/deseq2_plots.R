#!/usr/bin/env Rscript

# ----------------------------------------------
# Graphiques
# ----------------------------------------------
if (!requireNamespace("scales", quietly = TRUE)) {
    install.packages("scales", repos = "https://cloud.r-project.org")
}
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(scales)


# ==== Arguments Nextflow ====
args <- commandArgs(trailingOnly = TRUE)
deseq_results_path  <- args[1]
gene_pathway <- args[2]
link_path        <- args[3]

# -----------Chargement des données-----------
res <- read.csv(deseq_results_path,row.names = "gene_id")
rownames(res) <- sub("^gene-", "", rownames(res))

kegg_all <- read.table(
  gene_pathway,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

link <- read.table(
  link_path,
  header = TRUE,
  sep = "\t",
  quote = "",
  fill = TRUE,
  stringsAsFactors = FALSE,
  check.names = F
)


link[link$product == "peptidyl-tRNA hydrolase", ]$symbol="pth"


# Ajout du mask sur les données
res$log2BaseMean <- ifelse(res$baseMean > 0, log2(res$baseMean), NA)
res$is_sig <- !is.na(res$padj) & res$padj < 0.05
res$signif <- ifelse(res$is_sig, "Significant", "Non-Significant")

# -----------Premier plot-----------

limit <- 4


res_all <- res %>%
  mutate(
    # Clamp log2FC entre -limit et +limit
    log2FC_clamped = pmax(pmin(log2FoldChange, limit), -limit),

    shape_point = case_when(
      log2FoldChange < -limit ~ 25,   # triangle pointe vers le bas (dépassement négatif)
      log2FoldChange > limit  ~ 24,   # triangle pointe vers le haut (dépassement positif)
      TRUE                    ~ 16    # point rond normal
    )
  )


pdf("MA_plot_all_genes.pdf")

p <- ggplot(res_all, aes(
  x = baseMean,
  y = log2FC_clamped,
  color = signif,
  fill  = signif
)) +
  geom_point(
    aes(shape = factor(shape_point)),  
    size = 2,
    stroke = 0.8,
    alpha = 0.7
  ) +
  
  scale_shape_manual(values = c(
    "16" = 16,   # rond normal
    "24" = 24,   # triangle normal plein
    "25" = 25    # triangle inversé plein
  )) +

  scale_color_manual(values = c("Non-Significant" = "grey50", "Significant" = "red")) +
  scale_fill_manual(values = c("Non-Significant" = "grey50", "Significant" = "red")) +

  scale_y_continuous(limits = c(-limit, limit)) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +

  labs(
    x = "Mean Normalized Count",
    y = "Log2 Fold Change"
  ) +
  theme_minimal() +
  theme(legend.position = "none")   # légende désactivée

print(p)
dev.off()



# -----------Second plot-----------

# Préparation des données et identifications des gènes liés à la traduction

res$locus_tag <- rownames(res)
res_annot <- merge(res, link, by = "locus_tag", all.x = TRUE)

colnames(kegg_all)<- c("locus_tag","Pathway_ID")
res_annot <- merge(res_annot, kegg_all, by = "locus_tag", all.x = TRUE)

# Sélection de pathways 

# Gènes de traductions : 
# sao03011 : Ribosome,
#sao03009 : Ribosome biogenesis,
# sao03016 : Transfer RNA biogenesis,
# sao03012 : Translation factors

translation_genes = !is.na(res_annot$Pathway_ID) &
  (grepl("(^|;)sao03011(;|$)", res_annot$Pathway_ID) |
     grepl("(^|;)sao03009(;|$)", res_annot$Pathway_ID) |
     grepl("(^|;)sao03016(;|$)", res_annot$Pathway_ID) |
     grepl("(^|;)sao03012(;|$)", res_annot$Pathway_ID))


# Filtre uniquement pour les AA_tRNA_synthetases
AA_tRNA_synthetases = !is.na(res_annot$Pathway_ID) &
  grepl("(^|;)sao03016(;|$)", res_annot$Pathway_ID) &
  grepl("-tRNA synthetase", res_annot$product)

res_annot$is_AA_tRNA <- AA_tRNA_synthetases


tRNAsyn = subset( res_annot, symbol %in% c("infA", "infB", "frr", "infC", "tsf", "pth"))


# Plot

plot_df=res_annot[translation_genes, ]

pdf("MA_plot_of_genes_related_to_translation.pdf")
p <- ggplot(
  data = plot_df,
  aes(x = log2BaseMean, y = log2FoldChange, color = signif)
) +
  
  # couleurs (gris / rouge)
  geom_point(alpha = 0.8, size = 1.2) +
  
  geom_point(
    data = subset(plot_df, is_AA_tRNA),
    aes(shape = "AA_tRNA_synthetases"),
    color = "black", size = 1.3, stroke = 1.3
  ) +
  

  scale_color_manual(values = c("Non-Significant" = "grey50", "Significant" = "red")) +
  

  # Cercle vide pour les AA-tRNA
  scale_shape_manual(values = c("AA_tRNA_synthetases" = 1)) +
  
  # Axes du plot
  scale_x_continuous(
    limits = c(0, 20),
    breaks = seq(0, 20, by = 2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(-6, 5),
    breaks = seq(-6, 5, by = 1),
    expand = c(0, 0)
  )+
  

  geom_text_repel(
    data = tRNAsyn,
    aes(x = log2BaseMean, y = log2FoldChange, label = symbol),
    inherit.aes = FALSE,       
    fontface = "italic",
    size = 4,
    segment.color = "black",
    segment.size = 1,
    min.segment.length = 0,
    box.padding = 1.6,
    point.padding = 0
  )+

  # Dashed line
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  labs(
    x = expression(log[2]~"BaseMean"),
    y = expression(log[2]~"Fold Change"),
    color = NULL,
    shape = NULL,
    title = ""
  ) +
  
  theme_classic() +
  theme(
  panel.background = element_rect(fill = NA, color = NA),
  panel.border = element_rect(fill = NA, color = "black", size = 1)
) 


# Légende
p <- p + theme(
  legend.box = "horizontal"
)

p <- p + guides(
  color = guide_legend(order = 1),
  shape = guide_legend(order = 2)
)

p <- p + theme(
  legend.position = c(0.32, 0.08),
  legend.background = element_blank(),
  legend.box.background = element_blank(),
  legend.key = element_blank()
)
print(p)
dev.off()

