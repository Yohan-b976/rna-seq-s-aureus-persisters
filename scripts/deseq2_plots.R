#!/usr/bin/env Rscript

# ----------------------------------------------
# Graphiques
# ----------------------------------------------

library(ggplot2)
library(ggrepel)

# ==== Arguments Nextflow ====
args <- commandArgs(trailingOnly = TRUE)
deseq_results_path  <- args[1]
gene_pathway <- args[2]
mapping_path        <- args[3]

# Lire résultats deseq
res <- read.csv(deseq_results_path,row.names =1)

# Lire table kegg 
kegg_all <- read.table(
  gene_pathway,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Lire table de mapping
mapping <- read.table(
  mapping_path,
  header = TRUE,
  sep = "\t",
  quote = "",
  fill = TRUE,
  stringsAsFactors = FALSE,
  check.names = F
)

# On ajoute 'pth' à la main car il n'est pas présent dans la table d'origine
mapping[mapping$product == "peptidyl-tRNA hydrolase", ]$symbol="pth"


# === Préparation du dataframe pour le plot : ===
# Compute le log2base mean pour les abscisses (si diff de 0)
res$log2BaseMean <- ifelse(res$baseMean > 0, log2(res$baseMean), NA)
# colonne significatif / non significatif 
res$is_sig <- !is.na(res$padj) & res$padj < 0.05
res$signif <- ifelse(res$is_sig, "Significant", "Non-Significant")
# ===============================================


# == On merge les dataframes : ==
# On concatène les colonnes du mapping sur les résultats DESeq
res$locus_tag <- rownames(res)
res_annot <- merge(res, mapping,
                   by = "locus_tag",
                   all.x = TRUE)

# On rajoute les colonnes du kegg_file pour les pathways 
colnames(kegg_all)<-c("locus_tag","Pathway_ID")
res_annot <- merge(res_annot, kegg_all,
                   by = "locus_tag",
                   all.x = TRUE)
# ===============================



# === Sélection de pathways ===

# Les gènes en lien avec la traduction sont : 
# sao03011 : Ribosome, sao03009 : Ribosome biogenesis, sao03016 : Transfer RNA biogenesis, sao03012 : Translation factors
# On les trouve dans BRITE

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

# Typical genes 
typical_members = subset(
  res_annot,
  symbol %in% c("pth", "infA", "infB", "infC", "frr", "tsf")
)


# ===== MA-PLOT TRANSLATION GENES =====
plot_df=res_annot[translation_genes, ]

pdf("MA_plot_translationgenes.pdf")
p <- ggplot(
  data = plot_df,
  aes(x = log2BaseMean, y = log2FoldChange, color = signif)
) +
  
  # Tous les gènes de traduction (gris / rouge)
  geom_point(alpha = 0.8, size = 1.2) +
  
  # Cercle noir autour des AA-tRNA synthetases (restreint aux translation_genes)
  geom_point(
    data = subset(plot_df, is_AA_tRNA),
    aes(shape = "AA_tRNA_synthetases"),
    color = "black", size = 1.3, stroke = 1.3
  ) +
  
  # Couleurs gris / rouge
  scale_color_manual(values = c("Non-Significant" = "grey70",
                                "Significant" = "red")) +
  
  # Cercle vide pour les AA-tRNA
  scale_shape_manual(values = c("AA_tRNA_synthetases" = 1)) +
  
  # Axes
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
  
  # Labels des gènes typiques
  geom_text_repel(
    data = typical_members,
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

  # Ligne horizontale
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


# Gestion de la légende
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
#======================================

