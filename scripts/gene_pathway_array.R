#!/usr/bin/env Rscript

# ----------------------------------------------
# Création de la table de gène
# ----------------------------------------------

library(KEGGREST)

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]



# Associations gènes - pathways 
gene_pathway= keggLink(target = "brite", source = "sao") 

# Retirer les préfixes 'sao' et 'br'
Gene_ID = sub(pattern = "^sao:",  replacement = "", names(gene_pathway)) 
Pathway_ID = sub(pattern = "^br:",  replacement = "", gene_pathway) 


gene_pathway= data.frame(Gene_ID, Pathway_ID) 
gene_pathway = aggregate(. ~ Gene_ID, data = gene_pathway, FUN = function(x) paste(unique(x), collapse = ";"))

# Enregistrement de la table
write.table(gene_pathway, file = file.path(out_dir,"gene_pathway.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

