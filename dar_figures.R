library(openxlsx)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(annotatr)
library(org.Hs.eg.db)
library(plyranges)

# retrieve hallmark apoptosis genes
hallmark <- msigdbr(species = "Homo sapiens", category = "H") 
apoptosis <- hallmark[grepl("APOPTOSIS",hallmark$gs_name),]
apoptosis_genes <- apoptosis$gene_symbol
genelist <- c(apoptosis_genes, genes) %>% unique()

file <- "G:/diabneph/analysis/dkd/markers/dar.macs2.PCT_vs_PTVCAM1.markers.xlsx"
dar <- read.xlsx(file, rowNames = TRUE)

db <- c("hg38_genes_promoters","hg38_genes_enhancers")
anno <- build_annotations(genome = 'hg38', annotations=db)
dar.gr <- as_granges(dar)

# annotate peaks
dar <- annotate_regions(dar.gr, annotations=anno, ignore.strand=TRUE) %>%
  as.data.frame()

# intersect 
dar.genes <- dar[dar$gene %in% genes,]

dar.genes.filter <- dar.genes %>%
  dplyr::filter(p_val_adj < 0.05)

# visualize
dar.genes %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::mutate(fold_change = 2^avg_log2FC) %>%
  dplyr::mutate(label = paste0(gene)) %>%
  ggplot(aes(avg_log2FC, -log10(p_val_adj), label=label)) +
  geom_point() +
  geom_text_repel() +
  xlim(c(-0.25,0.25)) +
  xlab("Average log-fold change for PT_VCAM1 vs PCT") +
  ggtitle("Differentially accessible regions in PT_VCAM1 vs PCT", subtitle = "Adjusted p-value < 0.05") +
  theme_bw()

# intersect 
dar.genes <- dar[dar$gene %in% apoptosis_genes,]
dar.genes %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::mutate(fold_change = 2^avg_log2FC) %>%
  dplyr::mutate(label = paste0(gene)) %>%
  ggplot(aes(avg_log2FC, -log10(p_val_adj), label=label)) +
  geom_point() +
  geom_text_repel() +
  xlim(c(-0.25,0.25)) +
  xlab("Average log-fold change for PT_VCAM1 vs PCT") +
  ggtitle("Differentially accessible regions near apoptosis genes in PT_VCAM1 vs PCT", subtitle = "Adjusted p-value < 0.05") +
  theme_bw()
