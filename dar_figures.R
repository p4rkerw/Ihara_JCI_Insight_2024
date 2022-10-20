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
dar.gr <- read.xlsx(file, rowNames = TRUE) %>%
  rownames_to_column(var = "peak") %>%
  tidyr::separate(peak, sep = "-", into = c("seqnames","start","end")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

db <- c("hg38_genes_promoters","hg38_enhancers_fantom","hg38_genes_exons", "hg38_genes_introns", "hg38_genes_intergenic")
# db <- "hg38_basicgenes"
anno <- build_annotations(genome = 'hg38', annotations=db)

# annotate peaks
dar.anno <- annotate_regions(dar.gr, annotations=anno, ignore.strand=TRUE) %>%
  as.data.frame()
dar.anno <- dar.anno %>%
  dplyr::mutate(peak = paste0(seqnames,"-",start,"-",end)) %>%
  dplyr::distinct(peak, annot.type)

# prioritize annotations using promoters > enhancers > genes > intergenic
dar.anno$annot.type <- factor(dar.anno$annot.type, levels = db)
dar.anno <- dar.anno %>%
  group_by(peak) %>%
  arrange(annot.type) %>%
  slice(1)

# join annotation to df
dar <- as.data.frame(dar.gr) %>%
  dplyr::mutate(peak = paste0(seqnames,"-",start,"-",end)) %>%
  left_join(dar.anno, by = "peak")

# intersect 
dar.genes <- dar[dar$gene %in% genes,]

# visualize
dar.genes %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::mutate(logpval = -log10(p_val_adj)) %>%
  dplyr::mutate(logpval = ifelse(logpval > 100, 100, logpval)) %>%
  dplyr::mutate(fold_change = 2^avg_log2FC) %>%
  dplyr::mutate(label = paste0(gene)) %>%
  ggplot(aes(avg_log2FC, logpval, label=label, color=annot.type)) +
  geom_point() +
  geom_text_repel(show.legend = FALSE, max.overlaps=20) +
  xlim(c(-0.15,0.25)) +
  xlab("Average log-fold change for PT_VCAM1 vs PCT") +
  ylab("-log10(p_val_adj)") + 
  ggtitle("Differentially accessible regions in PT_VCAM1 vs PCT") +
  theme_bw() +
  theme(legend.title = element_blank())

# intersect 
dar.genes <- dar[dar$gene %in% apoptosis_genes,]
dar.genes %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::mutate(logpval = -log10(p_val_adj)) %>%
  dplyr::mutate(logpval = ifelse(logpval > 100, 100, logpval)) %>%
  dplyr::mutate(fold_change = 2^avg_log2FC) %>%
  dplyr::mutate(label = paste0(gene)) %>%
  ggplot(aes(avg_log2FC, logpval, label=label, color=annot.type)) +
  geom_point() +
  geom_text_repel(show.legend = FALSE, max.overlaps=20) +
  xlim(c(-0.2,0.25)) +
  xlab("Average log-fold change for PT_VCAM1 vs PCT") +
  ylab("-log10(p_val_adj)") + 
  ggtitle("Differentially accessible regions near apoptosis genes in PT_VCAM1 vs PCT") +
  theme_bw() +
  theme(legend.title = element_blank())
