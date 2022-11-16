library(openxlsx)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)

xl <- read.xlsx("G:/downloads/Joslin_New_46_prots.xlsx", sheet = "New_46_prots")
genes <- xl$Gene


file <- "G:/diabneph/analysis/dkd/markers/deg.celltype.diab_vs_ctrl.xlsx"
idents <- getSheetNames(file)
deg <- lapply(idents, function(ident){
  df <- read.xlsx(file, sheet = ident, rowNames = TRUE) %>%
    rownames_to_column(var = "gene")
  df$celltype <- ident
  return(df)
}) %>% bind_rows()
deg$celltype <- as.factor(deg$celltype)
levels(deg$celltype) <- unique(deg$celltype)


# intersect 
deg.genes <- deg[deg$gene %in% genes,]

# create a color variable and maintain cell type levels
deg.genes <- deg.genes %>%
  dplyr::mutate(color = factor(ifelse(p_val_adj < 0.05, as.character(celltype), "NS"))) %>%
  dplyr::mutate(fold_change = 2^avg_log2FC) %>%
  dplyr::mutate(label = ifelse(p_val_adj < 0.05, paste0(gene,"_",celltype), ""))
deg.genes$color <- factor(deg.genes$color, levels = c(levels(deg.genes$celltype),"NS"))

# plot the background NS degs
toplot1 <- deg.genes %>%
  dplyr::filter(color != "NS")
toplot2 <- deg.genes %>%
  dplyr::filter(color == "NS")

top_layer <- ggplot(data = toplot1, aes(x=avg_log2FC, y=-log10(p_val_adj), color=color, label=label)) + 
  geom_point() +
  geom_text_repel(show.legend = FALSE) +
  theme_bw() +
  xlab("Average log-fold change for DKD vs. Control") +
  ggtitle("Differentially expressed biomarker genes in DKD by cell type", subtitle = "Adjusted p-value < 0.05")

p1 <- top_layer + 
  geom_point(data = toplot2, aes(x=avg_log2FC, y=-log10(p_val_adj)), color = "gray") +
  theme(legend.title = element_blank())
  
p1$layers <- rev(p1$layers)
p1

######################################################################################################
library(openxlsx)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)

xl <- read.xlsx("G:/downloads/Joslin_New_46_prots.xlsx", sheet = "New_46_prots")
genes <- xl$Gene

hallmark <- msigdbr(species = "Homo sapiens", category = "H") 
apoptosis <- hallmark[grepl("APOPTOSIS",hallmark$gs_name),]
apoptosis_genes <- apoptosis$gene_symbol
genelist <- c(apoptosis_genes, genes) %>% unique()

# correlate expression of biomarkers with hallmark apoptosis genes
library(Seurat)
rnaAggr <- readRDS("G:/diabneph/analysis/dkd/rna_aggr_prep/step2_anno.rds")
DefaultAssay(rnaAggr) <- "RNA"
rnaAggr <- NormalizeData(rnaAggr)
counts <- GetAssayData(rnaAggr, slot = "data")
apoptosis_counts <- counts[rownames(counts) %in% apoptosis_genes, rnaAggr@meta.data$celltype %in% c("PTVCAM1")]
apoptosis_totals <- colSums(apoptosis_counts) %>% as.data.frame %>% t()
rownames(apoptosis_totals) <- "apoptosis"

counts <- counts[rownames(counts) %in% genes, rnaAggr@meta.data$celltype %in% c("PTVCAM1")]
allcounts <- rbind(counts, apoptosis_totals) %>%

library(Hmisc)
cor <- rcorr(t(as.matrix(allcounts)))

# pearson correlation coefficients
cor.df <- cor$r %>% as.data.frame()
cor.df$gene <- rownames(cor.df) 

# retrieve pval for correlation between gene and apoptosis counts
cor.df$pval <- cor$P[colnames(cor$P) == "apoptosis"]
cor.df <- cor.df[,c("apoptosis","gene","pval")]
cor.df <- cor.df %>% dplyr::arrange(gene)
cor.df$gene <- as.factor(cor.df$gene)

# prepare for plots
toplot <- melt(cor.df)
toplot %>% 
  dplyr::filter(variable %in% "apoptosis", gene %in% genes) %>%
  ggplot(aes(variable, gene, fill = value)) +
    geom_tile() + 
    scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0,
                       guide = "colorbar")


# avexp <- AverageExpression(rnaAggr)$RNA

# # subset for biomarkers and hallmark apoptosis genes
# library(reshape2)
# genelist <- c(apoptosis_genes, genes) %>% unique()
# mat <- avexp[rownames(avexp) %in% genelist,] %>% as.data.frame()
# # mat$gene <- rownames(mat)
# # mat <- melt(mat)
# # ggplot(mat, aes(gene, variable, fill = value)) + geom_tile()

# cor.exp <- cor(as.data.frame(t(mat))) %>% as.data.frame()

# # put the biomarkers on x-axis and the apoptosis genes on y-axis
# cor.exp <- cor.exp[rownames(cor.exp) %in% genes, colnames(cor.exp) %in% apoptosis_genes]

# cor.exp$gene <- rownames(cor.exp)
# cor.exp <- melt(cor.exp)
# cor.exp <- cor.exp %>% dplyr::arrange(desc(value))
# cor.exp$gene <- as.factor(cor.exp$gene)

# ggplot(cor.exp, aes(gene, variable, fill = value)) + geom_tile()

#####################################################################################################
file <- "G:/diabneph/analysis/dkd/markers/deg.PT_vs_PTVCAM1.markers.xlsx"
deg <- read.xlsx(file, rowNames = TRUE) %>%
    rownames_to_column(var = "gene")

# intersect 
deg.genes <- deg[deg$gene %in% genes,]

deg.genes.filter <- deg.genes %>%
  dplyr::filter(p_val_adj < 0.05)

# visualize
deg.genes %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::mutate(fold_change = 2^avg_log2FC) %>%
  dplyr::mutate(label = gene) %>%
  ggplot(aes(avg_log2FC, -log10(p_val_adj), label=label)) +
  geom_point() +
  geom_text_repel() +
  xlim(c(-1,1)) +
  xlab("Average log-fold change for PT_VCAM1 vs PT") +
  ggtitle("Differentially expressed genes in PT_VCAM1 vs PT", subtitle = "Adjusted p-value < 0.05") +
  theme_bw()

deg.genes %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::mutate(label = gene) %>%
  dplyr::mutate(fold_change = 2^avg_log2FC) %>%
  ggplot(aes(avg_log2FC, -log10(p_val_adj), label=label)) +
  geom_point() +
  geom_text_repel() +
  xlim(c(-2,2)) +
  xlab("Average log-fold change for PT_VCAM1 vs PT") +
  ggtitle("Differentially expressed genes in PT_VCAM1 vs PT", subtitle = "Unadjusted p-value < 0.05")

# intersect with hallmark apoptosis genes
deg.genes <- deg[deg$gene %in% apoptosis_genes,]

deg.genes.filter <- deg.genes %>%
  dplyr::filter(p_val_adj < 0.05)

# visualize
deg.genes %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::mutate(fold_change = 2^avg_log2FC) %>%
  dplyr::mutate(label = gene) %>%
  ggplot(aes(avg_log2FC, -log10(p_val_adj), label=label)) +
  geom_point() +
  geom_text_repel() +
  xlim(c(-1,1)) +
  xlab("Average log-fold change for PT_VCAM1 vs PT") +
  ggtitle("Differentially expressed genes in PT_VCAM1 vs PT", subtitle = "Adjusted p-value < 0.05") +
  theme_bw()


###################################################
file <- "G:/diabneph/analysis/dkd/markers/deg.celltype.markers.xlsx"
idents <- getSheetNames(file)
deg <- lapply(idents, function(ident){
  df <- read.xlsx(file, sheet = ident, rowNames = TRUE) %>%
    rownames_to_column(var = "gene")
  df$celltype <- ident
  return(df)
}) %>% bind_rows()
deg$celltype <- as.factor(deg$celltype)
levels(deg$celltype) <- unique(deg$celltype)


# intersect 
deg.genes <- deg[deg$gene %in% genes,]

deg.genes.filter <- deg.genes %>%
  dplyr::filter(p_val_adj < 0.05)

# visualize
deg.genes %>%
  dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0) %>%
  dplyr::mutate(fold_change = 2^avg_log2FC) %>%
  dplyr::mutate(label = paste0(gene,"_",celltype)) %>%
  ggplot(aes(avg_log2FC, -log10(p_val_adj), label=label, color=celltype)) +
  geom_point() +
  geom_text_repel() +
  xlim(c(0,2)) +
  xlab("Average log-fold change for DKD vs. Control") +
  ggtitle("Cell-specific genes by Cell Type", subtitle = "Adjusted p-value < 0.05") +
  theme_bw()

deg.genes %>%
  dplyr::filter(p_val < 0.05, avg_log2FC > 0) %>%
  dplyr::mutate(fold_change = 2^avg_log2FC) %>%
  dplyr::mutate(label = paste0(gene,"_",celltype)) %>%
  ggplot(aes(avg_log2FC, -log10(p_val_adj), label=label, color=celltype)) +
  geom_point() +
  geom_text_repel() +
  xlim(c(0,2)) +
  xlab("Average log-fold change for DKD vs. Control") +
  ggtitle("Differentially expressed genes in DKD by Cell Type", subtitle = "Unadjusted p-value < 0.05")

#################################################

file <- "G:/diabneph/analysis/dkd/markers/dar.macs2.PCT_vs_PTVCAM1.markers.xlsx"
dar <- read.xlsx(file, rowNames = TRUE)  

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



file <- "G:/diabneph/analysis/dkd/markers/dar.macs2.celltype.diab_vs_ctrl.xlsx"
idents <- getSheetNames(file)
dar <- lapply(idents, function(ident){
  df <- read.xlsx(file, sheet = ident, rowNames = TRUE)  
  df$celltype <- ident
  return(df)
}) %>% bind_rows()
dar$celltype <- as.factor(dar$celltype)
levels(dar$celltype) <- unique(dar$celltype)


# intersect 
dar.genes <- dar[dar$gene %in% genes,]

dar.genes.filter <- dar.genes %>%
  dplyr::filter(p_val_adj < 0.05)

# visualize
dar.genes %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::mutate(fold_change = 2^avg_log2FC) %>%
  dplyr::mutate(label = paste0(gene,"_",celltype)) %>%
  ggplot(aes(avg_log2FC, -log10(p_val_adj), label=label, color=celltype)) +
  geom_point() +
  geom_text_repel() +
  xlim(c(-0.25,0.25)) +
  xlab("Average log-fold change for DKD vs. Control") +
  ggtitle("Differentially accessible regions in DKD by Cell Type", subtitle = "Adjusted p-value < 0.05") +
  theme_bw()

deg.genes %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::mutate(fold_change = 2^avg_log2FC) %>%
  dplyr::mutate(label = paste0(gene,"_",celltype)) %>%
  ggplot(aes(avg_log2FC, -log10(p_val_adj), label=label, color=celltype)) +
  geom_point() +
  geom_text_repel() +
  xlim(c(-0.25,0.25)) +
  xlab("Average log-fold change for DKD vs. Control") +
  ggtitle("Differentially accessible regions in DKD by Cell Type", subtitle = "Unadjusted p-value < 0.05")


file <- "G:/diabneph/analysis/dkd/methylation/SD19_intersection_DMR_with_DAR_and_GR_cut_and_run.xlsx"
dmr <- read.xlsx(file, sheet = "ALL_DMR", rowNames = TRUE)  

# load databases
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(annotatr)
library(org.Hs.eg.db)
library(plyranges)
db <- c("hg38_genes_promoters")
anno <- build_annotations(genome = 'hg38', annotations=db)
dmr.gr <- as_granges(dmr)

# annotate peaks
dmr <- annotate_regions(dmr.gr, annotations=anno, ignore.strand=TRUE) %>%
  as.data.frame()

# intersect 
dmr.genes <- dmr[dmr$annot.gene_id %in% genes,]
