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

deg.genes.filter <- deg.genes %>%
  dplyr::filter(p_val_adj < 0.05)

# visualize
deg.genes %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::mutate(fold_change = 2^avg_log2FC) %>%
  dplyr::mutate(label = paste0(gene,"_",celltype)) %>%
  ggplot(aes(avg_log2FC, -log10(p_val_adj), label=label, color=celltype)) +
  geom_point() +
  geom_text_repel() +
  xlim(c(-1,1)) +
  xlab("Average log-fold change for DKD vs. Control") +
  ggtitle("Differentially expressed genes in DKD by Cell Type", subtitle = "Adjusted p-value < 0.05") +
  theme_bw()

deg.genes %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::mutate(label = paste0(gene,"_",celltype)) %>%
  dplyr::mutate(fold_change = 2^avg_log2FC) %>%
  ggplot(aes(avg_log2FC, -log10(p_val_adj), label=label, color=celltype)) +
  geom_point() +
  geom_text_repel() +
  xlim(c(-2,2)) +
  xlab("Average log-fold change for DKD vs. Control") +
  ggtitle("Differentially expressed genes in DKD by Cell Type", subtitle = "Unadjusted p-value < 0.05")

######################################################################################################
# retrieve hallmark apoptosis genes
hallmark <- msigdbr(species = "Homo sapiens", category = "H") 
apoptosis <- hallmark[grepl("APOPTOSIS",hallmark$gs_name),]
apoptosis_genes <- apoptosis$gene_symbol

# correlate expression of biomarkers with hallmark apoptosis genes
library(Seurat)
rnaAggr <- readRDS("G:/diabneph/analysis/dkd/rna_aggr_prep/step2_anno.rds")
DefaultAssay(rnaAggr) <- "RNA"
rnaAggr <- NormalizeData(rnaAggr)
counts <- GetAssayData(rnaAggr, slot = "data")
counts <- counts[rownames(counts) %in% genelist, rnaAggr@meta.data$celltype %in% c("PTVCAM1")]
cor.df <- cor(t(as.matrix(counts))) %>% as.data.frame()
cor.df$gene <- rownames(cor.df) 
cor.df <- melt(cor.df)
cor.df <- cor.df %>% dplyr::arrange(gene)
cor.df$gene <- as.factor(cor.df$gene)
cor.df <- dplyr::filter(cor.df, gene %in% apoptosis_genes, variable %in% genes)

ggplot(cor.df, aes(gene, variable, fill = value)) +
  geom_tile() + 
  scale_fill_gradient(low = "black",
                    high = "red",
                    breaks = c(0,0.25,0.5,1),
                    guide = "colorbar")


avexp <- AverageExpression(rnaAggr)$RNA

# subset for biomarkers and hallmark apoptosis genes
library(reshape2)
genelist <- c(apoptosis_genes, genes) %>% unique()
mat <- avexp[rownames(avexp) %in% genelist,] %>% as.data.frame()
# mat$gene <- rownames(mat)
# mat <- melt(mat)
# ggplot(mat, aes(gene, variable, fill = value)) + geom_tile()

cor.exp <- cor(as.data.frame(t(mat))) %>% as.data.frame()

# put the biomarkers on x-axis and the apoptosis genes on y-axis
cor.exp <- cor.exp[rownames(cor.exp) %in% genes, colnames(cor.exp) %in% apoptosis_genes]

cor.exp$gene <- rownames(cor.exp)
cor.exp <- melt(cor.exp)
cor.exp <- cor.exp %>% dplyr::arrange(desc(value))
cor.exp$gene <- as.factor(cor.exp$gene)

ggplot(cor.exp, aes(gene, variable, fill = value)) + geom_tile()

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
