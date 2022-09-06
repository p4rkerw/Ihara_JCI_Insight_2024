library(openxlsx)

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
