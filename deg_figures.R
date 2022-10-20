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
  ggtitle("Differentially expressed biomarker genes in DKD by cell type")

p1 <- top_layer + 
  geom_point(data = toplot2, aes(x=avg_log2FC, y=-log10(p_val_adj)), color = "gray") +
  theme(legend.title = element_blank())
  
p1$layers <- rev(p1$layers)
p1

##########################################################################
##########################################################################
##########################################################################

file <- "G:/diabneph/analysis/dkd/markers/deg.PT_vs_PTVCAM1.markers.xlsx"
deg <- read.xlsx(file, rowNames = TRUE) %>%
  rownames_to_column(var = "gene")

# intersect 
deg.genes <- deg[deg$gene %in% genes,]

# visualize
deg.genes %>%
  dplyr::mutate(color = factor(ifelse(p_val_adj < 0.05,"padj < 0.05","NS"), levels = c("padj < 0.05","NS"))) %>%
  dplyr::mutate(fold_change = 2^avg_log2FC) %>%
  dplyr::mutate(label = ifelse(p_val_adj < 0.05, gene, "")) %>%
  ggplot(aes(avg_log2FC, -log10(p_val_adj), label=label, color=color)) +
  geom_point() +
  geom_text_repel(show.legend = FALSE) +
  scale_color_manual(values = c("black","gray")) +
  xlim(c(-1,1)) +
  xlab("Average log-fold change for PT_VCAM1 vs PT") +
  ggtitle("Differentially expressed genes in PT_VCAM1 vs PT") +
  theme_bw() +
  theme(legend.title = element_blank())
  
  
  
  
