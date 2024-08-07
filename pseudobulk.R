library(openxlsx)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)
library(Hmisc)
library(reshape2)
library(Seurat)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(Seurat)
library(EnsDb.Hsapiens.v86)



######################
# subset for apoptosis genes
# correlation between hallmark apoptosis genes and biomarkers in PT_VCAM1
new_celltypes <- c("PT", "PT_VCAM1", "PEC", "iTAL", "TAL1", "TAL2", "DCT1", "DCT2", "PC", "ICA", "ICB", "PODO",
                   "ENDO","MES","FIB","LEUK")
new_levels <- c("PT", "PT_VCAM1", "PEC", "TAL1", "TAL2", "iTAL", "DCT1", "DCT2", "PC", "ICA", "ICB", "PODO",
                "ENDO","MES","FIB","LEUK")

xl <- read.xlsx("G:/krolewski/Joslin_New_46_prots.xlsx", sheet = "New_46_prots")
genes <- xl$Gene

# add grouping colors
xl2 <- read.xlsx("G:/krolewski/Proteins_list_46_for_Parker.xlsx") %>%
  dplyr::rename(Gene = "NAME_1") %>%
  dplyr::rename(color = Colors.in.Panel.A) %>%
  dplyr::mutate(color = ifelse(color == "Red", "tnf_and_apoptosis", "other")) %>%
  dplyr::select(Gene, Name_2, color)
xl2$color <- as.factor(xl2$color)
xl <- xl %>% left_join(xl2, by = "Gene")

idents <- getSheetNames("G:/diabneph/analysis/dkd/markers/deg.celltype.markers.xlsx")

df <- lapply(idents, function(ident) {
  res <- read.xlsx("G:/diabneph/analysis/dkd/markers/deg.celltype.markers.xlsx", sheet = ident)
  colnames(res)[1] <- "Gene"
  res$celltype <- ident
  return(res)
}) %>% bind_rows

celltype_update = data.frame(celltype = idents, new_celltype = new_celltypes)
df <- df %>% left_join(celltype_update, by = "celltype")

df$in_markers <- df$Gene %in% xl$Gene
df <- df %>%
  dplyr::mutate(gene_type = ifelse(Gene %in% xl$Gene[xl$color == "tnf_and_apoptosis"], "TNFR Signaling and \nApoptotic Process", "Other proteins"))

df <- df %>% 
  dplyr::filter(p_val_adj < 0.05)

toplot <- table(celltype=df$new_celltype, df$in_markers, gene_type = df$gene_type) %>%
  as.data.frame() %>%
  group_by(celltype) %>%
  mutate(prop = Freq / sum(Freq)) %>%
  dplyr::filter(Var2 == "TRUE") 
toplot[is.na(toplot)] <- 0

toplot$celltype <- factor(toplot$celltype, levels = new_levels)

panelA <- toplot %>%
  group_by(celltype) %>%
  dplyr::mutate(total_freq = sum(Freq)) %>%
  ggplot(aes(celltype, prop, fill=gene_type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c ("red", "blue")) + 
  theme_bw() +
  theme(text = element_text(size=8),
        plot.title = element_text(size=8, hjust = 0),
        axis.text = element_text(colour="black", size=8),
        axis.title=element_text(size=8),
        panel.border = element_rect(color="black",size=0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.justification="center",
        legend.text=element_text(size=8)) +
  guides(fill=guide_legend(nrow=1, ncol=2, reverse = TRUE)) +
  xlab("") +
  ylab("Proportion cell-specific genes") +
  geom_text(
    aes(label = total_freq, group = celltype), 
    stat = 'summary', fun = sum, vjust = -1
  ) +
  labs(fill = "Gene Type") +
  ylim(c(0,0.005)) +
  ggtitle("A)")
######################################################################################################
# deg plot PT vs PT_VCAM1
file <- "G:/diabneph/analysis/dkd/markers/deg.PT_vs_PTVCAM1.markers.xlsx"
deg <- read.xlsx(file, rowNames = TRUE) %>%
  rownames_to_column(var = "Gene") %>%
  left_join(xl, by = "Gene") %>%
  dplyr::filter(Gene %in% genes) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::mutate(fold_change = 2^avg_log2FC) %>%
  dplyr::mutate(label = paste0(Name_2)) %>%
  dplyr::mutate(color = ifelse(color %in% "tnf_and_apoptosis", "TNFR Signaling and \nApoptotic Process", "Other proteins"))

# visualize
panelB <- deg %>%
  ggplot(aes(avg_log2FC, -log10(p_val_adj), label=label)) +
  geom_point(aes(color=color), position="dodge", size=4, alpha=0.5) +
  scale_color_manual(values = c ("red", "blue")) + 
  geom_text_repel(show.legend = FALSE,
                  force=10,
                  max.overlaps = nrow(deg),
                  point.size = NA,
                  size=3) +
  xlab("Average log-fold change") +
  ggtitle("B)") +
  theme_bw() +
  ylim(c(0,150)) +
  xlim(c(-1,1)) +
  theme(text = element_text(size=8),
        plot.title = element_text(size=8, hjust = 0),
        axis.text = element_text(colour="black", size=8),
        axis.title=element_text(size=8),
        panel.border = element_rect(color="black",size=0.5),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.justification="center",
        legend.text=element_text(size=8)) +
  guides(color=guide_legend(nrow=1, ncol=2, reverse = TRUE))

######################################################################################################
# correlation between hallmark apoptosis genes and biomarkers in PT_VCAM1
xl <- read.xlsx("G:/krolewski/Joslin_New_46_prots.xlsx", sheet = "New_46_prots")
genes <- xl$Gene

# add grouping colors
xl2 <- read.xlsx("G:/krolewski/Proteins_list_46_for_Parker.xlsx") %>%
  dplyr::rename(Gene = "NAME_1") %>%
  dplyr::rename(color = Colors.in.Panel.A) %>%
  dplyr::mutate(color = ifelse(color == "Red", "TNFR Signaling and\nApoptotic Processes", "Other\nproteins")) %>%
  dplyr::select(Gene, Name_2, color)
xl2$color <- as.factor(xl2$color)
xl <- xl %>% left_join(xl2, by = "Gene")

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
apoptosis_counts <- counts[rownames(counts) %in% apoptosis_genes, rnaAggr@meta.data$celltype %in% c("PTVCAM1")] %>% 
  t() %>%
  as.data.frame()
apoptosis_counts$barcode <- rownames(apoptosis_counts)
anno <- data.frame(barcode = rownames(rnaAggr@meta.data), library_id = rnaAggr$orig.ident)
apoptosis_counts <- apoptosis_counts %>% left_join(anno, by = "barcode")
apoptosis_counts <- apoptosis_counts %>% group_by(library_id) %>% summarize(across(colnames(apoptosis_counts)[colnames(apoptosis_counts) %in% apoptosis_genes], sum))
apoptosis_counts <- apoptosis_counts[,colnames(apoptosis_counts)[colnames(apoptosis_counts) %in% apoptosis_genes]]
apoptosis_counts <- t(apoptosis_counts) 
colnames(apoptosis_counts) <- unique(rnaAggr$orig.ident)
apoptosis_totals <- colSums(apoptosis_counts)

counts <- counts[rownames(counts) %in% genes, rnaAggr@meta.data$celltype %in% c("PTVCAM1")]

# pseudobulk
counts.df <- as.data.frame(counts) %>%
  t() %>%
  as.data.frame()
counts.df$barcode <- rownames(counts.df)

anno <- data.frame(barcode = rownames(rnaAggr@meta.data), library_id = rnaAggr$orig.ident)
counts.df <- counts.df %>% left_join(anno, by = "barcode")
counts <- counts.df %>% group_by(library_id) %>% summarize(across(colnames(counts.df)[colnames(counts.df) %in% genes], sum))
counts <- counts %>% select(-library_id)

apoptosis_totals <- as.matrix(apoptosis_totals) %>% t()
rownames(apoptosis_totals) <- "apoptosis"

counts <- t(counts)

allcounts <- rbind(counts, apoptosis_totals)
cor <- rcorr(t(as.matrix(allcounts)))

# pearson correlation coefficients
cor.df <- cor$r %>% as.data.frame()
cor.df$gene <- rownames(cor.df) 

# retrieve pval for correlation between gene and apoptosis counts
cor.df$pval <- cor$P[colnames(cor$P) == "apoptosis"]
cor.df <- cor.df[,c("apoptosis","gene","pval")]

# join with group labels
cor.df <- cor.df %>% 
  dplyr::filter(gene %in% genes)  %>%
  dplyr::rename(Gene = gene) %>%
  left_join(xl, by = "Gene") 

# arrange by pval
cor.df <- cor.df %>%
  dplyr::arrange(desc(apoptosis)) %>%
  dplyr::mutate(star = ifelse(pval < 0.05, "*", "")) %>%
  dplyr::mutate(Gene = Name_2)

levels(cor.df$Gene) <- unique(cor.df$Gene)
cor.df$color <- as.factor(cor.df$color)

# Conditional statement to be used in plot
con <- ifelse(cor.df$color == "TNFR Signaling and\nApoptotic Processes", 'red', 'blue')

# prepare for plots
cor.df$variable <- ""
cor.df %>% 
  ggplot(aes(variable, Gene, fill = apoptosis, label = star)) +
  geom_tile() + 
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0,
                       guide = "colorbar",
                       limits = c(-0.1, 1),
                       breaks = c(-0.1, 1),
                       labels=c("-0.1","1")) +
  theme_bw() +
  geom_text(aes(label = star)) +
  ylim(rev(levels(cor.df$Gene))) +
  labs(fill = "", x = "", y = "") +
  ggtitle("C)") +
  xlab("Apoptosis\nIndex") +
  theme(text = element_text(size=8),
        plot.title = element_text(size=8, hjust = 0),
        axis.text.y = element_text(color=rev(con), size=8),
        legend.text=element_text(size=8),
        legend.title=element_text(size=8),
        axis.title=element_text(size=8),
        panel.border = element_rect(color="black",size=0.5),
        legend.pos = "bottom")
panel_levels = levels(cor.df$Gene)
