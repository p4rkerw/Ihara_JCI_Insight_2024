library(openxlsx)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)
library(Hmisc)
library(reshape2)

######################################################################################################
# correlation between hallmark apoptosis genes and biomarkers in PT_VCAM1
xl <- read.xlsx("G:/krolewski/Joslin_New_46_prots.xlsx", sheet = "New_46_prots")
genes <- xl$Gene

# add grouping colors
xl2 <- read.xlsx("G:/krolewski/Proteins_list_46_for_Parker.xlsx") %>%
  dplyr::rename(Gene = "NAME_1") %>%
  dplyr::rename(color = Colors.in.Panel.A) %>%
  dplyr::mutate(color = ifelse(color == "Red", "TNF Signaling", "Apoptotic Processes"))
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
apoptosis_counts <- counts[rownames(counts) %in% apoptosis_genes, rnaAggr@meta.data$celltype %in% c("PTVCAM1")]
apoptosis_totals <- colSums(apoptosis_counts) %>% as.data.frame %>% t()
rownames(apoptosis_totals) <- "apoptosis"

counts <- counts[rownames(counts) %in% genes, rnaAggr@meta.data$celltype %in% c("PTVCAM1")]
allcounts <- rbind(counts, apoptosis_totals)
cor <- rcorr(t(as.matrix(allcounts)))

# pearson correlation coefficients
cor.df <- cor$r %>% as.data.frame()
cor.df$gene <- rownames(cor.df) 

# retrieve pval for correlation between gene and apoptosis counts
cor.df$pval <- cor$P[colnames(cor$P) == "apoptosis"]
cor.df <- cor.df[,c("apoptosis","gene","pval")]
cor.df <- cor.df %>% dplyr::filter(gene %in% genes)
cor.df$label <- ifelse(cor.df$pval < 0.05, "*", "")
cor.df <- dplyr::arrange(cor.df, gene)
cor.df$gene <- as.factor(cor.df$gene)
levels(cor.df$gene) <- unique(cor.df$gene)

# prepare for plots
cor.df$variable <- ""
p3 <- cor.df %>% 
  ggplot(aes(variable, gene, fill = apoptosis, label = label)) +
  geom_tile() + 
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0,
                       guide = "colorbar") +
  theme_bw() +
  geom_text(aes(label = label)) +
  ylim(rev(levels(cor.df$gene))) +
  labs(fill = "Pearson r", x = "", y = "") +
  ggtitle("B) Correlation") +
  xlab("Apoptosis Index") +
  theme(plot.title = element_text(size=20, hjust = 0),
        axis.text = element_text(colour="black", size=8),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        axis.title=element_text(size=14))

#####################################################################################################
# deg plot PT vs PT_VCAM1
file <- "G:/diabneph/analysis/dkd/markers/deg.PT_vs_PTVCAM1.markers.xlsx"
deg <- read.xlsx(file, rowNames = TRUE) %>%
    rownames_to_column(var = "Gene") %>%
    left_join(xl, by = "Gene") %>%
    dplyr::filter(Gene %in% genes) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::mutate(fold_change = 2^avg_log2FC) %>%
    dplyr::mutate(label = paste0(Name_2))

# visualize
p1 <- deg %>%
  ggplot(aes(avg_log2FC, -log10(p_val_adj), label=label)) +
  geom_point(aes(color=color)) +
  scale_color_manual(values = c (`TNF Signaling` = "red", `Apoptotic Processes` = "blue")) + 
  geom_text_repel(size=4) +
  xlab("Average log-fold change") +
  ggtitle("A) DE Biomarker Genes") +
  theme_bw() +
  ylim(c(0,150)) +
  xlim(c(-1,1)) +
  theme(plot.title = element_text(size=20, hjust = 0),
        axis.text = element_text(colour="black", size=12),
        axis.title=element_text(size=14),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.justification="center",
        legend.text=element_text(size=12)) +
  guides(color=guide_legend(nrow=1, byrow=TRUE))

# # intersect with hallmark apoptosis genes
# deg.genes <- deg[deg$gene %in% apoptosis_genes,]

# # visualize
# pX <- deg.genes %>%
#   dplyr::filter(p_val_adj < 0.05) %>%
#   dplyr::mutate(fold_change = 2^avg_log2FC) %>%
#   dplyr::mutate(label = gene) %>%
#   ggplot(aes(avg_log2FC, -log10(p_val_adj), label=label)) +
#   geom_point() +
#   geom_text_repel(size=4) +
#   xlab("Average log-fold change") +
#   ggtitle("B) DE Apoptosis Genes") +
#   theme_bw() +
#   ylim(c(0,150)) +
#   xlim(c(-1,1)) +
#   theme(plot.title = element_text(size=20, hjust = 0),
#         axis.text = element_text(colour="black", size=12),
#         axis.title=element_text(size=14))

#################################################
# dar PT vs PT_VCAM1 volcano
file <- "G:/diabneph/analysis/dkd/markers/dar.macs2.PCT_vs_PTVCAM1.markers.xlsx"
dar <- read.xlsx(file, rowNames = TRUE)  

# intersect 
dar.genes <- dar[dar$gene %in% genes,]

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

db <- c("hg38_genes_promoters","hg38_genes_introns")
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

# recode the annot.type to shorten descriptions
dar <- dar %>%
  dplyr::mutate(annot.type = recode(annot.type,
                hg38_enhancers_fantom = "Enhancer",
                hg38_genes_introns = "Intron",
                hg38_genes_promoters = "Promoter",
                hg38_genes_intergenic = "Intergenic",
                hg38_genes_exons = "Exon"))

# intersect 
dar <- dar %>% 
  dplyr::filter(gene %in% genes) %>%
  dplyr::rename(Gene = gene) %>%
  left_join(xl, by = "Gene") %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::mutate(fold_change = 2^avg_log2FC) %>%
  dplyr::mutate(label = paste0(Name_2))

# visualize
p4 <- dar %>%
  na.omit() %>%
  dplyr::mutate(logpval = -log10(p_val_adj)) %>%
  dplyr::mutate(logpval = ifelse(logpval > 100, 100, logpval)) %>%
  ggplot(aes(avg_log2FC, logpval, label=label, shape=annot.type)) +
  geom_point(aes(color=color)) +
  scale_color_manual(values = c (`TNF Signaling` = "red", `Apoptotic Processes` = "blue")) + 
  geom_text_repel(show.legend = FALSE, max.overlaps=10) +
  xlim(c(-0.15,0.25)) +
  xlab("Average log-fold change") +
  ylab("-log10(p_val_adj)") + 
  ggtitle("C) DA Biomarker Genes") +
  theme_bw() +
  theme(legend.title = element_blank(),
        plot.title = element_text(size=20, hjust = 0),
        axis.text = element_text(colour="black", size=12),
        axis.title=element_text(size=14),
        legend.position="bottom",
        legend.justification="center",
        legend.text=element_text(size=12)) +
  guides(color=guide_legend(nrow=1, ncol=2))

# # intersect 
# dar.genes <- dar[dar$gene %in% apoptosis_genes,]
# pX <- dar.genes %>%
#   na.omit() %>%
#   dplyr::filter(p_val_adj < 0.05) %>%
#   dplyr::mutate(logpval = -log10(p_val_adj)) %>%
#   dplyr::mutate(logpval = ifelse(logpval > 100, 100, logpval)) %>%
#   dplyr::mutate(fold_change = 2^avg_log2FC) %>%
#   dplyr::mutate(label = paste0(gene)) %>%
#   ggplot(aes(avg_log2FC, logpval, label=label, color=annot.type)) +
#   geom_point() +
#   geom_text_repel(show.legend = FALSE, max.overlaps=10) +
#   xlim(c(-0.15,0.25)) +
#   xlab("Average log-fold change") +
#   ylab("-log10(p_val_adj)") + 
#   ggtitle("E) DA Apoptosis Genes") +
#   theme_bw() +
#   theme(legend.title = element_blank(),
#         plot.title = element_text(size=20, hjust = 0),
#         axis.text = element_text(colour="black", size=12),
#         axis.title=element_text(size=14),
#         legend.position="bottom",
#         legend.justification="center",
#         legend.text=element_text(size=12)) +
#   guides(color=guide_legend(nrow=1, byrow=TRUE))

# arrange
library(gridExtra)
pdf("G:/krolewski/figure.pdf",width=15, height=12)
margin = theme(plot.margin = unit(c(0.5,0.5,0.5,0.5,0.5), "cm"))
pl <- list(p1,p3,p4)
grid.arrange(grobs = lapply(pl, "+", margin), ncol=2)
dev.off()

# ################################################################################
# xl <- read.xlsx("G:/downloads/Joslin_New_46_prots.xlsx", sheet = "New_46_prots")
# genes <- xl$Gene

# file <- "G:/diabneph/analysis/dkd/markers/deg.celltype.diab_vs_ctrl.xlsx"
# idents <- getSheetNames(file)
# deg <- lapply(idents, function(ident){
#   df <- read.xlsx(file, sheet = ident, rowNames = TRUE) %>%
#     rownames_to_column(var = "gene")
#   df$celltype <- ident
#   return(df)
# }) %>% bind_rows()
# deg$celltype <- as.factor(deg$celltype)
# levels(deg$celltype) <- unique(deg$celltype)


# # intersect 
# deg.genes <- deg[deg$gene %in% genes,]

# # create a color variable and maintain cell type levels
# deg.genes <- deg.genes %>%
#   dplyr::mutate(color = factor(ifelse(p_val_adj < 0.05, as.character(celltype), "NS"))) %>%
#   dplyr::mutate(fold_change = 2^avg_log2FC) %>%
#   dplyr::mutate(label = ifelse(p_val_adj < 0.05, paste0(gene,"_",celltype), ""))
# deg.genes$color <- factor(deg.genes$color, levels = c(levels(deg.genes$celltype),"NS"))

# # plot the background NS degs
# toplot1 <- deg.genes %>%
#   dplyr::filter(color != "NS")
# toplot2 <- deg.genes %>%
#   dplyr::filter(color == "NS")

# top_layer <- ggplot(data = toplot1, aes(x=avg_log2FC, y=-log10(p_val_adj), color=color, label=label)) + 
#   geom_point() +
#   geom_text_repel(show.legend = FALSE) +
#   theme_bw() +
#   xlab("Average log-fold change for DKD vs. Control") +
#   ggtitle("Differentially expressed biomarker genes in DKD by cell type", subtitle = "Adjusted p-value < 0.05")

# p1 <- top_layer + 
#   geom_point(data = toplot2, aes(x=avg_log2FC, y=-log10(p_val_adj)), color = "gray") +
#   theme(legend.title = element_blank())
  
# p1$layers <- rev(p1$layers)
# p1

# ###################################################
# file <- "G:/diabneph/analysis/dkd/markers/deg.celltype.markers.xlsx"
# idents <- getSheetNames(file)
# deg <- lapply(idents, function(ident){
#   df <- read.xlsx(file, sheet = ident, rowNames = TRUE) %>%
#     rownames_to_column(var = "gene")
#   df$celltype <- ident
#   return(df)
# }) %>% bind_rows()
# deg$celltype <- as.factor(deg$celltype)
# levels(deg$celltype) <- unique(deg$celltype)


# # intersect 
# deg.genes <- deg[deg$gene %in% genes,]

# deg.genes.filter <- deg.genes %>%
#   dplyr::filter(p_val_adj < 0.05)

# # visualize
# deg.genes %>%
#   dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0) %>%
#   dplyr::mutate(fold_change = 2^avg_log2FC) %>%
#   dplyr::mutate(label = paste0(gene,"_",celltype)) %>%
#   ggplot(aes(avg_log2FC, -log10(p_val_adj), label=label, color=celltype)) +
#   geom_point() +
#   geom_text_repel() +
#   xlim(c(0,2)) +
#   xlab("Average log-fold change for DKD vs. Control") +
#   ggtitle("Cell-specific genes by Cell Type", subtitle = "Adjusted p-value < 0.05") +
#   theme_bw()

# deg.genes %>%
#   dplyr::filter(p_val < 0.05, avg_log2FC > 0) %>%
#   dplyr::mutate(fold_change = 2^avg_log2FC) %>%
#   dplyr::mutate(label = paste0(gene,"_",celltype)) %>%
#   ggplot(aes(avg_log2FC, -log10(p_val_adj), label=label, color=celltype)) +
#   geom_point() +
#   geom_text_repel() +
#   xlim(c(0,2)) +
#   xlab("Average log-fold change for DKD vs. Control") +
#   ggtitle("Differentially expressed genes in DKD by Cell Type", subtitle = "Unadjusted p-value < 0.05")
