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
  theme(plot.title = element_text(size=20, hjust = 0),
        axis.text = element_text(colour="black", size=12),
        axis.title=element_text(size=14),
        panel.border = element_rect(color="black",size=1),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.justification="center",
        legend.text=element_text(size=12)) +
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
                  point.size = NA) +
  xlab("Average log-fold change") +
  ggtitle("B)") +
  theme_bw() +
  ylim(c(0,150)) +
  xlim(c(-1,1)) +
  theme(plot.title = element_text(size=20, hjust = 0),
        axis.text = element_text(colour="black", size=12, face="bold"),
        axis.title=element_text(size=14),
        panel.border = element_rect(color="black",size=1),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.justification="center",
        legend.text=element_text(size=12)) +
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
panelC <- cor.df %>% 
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
  theme(plot.title = element_text(size=20, hjust = 0),
        axis.text.y = element_text(color=rev(con), size=10),
        legend.text=element_text(size=12, face="bold"),
        legend.title=element_text(size=12),
        axis.title=element_text(size=14),
        panel.border = element_rect(color="black",size=1),
        legend.pos = "bottom") 
panel_levels = levels(cor.df$Gene)

#################################################################################################
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
rnaAggr <- readRDS("G:/diabneph/analysis/dkd/rna_aggr_prep/step2_magic.rds")
DefaultAssay(rnaAggr) <- "MAGIC_RNA"
counts <- GetAssayData(rnaAggr, slot = "data", assay = "MAGIC_RNA")
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

# prepare for plots
cor.df$variable <- ""
panelD <- cor.df %>% 
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
  ylim(rev(panel_levels)) +
  labs(fill = "", x = "", y = "") +
  ggtitle("D)") +
  xlab("Imputed\nIndex") +
  theme(plot.title = element_text(size=20, hjust = 0),
        axis.text.y = element_text(color=rev(con), size=10),
        legend.text=element_text(size=12, face="bold"),
        legend.title=element_text(size=12),
        axis.title=element_text(size=14),
        panel.border = element_rect(color="black",size=1),
        legend.pos = "bottom") 

####################################################

library(gridExtra)
pdf("G:/krolewski/revised_figure.pdf",width=10, height=8.5)
margin = theme(plot.margin = unit(c(0.25,0.25,0.25,0.25,0.25), "cm"))
pl <- list(panelA, panelB, panelC, panelD)
grid.arrange(grobs = lapply(pl, "+", margin),
             ncol=3,
             layout_matrix = rbind(c(1,3,4),
                                   c(2,3,4)),
             widths = c(2,1,1))
dev.off()

#####################################################################################################
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
  dplyr::mutate(label = paste0(gene))

# visualize
supplemental_dar <- dar %>%
  na.omit() %>%
  dplyr::mutate(logpval = -log10(p_val_adj)) %>%
  dplyr::mutate(logpval = ifelse(logpval > 150, 150, logpval)) %>%
  ggplot(aes(avg_log2FC, logpval, label=label, shape=annot.type)) +
  geom_point(aes(color=color), position="dodge", size=4, alpha=0.5) +
  scale_color_manual(values = c (`TNFR Signaling and\nApoptotic Processes` = "red", `Other\nproteins` = "blue")) + 
  geom_text_repel(show.legend = FALSE,
                  force=10,
                  max.overlaps = nrow(dar)) +
  xlim(c(-0.15,0.25)) +
  xlab("Average log-fold change") +
  ylab("-log10(p_val_adj)") +
  ggtitle("C)") +
  theme_bw() +
  theme(legend.title = element_blank(),
        plot.title = element_text(size=20, hjust = 0),
        axis.text = element_text(colour="black", size=12, face="bold"),
        axis.title=element_text(size=14),
        panel.border = element_rect(color="black",size=1),
        legend.position="bottom",
        legend.justification="center",
        legend.text=element_text(size=12)) +
  guides(color=guide_legend(nrow=1, reverse = TRUE))
##############################################################
# lmer supplemental figures
# correlation between hallmark apoptosis genes and biomarkers in PT_VCAM1
xl <- read.xlsx("G:/krolewski/Joslin_New_46_prots.xlsx", sheet = "New_46_prots")
genes <- xl$Gene	
results.df <- read.csv("G:/diabneph/analysis/dkd/rna_aggr_prep/lmer_results.csv")

results.df <- results.df %>%
  dplyr::mutate(label = ifelse(gene %in% xl$Gene, gene, ""))

density_sel <- results.df %>%
  dplyr::filter(gene %in% xl$Gene) %>%
  ggplot(aes(x=estimate_exp, y=-log10(p.value_exp))) + 
  stat_density_2d(n=1000)

p_sel <- ggplot_build(density_sel)
p_sel <- p_sel$data[[1]]

density <- results.df %>%
  ggplot(aes(x=estimate_exp, y=-log10(p.value_exp))) + 
  stat_density_2d(n=1000)

p <- ggplot_build(density)
p <- p$data[[1]]
p$gene <- ""


p1 <- results.df %>%
  ggplot(aes(x=estimate_exp, y=-log10(p.value_exp))) + 
  geom_point(size=0.1, alpha=0.5) +
  stat_density_2d(linewidth=1, color="red", alpha=0.5) +
  geom_point(data = p_sel, aes(x, y), color = "blue", size = 0.75) +
  theme_bw() +
  coord_cartesian(xlim = c(-0.25,0.5), ylim = c(0,200)) + 
  theme(legend.pos = "none") +
  xlab("Beta coefficient apoptosis genes") +
  ylab("-log10(pval)") +
  ggtitle("A)")

p2 <- results.df %>%
  dplyr::filter(gene %in% xl$Gene) %>%
  ggplot(aes(x=estimate_exp, y=-log10(p.value_exp), label=gene)) + 
  geom_point() +
  stat_density_2d(linewidth=1, color="blue", alpha=0.5) +
  geom_point(data = p, aes(x, y), color = "red", size = 0.75) +
  theme_bw() +
  geom_text_repel() +
  coord_cartesian(xlim = c(-0.25,0.5), ylim = c(0,200)) + 
  theme(legend.pos = "none") +
  xlab("Beta coefficient apoptosis genes")+
  ylab("-log10(pval)")+
  ggtitle("B)")

grid.arrange(p1,p2)

