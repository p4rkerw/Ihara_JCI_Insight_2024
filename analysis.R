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
  dplyr::arrange(pval) %>%
  dplyr::mutate(star = ifelse(pval < 0.05, "*", "")) %>%
  dplyr::mutate(Gene = Name_2)

levels(cor.df$Gene) <- unique(cor.df$Gene)
cor.df$color <- as.factor(cor.df$color)

# Conditional statement to be used in plot
con <- ifelse(cor.df$color == "TNFR Signaling and\nApoptotic Processes", 'red', 'blue')

# prepare for plots
cor.df$variable <- ""
p3 <- cor.df %>% 
  ggplot(aes(variable, Gene, fill = apoptosis, label = star)) +
  geom_tile() + 
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0,
                       guide = "colorbar") +
  theme_bw() +
  geom_text(aes(label = star)) +
  ylim(rev(levels(cor.df$Gene))) +
  labs(fill = "Pearson r", x = "", y = "") +
  ggtitle("B)") +
  xlab("Apoptosis Index") +
  theme(plot.title = element_text(size=20, hjust = 0),
        axis.text.y = element_text(color=rev(con), size=10),
        legend.text=element_text(size=12, face="bold"),
        legend.title=element_text(size=12),
        axis.title=element_text(size=14),
        panel.border = element_rect(color="black",size=1))

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
  geom_point(aes(color=color), size=4) +
  scale_color_manual(values = c (`TNFR Signaling and\nApoptotic Processes` = "red", `Other\nproteins` = "blue")) + 
  geom_text_repel(show.legend = FALSE, max.overlaps=8) +
  xlab("Average log-fold change") +
  ggtitle("A)") +
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
  dplyr::mutate(logpval = ifelse(logpval > 150, 150, logpval)) %>%
  ggplot(aes(avg_log2FC, logpval, label=label, shape=annot.type)) +
  geom_point(aes(color=color), size=4) +
  scale_color_manual(values = c (`TNFR Signaling and\nApoptotic Processes` = "red", `Other\nproteins` = "blue")) + 
  geom_text_repel(show.legend = FALSE, max.overlaps=6, point.padding=0.5) +
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

# # arrange
# library(gridExtra)
# pdf("G:/krolewski/figure.pdf",width=8.5, height=11)
# margin = theme(plot.margin = unit(c(0.5,0.5,0.5,0.5,0.5), "cm"))
# pl <- list(p1,p3,p4)
# grid.arrange(grobs = lapply(pl, "+", margin),
#              ncol=2,
#              widths = c(0.9,0.7))
# dev.off()

library(gridExtra)
pdf("G:/krolewski/figure.pdf",width=8.5, height=11)
margin = theme(plot.margin = unit(c(0.5,0.5,0.5,0.5,0.5), "cm"))
pl <- list(p1,p3,p4)
grid.arrange(grobs = lapply(pl, "+", margin),
             ncol=2,
             layout_matrix = cbind(c(1,3), c(2,2)))
dev.off()


