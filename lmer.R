#!/usr/bin/env Rscript
# to run locally:
SCRATCH1=/mnt/g/scratch
docker run -it --rm \
--workdir $HOME \
-v /mnt/g:$HOME/project \
-v $HOME:$HOME \
-v $SCRATCH1:$SCRATCH1 \
-e SCRATCH1="/mnt/g/scratch" \
p4rkerw/allele_mod:latest R

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(Matrix)
library(future.apply)
library(tibble)
library(here)
library(reshape2)
library(broom.mixed)
library(openxlsx)
library(EnsDb.Hsapiens.v86)
library(Seurat)

library(lmerTest)
library(msigdbr)


# correlation between hallmark apoptosis genes and biomarkers in PT_VCAM1
xl <- read.xlsx("project/krolewski/Joslin_New_46_prots.xlsx", sheet = "New_46_prots")
genes <- xl$Gene

# add grouping colors
xl2 <- read.xlsx("project/krolewski/Proteins_list_46_for_Parker.xlsx") %>%
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

rnaAggr <- readRDS("project/diabneph/analysis/dkd/rna_aggr_prep/step2_magic.rds")
DefaultAssay(rnaAggr) <- "MAGIC_RNA"
counts <- GetAssayData(rnaAggr, slot = "data")
apoptosis_counts <- counts[rownames(counts) %in% apoptosis_genes, rnaAggr@meta.data$celltype %in% c("PTVCAM1")]
apoptosis_totals <- colSums(apoptosis_counts) %>% as.data.frame %>% t()
rownames(apoptosis_totals) <- "apoptosis"

counts <- counts[, rnaAggr@meta.data$celltype %in% c("PTVCAM1")]

meta <- data.frame(barcode = rownames(rnaAggr@meta.data), nCount_RNA = rnaAggr$nCount_RNA, sample=rnaAggr$orig.ident)

apoptosis.df <- t(apoptosis_totals) %>%
  as.data.frame() %>%
  rownames_to_column(var = "barcode")

results.ls <- mclapply(rownames(counts), function(gene) {
	  print(gene)
	  exp_mat <- counts[gene,] %>% as.data.frame()
	  colnames(exp_mat) <- "imprna"
	  exp_mat <- rownames_to_column(exp_mat, var="barcode") 

	  mat <- left_join(exp_mat, meta, by="barcode")
	  
	  # scale the values
	  mat <- left_join(mat, apoptosis.df, by = "barcode")
	  mat <- mat %>% dplyr::mutate(scaled_exp = scale(imprna),
	  							   scaled_nCount_RNA = scale(nCount_RNA),
	  							   scaled_apoptosis = scale(apoptosis))
	  
	  # run model
	  fit <- tryCatch({lmerTest::lmer(y ~ exp + nCount_RNA + (1|ident), data = data.frame(exp=mat$scaled_exp, y=mat$scaled_apoptosis, nCount_RNA=mat$scaled_nCount_RNA, ident=mat$sample))
	  }, error=function(e) NULL)

	  if(is.null(fit)) {
	  	# if model returns ERROR return null
		return(NULL)	
	  } else {
	    res <- tidy(fit, conf.int=TRUE) %>%
		  dplyr::select(term, estimate, std.error, p.value, conf.low, conf.high) %>%
		  as.data.frame() %>%
		  melt(id.var = "term") %>%
		  dcast(1~variable+term)
	    
            # reorder fixed and random effect columns and format results
	    exp_cols <- colnames(res)[str_ends(colnames(res), "_exp")]
	    ranef_cols <- colnames(res)[str_detect(colnames(res), "sd__")][1] # only take sd of ranef

	    # designate column output order
        output_cols <- unlist(list(exp_cols, ranef_cols))
	    res <- res[,output_cols]
		  
	    res <- cbind(gene=gene, res)
	    }
   }, mc.cores=16)
results.df <- bind_rows(results.ls) %>% as.data.frame()

write.csv(results.df, file="project/diabneph/analysis/dkd/rna_aggr_prep/lmer_results.csv")

results.df <- results.df %>%
  dplyr::mutate(label = ifelse(gene %in% xl$Gene, gene, ""))

p1 <- results.df %>%
  ggplot(aes(x=estimate_exp, y=-log10(p.value_exp))) + 
  geom_point(size=0.1, alpha=0.5) +
  stat_density_2d(linewidth=1, color="blue", alpha=0.5) +
  theme_bw() +
  coord_cartesian(xlim = c(-0.25,0.5), ylim = c(0,200)) + 
  theme(legend.pos = "none")

p2 <- results.df %>%
  dplyr::filter(gene %in% xl$Gene) %>%
  ggplot(aes(x=estimate_exp, y=-log10(p.value_exp), label=gene)) + 
  geom_point() +
  stat_density_2d(linewidth=1, color="blue", alpha=0.5) +
  theme_bw() +
  geom_text_repel() +
  coord_cartesian(xlim = c(-0.25,0.5), ylim = c(0,200)) + 
  theme(legend.pos = "none")

grid.arrange(p1,p2)

#################

counts <- GetAssayData(rnaAggr, slot = "data")
apoptosis_counts <- counts[rownames(counts) %in% apoptosis_genes, rnaAggr@meta.data$celltype %in% c("PTVCAM1")]
apoptosis_totals <- colSums(apoptosis_counts) %>% as.data.frame %>% t()
rownames(apoptosis_totals) <- "apoptosis"

# random 46 genes
genes <- genes(EnsDb.Hsapiens.v86)
genes <- genes[genes$gene_biotype == "protein_coding"]
genes <- sample(genes$symbol, 46)

counts <- counts[rownames(counts) %in% genes, rnaAggr@meta.data$celltype %in% c("PTVCAM1")]

meta <- data.frame(barcode = rownames(rnaAggr@meta.data), nCount_RNA = rnaAggr$nCount_RNA, sample=rnaAggr$orig.ident)

apoptosis.df <- t(apoptosis_totals) %>%
  as.data.frame() %>%
  rownames_to_column(var = "barcode")

results.ls <- lapply(rownames(counts), function(gene) {
	  print(gene)
	  exp_mat <- counts[gene,] %>% as.data.frame()
	  colnames(exp_mat) <- "imprna"
	  exp_mat <- rownames_to_column(exp_mat, var="barcode") 

	  mat <- left_join(exp_mat, meta, by="barcode")
	  
	  # scale the values
	  mat <- left_join(mat, apoptosis.df, by = "barcode")
	  mat <- mat %>% dplyr::mutate(scaled_exp = scale(imprna),
	  							   scaled_nCount_RNA = scale(nCount_RNA),
	  							   scaled_apoptosis = scale(apoptosis))
	  
	  # run model
	  fit <- tryCatch({lmerTest::lmer(y ~ exp + nCount_RNA + (1|ident), data = data.frame(exp=mat$scaled_exp, y=mat$scaled_apoptosis, nCount_RNA=mat$scaled_nCount_RNA, ident=mat$sample))
	  }, error=function(e) NULL)

	  if(is.null(fit)) {
	  	# if model returns ERROR return null
		return(NULL)	
	  } else {
	    res <- tidy(fit, conf.int=TRUE) %>%
		  dplyr::select(term, estimate, std.error, p.value, conf.low, conf.high) %>%
		  as.data.frame() %>%
		  melt(id.var = "term") %>%
		  dcast(1~variable+term)
	    
            # reorder fixed and random effect columns and format results
	    exp_cols <- colnames(res)[str_ends(colnames(res), "_exp")]
	    ranef_cols <- colnames(res)[str_detect(colnames(res), "sd__")][1] # only take sd of ranef

	    # designate column output order
        output_cols <- unlist(list(exp_cols, ranef_cols))
	    res <- res[,output_cols]
		  
	    res <- cbind(gene=gene, res)
	    }
   })
results.df <- bind_rows(results.ls) %>% as.data.frame()


p2 <- results.df %>%
  ggplot(aes(x=estimate_exp, y=-log10(p.value_exp), label=gene)) + 
  geom_point() +
  theme_bw() +
  geom_text_repel()   
