#!/usr/bin/env Rscript
# to run locally:
# SCRATCH1=/mnt/g/scratch
# docker run -it \
# --workdir $HOME \
# -v /mnt/g/diabneph:$HOME/project \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# p4rkerw/sctools:R4.1.0

library(Rmagic)
library(Seurat)
rna <- readRDS("project/analysis/dkd/rna_aggr_prep/step2_anno.rds")
DefaultAssay(rna) <- "RNA"
magic <- magic(rna, solver="approximate")

saveRDS(magic, "project/analysis/dkd/rna_aggr_prep/step2_magic.rds", compress=FALSE)
