
########################################
## Script to do the mapping using MNN ##
########################################

# Note: the atlas data set is not provided alongside this data set. 
# You can find the output of the mapping in the file "rna/mapping.rds".

library(SingleCellExperiment)
library(data.table)
library(purrr)
library(ggplot2)

source("/Users/C02RF23NFVH8/scnmt_eb/rna/mapping/10x/Mapping2gastrulationAtlas.R")

path2atlas <- "/Users/ricard/data/gastrulation10x"
path2scNMT <- "/Users/ricard/data/scnmt_eb"

####################
## Load 10x atlas ##
####################

sce_atlas  <- readRDS(paste0(path2atlas, "/processed/SingleCellExperiment.rds"))
meta_atlas <- readRDS(paste0(path2atlas, "/processed/sample_metadata.rds"))

##########################
## Load scNMT-seq query ##
##########################

sce_nmt  <- readRDS(paste0(path2scNMT, "/rna/SingleCellExperiment.rds"))
meta_nmt <- read.table(paste0(path2scNMT, "/sample_metadata.txt"), header = TRUE, sep = "\t", stringsAsFactors = F)
meta_nmt$stage <- meta_nmt$day

# Filter
meta_nmt <- meta_nmt[meta_nmt$pass_rnaQC==T,]
sce_nmt <- sce_nmt[,meta_nmt$id_rna] 

#############
## Prepare ## 
#############

# Data structure required for the mapping...
meta_scnmt <- list()
meta_scnmt$cell <- meta_nmt$id_rna[match(colnames(sce_nmt), meta_nmt$id_rna)]
meta_scnmt$cells <- meta_nmt$id_rna[match(colnames(sce_nmt), meta_nmt$id_rna)]
meta_scnmt$stage <- meta_nmt$day[match(colnames(sce_nmt), meta_nmt$id_rna)]

# Subset genes that are present in both data sets
genes <- intersect(rownames(sce_nmt), rownames(sce_atlas))
sce_nmt  <- sce_nmt[genes,]
sce_atlas <- sce_atlas[genes,]

#########
## Map ##
#########

mapping  <- mapWrap(
  atlas_sce = sce_atlas, atlas_meta = meta_atlas,
  map_sce = sce_nmt, map_meta = meta_scnmt, 
  k = 25
)

##########
## Save ##
##########

saveRDS(mapping, "/Users/ricard/data/scnmt_eb/rna/mapping/mapping.rds")

