###################################################################
## Plot dimensionality reduction of EB cells mapped to the atlas ##
###################################################################

# This script requires the cell metadata from the atlas, which contains the precomputed UMAP coordinates

library(SingleCellExperiment)
library(data.table)
library(purrr)
library(ggplot2)

source("/Users/C02RF23NFVH8/scnmt_eb/rna/mapping/plot/plot_settings.R")

#########
## I/O ##
#########

io <- list()
io$path2atlas <- "/Users/C02RF23NFVH8/data/gastrulation10x"
io$path2query <- "/Users/C02RF23NFVH8/data/scnmt_eb"
io$mapping <- "/Users/C02RF23NFVH8/data/scnmt_eb/rna/results/mapping/mapping.rds"
io$outdir <- "/Users/C02RF23NFVH8/data/scnmt_eb/rna/results/mapping/pdf"

#############
## Options ##
#############

opts <- list()

# Dot size
opts$size.mapped <- 0.6
opts$size.nomapped <- 0.1

# Transparency
opts$alpha.mapped <- 1.0
opts$alpha.nomapped <- 0.5

####################
## Load 10x atlas ##
####################

# Load atlas cell metadata
meta_atlas <- fread(paste0(io$path2atlas, "/sample_metadata.txt"))

# Extract precomputed dimensionality reduction coordinates
umap <- meta_atlas[,c("cell","umapX","umapY")] %>%
  setnames(c("umapX","umapY"),c("V1","V2"))

#####################
## Load query data ##
#####################

# Load query cell metadata
meta_query <- fread(paste0(io$path2query, "/sample_metadata.txt")) %>%
  .[!lineage10x_2 %in% c("PGC","Allantois","NOIDEA",NA)] # remove weird lineages 

# Load precomputed mapping
mapping <- readRDS(io$mapping)$mapping
mapping.dt <- data.table(
  # id_rna            = mapping$cell, 
  id_rna            = stringr::str_replace_all(mapping$cell,"map_",""), 
  celltype.mapped = mapping$celltype.mapped,
  stage.mapped    = mapping$stage.mapped,
  closest.cell    = as.character(mapping$closest.cell)
)

###################################
## Plot dimensionality reduction ##
###################################

# Prepare query data.frame to plot
plot_df_query = mapping.dt %>% merge(meta_query, by=c("id_rna"))

# Prepare atlas data.frame to plot
plot_df_atlas = umap %>% merge(meta_atlas, by=c("cell"))
plot_df_atlas <- plot_df_atlas[celltype!="ExE ectoderm"]

## Day 2 ##
plot_df_atlas[,index.wt:=match(plot_df_atlas$cell, plot_df_query[day=="Day2" & genotype=="WT",closest.cell] )]
plot_df_atlas[,index.ko:=match(plot_df_atlas$cell, plot_df_query[day=="Day2" & genotype=="KO",closest.cell] )]
plot_df_atlas[,mapped.wt:=c(0,-10)[as.numeric(as.factor(!is.na(index.wt)))]]
plot_df_atlas[,mapped.ko:=c(0,10)[as.numeric(as.factor(!is.na(index.ko)))]]
plot_df_atlas[,mapped:=factor(mapped.wt + mapped.ko, levels=c("0","-10","10"))] %>% setorder(mapped)

p <- plot.dimred.wtko(plot_df_atlas)

pdf(paste0(io$outdir,"/umap_mapped_day2.pdf"), width=4, height=4)
print(p)
dev.off()

## Day 4-5 ##
plot_df_atlas[,index.wt:=match(plot_df_atlas$cell, plot_df_query[day2=="Day4-5" & genotype=="WT",closest.cell] )]
plot_df_atlas[,index.ko:=match(plot_df_atlas$cell, plot_df_query[day2=="Day4-5" & genotype=="KO",closest.cell] )]
plot_df_atlas[,mapped.wt:=c(0,-10)[as.numeric(as.factor(!is.na(index.wt)))]]
plot_df_atlas[,mapped.ko:=c(0,10)[as.numeric(as.factor(!is.na(index.ko)))]]
plot_df_atlas[,mapped:=factor(mapped.wt + mapped.ko, levels=c("0","-10","10"))] %>% setorder(mapped)

p <- plot.dimred.wtko(plot_df_atlas)

pdf(paste0(io$outdir,"/umap_mapped_day45.pdf"), width=4, height=4)
print(p)
dev.off()



## Day 6-7 ##
plot_df_atlas[,index.wt:=match(plot_df_atlas$cell, plot_df_query[day=="Day6-7" & genotype=="WT",closest.cell] )]
plot_df_atlas[,index.ko:=match(plot_df_atlas$cell, plot_df_query[day=="Day6-7" & genotype=="KO",closest.cell] )]
plot_df_atlas[,mapped.wt:=c(0,-10)[as.numeric(as.factor(!is.na(index.wt)))]]
plot_df_atlas[,mapped.ko:=c(0,10)[as.numeric(as.factor(!is.na(index.ko)))]]
plot_df_atlas[,mapped:=factor(mapped.wt + mapped.ko, levels=c("0","-10","10"))] %>% setorder(mapped)

p <- plot.dimred.wtko(plot_df_atlas)

pdf(paste0(io$outdir,"/umap_mapped_day67.pdf"), width=4, height=4)
print(p)
dev.off()

## All days ##
# plot_df_atlas[,index.wt:=match(plot_df_atlas$cell, plot_df_query[genotype=="WT",closest.cell] )]
# plot_df_atlas[,index.ko:=match(plot_df_atlas$cell, plot_df_query[genotype=="KO",closest.cell] )]
# plot_df_atlas[,mapped.wt:=c(0,-10)[as.numeric(as.factor(!is.na(index.wt)))]]
# plot_df_atlas[,mapped.ko:=c(0,10)[as.numeric(as.factor(!is.na(index.ko)))]]
# plot_df_atlas[,mapped:=factor(mapped.wt + mapped.ko, levels=c("0","-10","10"))] %>% setorder(mapped)
# 
# p <- plot.dimred.wtko(plot_df_atlas)
# 
# pdf(paste0(io$outdir,"/umap_mapped_alldays.pdf"), width=8, height=6.5)
# print(p)
# dev.off()
