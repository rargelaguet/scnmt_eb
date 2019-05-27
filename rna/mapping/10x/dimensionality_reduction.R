
####################################################################
## Plot dimensionality reduction of EB cells mapped to the atlas ##
####################################################################

library(SingleCellExperiment)
library(data.table)
library(purrr)
library(ggplot2)

source("/Users/ricard/NMT-seq_EB+ESC/rna/mapping/settings.R")

#####################
## I/O and options ##
#####################

io <- list()
io$path2atlas <- "/Users/ricard/data/gastrulation10x"
io$path2query <- "/Users/ricard/data/NMT-seq_EB+ESC"
io$outdir <- "/Users/ricard/data/NMT-seq_EB+ESC/rna/mapping"

opts <- list()

# Dot size
opts$size.mapped <- 0.6
opts$size.nomapped <- 0.1

# Transparency
opts$alpha.mapped <- 1.0
opts$alpha.nomapped <- 0.5

##########################
## Load scNMT-seq atlas ##
##########################

# Load atlas cell metadata
# meta_atlas <- fread(paste0(io$path2atlas, "/sample_metadata.txt")) %>%
#   .[pass_rnaQC==T] %>%
#   .[,c("sample","id_rna","stage","lineage10x_2")]
# # Load precomputed dimensionality reduction coordinates
# atlas.coordinates <- fread("/Users/ricard/data/gastrulation_norsync_stuff/metaccrna/mofa/full/umap_coordinates.txt")

####################
## Load 10x atlas ##
####################

# Load atlas cell metadata
# meta_atlas <- readRDS(paste0(io$path2atlas, "/sample_metadata_atlas.rds")) %>% as.data.table 
meta_atlas <- fread(paste0(io$path2atlas, "/sample_metadata.txt"))# %>%
  # .[,aggregated_celltype:=stringr::str_replace_all(celltype,aggregated_celltypes)]

# Load precomputed dimensionality reduction coordinates
# pca <- readRDS(paste0(io$path2atlas, "/dimensionality_reduction/corrected_pcas.rds"))$all
# umap <- readRDS(paste0(io$path2atlas,"/umap_coordinates.rds"))
umap <- meta_atlas[,c("cell","umapX","umapY")] %>%
  setnames(c("umapX","umapY"),c("V1","V2"))

#####################
## Load query data ##
#####################

# Load query cell metadata
meta_query <- fread(paste0(io$path2query, "/sample_metadata.txt")) %>%
  .[,day:=ifelse(day%in%c("Day4","Day5"),"Day45",day)] %>%
  .[,day:=ifelse(day%in%c("Day6","Day7"),"Day67",day)] %>%
  .[!lineage10x_2 %in% c("PGC","Allantois","NOIDEA",NA)]

# Load precomputed mapping
mapping <- readRDS("/Users/ricard/data/NMT-seq_EB+ESC/rna/mapping/mapping.rds")$mapping
mapping.dt <- data.table(
  # id_rna            = mapping$cell, 
  id_rna            = stringr::str_replace_all(mapping$cell,"map_",""), 
  celltype.mapped = mapping$celltype.mapped,
  stage.mapped    = mapping$stage.mapped,
  closest.cell    = as.character(mapping$closest.cell)
)

################
## Parse data ##
################

# plot_df = as.data.frame(umap)
# plot_df$cell = meta_atlas$cell
# plot_df$sample = meta_atlas$sample
# plot_df$stage = meta_atlas$stage
# plot_df$cluster = meta_atlas$cluster
# plot_df$celltype = meta_atlas$celltype
# plot_df = plot_df[sample(nrow(plot_df), nrow(plot_df), replace = FALSE), ]

# plot_df <- plot_df %>% .[complete.cases(.),]

# Remove some lineages
# plot_df <- plot_df[plot_df$celltype!="ExE ectoderm",]

###################################
## Plot dimensionality reduction ##
###################################

# Prepare query data.frame to plot
plot_df_query = mapping.dt %>% merge(meta_query, by=c("id_rna"))

# Prepare atlas data.frame to plot
plot_df_atlas = umap %>% merge(meta_atlas, by=c("cell"))
plot_df_atlas <- plot_df_atlas[celltype!="ExE ectoderm"]

## Day 1 ##
# plot_df_atlas[,index:=match(plot_df_atlas$cell, plot_df_query[day=="Day1",closest.cell] )]
# plot_df_atlas[,mapped:=as.factor(!is.na(index))] %>% setorder(mapped)
# p <- plot.dimred(plot_df_atlas)
# pdf(paste0(io$outdir,"/pdf/umap_mapped_day1.pdf"), width=4, height=4)
# print(p)
# dev.off()

## Day 2 ##
plot_df_atlas[,index.wt:=match(plot_df_atlas$cell, plot_df_query[day=="Day2" & genotype=="WT",closest.cell] )]
plot_df_atlas[,index.ko:=match(plot_df_atlas$cell, plot_df_query[day=="Day2" & genotype=="KO",closest.cell] )]
plot_df_atlas[,mapped.wt:=c(0,-10)[as.numeric(as.factor(!is.na(index.wt)))]]
plot_df_atlas[,mapped.ko:=c(0,10)[as.numeric(as.factor(!is.na(index.ko)))]]
plot_df_atlas[,mapped:=factor(mapped.wt + mapped.ko, levels=c("0","-10","10"))] %>% setorder(mapped)

p <- plot.dimred.wtko(plot_df_atlas)

pdf(paste0(io$outdir,"/pdf/umap_mapped_day2.pdf"), width=4, height=4)
print(p)
dev.off()

## Day 3 ##
# plot_df_atlas[,index:=match(plot_df_atlas$cell, plot_df_query[day=="Day3",closest.cell] )]
# plot_df_atlas[,mapped:=as.factor(!is.na(index))] %>% setorder(mapped)
# 
# p <- plot.dimred(plot_df_atlas)
# 
# pdf(paste0(io$outdir,"/pdf/umap_mapped_day3.pdf"), width=4, height=4)
# print(p)
# dev.off()



## Day 4/5 ##
plot_df_atlas[,index.wt:=match(plot_df_atlas$cell, plot_df_query[day=="Day45" & genotype=="WT",closest.cell] )]
plot_df_atlas[,index.ko:=match(plot_df_atlas$cell, plot_df_query[day=="Day45" & genotype=="KO",closest.cell] )]
plot_df_atlas[,mapped.wt:=c(0,-10)[as.numeric(as.factor(!is.na(index.wt)))]]
plot_df_atlas[,mapped.ko:=c(0,10)[as.numeric(as.factor(!is.na(index.ko)))]]
plot_df_atlas[,mapped:=factor(mapped.wt + mapped.ko, levels=c("0","-10","10"))] %>% setorder(mapped)

p <- plot.dimred.wtko(plot_df_atlas)

pdf(paste0(io$outdir,"/pdf/umap_mapped_day45.pdf"), width=4, height=4)
print(p)
dev.off()



## Day 6/7 ##
plot_df_atlas[,index.wt:=match(plot_df_atlas$cell, plot_df_query[day=="Day67" & genotype=="WT",closest.cell] )]
plot_df_atlas[,index.ko:=match(plot_df_atlas$cell, plot_df_query[day=="Day67" & genotype=="KO",closest.cell] )]
plot_df_atlas[,mapped.wt:=c(0,-10)[as.numeric(as.factor(!is.na(index.wt)))]]
plot_df_atlas[,mapped.ko:=c(0,10)[as.numeric(as.factor(!is.na(index.ko)))]]
plot_df_atlas[,mapped:=factor(mapped.wt + mapped.ko, levels=c("0","-10","10"))] %>% setorder(mapped)

# overlaps <- which(plot_df_atlas$mapped.wt==(-10) & plot_df_atlas$mapped.ko==10)
# plot_df_atlas$mapped.wt[overlaps] <- sample(x=c(0,-10), size=length(overlaps), replace=T)
# plot_df_atlas$mapped.ko[overlaps] <- sample(x=c(0,10), size=length(overlaps), replace=T)

p <- plot.dimred.wtko(plot_df_atlas)

pdf(paste0(io$outdir,"/pdf/umap_mapped_day67.pdf"), width=4, height=4)
print(p)
dev.off()
