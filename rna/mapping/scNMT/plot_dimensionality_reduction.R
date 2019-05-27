
####################################################################
## Plot dimensionality reduction of scNMT-seq mapped to the atlas ##
####################################################################

library(SingleCellExperiment)
library(data.table)
library(purrr)
library(ggplot2)

# source("/Users/ricard/NMT-seq_EB_2/rna/mapping/scNMT/Mapping2gastrulationAtlas.R")

#####################
## I/O and options ##
#####################

io <- list()
io$path2atlas <- "/Users/ricard/data/gastrulation"
io$path2query <- "/Users/ricard/data/NMT-seq_EB+ESC"
io$outdir <- "/Users/ricard/data/NMT-seq_EB+ESC/mapping"

opts <- list()

################
## Load atlas ##
################

# Load atlas cell metadata
meta_atlas <- fread(paste0(io$path2atlas, "/sample_metadata.txt")) %>%
  .[pass_rnaQC==T] %>%
  .[,c("sample","id_rna","stage","lineage10x_2")]


# Load precomputed dimensionality reduction coordinates
atlas.coordinates <- fread("/Users/ricard/data/gastrulation_norsync_stuff/metaccrna/mofa/full/umap_coordinates.txt")

#####################
## Load query data ##
####################

# Load query cell metadata
meta_query <- fread(paste0(io$path2query, "/sample_metadata.txt")) %>%
  # .[,culture:=ifelse(culture=="ESC_2i","EB_Day0/1",culture)] %>%
  .[,culture:=ifelse(culture=="ESC_Serum","EB_Day1",culture)] %>%
  .[,culture:=ifelse(culture=="EB_Day3","EB_Day3",culture)] %>%
  # .[,culture:=ifelse(culture=="EB_Day5","EB_Day6/7",culture)] %>%
  .[,culture:=ifelse(culture=="EB_Day6","EB_Day6/7",culture)] %>%
  .[,culture:=ifelse(culture=="EB_Day7","EB_Day6/7",culture)]

# Load precomputed mapping
mapping <- readRDS("/Users/ricard/data/NMT-seq_EB+ESC/mapping/mapping_scNMT_2.rds")$mapping
mapping.dt <- data.table(
  id_rna            = mapping$cell, 
  celltype.mapped = mapping$celltype.mapped,
  stage.mapped    = mapping$stage.mapped,
  closest.cell    = as.character(mapping$closest.cell)
)

###################################
## Plot dimensionality reduction ##
###################################

# Prepare data.frame to plot
plot_df_atlas = atlas.coordinates %>% merge(meta_atlas, by=c("sample"))
plot_df_query = mapping.dt %>% merge(meta_query, by=c("id_rna"))

# index <- match(plot_df$id_rna, as.character(mapping$closest.cell) )
# plot_df$mapped <- as.factor(!is.na(index))


## Day 0 ##
# index <- match(plot_df_atlas$id_rna, plot_df_query[culture=="EB_Day0",closest.cell] )
# plot_df_atlas$mapped_day0 <- as.factor(!is.na(index))
# 
# p <- ggplot(data=plot_df_atlas, mapping=aes(x=V1, y=V2, colour=mapped_day0)) +
#   ggrastr::geom_point_rast(aes(size=mapped_day0), alpha=opts$dot_alpha) +
#   scale_size_manual(values = c("TRUE"=1.0, "FALSE"=0.5)) +
#   labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
#   scale_colour_manual(values = c("TRUE"="red", "FALSE"="grey")) +
#   guides(colour = guide_legend(override.aes = list(size=6))) +
#   theme_classic() +
#   theme(
#     legend.position = "none",
#     axis.text = element_blank(),
#     axis.ticks = element_blank()
#   )
# pdf(paste0(io$outdir,"/pdf/umap_mapped_day0.pdf"), width=4, height=4)
# print(p)
# dev.off()

## Day 1 ##
index <- match(plot_df_atlas$id_rna, plot_df_query[culture=="EB_Day1",closest.cell] )
plot_df_atlas$mapped <- as.factor(!is.na(index))

p <- ggplot(data=plot_df_atlas, mapping=aes(x=V1, y=V2)) +
  ggrastr::geom_point_rast(aes(size=mapped, alpha=mapped, colour=mapped)) +
  scale_size_manual(values = c("TRUE"=0.9, "FALSE"=0.4)) +
  scale_alpha_manual(values = c("TRUE"=1.0, "FALSE"=0.5)) +
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  scale_colour_manual(values = c("TRUE"="red", "FALSE"="grey")) +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

pdf(paste0(io$outdir,"/pdf/umap_mapped_day1.pdf"), width=4, height=4)
print(p)
dev.off()


## Day 3 ##
index <- match(plot_df_atlas$id_rna, plot_df_query[culture=="EB_Day3",closest.cell] )
plot_df_atlas$mapped <- as.factor(!is.na(index))

p <- ggplot(data=plot_df_atlas, mapping=aes(x=V1, y=V2)) +
  ggrastr::geom_point_rast(aes(size=mapped, alpha=mapped, colour=mapped)) +
  scale_size_manual(values = c("TRUE"=0.9, "FALSE"=0.4)) +
  scale_alpha_manual(values = c("TRUE"=1.0, "FALSE"=0.5)) +
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  scale_colour_manual(values = c("TRUE"="red", "FALSE"="grey")) +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

pdf(paste0(io$outdir,"/pdf/umap_mapped_day3.pdf"), width=4, height=4)
print(p)
dev.off()


## Day 5 ##
# index.wt <- match(plot_df_atlas$id_rna, plot_df_query[culture=="EB_Day5" & phenotype=="Tet_WT",closest.cell] )
# index.ko <- match(plot_df_atlas$id_rna, plot_df_query[culture=="EB_Day5" & phenotype=="Tet_KO",closest.cell] )
# plot_df_atlas$mapped_day5.wt <- c(0,-10)[as.numeric(as.factor(!is.na(index.wt)))]
# plot_df_atlas$mapped_day5.ko <- c(0,10)[as.numeric(as.factor(!is.na(index.ko)))]
# plot_df_atlas$mapped_day5 <- as.factor( plot_df_atlas$mapped_day5.wt + plot_df_atlas$mapped_day5.ko )
#   
# p <- ggplot(data=plot_df_atlas, mapping=aes(x=V1, y=V2, colour=mapped_day5)) +
#   ggrastr::geom_point_rast(aes(size=mapped_day5), alpha=opts$dot_alpha) +
#   scale_size_manual(values = c("-10"=1.0, "10"=1.0, "0"=0.5)) +
#   labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
#   scale_colour_manual(values = c("-10"="red", "10"="blue", "0"="grey")) +
#   guides(colour = guide_legend(override.aes = list(size=6))) +
#   theme_classic() +
#   theme(
#     legend.position = "none",
#     axis.text = element_blank(),
#     axis.ticks = element_blank()
#   )
# print(p)
# 
# pdf(paste0(io$outdir,"/pdf/umap_mapped_day5.pdf"), width=4, height=4)
# print(p)
# dev.off()


## Day 0/1/3 ##
# index <- match(plot_df_atlas$id_rna, plot_df_query[culture=="EB_Day0/1/3",closest.cell] )
# plot_df_atlas$mapped <- as.factor(!is.na(index))
# 
# p <- ggplot(data=plot_df_atlas, mapping=aes(x=V1, y=V2)) +
#   ggrastr::geom_point_rast(aes(size=mapped, alpha=mapped, colour=mapped)) +
#   scale_size_manual(values = c("TRUE"=0.9, "FALSE"=0.4)) +
#   scale_alpha_manual(values = c("TRUE"=1.0, "FALSE"=0.5)) +
#   labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
#   scale_colour_manual(values = c("TRUE"="red", "FALSE"="grey")) +
#   guides(colour = guide_legend(override.aes = list(size=6))) +
#   theme_classic() +
#   theme(
#     legend.position = "none",
#     axis.text = element_blank(),
#     axis.ticks = element_blank()
#   )
# 
# pdf(paste0(io$outdir,"/pdf/umap_mapped_day013.pdf"), width=4, height=4)
# print(p)
# dev.off()

## Day 6/7 ##

index.wt <- match(plot_df_atlas$id_rna, plot_df_query[culture=="EB_Day6/7" & phenotype=="Tet_WT",closest.cell] )
index.ko <- match(plot_df_atlas$id_rna, plot_df_query[culture=="EB_Day6/7" & phenotype=="Tet_KO",closest.cell] )
plot_df_atlas$mapped.wt <- c(0,-10)[as.numeric(as.factor(!is.na(index.wt)))]
plot_df_atlas$mapped.ko <- c(0,10)[as.numeric(as.factor(!is.na(index.ko)))]

# overlaps <- which(plot_df_atlas$mapped.wt==(-10) & plot_df_atlas$mapped.ko==10)
# plot_df_atlas$mapped.wt[overlaps] <- sample(x=c(0,-10), size=length(overlaps), replace=T)
# plot_df_atlas$mapped.ko[overlaps] <- sample(x=c(0,10), size=length(overlaps), replace=T)

plot_df_atlas$mapped <- as.factor( plot_df_atlas$mapped.wt + plot_df_atlas$mapped.ko )

p <- ggplot(data=plot_df_atlas, mapping=aes(x=V1, y=V2)) +
  ggrastr::geom_point_rast(aes(size=mapped, alpha=mapped, colour=mapped)) +
  scale_size_manual(values = c("-10"=0.9, "10"=0.9, "0"=0.4)) +
  scale_alpha_manual(values = c("-10"=1.0, "10"=1.0, "0"=0.5)) +
  scale_colour_manual(values = c("-10"="red", "10"="blue", "0"="grey")) +
  labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

pdf(paste0(io$outdir,"/pdf/umap_mapped_day67.pdf"), width=4, height=4)
print(p)
dev.off()