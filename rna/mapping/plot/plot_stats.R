
#######################################
## Script to plot mapping statistics ##
#######################################

library(data.table)
library(purrr)
library(ggplot2)
library(ggpubr)

source("/Users/ricard/scnmt_eb/rna/mapping/plot/plot_settings.R")

################
## Define I/O ##
################

io <- list()
io$sample_metadata <- "/Users/ricard/data/scnmt_eb/sample_metadata.txt"
io$outdir <- "/Users/ricard/data/scnmt_eb/rna/results/mapping/pdf"

####################
## Define options ##
####################

opts <- list()

# Figure dimensions
opts$width <- 3
opts$height <- 4

# colors
colors <- aggregated_celltype_colours

###############
## Load data ##
###############

sample_metadata <- fread(io$sample_metadata) %>%
  .[!lineage10x_2 %in% c("PGC","Allantois","NOIDEA",NA)]

#############################################
## Plot number of cells for each cell type ##
#############################################

foo <- sample_metadata %>%
  .[, Ncells:=.N, by=c("day","genotype")] %>%
  .[, .(prop=.N/unique(Ncells)), by=c("lineage10x_2","day2","genotype")] %>%
  .[day2=="Day4-5"]

xlim.max <- max(sample_metadata[,.N,by=c("lineage10x_2","day2","genotype")] %>% .[,N])


## Across all days ##

# to.plot <- sample_metadata %>%
#   .[day2 %in% c("Day2","Day4-5","Day6-7")] %>%
#   .[,.N,by=c("lineage10x_2","genotype","day2")]
# to.plot[,lineage10x_2:=factor(lineage10x_2,levels=rev(names(colors)))]
# to.plot[,genotype:=factor(genotype,levels=c("WT","KO"))]
# 
# p <- barplot.pub(to.plot, x="lineage10x_2", colors=colors, xlim.max=xlim.max) +
#   facet_wrap(~day+genotype, nrow=1, scales="free_x") +
#   theme(
#     strip.background = element_blank(),
#     strip.text = element_blank()
#   )
# 
# pdf(paste0(io$outdir,"/mapping_stats_all.pdf"), width=opts$width*2.5, height=opts$height)
# print(p)
# dev.off()

## Day 2, WT vs KO ##

to.plot <- sample_metadata %>%
  .[day2 %in% c("Day2")] %>%
  .[,.N,by=c("lineage10x_2","genotype")]
colors.tmp <- colors[names(colors) %in% to.plot$lineage10x_2]
to.plot[,lineage10x_2:=factor(lineage10x_2,levels=rev(names(colors)))]
to.plot[,genotype:=factor(genotype,levels=c("WT","KO"))]

p <- barplot.pub(to.plot, x="lineage10x_2", colors=colors, xlim.max=xlim.max) +
  facet_wrap(~genotype, nrow=1, scales="free_x") +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )

pdf(paste0(io$outdir,"/mapping_stats_day2.pdf"), width=opts$width, height=opts$height)
print(p)
dev.off()

## Day 4-5, WT vs KO ##

to.plot <- sample_metadata %>%
  .[day2 %in% c("Day4-5")] %>%
  .[,.N,by=c("lineage10x_2","genotype")]
colors.tmp <- colors[names(colors) %in% to.plot$lineage10x_2]
to.plot[,lineage10x_2:=factor(lineage10x_2,levels=rev(names(colors)))]
to.plot[,genotype:=factor(genotype,levels=c("WT","KO"))]

p <- barplot.pub(to.plot, x="lineage10x_2", colors=colors, xlim.max=xlim.max) +
  facet_wrap(~genotype, nrow=1, scales="free_x") +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )

pdf(paste0(io$outdir,"/mapping_stats_day45.pdf"), width=opts$width, height=opts$height)
print(p)
dev.off()


## Day 6-7, WT vs KO ##

to.plot <- sample_metadata %>%
  .[day %in% c("Day67")] %>%
  .[,.N,by=c("lineage10x_2","genotype")]
colors.tmp <- colors[names(colors) %in% to.plot$lineage10x_2]
to.plot[,lineage10x_2:=factor(lineage10x_2,levels=rev(names(colors)))]
to.plot[,genotype:=factor(genotype,levels=c("WT","KO"))]

p <- barplot.pub(to.plot, x="lineage10x_2", colors=colors, xlim.max=xlim.max) +
  facet_wrap(~genotype, nrow=1, scales="free_x") +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )

pdf(paste0(io$outdir,"/mapping_stats_day67.pdf"), width=opts$width, height=opts$height)
print(p)
dev.off()
