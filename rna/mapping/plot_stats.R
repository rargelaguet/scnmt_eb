library(data.table)
library(purrr)
library(ggplot2)
library(ggpubr)

source("/Users/ricard/NMT-seq_EB+ESC/rna/mapping/settings.R")

####################
## Define options ##
####################

io <- list()
io$sample_metadata <- "/Users/ricard/data/NMT-seq_EB+ESC/sample_metadata.txt"
io$outdir <- "/Users/ricard/data/NMT-seq_EB+ESC/rna/mapping"

opts <- list()
opts$colors <- c(
  "Epiblast" = "grey70",
  "Ectoderm" = "steelblue",
  "Primitive Streak" = "sandybrown",
  "Mesoderm" = "violetred",
  "Endoderm" = "#43CD80",
  "Blood" = "darkred"
)

opts$width <- 3
opts$height <- 4

###############
## Load data ##
###############

sample_metadata <- fread(io$sample_metadata) %>%
  .[,day:=ifelse(day%in%c("Day4","Day5"),"Day45",day)] %>%
  .[,day:=ifelse(day%in%c("Day6","Day7"),"Day67",day)] %>%
  .[!lineage10x_2 %in% c("PGC","Allantois","NOIDEA",NA)]

#############################################
## Plot number of cells for each cell type ##
#############################################

foo <- sample_metadata %>%
  .[, Ncells:=.N, by=c("day","genotype")] %>%
  .[, .(prop=.N/unique(Ncells)), by=c("lineage10x_2","day","genotype")] %>%
  .[day=="Day45"]

xlim.max <- max(sample_metadata[,.N,by=c("lineage10x_2","day","genotype")] %>% .[,N])

## ALL ##

to.plot <- sample_metadata %>%
  .[day %in% c("Day2","Day45","Day67")] %>%
  .[,.N,by=c("lineage10x_2","genotype","day")]
to.plot[,lineage10x_2:=factor(lineage10x_2,levels=rev(names(opts$colors)))]
to.plot[,genotype:=factor(genotype,levels=c("WT","KO"))]

p <- barplot.pub(to.plot, x="lineage10x_2", colors=opts$colors, xlim.max=xlim.max) +
  facet_wrap(~day+genotype, nrow=1, scales="free_x") +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )

pdf(paste0(io$outdir,"/pdf/mapping_stats_all.pdf"), width=opts$width*2.5, height=opts$height)
print(p)
dev.off()

# pieplot.pub(to.plot, x="lineage10x_2", colors=opts$colors) +
#   facet_wrap(~genotype+day, nrow=2, scales="free") +
#   theme(
#     legend.position = "none"
#   )


## Day 1 ##

# to.plot <- sample_metadata %>%
#   .[day %in% c("Day1")] %>%
#   .[,.N,by=c("day","lineage10x_2")]
# colors.tmp <- opts$colors[names(opts$colors) %in% to.plot$lineage10x_2]
# to.plot[,lineage10x_2:=factor(lineage10x_2,levels=rev(names(opts$colors)))]
# 
# p <- barplot.pub(to.plot, x="lineage10x_2", color=opts$colors)
# 
# pdf(paste0(io$outdir,"/mapping_stats_day1.pdf"), width=4, height=opts$height)
# print(p)
# dev.off()

## Day 3 ##

# to.plot <- sample_metadata %>%
#   .[day %in% c("Day3")] %>%
#   .[,.N,by=c("day","lineage10x_2")]
# colors.tmp <- opts$colors[names(opts$colors) %in% to.plot$lineage10x_2]
# to.plot[,lineage10x_2:=factor(lineage10x_2,levels=rev(names(opts$colors)))]
# 
# p <- barplot.pub(to.plot, x="lineage10x_2", color=opts$colors)
# 
# pdf(paste0(io$outdir,"/pdf/mapping_stats_day3.pdf"), width=opts$width, height=opts$height)
# print(p)
# dev.off()

## Day 2 WT vs KO ##

to.plot <- sample_metadata %>%
  .[day %in% c("Day2")] %>%
  .[,.N,by=c("lineage10x_2","genotype")]
colors.tmp <- opts$colors[names(opts$colors) %in% to.plot$lineage10x_2]
to.plot[,lineage10x_2:=factor(lineage10x_2,levels=rev(names(opts$colors)))]
to.plot[,genotype:=factor(genotype,levels=c("WT","KO"))]

p <- barplot.pub(to.plot, x="lineage10x_2", colors=opts$colors, xlim.max=xlim.max) +
  facet_wrap(~genotype, nrow=1, scales="free_x") +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )

pdf(paste0(io$outdir,"/pdf/mapping_stats_day2.pdf"), width=opts$width, height=opts$height)
print(p)
dev.off()

## Day 4/5 WT vs KO ##

to.plot <- sample_metadata %>%
  .[day %in% c("Day45")] %>%
  .[,.N,by=c("lineage10x_2","genotype")]
colors.tmp <- opts$colors[names(opts$colors) %in% to.plot$lineage10x_2]
to.plot[,lineage10x_2:=factor(lineage10x_2,levels=rev(names(opts$colors)))]
to.plot[,genotype:=factor(genotype,levels=c("WT","KO"))]

p <- barplot.pub(to.plot, x="lineage10x_2", colors=opts$colors, xlim.max=xlim.max) +
  facet_wrap(~genotype, nrow=1, scales="free_x") +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )

pdf(paste0(io$outdir,"/pdf/mapping_stats_day45.pdf"), width=opts$width, height=opts$height)
print(p)
dev.off()


## Day 6 and 7 ##

to.plot <- sample_metadata %>%
  .[day %in% c("Day67")] %>%
  .[,.N,by=c("lineage10x_2","genotype")]
colors.tmp <- opts$colors[names(opts$colors) %in% to.plot$lineage10x_2]
to.plot[,lineage10x_2:=factor(lineage10x_2,levels=rev(names(opts$colors)))]
to.plot[,genotype:=factor(genotype,levels=c("WT","KO"))]

p <- barplot.pub(to.plot, x="lineage10x_2", colors=opts$colors, xlim.max=xlim.max) +
  facet_wrap(~genotype, nrow=1, scales="free_x") +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )

pdf(paste0(io$outdir,"/pdf/mapping_stats_day67.pdf"), width=opts$width, height=opts$height)
print(p)
dev.off()



#########################################
## Plot number of cells for each stage ##
#########################################

# to.plot <- sample_metadata %>%
#   .[,stage.mapped:=factor(stage.mapped)] %>%
#   .[,.N,by=c("stage.mapped","genotype","day")] %>%
#   .[complete.cases(.)]
# 
# 
# p <- barplot.pub(to.plot[day%in%c("Day1")], x="stage.mapped", colors=NULL) +
#   facet_wrap(~genotype, nrow=1, scales="free_x") +
#   theme(
#     strip.background = element_blank(),
#     strip.text = element_blank()
#   )
# pdf(paste0(io$outdir,"/pdf/mapping_stages_day1.pdf"), width=3, height=6)
# print(p)
# dev.off()
# 
# p <- barplot.pub(to.plot[day%in%c("Day3")], x="stage.mapped", colors=NULL) +
#   facet_wrap(~genotype, nrow=1, scales="free_x") +
#   theme(
#     strip.background = element_blank(),
#     strip.text = element_blank()
#   )
# pdf(paste0(io$outdir,"/pdf/mapping_stages_day3.pdf"), width=3, height=6)
# print(p)
# dev.off()
# 
# p <- barplot.pub(to.plot[day%in%c("Day6/7")], x="stage.mapped", colors=NULL) +
#   facet_wrap(~genotype, nrow=1, scales="free_x") +
#   theme(
#     strip.background = element_blank(),
#     strip.text = element_blank()
#   )
# pdf(paste0(io$outdir,"/pdf/mapping_stages_day67.pdf"), width=3, height=6)
# print(p)
# dev.off()