
#####################
## Define settings ##
#####################

## Define I/O ##
io <- list()
io$basedir <- "/Users/ricard/data/NMT-seq_EB+ESC"
io$sample.metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$met.dir <- paste0(io$basedir,"/met/cpg_level")
io$acc.dir <- paste0(io$basedir,"/acc/gpc_level")
io$annos_dir <- paste0(io$basedir,"/features/filt")
io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
io$outdir <- "/Users/ricard/NMT-seq_EB+ESC/metacc/profiles/out"

# Folders with the differential analysis results
# io$diff.met <- paste0(io$basedir,"/met/differential/feature_level")
# io$diff.acc <- paste0(io$basedir,"/acc/differential/feature_level")

# Folders with the global statistics per cell
io$met.stats <- paste0(io$basedir,"/met/stats/samples/sample_stats.txt")
io$acc.stats <- paste0(io$basedir,"/acc/stats/samples/sample_stats.txt")

## Define options ##
opts <- list()


# Define genomic contexts for methylation
opts$annos <- c(
  "H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers",
  "H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers"
)

# Define window positions
opts$positions <- c(
  "H3K27ac_distal_E7.5_Ect_intersect12"="center",
  "H3K27ac_distal_E7.5_Mes_intersect12"="center",
  "H3K27ac_distal_E7.5_End_intersect12"="center"
)
opts$window_size <- 2000
opts$met.tile <- 200
opts$acc.tile <- 150



# How to select differential hits?
#   Option 1 (more liberal): (lineage_A) vs (lineage_B,lineage_C)
#   Option 2 (more conservative): (lineage_A vs lineage_B) AND (lineageA vs lineage_C)
# opts$diff.type <- 2
# opts$min.fdr <- 0.10
# opts$min.met.diff <- 5
# opts$min.acc.diff <- 5


# Define which cells to use
opts$day_lineage <- c(
  # Day 2
  "Day2_Epiblast",
  # "Day2_Primitive Streak",
  
  # Day 4/5
  "Day4_Epiblast",
  "Day4_Primitive Streak",
  "Day4_Mesoderm",
  "Day5_Epiblast",
  "Day5_Primitive Streak",
  "Day5_Mesoderm",
  
  # Day 6/7
  "Day6_Primitive Streak",
  "Day6_Mesoderm",
  "Day7_Primitive Streak",
  "Day7_Mesoderm"
)

opts$genotype <- c(
  "WT"
  # "KO"
)

# Define colors
opts$colors <- c(
  Epiblast="#63B8FF",
  EB_Day0_Epiblast="#63B8FF",
  EB_Day1_Epiblast="#63B8FF",
  EB_Day3_Epiblast="#63B8FF",
  EB_Day5_Epiblast="#63B8FF",
  Primitive_Streak="sandybrown",
  Mesoderm="#CD3278",
  Endoderm="#43CD80"
  # Ectoderm="steelblue"
)

tmp <- fread(io$sample.metadata) %>%
  .[,day_lineage:=paste(day,lineage10x_2, sep="_")] %>%
  .[day_lineage%in%opts$day_lineage] %>%
  .[genotype%in%opts$genotype]
opts$met.cells <- tmp %>% .[pass_metQC==T, id_met]
opts$acc.cells <- tmp %>% .[pass_accQC==T, id_acc]


theme_boxplot <- function() {
  theme(
    plot.margin = unit(c(t=1,r=1,b=1,l=1), "cm"),
    plot.title = element_text(size=25,hjust=0.5),
    axis.text=element_text(size=15, colour="black"),
    axis.title.x=element_text(size=17, margin=margin(10,0,0,0)),
    axis.title.y=element_text(size=17, margin=margin(0,10,0,0)),
    axis.line = element_line(size=rel(1.0)),
    axis.ticks = element_line(size=rel(1.3), color="black"),
    legend.key = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    # legend.key.size= unit(0.5, "cm"),
    legend.key.width=unit(1.2,"line"),
    legend.key.height=unit(1.0,"line"),
    legend.margin = margin(t=10, r=0, b=0, l=0, unit="pt"),
    legend.title = element_blank(),
    legend.text = element_text(size=15),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank()
  )
}

f <- function(x) { return(data.frame(y=mean(x), ymin=mean(x)-sd(x), ymax=mean(x)+sd(x))) }