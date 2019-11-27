library(data.table)
library(purrr)

#####################
## Define settings ##
#####################

## Define I/O ##
io <- list()
io$basedir <- "/Users/C02RF23NFVH8/data/scnmt_eb"
io$sample.metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$met.dir <- paste0(io$basedir,"/met/feature_level")
io$acc.dir <- paste0(io$basedir,"/acc/feature_level")
io$annos_dir <- paste0(io$basedir,"/features/genomic_contexts")
io$gene_metadata <- paste0(io$basedir,"/features/genes/Mmusculus_genes_BioMart.87.txt")
io$outdir <- paste0(io$basedir,"/metacc/boxplots")

# Folders with the global statistics per cell
io$met.stats <- paste0(io$basedir,"/met/results/stats/samples/sample_stats.txt")
io$acc.stats <- paste0(io$basedir,"/acc/results/stats/samples/sample_stats.txt")

# Folders with the differential analysis results
io$basedir2 <- "/Users/C02RF23NFVH8/data/gastrulation"
io$diff.met <- paste0(io$basedir2,"/met/results/differential/feature_level")
io$diff.acc <- paste0(io$basedir2,"/acc/results/differential/feature_level")



## Define options ##
opts <- list()

# Define genomic contexts for methylation
opts$met.annos <- c(
  "H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers",
  "H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers"
  # "prom_2000_2000"="Promoters",
  # "genebody"="Gene bodies"
)

# Define genomic contexts for accessibility
opts$acc.annos <- c(
  "H3K27ac_distal_E7.5_Mes_intersect12"="Mesoderm enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="Endoderm enhancers",
  "H3K27ac_distal_E7.5_Ect_intersect12"="Ectoderm enhancers"
  # "prom_2000_2000"="Promoters",
  # "genebody"="Gene bodies"
)

# How to select differential hits?
#   Option 1 (more liberal): (lineage_A) vs (lineage_B,lineage_C)
#   Option 2 (more conservative): (lineage_A vs lineage_B) AND (lineageA vs lineage_C)
opts$diff.type <- 2
opts$min.fdr <- 0.10
opts$min.met.diff <- 5
opts$min.acc.diff <- 5


# Define which cells to use
opts$day_lineage <- c(
  # Day 2
  "Day2_Epiblast",
  "Day2_Primitive Streak",

  # Day 6/7
  "Day6-7_Mesoderm"
  # "Day6-7_Blood"
)

opts$genotype <- c(
  "WT",
  "KO"
)

# Define genomic context colors
opts$enhancer.colors <- c(
  "Mesoderm enhancers"="#CD3278",
  "Endoderm enhancers"="#43CD80",
  "Ectoderm enhancers"="steelblue",
  "Promoters"="gray70",
  "Gene bodies"="gray70"
)

# Define which cells to use
tmp <- fread(io$sample.metadata) %>%
  .[,day_lineage:=paste(day2,lineage10x_2, sep="_")] %>%
  .[day_lineage%in%opts$day_lineage] %>%
  .[genotype%in%opts$genotype]
opts$met_cells <- tmp %>% .[pass_metQC==T, id_met]
opts$acc_cells <- tmp %>% .[pass_accQC==T, id_acc]

# Load sample metadata
sample_metadata <- fread(io$sample.metadata) %>%
  .[,day_lineage:=paste(day,lineage10x_2, sep="_")] %>%
  .[,genotype:=factor(genotype, levels=c("WT","KO"))] %>%
  .[id_met%in%opts$met_cells | id_acc %in% opts$acc_cells ]

