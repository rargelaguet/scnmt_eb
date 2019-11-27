
#####################
## Define settings ##
#####################

## Define I/O ##
io <- list()
io$basedir <- "/Users/ricard/data/scnmt_eb"
io$sample.metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$met.dir <- paste0(io$basedir,"/met/cpg_level")
io$acc.dir <- paste0(io$basedir,"/acc/gpc_level")
io$annos_dir <- paste0(io$basedir,"/features/genomic_contexts")
io$gene_metadata <- paste0(io$basedir,"/features/genes/Mmusculus_genes_BioMart.87.txt")
io$outdir <- paste0(io$basedir,"/metacc/pseudobulk_profiles")

# Folders with the global statistics per cell
io$met.stats <- paste0(io$basedir,"/met/results/stats/samples/sample_stats.txt")
io$acc.stats <- paste0(io$basedir,"/acc/results/stats/samples/sample_stats.txt")

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

# Define window size 
opts$window_size <- 2000
opts$met.tile <- 200
opts$acc.tile <- 150

# Define which days and lineages to use
opts$day_lineage <- c(
  # Day 2
  "Day2_Epiblast",
  
  # Day 4/5
  "Day4-5_Epiblast",
  "Day4-5_Primitive Streak",
  # "Day4-5_Mesoderm",
  
  # Day 6/7
  # "Day6-7_Primitive Streak",
  "Day6-7_Mesoderm",
  # "Day6-7_Blood"
)

# Define which genotypes to use
opts$genotype <- c(
  "WT"
  # "KO"
)

# Define which cells to use
tmp <- fread(io$sample.metadata) %>%
  .[,day_lineage:=paste(day2,lineage10x_2, sep="_")] %>%
  .[day_lineage%in%opts$day_lineage] %>%
  .[genotype%in%opts$genotype]
opts$met.cells <- tmp %>% .[pass_metQC==T, id_met]
opts$acc.cells <- tmp %>% .[pass_accQC==T, id_acc]

f <- function(x) { return(data.frame(y=mean(x), ymin=mean(x)-sd(x), ymax=mean(x)+sd(x))) }