#####################
## Define settings ##
#####################

## Define I/O ##
io <- list()
io$basedir <- "/Users/ricard/data/NMT-seq_EB+ESC"
io$sample.metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$met.dir <- paste0(io$basedir,"/met/feature_level")
io$acc.dir <- paste0(io$basedir,"/acc/feature_level")
io$annos_dir <- paste0(io$basedir,"/features/filt")
io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
io$outdir <- "/Users/ricard/NMT-seq_EB+ESC/metacc/boxplots/out"

# Folders with the differential analysis results
io$basedir2 <- "/Users/ricard/data/gastrulation"
io$diff.met <- paste0(io$basedir2,"/met/differential/feature_level")
io$diff.acc <- paste0(io$basedir2,"/acc/differential/feature_level")


# Folders with the global statistics per cell
io$met.stats <- paste0(io$basedir,"/met/stats/samples/sample_stats.txt")
io$acc.stats <- paste0(io$basedir,"/acc/stats/samples/sample_stats.txt")

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
  
  # Day 4/5
  # "Day4_Epiblast",
  # "Day4_Primitive Streak",
  # "Day4_Mesoderm",
  # "Day5_Epiblast",
  # "Day5_Primitive Streak",
  # "Day5_Mesoderm",
  
  # Day 6/7
  "Day6_Mesoderm",
  "Day7_Mesoderm"
)

opts$genotype <- c(
  "WT",
  "KO"
)

# Define colors
opts$lineage.colors <- c(
  Epiblast="#63B8FF",
  EB_Day0_Epiblast="#63B8FF",
  EB_Day3_Epiblast="#63B8FF",
  EB_Day5_Epiblast="#63B8FF",
  EB_Day6_Epiblast="#63B8FF",
  EB_Day7_Epiblast="#63B8FF",
  "Primitive Streak"="sandybrown",
  Mesoderm="#CD3278",
  Blood="#CD3278",
  Endoderm="#43CD80",
  Ectoderm="steelblue",
  "Epi/PS"="steelblue"
)

opts$enhancer.colors <- c(
  "Mesoderm enhancers"="#CD3278",
  "Endoderm enhancers"="#43CD80",
  "Ectoderm enhancers"="steelblue",
  "Promoters"="gray70",
  "Gene bodies"="gray70"
)

tmp <- fread(io$sample.metadata) %>%
  .[,day_lineage:=paste(day,lineage10x_2, sep="_")] %>%
  .[day_lineage%in%opts$day_lineage] %>%
  .[genotype%in%opts$genotype]
opts$met_cells <- tmp %>% .[pass_metQC==T, id_met]
opts$acc_cells <- tmp %>% .[pass_accQC==T, id_acc]

sample_metadata <- fread(io$sample.metadata) %>%
  .[id_met%in%opts$met_cells | id_acc %in% opts$acc_cells ] %>%
  .[,c("sample","id_rna","id_met","id_acc","lineage10x","lineage10x_2","day","genotype")] %>%
  .[,day_lineage:=paste(day,lineage10x_2, sep="_")] %>%
  droplevels()

