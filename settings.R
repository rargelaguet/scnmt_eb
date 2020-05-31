suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))

#########
## I/O ##
#########

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/scnmt_eb"
  io$gene.metadata <- "/Users/ricard/data/ensembl/"
} else if (grepl("ebi",Sys.info()['nodename'])) {
  io$basedir <- "/hps/nobackup2/research/stegle/users/ricard/scnmt_eb"
  io$gene.metadata <- "/hps/nobackup2/research/stegle/users/ricard/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
} else {
  stop("Computer not recognised")
}

io$metadata <- paste0(io$basedir,"/sample_metadata.txt")

io$met_data_raw <- paste0(io$basedir,"/met/cpg_level")
io$met_data_parsed <- paste0(io$basedir,"/met/feature_level")
io$met.stats <- paste0(io$basedir,"/met/stats/sample_stats.txt")

io$acc_data_raw <- paste0(io$basedir,"/acc/gpc_level")
io$acc_data_parsed <- paste0(io$basedir,"/acc/feature_level")
io$acc.stats <- paste0(io$basedir,"/acc/stats/sample_stats.txt")

io$rna <- paste0(io$basedir,"/rna/SingleCellExperiment.rds")

io$features.dir <- paste0(io$basedir,"/features/genomic_contexts")
# io$cpg.density <- paste0(io$basedir,"/met/stats/features/cpg_density_perfeature.txt.gz")


#############
## Options ##
#############

opts <- list()

# Define which cells to use
opts$days <- c(
  # "Day2",
  # "Day4",
  # "Day5",
  # "Day6",
  # "Day7"
  "Day2",
  "Day4-5",
  "Day6-7"
)

opts$day_lineage <- c(
  
  # Day 2
  "Day2_Epiblast",
  "Day2_Primitive Streak",
  
  # Day 4/5
  "Day4-5_Epiblast",
  "Day4-5_Primitive Streak",
  "Day4-5_Mesoderm",
  
  # Day 6/7
  "Day6-7_Mesoderm",
  "Day6-7_Endoderm",
  "Day6-7_Blood"
)

# Use both WT and KO cells
opts$genotype <- c(
  "WT",
  "KO"
)

# opts$celltype.colors <- c(
#   "Epiblast"="grey70",
#   "Mesoderm"="#CD3278",
#   "Primitive_Streak"="sandybrown",
#   "Endoderm"="#43CD80",
#   "Ectoderm"="steelblue",
#   "Epiblast/Ectoderm"="steelblue",
#   "Visceral_endoderm"="darkgreen"
# )

opts$celltype.colors = c(
  "Epiblast" = "#635547",
  "Primitive Streak" = "#DABE99",
  "Caudal epiblast" = "#9e6762",
  "PGC" = "#FACB12",
  "Anterior Primitive Streak" = "#c19f70",
  "Notochord" = "#0F4A9C",
  "Def. endoderm" = "#F397C0",
  "Gut" = "#EF5A9D",
  "Nascent mesoderm" = "#C594BF",
  "Mixed mesoderm" = "#DFCDE4",
  "Intermediate mesoderm" = "#139992",
  "Caudal Mesoderm" = "#3F84AA",
  "Paraxial mesoderm" = "#8DB5CE",
  "Somitic mesoderm" = "#005579",
  "Pharyngeal mesoderm" = "#C9EBFB",
  "Cardiomyocytes" = "#B51D8D",
  "Allantois" = "#532C8A",
  "ExE mesoderm" = "#8870ad",
  "Mesenchyme" = "#cc7818",
  "Haematoendothelial progenitors" = "#FBBE92",
  "Endothelium" = "#ff891c",
  "Blood progenitors 1" = "#f9decf",
  "Blood progenitors 2" = "#c9a997",
  "Erythroid1" = "#C72228",
  "Erythroid2" = "#f79083",
  "Erythroid3" = "#EF4E22",
  "NMP" = "#8EC792",
  "Rostral neurectoderm" = "#65A83E",
  "Caudal neurectoderm" = "#354E23",
  "Neural crest" = "#C3C388",
  "Forebrain/Midbrain/Hindbrain" = "#647a4f",
  "Spinal cord" = "#CDE088",
  "Surface ectoderm" = "#f7f79e",
  "Visceral endoderm" = "#F6BFCB",
  "ExE endoderm" = "#7F6874",
  "ExE ectoderm" = "#989898",
  "Parietal endoderm" = "#1A1A1A"
)

opts$celltype2.colors <- c(
  "Epiblast" = "#63B8FF",
  "Mesoderm" = "#CD3278",
  "Primitive Streak"="sandybrown",
  "Endoderm" = "#43CD80",
  "Ectoderm" = "steelblue",
  "Blood" = "darkred"
)


##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>% 
  .[,day_lineage:=paste(day2,lineage10x_2, sep="_")]
  