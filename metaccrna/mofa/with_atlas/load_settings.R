library(data.table)
library(purrr)

## Define I/O ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir.1 <- "/Users/ricard/data/NMT-seq_EB+ESC"
  io$basedir.2 <- "/Users/ricard/data/gastrulation"
  io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  io$outdir <- "/Users/ricard/data/NMT-seq_EB+ESC/mofa"
} else {
  io$basedir.1 <- "/hps/nobackup/stegle/users/ricard/NMT-seq_EB+ESC"
  io$basedir.2 <- "/hps/nobackup/stegle/users/ricard/gastrulation"
  io$gene_metadata <- "/hps/nobackup/stegle/users/ricard/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  io$outdir <- "/homes/ricard/NMT-seq_EB+ESC/metaccrna/mofa"
}
io$sample.metadata.1 <- paste0(io$basedir.1,"/sample_metadata.txt")
io$sample.metadata.2 <- paste0(io$basedir.2,"/sample_metadata.txt")
io$met.dir.1 <- paste0(io$basedir.1,"/met/feature_level")
io$met.dir.2 <- paste0(io$basedir.2,"/met/parsed")
io$acc.dir.1 <- paste0(io$basedir.1,"/acc/feature_level")
io$acc.dir.2 <- paste0(io$basedir.2,"/acc/parsed")
io$rna.file.1 <- paste0(io$basedir.1,"/rna/SingleCellExperiment.rds")
io$rna.file.2 <- paste0(io$basedir.2,"/rna/parsed/SingleCellExperiment.rds")

io$annos_dir  <- paste0(io$basedir, "/features/filt")

## Define options ##
opts <- list()

# Define which annotations to look at
opts$met.annos <- c(
  # "genebody",
  # "prom_2000_2000",
  "H3K27ac_distal_E7.5_Mes_intersect12",
  "H3K27ac_distal_E7.5_Ect_intersect12",
  "H3K27ac_distal_E7.5_End_intersect12"
  # "H3K4me3_E7.5_Mes",
  # "H3K4me3_E7.5_End",
  # "H3K4me3_E7.5_Ect",
)

opts$acc.annos <- c(
  # "genebody",
  # "prom_200_200",
  "H3K27ac_distal_E7.5_Mes_intersect12",
  "H3K27ac_distal_E7.5_Ect_intersect12",
  "H3K27ac_distal_E7.5_End_intersect12"
  # "H3K4me3_E7.5_Mes",
  # "H3K4me3_E7.5_End",
  # "H3K4me3_E7.5_Ect",
)


## Define which cells to use for EB data ##
opts$culture <- c(
  "EB_Day5",
  "EB_Day3"
)
# opts$celltype.mapped <- c(
#   "Epiblast",
#   "Primitive_Streak",
#   "Mature_mesoderm",
#   "Nascent_mesoderm",
#   "Embryonic_endoderm"
# )

# opts$stage.mapped <- c(
#   "E4.5",
#   "E5.5",
#   "E6.5",
#   "E7.5"
# )

## Define which cells to use for in vivo data ##
opts$stage_lineage <- c(
  
  "E4.5_Epiblast",
  "E4.5_Primitive_endoderm",
  
  "E5.5_Epiblast",
  "E5.5_Visceral_endoderm",
  
  "E6.5_Epiblast",
  "E6.5_Primitive_Streak",
  "E6.5_Nascent_mesoderm",
  "E6.5_Visceral_endoderm",
  
  "E7.5_Epiblast",
  "E7.5_Ectoderm",
  "E7.5_Nascent_mesoderm",
  "E7.5_Mature_mesoderm",
  "E7.5_Primitive_Streak",
  "E7.5_Embryonic_endoderm",
  "E7.5_Visceral_endoderm"
  
)

# Filtering options for methylation
opts$met_min.CpGs <- 1        # minimum number of CpG sites per feature
opts$met_min.cells <- 50      # minimum number of cells per feature 
opts$met_nfeatures <- 1000    # maximum number of features per view (filter based on variance)

# Filtering options for accessibility
opts$acc_min.GpCs <- 5        # minimum number of GpC sites per feature
opts$acc_min.cells <- 50      # minimum number of cells per feature
opts$acc_nfeatures <- 1000    # maximum number of features per view (filter based on variance)

# Filtering options for RNA
opts$rna_min.cdr <- 0.25      # Remove genes with cellular detection rate smaller than opts$min.cdr
opts$rna_ngenes <- 2500       # maximum number of genes (filter based on variance)

# Define colors
opts$colors <- c(
  Epiblast="#63B8FF",
  Nascent_mesoderm="#FF82AB",
  Mature_mesoderm="#CD3278",
  Primitive_Streak="sandybrown",
  Embryonic_endoderm="#43CD80",
  Endoderm="#43CD80",
  Visceral_endoderm="#2E8B57",
  Ectoderm="steelblue"
)


# Define which cells to use for EB data
tmp <- fread(io$sample.metadata.1) %>%
  .[culture%in%opts$culture]# %>%
  # .[celltype.mapped%in%opts$celltype.mapped] %>%
  # .[stage.mapped%in%opts$stage.mapped]
opts$met_cells.1 <- tmp %>% .[pass_metQC==T, id_met]
opts$rna_cells.1 <- tmp %>% .[pass_rnaQC==T, id_rna]
opts$acc_cells.1 <- tmp %>% .[pass_accQC==T, id_acc]


# Define which cells to use for in vivo data
tmp <- fread(io$sample.metadata.2) %>% 
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")]  %>%
  .[stage_lineage%in%opts$stage_lineage] 
opts$met_cells.2 <- tmp %>% .[pass_metQC==T, id_met]
opts$rna_cells.2 <- tmp %>% .[pass_rnaQC==T, id_rna]
opts$acc_cells.2 <- tmp %>% .[pass_accQC==T, id_acc]

##########################
## Load sample metadata ##
##########################

sample_metadata.1 <- fread(io$sample.metadata.1) %>%
  .[id_met%in%opts$met_cells.1 | id_rna %in% opts$rna_cells.1 | id_acc %in% opts$acc_cells.1 ] %>%
  .[,c("sample","id_met","id_acc","id_rna","stage.mapped","celltype.mapped")] %>%
  setnames(c("stage.mapped","celltype.mapped"),c("stage","lineage10x_2")) %>%
  .[,assay:="EB"] %>%
  droplevels()

sample_metadata.2 <- fread(io$sample.metadata.2) %>%
  .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")]  %>%
  .[id_met%in%opts$met_cells.2 | id_rna %in% opts$rna_cells.2 | id_acc %in% opts$acc_cells.2 ] %>%
  .[,c("sample","id_met","id_acc","id_rna","stage","lineage10x_2")] %>%
  .[,assay:="in vivo"] %>%
  droplevels()

sample_metadata <- rbind(sample_metadata.1,sample_metadata.2)
