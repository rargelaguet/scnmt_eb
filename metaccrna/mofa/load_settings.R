library(data.table)
library(purrr)

## Define I/O ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/NMT-seq_EB+ESC"
  io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  io$outdir <- "/Users/ricard/data/NMT-seq_EB+ESC/mofa"
} else {
  io$basedir <- "/hps/nobackup/stegle/users/ricard/NMT-seq_EB+ESC"
  io$gene_metadata <- "/hps/nobackup/stegle/users/ricard/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  # io$outdir <- "/homes/ricard/NMT-seq_EB+ESC/metaccrna/mofa"
}
io$annos_dir  <- paste0(io$basedir, "/features/filt")
io$sample.metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$met.dir <- paste0(io$basedir,"/met/feature_level")
io$acc.dir <- paste0(io$basedir,"/acc/feature_level")
io$rna.file <- paste0(io$basedir,"/rna/SingleCellExperiment.rds")

io$met.stats <- paste0(io$basedir,"/met/stats/samples/sample_stats.txt")
io$acc.stats <- paste0(io$basedir,"/acc/stats/samples/sample_stats.txt")


## Define options ##
opts <- list()

# Define which annotations to look at
opts$met.annos <- c(
  # "genebody",
  "prom_2000_2000",
  "H3K27ac_distal_E7.5_Mes_intersect12_500",
  "H3K27ac_distal_E7.5_Ect_intersect12_500",
  "H3K27ac_distal_E7.5_End_intersect12_500"
  # "H3K4me3_E7.5_Mes",
  # "H3K4me3_E7.5_End",
  # "H3K4me3_E7.5_Ect",
)

opts$acc.annos <- c(
  # "genebody",
  "prom_2000_2000",
  "H3K27ac_distal_E7.5_Mes_intersect12",
  "H3K27ac_distal_E7.5_Ect_intersect12",
  "H3K27ac_distal_E7.5_End_intersect12"
  # "H3K4me3_E7.5_Mes",
  # "H3K4me3_E7.5_End",
  # "H3K4me3_E7.5_Ect",
)


# Define which cells to use
opts$culture <- c(
  # "EB_Day5",
  # "EB_Day3"
  "EB_Day6",
  "EB_Day7"
)
opts$celltype.mapped <- c(
  "Epiblast",
  "Primitive_Streak",
  "Mesoderm",
  "Endoderm",
  "Ectoderm"
)

# opts$stage.mapped <- c(
#   "E6.5",
#   "E7.5"
# )

# Filtering options for methylation
opts$met_min.CpGs <- 1        # minimum number of CpG sites per feature
opts$met_min.cells <- 10      # minimum number of cells per feature
opts$met_nfeatures <- 1000    # maximum number of features per view (filter based on variance)

# Filtering options for accessibility
opts$acc_min.GpCs <- 5        # minimum number of GpC sites per feature
opts$acc_min.cells <- 10      # minimum number of cells per feature
opts$acc_nfeatures <- 1000    # maximum number of features per view (filter based on variance)

# Filtering options for RNA
opts$rna_min.cdr <- 0.25      # Remove genes with cellular detection rate smaller than opts$min.cdr
opts$rna_ngenes <- 2500       # maximum number of genes (filter based on variance)

# Define colors
opts$colors <- c(
  Epiblast="#63B8FF",
  Mesoderm="#CD3278",
  Primitive_Streak="sandybrown",
  Endoderm="#43CD80",
  Ectoderm="steelblue"
)


# window length for the overlap between genes and features
opts$overlapGenes  <- FALSE
opts$gene_window  <- 5e4

# Define which cells to use
tmp <- fread(io$sample.metadata) %>%
  .[culture%in%opts$culture] %>%
  .[lineage10x_2%in%opts$celltype.mapped]
  # .[stage.mapped%in%opts$stage.mapped]
opts$met_cells <- tmp %>% .[pass_metQC==T, id_met]
opts$rna_cells <- tmp %>% .[pass_rnaQC==T, id_rna]
opts$acc_cells <- tmp %>% .[pass_accQC==T, id_acc]

sample_metadata <- fread(io$sample.metadata,stringsAsFactors=T) %>%
  # .[,stage_lineage:=as.factor(paste(stage,lineage10x_2,sep="_"))] %>%
  .[id_met%in%opts$met_cells | id_rna %in% opts$rna_cells | id_acc %in% opts$acc_cells ] %>%
  droplevels()
