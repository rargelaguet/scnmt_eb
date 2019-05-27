suppressMessages(library(scater))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(ggplot2))
suppressMessages(library(edgeR))
suppressMessages(library(argparse))

## Initialize argument parser ##
p <- ArgumentParser(description='')
p$add_argument('-s1', '--lineage1', type="character",  nargs='+',  help='lineage 1 (E4.5_EPI, E5.5_VE,...)')
p$add_argument('-s2', '--lineage2', type="character",  nargs='+',  help='lineage 2 (E4.5_EPI, E5.5_VE,...)')
p$add_argument('-o',  '--outfile',        type="character",              help='Output file')
args <- p$parse_args(commandArgs(TRUE))


## I/O ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/NMT-seq_EB+ESC"
  io$gene.metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  source("/Users/ricard/gastrulation/rna/differential/utils.R")
} else {
  io$basedir <- "/hps/nobackup/stegle/users/ricard/NMT-seq_EB+ESC"
  io$gene.metadata <- "/hps/nobackup/stegle/users/ricard/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  source("/homes/ricard/gastrulation/rna/differential/utils.R")
}
io$sample_metadata <- paste0(io$basedir,"/sample_metadata2.txt")
io$rna.infile <- paste(io$basedir,"rna/SingleCellExperiment.rds",sep="/")
io$outfile <- args$outfile

## Options ##
opts <- list()

# Define stage and lineage
opts$groupA <- args$lineage1
opts$groupB <- args$lineage2

opts$groupA <- c("Epiblast_Tet_WT")
opts$groupB <- c("Epiblast_Tet_KO")

# Define FDR threshold
opts$threshold_fdr <- 0.01

# Define minimum logFC for significance
opts$min.logFC <- 1.0

# Define which cells to use
opts$cells <- fread(io$sample_metadata) %>% 
  .[,lineage:=lineage10x_2] %>% 
  .[,lineage_phenotype:=paste(lineage,phenotype,sep="_")] %>%
  .[culture%in%c("EB_Day2")] %>%
  # .[pass_rnaQC==T & lineage%in%c(opts$groupA,opts$groupB),id_rna]
  .[pass_rnaQC==T & lineage_phenotype%in%c(opts$groupA,opts$groupB),id_rna]

###############
## Load data ##
###############

# Load sample metadata
sample_metadata <- fread(io$sample_metadata) %>% 
  .[id_rna %in% opts$cells] %>% 
  .[,lineage:=lineage10x_2] %>%
  .[,lineage_phenotype:=paste(lineage,phenotype,sep="_")] 

# Load SingleCellExperiment object
sce <- readRDS(io$rna.infile)[,opts$cells]
sce$lineage <- sample_metadata$lineage
sce$lineage_phenotype <- sample_metadata$lineage_phenotype

# Define the two exclusive groups
# sample_metadata[,group:=as.factor(as.numeric(lineage%in%opts$groupB))]
# sce$group <- as.factor(as.numeric(sce$lineage%in%opts$groupB))
sample_metadata[,group:=as.factor(as.numeric(lineage_phenotype%in%opts$groupB))]
sce$group <- as.factor(as.numeric(sce$lineage_phenotype%in%opts$groupB))

# Load gene metadata
gene_metadata <- rowData(sce) %>% as.data.frame(row.names=rownames(sce)) %>% 
  tibble::rownames_to_column("ens_id") %>% .[,c("symbol","ens_id")] %>% 
  as.data.table %>% setnames("ens_id","id")
gene_metadata[,c("symbol","id"):=list(as.factor(symbol),as.factor(id))]


################################################
## Differential expression testing with edgeR ##
################################################

out <- doDiffExpr(sce, sample_metadata) %>%
  merge(gene_metadata, all.y=T, by=c("id")) %>% setorderv("padj_fdr", na.last=T)

##################
## Save results ##
##################

fwrite(out, file=io$outfile, sep="\t", na="NA", quote=F)
