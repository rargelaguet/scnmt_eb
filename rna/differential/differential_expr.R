
#########################################################
## Script to do differential RNA expression with edgeR ##
#########################################################

####################
## Load libraries ##
####################

suppressMessages(library(scater))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(ggplot2))
suppressMessages(library(edgeR))
suppressMessages(library(argparse))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('-s1', '--lineage_genotype1', type="character",  nargs='+',  help='lineage_genotype A (i.e. Epiblast_WT)')
p$add_argument('-s2', '--lineage_genotype2', type="character",  nargs='+',  help='lineage_genotype B (i.e. Mesoderm_KO)')
p$add_argument('-o',  '--outfile',           type="character",              help='Output file')
args <- p$parse_args(commandArgs(TRUE))

#########
## I/O ##
#########

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/scnmt_eb"
  io$gene.metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  source("/Users/ricard/gastrulation/rna/differential/utils.R")
} else {
  io$basedir <- "/hps/nobackup/stegle/users/ricard/scnmt_eb"
  io$gene.metadata <- "/hps/nobackup/stegle/users/ricard/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
  source("/homes/ricard/gastrulation/rna/differential/utils.R")
}
io$sample_metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$rna.infile <- paste(io$basedir,"rna/SingleCellExperiment.rds",sep="/")
io$outfile <- args$outfile

#############
## Options ##
#############

opts <- list()

# Define stage and lineage
opts$groupA <- args$lineage_genotype1
opts$groupB <- args$lineage_genotype2

# Define FDR threshold
opts$threshold_fdr <- 0.01

# Define minimum logFC for significance
opts$min.logFC <- 1.0

# Define which cells to use
opts$cells <- fread(io$sample_metadata) %>% 
  .[,lineage_genotype:=paste(lineage10x_2,genotype,sep="_")] %>%s
  .[pass_rnaQC==T & lineage_genotype%in%c(opts$groupA,opts$groupB),id_rna]

###############
## Load data ##
###############

# Load sample metadata
sample_metadata <- fread(io$sample_metadata) %>% 
  .[id_rna %in% opts$cells] %>% 
  .[,lineage:=lineage10x_2] %>%
  .[,lineage_genotype:=paste(lineage,genotype,sep="_")] 

# Load SingleCellExperiment object
sce <- readRDS(io$rna.infile)[,opts$cells]
sce$lineage <- sample_metadata$lineage
sce$lineage_genotype <- sample_metadata$lineage_genotype

# Define the two exclusive groups
sample_metadata[,group:=as.factor(as.numeric(lineage_genotype%in%opts$groupB))]
sce$group <- as.factor(as.numeric(sce$lineage_genotype%in%opts$groupB))

# Load gene metadata
gene_metadata <- rowData(sce) %>% as.data.frame(row.names=rownames(sce)) %>% 
  tibble::rownames_to_column("ens_id") %>% .[,c("symbol","ens_id")] %>% 
  as.data.table %>% setnames("ens_id","id")
gene_metadata[,c("symbol","id"):=list(as.factor(symbol),as.factor(id))]


################################################
## Differential expression testing with edgeR ##
################################################

out <- doDiffExpr(sce, sample_metadata) %>%
  merge(gene_metadata, all.y=T, by="id") %>% setorderv("padj_fdr", na.last=T)

##################
## Save results ##
##################

fwrite(out, file=io$outfile, sep="\t", na="NA", quote=F)
