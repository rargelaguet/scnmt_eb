library(SingleCellExperiment)
library(data.table)
library(purrr)
library(ggplot2)

source("/Users/ricard/NMT-seq_EB_2/rna/mapping/scNMT/Mapping2gastrulationAtlas.R")

####################
## Define options ##
####################

io <- list()
io$path2atlas <- "/Users/ricard/data/gastrulation"
io$path2target <- "/Users/ricard/data/NMT-seq_EB+ESC"

#######################
## Load target atlas ##
#######################

sce_atlas  <- readRDS(paste0(io$path2atlas, "/rna/parsed/SingleCellExperiment.rds"))
meta_atlas <- read.table(file = paste0(io$path2atlas, "/sample_metadata.txt"), header = TRUE, sep = "\t", stringsAsFactors = F)
meta_atlas$dataset <- "atlas"

# Filter
meta_atlas <- meta_atlas[meta_atlas$pass_rnaQC==T,]
sce_atlas <- sce_atlas[,meta_atlas$id_rna] 

#####################
## Load query data ##
####################

sce_target  <- readRDS(paste0(io$path2target, "/rna/SingleCellExperiment.rds"))

meta_target <- read.table(file = paste0(io$path2target, "/sample_metadata.txt"), header = TRUE, sep = "\t", stringsAsFactors = F)
meta_target$dataset <- "target"
# meta_target$stage <- "Day5"
meta_target$stage <- meta_target$culture

# Filter
meta_target <- meta_target[meta_target$pass_rnaQC==T,]
sce_target <- sce_target[,meta_target$id_rna] 

meta_scnmt <- list()
meta_scnmt$cell <- meta_target$id_rna[match(colnames(sce_target), meta_target$id_rna)]
meta_scnmt$cells <- meta_target$id_rna[match(colnames(sce_target), meta_target$id_rna)]
# meta_scnmt$stage <- "Day5"
meta_scnmt$stage <- meta_target$stage

#############
## Prepare ## 
#############

genes <- intersect(rownames(sce_target),rownames(sce_atlas))
sce_target  <- sce_target[genes, ]
sce_atlas <- sce_atlas[genes, ]

#########
## Map ##
#########

mapping  <- mapWrap(
  atlas_sce = sce_atlas, 
  atlas_meta = meta_atlas,
  map_sce = sce_target, 
  map_meta = meta_scnmt, 
  k = 25,
  nPCs = 50
)

saveRDS(mapping, "/Users/ricard/data/NMT-seq_EB+ESC/mapping/mapping_scNMT_2.rds")

################
## START TEST ##
################

# big_sce <- SingleCellExperiment::SingleCellExperiment(
#   list(counts=Matrix::Matrix(cbind(counts(sce_atlas),counts(sce_target)),sparse=TRUE))
#   ) %>% scater::normalize()

# foo <- meta_target[,c("id_rna","stage","dataset")]
# colnames(foo) <- c("id_rna","stage","dataset")
# bar <- meta_atlas[,c("id_rna","stage","dataset")]
# colnames(bar) <- c("id_rna","stage","dataset")
# meta_merged <- rbind(foo,bar)
# colData(big_sce) <- DataFrame(meta_merged)

# big_sce@reducedDims$PCA <- irlba::prcomp_irlba(t(logcounts(big_sce)), n = 2)$x
# plotPCA(big_sce, colour_by="dataset", ncomponents=c(1,2))

##############
## END TEST ##
##############

################################
## Merge with sample metadata ##
################################

foo <- fread(paste0(io$path2target,"/sample_metadata.txt"))# %>%
  # .[,c("celltype.mapped","stage.mapped","celltype.multinomial.prob"):=NULL]

bar <- mapping$mapping %>%
  .[,c("cell","celltype.mapped","stage.mapped")] %>%
  as.data.table %>%
  setnames("cell","id_rna")

sample_metadata_updated <- merge(foo,bar, by="id_rna", all.x=T)

# fwrite(sample_metadata_updated, paste0(io$path2target,"/sample_metadata2.txt"), sep="\t", col.names=T, row.names=F, quote=F, na="NA")


sample_metadata_updated[,.N,by=c("culture","stage.mapped.y")] %>% View
