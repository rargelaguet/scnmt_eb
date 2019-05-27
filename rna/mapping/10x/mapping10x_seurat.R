library(SingleCellExperiment)
library(data.table)
library(purrr)
library(ggplot2)
library(Seurat)

source("/Users/ricard/NMT-seq_EB+ESC/rna/mapping/10x/Mapping2gastrulationAtlas.R")

path2atlas <- "/Users/ricard/data/gastrulation10x/processed"
path2scNMT <- "/Users/ricard/data/NMT-seq_EB+ESC"

####################
## Load 10x atlas ##
####################

sce_atlas  <- readRDS(paste0(path2atlas, "/SingleCellExperimentAtlas.rds"))
meta_atlas <- readRDS(paste0(path2atlas, "/sample_metadata_atlas.rds"))

####################
## Load scNMT-seq ##
####################

sce_query  <- readRDS(paste0(path2scNMT, "/rna/SingleCellExperiment.rds"))
meta_query <- read.table(file = paste0(path2scNMT, "/sample_metadata.txt"), header = TRUE, sep = "\t", stringsAsFactors = F)
meta_query$stage <- meta_query$population

# Filter
meta_query <- meta_query[meta_query$pass_rnaQC==T,]
sce_query <- sce_query[,meta_query$id_rna] 

meta_scnmt <- list()
meta_scnmt$cell <- meta_query$id_rna[match(colnames(sce_query), meta_query$id_rna)]
meta_scnmt$cells <- meta_query$id_rna[match(colnames(sce_query), meta_query$id_rna)]
meta_scnmt$stage <- meta_query$culture[match(colnames(sce_query), meta_query$id_rna)]

#############
## Prepare ## 
#############

genes <- intersect(rownames(sce_query), rownames(sce_atlas))
sce_query  <- sce_query[genes,]
sce_atlas <- sce_atlas[genes,]

###########################
## Create Seurat objects ##
###########################

seurat_atlas   <- CreateSeuratObject(as.matrix(counts(sce_atlas)), project = "ATLAS", min.cells = 5)
seurat_atlas@meta.data$map <- "ATLAS"
seurat_atlas   <- NormalizeData(seurat_atlas)
rownames(meta_atlas) <- meta_atlas$cell
seurat_atlas   <- AddMetaData(object = seurat_atlas, metadata = meta_atlas)
# seurat_atlas   <- ScaleData(seurat_atlas, display.progress = TRUE, vars.to.regress = "sample")
seurat_atlas   <- ScaleData(seurat_atlas)

seurat_query <- CreateSeuratObject(as.matrix(counts(sce_query)), project = "QUERY", min.cells = 5)
seurat_query@meta.data$map <- "QUERY"
seurat_query <- NormalizeData(seurat_query)
seurat_query <- ScaleData(seurat_query)

#########
## Map ##
#########

# Gene selection for input to CCA
seurat_atlas <- FindVariableFeatures(seurat_atlas,selection.method = "vst", do.plot = F)
seurat_query <- FindVariableFeatures(seurat_query,selection.method = "vst", do.plot = F)

# Mapping
anchors <- FindTransferAnchors(reference = seurat_atlas, query = seurat_query, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = seurat_atlas$celltype, dims = 1:30)
seurat_query <- AddMetaData(object = seurat_query, metadata = predictions)

table(seurat_query$predicted.id)


################################
## Merge with sample metadata ##
################################


sample_metadata.file <- "/Users/ricard/data/NMT-seq_EB+ESC/sample_metadata.txt"
foo <- fread(sample_metadata.file)# %>%
  # .[,c("lineage10x","stage.mapped","celltype.multinomial.prob"):=NULL]

bar <- mapping$mapping %>% 
  .[,c("cell","celltype.mapped","stage.mapped","celltype.multinomial.prob")] %>%
  as.data.table %>%
  .[,cell:=stringr::str_replace_all(cell,"map_","")] %>%
  setnames("cell","id_rna")

sample_metadata_updated <- merge(foo,bar, by="id_rna")

# foobar[culture=="EB_Day5",.N,by=c("phenotype","celltype.mapped.y")] %>% View

##########
## Save ##
##########

saveRDS(mapping, "/Users/ricard/data/NMT-seq_EB+ESC/mapping/mapping.rds")

fwrite(sample_metadata_updated, sample_metadata.file, sep="\t", col.names=T, row.names=F, quote=F, na="NA")

