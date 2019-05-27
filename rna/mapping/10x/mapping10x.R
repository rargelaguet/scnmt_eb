library(SingleCellExperiment)
library(data.table)
library(purrr)
library(ggplot2)

source("/Users/ricard/NMT-seq_EB+ESC/rna/mapping/10x/Mapping2gastrulationAtlas.R")

path2atlas <- "/Users/ricard/data/gastrulation10x"
path2scNMT <- "/Users/ricard/data/NMT-seq_EB+ESC"

####################
## Load 10x atlas ##
####################

sce_atlas  <- readRDS(paste0(path2atlas, "/processed/SingleCellExperiment.rds"))
meta_atlas <- readRDS(paste0(path2atlas, "/processed/sample_metadata.rds"))
# meta_atlas <- read.table(paste0(path2atlas,"/sample_metadata.txt"), header=T, sep="\t", stringsAsFactors=F)

####################
## Load scNMT-seq ##
####################

sce_nmt  <- readRDS(paste0(path2scNMT, "/rna/SingleCellExperiment.rds"))
meta_nmt <- read.table(file = paste0(path2scNMT, "/sample_metadata2.txt"), header = TRUE, sep = "\t", stringsAsFactors = F)
meta_nmt$stage <- meta_nmt$population

# Filter
meta_nmt <- meta_nmt[meta_nmt$pass_rnaQC==T,]
# meta_nmt <- meta_nmt[meta_nmt$culture %in% c("EB_Day5","EB_Day6","EB_Day7") & meta_nmt$pass_rnaQC==T,]
sce_nmt <- sce_nmt[,meta_nmt$id_rna] 

meta_scnmt <- list()
meta_scnmt$cell <- meta_nmt$id_rna[match(colnames(sce_nmt), meta_nmt$id_rna)]
meta_scnmt$cells <- meta_nmt$id_rna[match(colnames(sce_nmt), meta_nmt$id_rna)]
meta_scnmt$stage <- meta_nmt$culture[match(colnames(sce_nmt), meta_nmt$id_rna)]

#############
## Prepare ## 
#############

genes <- intersect(rownames(sce_nmt), rownames(sce_atlas))
sce_nmt  <- sce_nmt[genes,]
sce_atlas <- sce_atlas[genes,]

#########
## Map ##
#########

mapping  <- mapWrap(
  atlas_sce = sce_atlas, 
  atlas_meta = meta_atlas,
  map_sce = sce_nmt, 
  map_meta = meta_scnmt, 
  k = 25
)


################################
## Merge with sample metadata ##
################################

sample_metadata.file <- "/Users/ricard/data/NMT-seq_EB+ESC/sample_metadata2.txt"

## START TEST ##
# foo <- fread(sample_metadata.file) %>% .[,c("id_rna","lineage10x","stage.mapped")]
# bar <- mapping$mapping %>% as.data.table %>%
#   .[,c("cell","celltype.mapped","stage.mapped","celltype.multinomial.prob")] %>%
#   .[,cell:=stringr::str_replace_all(cell,"map_","")] %>%
#   setnames("cell","id_rna")
# sample_metadata_updated <- merge(foo,bar, by="id_rna")
## END TEST ##

foo <- fread(sample_metadata.file) %>%
  .[,c("lineage10x","lineage10x_2","stage.mapped","celltype.multinomial.prob"):=NULL]

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

saveRDS(mapping, "/Users/ricard/data/NMT-seq_EB+ESC/mapping/mapping_v3.rds")

fwrite(sample_metadata_updated, sample_metadata.file, sep="\t", col.names=T, row.names=F, quote=F, na="NA")

