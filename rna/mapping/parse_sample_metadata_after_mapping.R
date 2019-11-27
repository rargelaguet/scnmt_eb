
############################################################
## Script to add mapping information to the cell metadata ##
############################################################

####################
## Load libraries ##
####################

library(data.table)
library(purrr)

#########
## I/O ##
#########

io$sample_metadata <- "/Users/C02RF23NFVH8/data/scnmt_eb/sample_metadata.txt"
io$mapping <- "/Users/C02RF23NFVH8/data/scnmt_eb/rna/results/mapping/mapping.rds"

###############
## Load data ##
###############

# Load mapping information
mapping <- readRDS(io$mapping)$mapping

# Load cell metadata and remove any previous mapping information
sample.metadata <- fread(sample_metadata.file) %>%
  .[,c("lineage10x","lineage10x_2","stage.mapped","celltype.multinomial.prob"):=NULL]

mapping.dt <- mapping %>% as.data.table %>%
  .[,c("cell","celltype.mapped","stage.mapped","celltype.multinomial.prob")] %>%
  .[,cell:=stringr::str_replace_all(cell,"map_","")] %>%
  setnames("cell","id_rna")

##################################################
## Merge cell metadata with mapping information ##
##################################################

sample.metadata <- merge(sample.metadata, mapping.dt, by="id_rna") %>%
  setnames("celltype.mapped","lineage10x")

################################
## Aggregate similar lineages ##
################################

sample.metadata %>%
  .[,lineage10x_2:=lineage10x] %>%
  # No ectoderm cells are expected at this stage, they are just epiblast (common error with the mapping)
  .[lineage10x=="Caudal neurectoderm",c("lineage10x_2","lineage10x"):="Epiblast"] %>%
  # Mesoderm
  .[lineage10x%in%c("Pharyngeal mesoderm","Paraxial mesoderm","ExE mesoderm","Mesenchyme","Intermediate mesoderm", "Mixed mesoderm", "Nascent mesoderm","Allantois"), lineage10x_2:="Mesoderm"] %>%
  # Blood
  .[lineage10x%in%c("Blood progenitors 2", "Blood progenitors 1", "Endothelium", "Erythroid2", "Erythroid3", "Haematoendothelial progenitors","Erythroid1","Cardiomyocytes"), lineage10x_2:="Blood"] %>%
  # Endoderm
  .[lineage10x%in%c("Gut","Def. endoderm","Notochord","Parietal endoderm","ExE endoderm","Visceral endoderm"), lineage10x_2:="Endoderm"] %>%
  # Primitive streak
  .[lineage10x%in%c("Caudal epiblast","Anterior Primitive Streak"), lineage10x_2:="Primitive Streak"] %>%
  # Ectoderm
  .[lineage10x%in%c("Rostral neurectoderm","Surface ectoderm"), lineage10x_2:="Ectoderm"]
  # PGC are not expected, they are just epiblast cells (common error with the mapping)
  .[lineage10x%in%c("PGC"), lineage10x_2:="Epiblast"] %>%

#################
## Save output ##
#################

fwrite(sample.metadata, file=io$sample_metadata, sep="\t", col.names=T, row.names=F, na="NA", quote=F)


