library(data.table)
library(purrr)

sample_metadata <- fread("/Users/ricard/data/NMT-seq_EB+ESC/sample_metadata2.txt") %>%
  setnames("celltype.mapped","lineage10x")

# sample_metadata[culture=="EB_Day5",.N,by=c("lineage10x_2","culture")]


sample_metadata %>%
  .[,lineage10x_2:=lineage10x] %>%
  # Epiblast
  .[culture%in%c("EB_Day2","ESC_Serum","ESC_2i") & lineage10x=="Caudal neurectoderm",c("lineage10x_2","lineage10x"):="Epiblast"] %>%
  # Mesoderm
  .[lineage10x%in%c("Pharyngeal mesoderm","Paraxial mesoderm","ExE mesoderm","Mesenchyme","Intermediate mesoderm", "Mixed mesoderm", "Nascent mesoderm","Allantois"), lineage10x_2:="Mesoderm"] %>%
  .[lineage10x%in%c("Blood progenitors 2", "Blood progenitors 1", "Endothelium", "Erythroid2", "Erythroid3", "Haematoendothelial progenitors","Erythroid1","Cardiomyocytes"), lineage10x_2:="Blood"] %>%
  # Endoderm
  .[lineage10x%in%c("Gut","Def. endoderm","Notochord","Parietal endoderm","ExE endoderm","Visceral endoderm"), lineage10x_2:="Endoderm"] %>%
  # Primitive streak
  .[lineage10x%in%c("Caudal epiblast","Anterior Primitive Streak"), lineage10x_2:="Primitive Streak"] %>%
  # Ectoderm
  .[lineage10x%in%c("Rostral neurectoderm","Surface ectoderm"), lineage10x_2:="Ectoderm"]
  # # PGC are not expected
  # .[lineage10x%in%c("PGC"), lineage10x_2:=NA] %>%
  # .[stage=="E6.5" & lineage10x%in%c("Rostral_neurectoderm","Surface_ectoderm"), lineage10x_2:="Epiblast"]

unique(sample_metadata$lineage10x_2)
## START TEST ##
# sample_metadata %>%
#   .[lineage10x_2%in%c("Mature_mesoderm","Nascent_mesoderm"), lineage10x_2:="Mesoderm"] %>%
#   .[stage=="E7.5" & lineage10x_2%in%c("Embryonic_endoderm","Visceral_endoderm"), lineage10x_2:="Endoderm"]  %>%
#   .[stage=="E5.5" & lineage10x_2%in%c("Endoderm"), lineage10x_2:="Visceral_endoderm"] %>%
#   .[stage=="E6.5" & lineage10x_2%in%c("Endoderm"), lineage10x_2:="Visceral_endoderm"] 
## END TEST ##

# table(sample_metadata[pass_metQC==T & stage=="E7.5",lineage10x_2])

fwrite(sample_metadata, "/Users/ricard/data/NMT-seq_EB+ESC/sample_metadata2.txt", sep="\t", col.names=T, row.names=F, na="NA", quote=F)


