
## Define I/O ##
io <- list()
io$basedir <- "/Users/C02RF23NFVH8/data/scnmt_eb"
io$in.metadata <- paste0(io$basedir,"/sample_metadata.txt")
io$in.data <- paste0(io$basedir,"/acc/feature_level")
io$in.features  <- "/Users/C02RF23NFVH8/data/scnmt_eb/features/genomic_contexts"
io$outdir <- "/Users/C02RF23NFVH8/data/scnmt_eb/acc/stats/features"

## Define options ##
opts <- list()

# Define which annotations to look at
opts$annos <- c(
  "CGI" = "CpG islands",
  "prom_2000_2000_cgi" = "CGI promoters",
  "prom_2000_2000_noncgi" = "non-CGI promoters",
  "genebody" = "Gene body",
  "LINE" = "LINE",
  "LTR" = "LTR"
)


# Define which cells to use
opts$day_lineage <- c(
  # Day 2
  "Day2_Epiblast",
  "Day2_Primitive Streak",
  
  # Day 4/5
  "Day4_Epiblast",
  "Day4_Primitive Streak",
  "Day4_Mesoderm",
  "Day4_Endoderm",
  "Day5_Epiblast",
  "Day5_Primitive Streak",
  "Day5_Mesoderm",
  "Day5_Endoderm",
  
  # Day 6/7
  "Day6_Endoderm",
  "Day6_Mesoderm",
  "Day6_Blood",
  "Day7_Endoderm",
  "Day7_Mesoderm",
  "Day7_Blood"
)

opts$genotype <- c(
  "WT"
  # "KO"
)

opts$cells <- fread(io$in.metadata) %>%
  .[,day_lineage:=paste(day,lineage10x_2, sep="_")] %>%
  .[day_lineage%in%opts$day_lineage] %>%
  .[genotype%in%opts$genotype] %>%
  .[pass_accQC==T,id_acc]