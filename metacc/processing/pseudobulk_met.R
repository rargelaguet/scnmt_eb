library(furrr)

######################
## Define functions ##
######################

merge_and_sum <- function(dt1, dt2){
  merge(dt1, dt2, by=c("chr","pos"), all = TRUE) %>%
    .[is.na(met_reads.x), met_reads.x := 0L] %>%
    .[is.na(met_reads.y), met_reads.y := 0L] %>%
    .[is.na(nonmet_reads.x), nonmet_reads.x := 0L] %>%
    .[is.na(nonmet_reads.y), nonmet_reads.y := 0L] %>%
    .[,.(chr=chr, pos=pos, met_reads=met_reads.x+met_reads.y, nonmet_reads=nonmet_reads.x+nonmet_reads.y)]
}

fread_and_merge <- function(dt, file){
  fread(file, colClasses=list(factor=1L)) %>% 
    setnames(c("chr","pos","met_reads","nonmet_reads","rate")) %>%
    .[,rate:=NULL] %>%
    merge_and_sum(dt)
}

################
## Define I/O ##
################

source("/homes/ricard/scnmt_eb/settings.R")
io$outdir <- paste0(io$basedir,"/met/cpg_level/pseudobulk")

####################
## Define options ##
####################

# opts$day_lineage <- c(
#   # "Day2_Epiblast",
#   # "Day2_Primitive Streak",
#   # "Day4-5_Epiblast",
#   # "Day4-5_Primitive Streak",
#   "Day4-5_Mesoderm",
#   "Day6-7_Mesoderm",
#   "Day6-7_Endoderm",
#   "Day6-7_Blood"
# )

# opts$genotype <- c(
#   "WT",
#   "KO"
# )

opts$days <- c("Day4-5","Day2")

# opts$group <- c(
#   # "Day2_Epiblast",
#   # "Day2_Primitive Streak",
#   # "Day4-5_Epiblast",
#   # "Day4-5_Primitive Streak",
#   # "WT_Epiblast",
#   # "KO_Epiblast",
#   # "WT_Mesoderm",
#   # "KO_Mesoderm",
#   # "WT_Endoderm",
#   # "KO_Endoderm"
#   "WT_Blood"
# )
opts$group <- c("KO","WT")

# Parallel processing
opts$parallel <- TRUE    # do parallel processing?
opts$ncores <- 2         # number of cores
opts$chunk_size <- 10    # chunk_size: the higher the less memory it is required????

# Update sample metadata 
sample_metadata <- sample_metadata %>% 
  .[pass_metQC==T & day2%in%opts$days] %>%
  # .[,group:=paste(genotype,lineage10x_2,sep="_")]
  .[,group:=genotype]
table(sample_metadata$day2)
table(sample_metadata$group)

##############################
## Load data and pseudobulk ##
##############################

# Parallel processing options
if (opts$parallel){
  plan(multiprocess, workers=opts$ncores)
} else {
  plan(sequential)
}

for (i in opts$group) {
  print(i)
  
  # Define input files 
  cells <- sample_metadata[group%in%i,id_met]
  
  # cells <- head(cells,n=2)
  files <- paste0(io$met_data_raw, "/", cells, ".tsv.gz")
  
  # split into chunks for parallel processing
  if (opts$parallel) {
    chunks <- ceiling(seq_along(files)/opts$chunk_size)
    file_list <- split(files, chunks)
  } else {
    file_list <- list(files)
  }
  
  # pseudobulk
  init <- data.table(chr=as.factor(NA), pos=as.integer(NA), met_reads=as.integer(NA), nonmet_reads=as.integer(NA))
  data <- future_map(file_list, purrr::reduce, fread_and_merge, .init=init, .progress=F) %>%
  # data <- map(file_list, purrr::reduce, fread_and_merge, .init=init) %>%
    purrr::reduce(merge_and_sum) %>%
    .[,rate:=round(100*met_reads/(met_reads+nonmet_reads))] %>%
    .[,.(chr,pos,met_reads,nonmet_reads,rate)]
  
  # filter
  data <- data[!is.na(rate)]
  
  # Save
  outfile = sprintf("%s/%s.tsv.gz",io$outdir,i)
  fwrite(data, outfile, quote=F, col.names=T, sep="\t")
}
