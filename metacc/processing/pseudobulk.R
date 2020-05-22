library(data.table)
library(purrr)
library(furrr)

######################
## Define functions ##
######################

merge_and_sum <- function(dt1, dt2){
  merge(dt1, dt2, by=c("chr","pos"), all = TRUE) %>%
    .[is.na(met_cpgs.x), met_cpgs.x := 0L] %>%
    .[is.na(met_cpgs.y), met_cpgs.y := 0L] %>%
    .[is.na(nonmet_cpgs.x), nonmet_cpgs.x := 0L] %>%
    .[is.na(nonmet_cpgs.y), nonmet_cpgs.y := 0L] %>%
    .[,.(chr=chr, pos=pos, met_cpgs=met_cpgs.x+met_cpgs.y, nonmet_cpgs=nonmet_cpgs.x+nonmet_cpgs.y)]
}

fread_and_merge <- function(dt, file){
  fread(file, colClasses=list(factor=1L)) %>% 
    setnames(c("chr","pos","met_cpgs","nonmet_cpgs","rate")) %>%
    .[,rate:=NULL] %>%
    merge_and_sum(dt)
}

################
## Define I/O ##
################

source("/Users/ricard/scnmt_eb/settings.R")
io$outdir <- paste0(io$basedir,"/met/cpg_level/pseudobulk")


####################
## Define options ##
####################

opts$day_lineage <- c(
  "Day2_Epiblast",
  "Day2_Primitive Streak",
  "Day4-5_Epiblast",
  "Day4-5_Primitive Streak",
  "Day4-5_Mesoderm",
  "Day6-7_Mesoderm",
  "Day6-7_Endoderm",
  "Day6-7_Blood"
)

opts$genotype <- c(
  "WT",
  "KO"
)

# Parallel processing
opts$parallel <- TRUE    # do parallel processing?
opts$ncores <- 2         # number of cores
opts$chunk_size <- 10    # chunk_size: the higher the less memory it is required????

# Update sample metadata 
sample_metadata <- sample_metadata %>% .[pass_metQC==T & day_lineage%in%opts$day_lineage]
  

##############################
## Load data and pseudobulk ##
##############################

# Parallel processing options
if (opts$parallel){
  plan(multiprocess, workers=opts$ncores)
} else {
  plan(sequential)
}

for (i in opts$day_lineage) {
  print(i)
  
  # Define input files 
  cells <- sample_metadata[day_lineage%in%i,id_met]
  
  # cells <- head(cells,n=2)
  files <- paste0(io$data, "/", cells, ".tsv.gz")
  
  # split into chunks for parallel processing
  if (opts$parallel) {
    chunks <- ceiling(seq_along(files)/opts$chunk_size)
    file_list <- split(files, chunks)
  } else {
    file_list <- list(files)
  }
  
  # pseudobulk
  init <- data.table(chr=as.factor(NA), pos=as.integer(NA), met_cpgs=as.integer(NA), nonmet_cpgs=as.integer(NA))
  data <- future_map(file_list, purrr::reduce, fread_and_merge, .init=init, .progress=F) %>%
  # data <- map(file_list, purrr::reduce, fread_and_merge, .init=init) %>%
    purrr::reduce(merge_and_sum) %>%
    .[,rate:=round(100*met_cpgs/(met_cpgs+nonmet_cpgs))] %>%
    .[,.(chr,pos,met_cpgs,nonmet_cpgs,rate)]
  
  # filter
  data <- data[complete.cases(.)]
  
  # Save
  outfile = sprintf("%s/%s.tsv.gz",io$outdir,i)
  fwrite(data, outfile, quote=F, col.names=T, sep="\t")
}
