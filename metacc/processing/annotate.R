##  Script to overlap bismark output files with genomic features ##
###################################################################

options(warn=-1)
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))
suppressMessages(library(argparse))

# Initialize argument parser
p <- ArgumentParser(description='')
p$add_argument('-c','--context', type="character",help='cg/CG or gc/GC')
p$add_argument('-n','--cores', type="integer",help='Number of cores')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define options ##
####################

# args <- list()
# args$context <- "GC"
# args$cores <- 1

io <- list()
opts <- args

## I/O ##

# Define what context to look at: CG (MET) or GC (ACC)
opts$context <- toupper(opts$context)
stopifnot(opts$context %in% c("CG","GC"))

## Own computer ##
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/scnmt_eb"
  
  # GC
  if (opts$context == "GC") {
    io$in.folder <- paste0(io$basedir,"/acc/gpc_level")
    io$out.folder <- paste0(io$basedir,"/acc/feature_level")
  # CG
  } else if (opts$context=="CG") {
    io$in.folder <- paste0(io$basedir,"/met/cpg_level")
    io$out.folder <- paste0(io$basedir,"/met/feature_level")
  }
  
## Cluster ##
} else {
  io$basedir <- "/hps/nobackup/stegle/users/ricard/scnmt_eb"

  # GC
  if (opts$context == "GC") {
    io$in.folder <- paste0(io$basedir,"/acc/gpc_level")
    io$out.folder <- paste0(io$basedir,"/acc/feature_level")
  # CG
  } else if (opts$context=="CG") {
    io$in.folder <- paste0(io$basedir,"/met/cpg_level")
    io$out.folder <- paste0(io$basedir,"/met/feature_level")
  }
}

io$anno.folder <- paste0(io$basedir,"/features/genomic_contexts")
io$in.sample_metadata <- paste0(io$basedir,"/sample_metadata.txt")

## Options ##

# Annotations to analyse
# opts$annos <- "all"
opts$annos <- c(
  # "CGI",
  # "H3K27ac_distal_E7.5_Ect_intersect12",
  # "H3K27ac_distal_E7.5_End_intersect12",
  # "H3K27ac_distal_E7.5_Mes_intersect12",
  # "H3K27ac_distal_E7.5_Ect_intersect12_500",
  # "H3K27ac_distal_E7.5_End_intersect12_500",
  # "H3K27ac_distal_E7.5_Mes_intersect12_500",
  # "H3K4me3_E7.5_Ect",
  # "H3K4me3_E7.5_End",
  # "H3K4me3_E7.5_Mes",
  # "genebody",
  # "prom_2000_2000"
  # "prom_2000_2000_cgi",
  # "prom_2000_2000_noncgi",
  # "LINE",
  # "LTR"
  # "window2000_step1000"
)



if (opts$annos == "all")
  opts$annos <- sapply(str_split(list.files(io$anno.folder, pattern = "\\.bed$"),"\\.bed"),"[[", 1)

cat("\nProcessing methylation samples with the following options:\n")
cat(sprintf("- Input folder for annotation: %s\n",io$anno.folder))
cat(sprintf("- Input folder for bismark files: %s\n",io$in.folder))
cat(sprintf("- Output folder: %s\n",io$out.folder))
cat(sprintf("- Annotations: %s\n", paste(opts$annos, collapse=" ")))
cat(sprintf("- Number of cores: %d\n",opts$cores))
cat(sprintf("- Annotating CG or GC?: %s\n",opts$context))
cat("\n")

###############
## Load data ##
###############

# Load samples to be kept
if (opts$context=="CG") {
  samples_keep <- fread(io$in.sample_metadata, header=T) %>% .[!is.na(id_met),id_met]
} else if (opts$context=="GC") {
  samples_keep <- fread(io$in.sample_metadata, header=T) %>% .[!is.na(id_acc),id_acc]
} else{
  stop()
}
stopifnot(all(!duplicated(samples_keep)))


############################
## Preprocess annotations ##
############################

anno_list <- list()
for (i in 1:length(opts$annos)) {
  
  # Read annotation file
  anno.file <- sprintf("%s/%s.bed",io$anno.folder,opts$annos[i])
  dat_anno <- fread(anno.file ,sep="\t", header=F, select=c(1,2,3,4,5), verbose=F) %>% 
    setnames(c("chr","start","end","strand","id"))
  
  # Check that there are no weird chromosomes
  anno_list[[i]] <- dat_anno %>% setkey(chr,start,end)
}
names(anno_list) <- opts$annos



#########################################
## Preprocess and annotate the samples ##
#########################################

# Create ouput temporary folder
dir.create(sprintf("%s/tmp",io$out.folder), recursive=T)

for (i in 1:length(samples_keep)) {
  sample=samples_keep[i]
  files_processed <- list.files(sprintf("%s/tmp",io$out.folder))
  if (all(sprintf("%s_%s.gz",sample,opts$annos) %in% files_processed)) {
    cat(sprintf("Sample %s already processed for all required annotations...\n",sample)) 
  } else {
    cat(sprintf("Sample %s has not been processed, annotating...\n",sample))  
    
    # Read and parse raw methylation data
    filename <- sprintf("%s/%s.tsv.gz",io$in.folder,sample)
    print(filename)
    if (file.exists(filename)) {
      dat_sample <- fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$in.folder,sample), sep="\t", verbose=F, showProgress=F) 

      if (ncol(dat_sample)==3) {
        dat_sample <- dat_sample %>%
          setnames(c("chr","pos","rate"))
      } else {
        dat_sample <- dat_sample %>%
          .[,c("chr","pos","rate")] %>%
          .[,rate:=rate*100]
      }
      dat_sample <- dat_sample %>%
        .[,c("start","end") := list(pos,pos)] %>%
        .[,c("chr","pos"):=list(as.factor(chr),NULL)] %>% 
        setkey(chr,start,end)

      stopifnot(all(dat_sample$rate %in% c(0,100)))
      
      # Overlap data with annotations
      for (anno in opts$annos) {
        # fname.out <- sprintf("%s/tmp/%s_%s.tsv",io$out.folder,sample,anno)
        fname.out <- sprintf("%s/tmp/%s_%s",io$out.folder,sample,anno)
        if (file.exists(paste0(fname.out,".gz"))) {
          cat(sprintf("Annotation for %s with %s already found, loading...\n",sample,anno))
        } else {
          cat(sprintf("Annotating %s with %s annotations...\n",sample,anno))
          
          # Overlap and calculate methylation status for each region in the annotation by summarising over all sites
          ov <- foverlaps(dat_sample, anno_list[[anno]], nomatch=0) %>% 
            .[,"i.end":=NULL] %>% setnames("i.start","pos") %>%
            .[,c("sample","anno") := list(sample,anno)] %>%
            # Compute number of methylated CpGs and the corresponding methylation rates
            .[,.(rate=round(mean(rate)), Nmet=sum(rate==100), N=.N), keyby=.(sample,id,anno)] %>%
            # Reorder columns
            .[,c("sample","id","anno","Nmet","N","rate")]
          
          # Store and save results
          fwrite(ov, fname.out, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
          system(sprintf("gzip -f %s",fname.out))
          
          rm(ov)
        }
      }
      rm(dat_sample)
    }
  }
}



# Concatenate everything and save it
for (i in opts$annos) {
  outfile <- sprintf("%s/%s.tsv.gz", io$out.folder, i)
  if(file.exists(outfile)) {
    cat(sprintf("File %s already exists, ignoring...\n", outfile))
  } else {
    # files <- list.files(paste(io$out.folder,"/tmp"), pattern=sprintf(".*%s.tsv.gz",i), full.names = T)
    system(sprintf("cat %s/tmp/*_%s.gz | zgrep -E -v 'sample' | gzip > %s", io$out.folder, i, outfile))
    cat("\n")
  }
}
