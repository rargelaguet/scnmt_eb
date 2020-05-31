suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))

################
## Define I/O ##
################

source("/Users/ricard/scnmt_eb/settings.R")
io$data <- paste0(io$basedir,"/met/cpg_level/pseudobulk")
io$outdir <- paste0(io$basedir,"/met/results/hemimethylation")

####################
## Define options ##
####################

# Define which samples to use
# opts$samples <- sample_metadata %>% .[!is.na(id_met) & pass_metQC==TRUE,id_met]
opts$samples <- c("WT","KO")

###############
## Load data ##
###############

# Load methylation data and calculate hemimethylation 
for (i in opts$samples) {
  file <- sprintf("%s/%s.tsv.gz",io$data,i)

  # Load sample methylation data
  data <- fread(file, colClasses=c("chr"="factor", "pos"="integer", "met_reads"="integer", "nonmet_reads"="integer", "rate"="numeric"))
  
  data %>% .[,total_reads:=met_reads+nonmet_reads]
  
  # Get genomic sequence
  seq <- unname( as.character(getSeq(Mmusculus, paste0("chr",data$chr), data$pos-1, data$pos+1)) )
  data[,c("base_up","base","base_down") := list(substr(seq,1,1),substr(seq,2,2),substr(seq,3,3))]
  
  print(table(data$base))
  
  data <- data %>%
    .[,strand:=ifelse(base=="C","+","-")] %>%   # Add strand information
    .[,pos:=ifelse(strand=="-",pos,pos+1)] %>%  # Convert all positions to the positive strand
    .[,N:=.N,by=c("chr","pos")] %>%
    .[N==2,.(rate=mean(rate)),c("chr","pos")]
  
  print(table(data$rate))

  stats_full[[i]] <- data %>%
    .[,.(number_reads=sum(total_reads), hemimethylation=mean(!rate%in%c(0,1))), by=c("sample","chr","pos")]
}

fwrite(rbindlist(stats_short), paste0(io$outdir,"/hemimethylation_short.txt.gz"), sep="\t")
fwrite(rbindlist(stats_full), paste0(io$outdir,"/hemimethylation_full.txt.gz"), sep="\t")
  