
#####################
## Define settings ##
#####################

source("/Users/ricard/scnmt_eb/settings.R")
io$pseudobulk.met <- paste0(io$basedir,"/met/cpg_level/pseudobulk")
io$outdir <- paste0(io$basedir,"/met/results/differential/pseudobulk/cpg_level")

opts$groups <- c("KO","WT")

###############
## Load data ##
###############

data <- opts$groups %>% map(function(i) {
  file <- sprintf("%s/%s.tsv.gz",io$pseudobulk.met,i)
  fread(file, colClasses=c("chr"="factor", "pos"="integer", "met_reads"="integer", "nonmet_reads"="integer", "rate"="numeric")) %>%
  .[,group:=factor(i,levels=c("WT","KO"))]
}) %>% rbindlist

data %>% .[,total_reads:=met_reads+nonmet_reads]


###########################
## Differential analysis ##
###########################

foo <- data[chr==1] %>% 
  dcast(chr+pos~group, value.var=c("rate","total_reads")) %>%
  .[,diff:=rate_KO-rate_WT]

##########
## Save ##
##########

outfile <- paste0(io$outdir,"/differential_pseudobulk_cpglevel.tsv.gz")
fwrite(foo, outfile, na="NA", quote=F, sep="\t")
