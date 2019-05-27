library(scater)
library(data.table)
library(purrr)
library(ggplot2)
library(ggpubr)

# source("/Users/ricard/gastrulation/rna/differential/utils.R")

################
## Define I/O ##
################

io <- list()

io$basedir.1 <- "/Users/ricard/data/gastrulation"
io$sample_metadata.1 <- paste0(io$basedir.1,"/sample_metadata.txt")
io$rna.file.1 <- paste(io$basedir.1,"rna/parsed/SingleCellExperiment.rds",sep="/")
io$diff.1 <- paste(io$basedir.1,"rna/differential",sep="/")

io$basedir.2 <- "/Users/ricard/data/NMT-seq_EB+ESC"
io$sample_metadata.2 <- paste0(io$basedir.2,"/sample_metadata.txt")
io$rna.file.2 <- paste(io$basedir.2,"rna/SingleCellExperiment.rds",sep="/")
io$diff.2 <- paste(io$basedir.2,"rna/differential",sep="/")

io$outdir <- "/Users/ricard/data/NMT-seq_EB+ESC/rna/differential/pdf"

####################
## Define options ##
####################

opts <- list()
opts$comparisons.1 <- c(
  "E7.5Endoderm_vs_E7.5Mesoderm" = "Endoderm_vs_Mesoderm",
  
  "E4.5Epiblast_vs_E7.5Mesoderm" = "E4.5Epiblast_vs_Mesoderm",
  "E5.5Epiblast_vs_E7.5Mesoderm" = "E5.5Epiblast_vs_Mesoderm",
  "E6.5Epiblast_vs_E7.5Mesoderm" = "E6.5Epiblast_vs_Mesoderm",
  "E7.5Epiblast_vs_E7.5Mesoderm" = "E7.5Epiblast_vs_Mesoderm",
  "E7.5Ectoderm_vs_E7.5Mesoderm" = "E7.5Ectoderm_vs_Mesoderm",
  
  "E4.5Epiblast_vs_E7.5Endoderm" = "E4.5Epiblast_vs_Endoderm",
  "E5.5Epiblast_vs_E7.5Endoderm" = "E5.5Epiblast_vs_Endoderm",
  "E6.5Epiblast_vs_E7.5Endoderm" = "E6.5Epiblast_vs_Endoderm",
  "E7.5Epiblast_vs_E7.5Endoderm" = "E7.5Epiblast_vs_Endoderm",
  "E7.5Ectoderm_vs_E7.5Endoderm" = "E7.5Ectoderm_vs_Endoderm"
)

opts$comparisons.2 <- c(
  "Endoderm_vs_Mesoderm" = "Endoderm_vs_Mesoderm",
  
  "Epiblast_vs_Mesoderm" = "E4.5Epiblast_vs_Mesoderm",
  "Epiblast_vs_Mesoderm" = "E5.5Epiblast_vs_Mesoderm",
  "Epiblast_vs_Mesoderm" = "E6.5Epiblast_vs_Mesoderm",
  "Epiblast_vs_Mesoderm" = "E7.5Epiblast_vs_Mesoderm",
  "Epiblast_vs_Mesoderm" = "E7.5Ectoderm_vs_Mesoderm",
  
  "Epiblast_vs_Endoderm" = "E4.5Epiblast_vs_Endoderm",
  "Epiblast_vs_Endoderm" = "E5.5Epiblast_vs_Endoderm",
  "Epiblast_vs_Endoderm" = "E6.5Epiblast_vs_Endoderm",
  "Epiblast_vs_Endoderm" = "E7.5Epiblast_vs_Endoderm",
  "Epiblast_vs_Endoderm" = "E7.5Ectoderm_vs_Endoderm"
)

##############################################
## Load differential RNA expression results ##
##############################################

diff.1 <- names(opts$comparisons.1) %>% 
  map(function(x) fread(cmd=sprintf("zcat < %s/%s.txt.gz",io$diff.1,x)) %>% .[,comparison:=opts$comparisons.1[x]]) %>%
  rbindlist() %>% .[,assay:="in_vivo"]

diff.2 <- 1:length(opts$comparisons.2) %>% 
  map(function(x) fread(cmd=sprintf("zcat < %s/%s.txt.gz",io$diff.2,names(opts$comparisons.2)[x])) %>% .[,comparison:=opts$comparisons.2[x]]) %>%
  rbindlist() %>% .[,assay:="EB"]


diff <- rbind(diff.1,diff.2) %>% .[,c("symbol","logFC","comparison","assay")] %>%
  data.table::dcast(symbol+comparison~assay, value.var=c("logFC")) %>%
  .[complete.cases(.)]

##########
## Plot ##
##########

for (i in unique(diff$comparison)) {
  
  to.plot <- diff[comparison==i]
  
  # Highlight top N markers
  marker.genes <- to.plot %>% copy %>%
    .[,sum:=abs(EB+in_vivo)] %>% setorder(-sum) %>% .[,sign:=sign(EB+in_vivo)] %>%
    .[,.(genes=head(.SD[,symbol],n=5)),by=c("sign")] %>% .$genes
  to.plot[,marker:=symbol%in%marker.genes]

  p <- ggscatter(to.plot, x="EB", y="in_vivo", size="marker", color="marker", 
            add="reg.line", add.params = list(color="black", fill="lightgray"), conf.int=TRUE) +
    stat_cor(method = "pearson") +
    ggrepel::geom_text_repel(data=to.plot[marker==T], aes(x=EB, y=in_vivo, label=symbol), size=5, color="red") +
    scale_size_manual(values=c("TRUE"=1.5, "FALSE"=0.5)) +
    scale_color_manual(values=c("TRUE"="red", "FALSE"="grey70")) +
    labs(x="differential RNA expression (Embryoid body)", y="differential RNA expression (in vivo)") + 
    theme(legend.position = "none")
  # print(p)

  pdf(sprintf("%s/%s.pdf",io$outdir,i), useDingbats = F, width=5, height=5)
  print(p)
  dev.off()
  
}
