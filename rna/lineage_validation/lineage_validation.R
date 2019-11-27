
###################################################################################################
## Script to plot differential RNA expression between lineages, comparing EB and in vivo results ##
###################################################################################################

library(data.table)
library(purrr)
library(ggplot2)
library(ggpubr)

# This 
################
## Define I/O ##
################

io <- list()

# In vivo gastrulation data
io$basedir.1 <- "/Users/ricard/data/gastrulation"
io$sample_metadata.1 <- paste0(io$basedir.1,"/sample_metadata.txt")
io$rna.file.1 <- paste(io$basedir.1,"rna/SingleCellExperiment.rds",sep="/")
io$diff.1 <- paste(io$basedir.1,"rna/results/differential",sep="/")

# In vitro EB data
io$basedir.2 <- "/Users/ricard/data/scnmt_eb"
io$sample_metadata.2 <- paste0(io$basedir.2,"/sample_metadata.txt")
io$rna.file.2 <- paste(io$basedir.2,"rna/SingleCellExperiment.rds",sep="/")
io$diff.2 <- paste(io$basedir.2,"rna/results/differential",sep="/")

io$outdir <- "/Users/ricard/data/scnmt_eb/rna/results/differential/pdf"

####################
## Define options ##
####################

opts <- list()

# Differential RNA expression results from the gastrulation in vivo data
opts$comparisons.1 <- c(
  "E6.5Epiblast_vs_E7.5Mesoderm" = "Epiblast_vs_Mesoderm",
  "E6.5Epiblast_vs_E7.5Endoderm" = "Epiblast_vs_Endoderm",
  "E6.5E7.5Epiblast_vs_E6.5E7.5Primitive_Streak" = "Epiblast_vs_PrimitiveStreak"
)

# Differential RNA expression results from the EB in vitro data
opts$comparisons.2 <- c(
  "Epiblast_vs_Mesoderm" = "Epiblast_vs_Mesoderm",
  "Epiblast_vs_Endoderm" = "Epiblast_vs_Endoderm",
  "Epiblast_vs_PrimitiveStreak" = "Epiblast_vs_PrimitiveStreak"
)

##############################################
## Load differential RNA expression results ##
##############################################

# Load in vivo results
diff.1 <- names(opts$comparisons.1) %>% 
  map(function(x) fread(sprintf("%s/%s.txt.gz",io$diff.1,x)) %>% .[,comparison:=opts$comparisons.1[x]]) %>%
  rbindlist() %>% .[,assay:="in_vivo"]

# Load in vitro results
diff.2 <- 1:length(opts$comparisons.2) %>% 
  map(function(x) fread(sprintf("%s/%s.txt.gz",io$diff.2,names(opts$comparisons.2)[x])) %>% .[,comparison:=opts$comparisons.2[x]]) %>%
  rbindlist() %>% .[,assay:="EB"]

# Concatenate
diff <- rbind(diff.1,diff.2) %>% .[,c("symbol","logFC","comparison","assay")] %>%
  data.table::dcast(symbol+comparison~assay, value.var=c("logFC")) %>%
  .[complete.cases(.)]

##################################################################################################
## Scatter plots of differential RNA expression in vitro vs differential RNA expression in vivo ##
##################################################################################################

for (i in unique(diff$comparison)) {
  
  # Subset comparison  
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
    labs(x="differential RNA expression (EBs)", y="differential RNA expression (in vivo)") + 
    theme(legend.position = "none")

  # pdf(sprintf("%s/%s.pdf",io$outdir,i), useDingbats = F, width=4, height=4)
  print(p)
  # dev.off()
  
}
