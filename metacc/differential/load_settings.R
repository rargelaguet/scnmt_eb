#####################
## Define settings ##
#####################

## Define I/O ##
io <- list()
io$sample.metadata <- "/Users/ricard/data/NMT-seq_EB+ESC/sample_metadata.txt"
io$met.dir <- "/Users/ricard/data/NMT-seq_EB+ESC/met/feature_level"
io$acc.dir <- "/Users/ricard/data/NMT-seq_EB+ESC/acc/feature_level"
io$outdir <- "/Users/ricard/NMT-seq_EB+ESC/metacc/differential/out"
io$annos_dir <- "/Users/ricard/data/NMT-seq_EB+ESC/features/filt"
io$gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"

# Folders with the differential analysis results
io$diff.met <- "/Users/ricard/data/NMT-seq_EB+ESC/met/differential/feature_level"
io$diff.acc <- "/Users/ricard/data/NMT-seq_EB+ESC/acc/differential/feature_level"
# io$diff.rna <- "/Users/ricard/data/NMT-seq_EB+ESC/rna/differential"

## Define options ##
opts <- list()

# Define genomic contexts for methylation
opts$met.annos <- c(
  # "genebody"="Gene body",
  "prom_2000_2000"="Promoters",
  # "prom_2000_2000_cgi"="CGI promoters",
  # "prom_2000_2000_noncgi"="non-CGI promoters",
  "H3K27ac_distal_E7.5_Mes_intersect12"="Mes- enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="End- enhancers",
  "H3K27ac_distal_E7.5_Ect_intersect12"="Ect- enhancers",
  "H3K4me3_E7.5_Mes"="Mes- H3K4me3",
  "H3K4me3_E7.5_End"="End- H3K4me3",
  "H3K4me3_E7.5_Ect"="Ect- H3K4me3"
)

# Define genomic contexts for accessibility
opts$acc.annos <- c(
  # "genebody"="Gene body",
  # "prom_2000_2000"="Promoters",
  # "prom_2000_2000_cgi"="CGI promoters",
  # "prom_2000_2000_noncgi"="non-CGI promoters",
  "H3K27ac_distal_E7.5_Mes_intersect12"="Mes- enhancers",
  "H3K27ac_distal_E7.5_End_intersect12"="End- enhancers",
  "H3K27ac_distal_E7.5_Ect_intersect12"="Ect- enhancers"
  # "H3K4me3_E7.5_Mes"="Mes- H3K4me3",
  # "H3K4me3_E7.5_End"="End- H3K4me3",
  # "H3K4me3_E7.5_Ect"="Ect- H3K4me3"
)

# How to select differential hits?
#   Option 1 (more liberal): (lineage_A) vs (lineage_B,lineage_C)
#   Option 2 (more conservative): (lineage_A vs lineage_B) AND (lineageA vs lineage_C)
opts$diff.type <- 2
opts$min.fdr <- 0.10
opts$min.met.diff <- 10
opts$min.acc.diff <- 10

######################
## Define functions ##
######################

gg_barplot <- function(tmp, title = "", ylim=NULL) {
  
  if (is.null(ylim)) {
    ylim <- c(min(tmp$value, na.rm=T), max(tmp$value, na.rm=T))
  }
  
  p <- ggplot(tmp, aes(x=anno, y=value)) +
    geom_bar(aes(fill=assay), color="black", stat="identity", position="dodge", size=0.25) +
    scale_fill_manual(values=c("met"="#F37A71", "acc"="#00BFC4")) +
    geom_hline(yintercept=0, color="black") +
    scale_y_continuous(limits=c(ylim[1],ylim[2])) +
    labs(title=i, x="", y="Number of hits") +
    theme_bw() +
    theme(
      plot.title = element_text(size=11, face='bold', hjust=0.5),
      axis.text = element_text(size=rel(1.0), color='black'),
      axis.text.x = element_text(size=rel(1.0), angle=60, hjust=1, vjust=1, color="black"),
      axis.ticks.x = element_blank(),
      axis.title = element_text(size=rel(1.0), color='black'),
      axis.line = element_line(color="black"),
      legend.position="none"
    )
  
  return(p)
}

theme_pub <- function() {
  theme_bw() +
  theme(
    axis.text.x = element_text(size=rel(1.2), angle=60, hjust=1, vjust=1, color="black"),
    axis.text.y = element_text(size=rel(1.2), color="black"),
    axis.title.y = element_text(size=rel(1.2), color="black"),
    legend.position = "right"
    )
}
