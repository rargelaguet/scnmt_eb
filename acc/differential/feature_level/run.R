###################################################################################
## Script to compute (in parallel) differential accessibility at the feature level ##
###################################################################################

## I/O ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$script <- "/Users/ricard/NMT-seq_EB+ESC/acc/differential/feature_level/diffacc.R"
  io$outdir <- "/Users/ricard/data/NMT-seq_EB+ESC/acc/differential/feature_level"
} else {
  io$script <- "/homes/ricard/NMT-seq_EB+ESC/acc/differential/feature_level/diffacc.R"
  io$outdir <- "/hps/nobackup/stegle/users/ricard/NMT-seq_EB+ESC/acc/differential/feature_level"
  io$tmpdir <- "/hps/nobackup/stegle/users/ricard/NMT-seq_EB+ESC/acc/differential/feature_level/tmp"
}
dir.create(io$outdir, showWarnings=F)


## Options ##
opts <- list()

opts$groups <- list(
  
  # Day 2  
  "Epiblast_Day2_WT_vs_Epiblast_Day2_KO" = list(c("Epiblast_Day2_WT"), c("Epiblast_Day2_KO")),
  "All_Day2_WT_vs_All_Day2_KO" = list(c("Epiblast_Day2_WT","Primitive Streak_Day2_WT"), c("Epiblast_Day2_KO","Primitive Streak_Day2_KO")),

  # Day 5
  "All_Day5_WT_vs_All_Day5_KO" = list(c("Epiblast_Day5_WT","Primitive Streak_Day5_WT"), c("Epiblast_Day5_KO","Primitive Streak_Day5_KO","Mesoderm_Day5_KO")),
  
  # Day 7
  "Mesoderm_Day7_WT_vs_Mesoderm_Day7_KO" = list(c("Mesoderm_Day7_WT"), c("Mesoderm_Day7_KO"))

  
)


# Minimum number of observed cells per group, for each feature tested
opts$min.cells <- 5

# Genomic contexts
opts$anno <- c(
  # "CGI",
  "H3K27ac_distal_E7.5_Ect_intersect12",
  "H3K27ac_distal_E7.5_End_intersect12",
  "H3K27ac_distal_E7.5_Mes_intersect12",
  # "H3K27ac_distal_E7.5_Ect_intersect12_500",
  # "H3K27ac_distal_E7.5_End_intersect12_500",
  # "H3K27ac_distal_E7.5_Mes_intersect12_500"
  "prom_2000_2000"
)

# opts$anno <- c("H3K27ac_distal_E7.5_Ect_intersect12")

for (group in names(opts$groups)) {
  groupA <- opts$groups[[group]][[1]]
  groupB <- opts$groups[[group]][[2]]
  for (anno in opts$anno) {
    outfile <- sprintf("%s/%s_%s.txt", io$outdir, group, anno)
    lsf <- sprintf("bsub -M 2048 -n 1 -q standard -o %s/%s_%s.txt", io$tmpdir, group, anno)
    # lsf <- ""
    cmd <- sprintf("%s Rscript %s --anno %s --groupA %s --groupB %s --min.cells %d --outfile %s", 
                   lsf, io$script, anno, paste(groupA, collapse=" "), paste(groupB, collapse=" "), opts$min.cells, outfile)
    system(cmd)
  }
}
