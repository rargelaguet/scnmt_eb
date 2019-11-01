#################################################################
## Script to compute (in parallel) differential RNA expression ##
#################################################################

## I/O ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$script <- "/Users/ricard/scnmt_eb/rna/differential/differential_expr.R"
  io$outdir <- "/Users/ricard/data/NMT-seq_EB+ESC/rna/differential"; dir.create(io$outdir)
} else {
  io$script <- "/homes/ricard/scnmt_eb/rna/differential/differential_expr.R"
  io$outdir <- "/hps/nobackup/stegle/users/ricard/NMT-seq_EB+ESC/rna/differential"; dir.create(io$outdir, showWarnings=F)
  io$tmpdir <- "/hps/nobackup/stegle/users/ricard/NMT-seq_EB+ESC/rna/differential/tmp"; dir.create(io$tmpdir, showWarnings=F)
}


## Options ##
opts <- list()

opts$groups <- list(
  
  # "Epiblast_vs_Mesoderm" = list(c("Epiblast_WT"), c("Mesoderm_WT")),
  # "Epiblast_vs_Endoderm" = list(c("Epiblast_WT"), c("Endoderm_WT")),
  # "Epiblast_vs_PrimitiveStreak" = list(c("Epiblast_WT"), c("Primitive Streak_WT"))
  # "Epiblast_vs_MesodermEndoderm" = list(c("Epiblast_WT"), c("Mesoderm_WT","Endoderm_WT")),

  # "Endoderm_vs_MesodermEpiblast" = list(c("Endoderm_WT"), c("Epiblast_WT","Mesoderm_WT")),
  # "Endoderm_vs_Mesoderm" = list(c("Endoderm_WT"), c("Mesoderm_WT")),
  # "Endoderm_vs_Epiblast" = list(c("Endoderm_WT"), c("Epiblast_WT")),

  # "Mesoderm_vs_EndodermEpiblast" = list(c("Mesoderm_WT"), c("Epiblast_WT","Endoderm_WT")),
  # "Mesoderm_vs_Epiblast" = list(c("Mesoderm_WT"), c("Epiblast_WT")),
  # "Mesoderm_vs_Endoderm" = list(c("Mesoderm_WT"), c("Endoderm_WT"))
  
)

for (group in names(opts$groups)) {
  lineage_genotype1 <- opts$groups[[group]][[1]]
  lineage_genotype2 <- opts$groups[[group]][[2]]
  outfile <- sprintf("%s/%s.txt", io$outdir, group)
  # lsf <- sprintf("bsub -M 8192 -n 1 -q standard -o %s/%s.txt", io$tmpdir, group)
  lsf <- ""
  cmd <- sprintf("%s Rscript %s --lineage_genotype1 %s --lineage_genotype2 %s --outfile %s", 
                 lsf, io$script, paste(lineage_genotype1, collapse=" "), paste(lineage_genotype2, collapse=" "), outfile)
  system(cmd)
}
