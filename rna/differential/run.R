#################################################################
## Script to compute (in parallel) differential RNA expression ##
#################################################################

## I/O ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$script <- "/Users/ricard/NMT-seq_EB+ESC/rna/differential/differential_expr.R"
  io$outdir <- "/Users/ricard/data/NMT-seq_EB+ESC/rna/differential"; dir.create(io$outdir)
} else {
  io$script <- "/homes/ricard/NMT-seq_EB+ESC/rna/differential/differential_expr.R"
  io$outdir <- "/hps/nobackup/stegle/users/ricard/NMT-seq_EB+ESC/rna/differential"; dir.create(io$outdir, showWarnings=F)
  io$tmpdir <- "/hps/nobackup/stegle/users/ricard/NMT-seq_EB+ESC/rna/differential/tmp"; dir.create(io$tmpdir, showWarnings=F)
}



## Options ##
opts <- list()

opts$groups <- list(
  
  "Epiblast_vs_Mesoderm" = list(c("Epiblast"), c("Mesoderm")),
  "Epiblast_vs_Endoderm" = list(c("Epiblast"), c("Endoderm")),
  "Epiblast_vs_MesodermEndoderm" = list(c("Epiblast"), c("Mesoderm","Endoderm")),

  "Endoderm_vs_MesodermEpiblast" = list(c("Endoderm"), c("Epiblast","Mesoderm")),
  "Endoderm_vs_Mesoderm" = list(c("Endoderm"), c("Mesoderm")),
  "Endoderm_vs_Epiblast" = list(c("Endoderm"), c("Epiblast")),
  
  "Mesoderm_vs_EndodermEpiblast" = list(c("Mesoderm"), c("Epiblast","Endoderm")),
  "Mesoderm_vs_Epiblast" = list(c("Mesoderm"), c("Epiblast")),
  "Mesoderm_vs_Endoderm" = list(c("Mesoderm"), c("Endoderm"))
  
)

for (group in names(opts$groups)) {
  lineage1 <- opts$groups[[group]][[1]]
  lineage2 <- opts$groups[[group]][[2]]
  outfile <- sprintf("%s/%s.txt", io$outdir, group)
  # lsf <- sprintf("bsub -M 8192 -n 1 -q standard -o %s/%s.txt", io$tmpdir, group)
  lsf <- ""
  cmd <- sprintf("%s Rscript %s --lineage1 %s --lineage2 %s --outfile %s", 
                 lsf, io$script, paste(lineage1, collapse=" "), paste(lineage2, collapse=" "), outfile)
  system(cmd)
}
