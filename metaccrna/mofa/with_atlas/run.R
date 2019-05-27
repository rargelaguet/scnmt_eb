suppressMessages(library(MOFA))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(ggplot2))
suppressMessages(library(scater))
suppressMessages(library(reticulate))
suppressMessages(library(argparse))

# print(py_config())
# print(.libPaths())

p <- ArgumentParser(description='')
p$add_argument('-i','--trial', type="character", help='Trial number')
args <- p$parse_args(commandArgs(TRUE))
if (is.null(args$trial)) args$trial <- ""

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/NMT-seq_EB+ESC/metaccrna/mofa/with_atlas/load_settings.R")
  use_python("/Users/ricard/anaconda3/bin/python")
} else {
  source("/homes/ricard/NMT-seq_EB+ESC/metaccrna/mofa/with_atlas/load_settings.R")
  use_python("/nfs/software/stegle/users/ricard/conda-envs/py3/bin/python")
}


if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/NMT-seq_EB+ESC/metaccrna/mofa/with_atlas/load_data.R")
} else {
  source("/homes/ricard/NMT-seq_EB+ESC/metaccrna/mofa/with_atlas/load_data.R")
}


# Create MOFAobject
MOFAobject <- createMOFAobject(all_matrix_list)

# Data processing options
DataOptions <- getDefaultDataOptions()

# Model options
ModelOptions <- getDefaultModelOptions(MOFAobject)
ModelOptions$numFactors <- 10

# Training options
TrainOptions <- getDefaultTrainOptions()
TrainOptions$maxiter <- 5000
TrainOptions$tolerance <- 0.50
TrainOptions$DropFactorThreshold <- 0
TrainOptions$seed <- 1

# Prepare MOFAobject for training
MOFAmodel <- prepareMOFA(MOFAobject, 
                         DataOptions = DataOptions, 
                         ModelOptions = ModelOptions, 
                         TrainOptions = TrainOptions
)

# Train the model
outfile <- sprintf("%s/hdf5/test2.hdf5",io$outdir)
model <- runMOFA(MOFAmodel, outfile)
