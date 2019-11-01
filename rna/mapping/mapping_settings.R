##########################################################################################
## Functions for mapping scRNAseq experiments to the 10X mouse gastrulation atlas 
##
## Code adapted from: 
## https://github.com/MarioniLab/EmbryoTimecourse2018/blob/master/analysis_scripts/atlas/core_functions.R
##
## Author: Ivan Imaz
##
## Packages to be installed:
## scran, scater, biomaRt, Matrix, irlba, BiocNeighbors, SingleCellExperiment
##########################################################################################

#ensure counts has columns names for the cells
#match timepoints,samples to the count table
#timepoint_order, sample_order should contain each sample/timepoint ONCE, in correct order for correction
doBatchCorrect <- function(counts, timepoints, samples, timepoint_order, sample_order, npc = 50, pc_override = NULL, BPPARAM = SerialParam()){
  require(BiocParallel)
  
  if(!is.null(pc_override)){
    pca = pc_override
  } else {
    pca = irlba::prcomp_irlba(t(counts), n = npc)$x
    rownames(pca) = colnames(counts)
  }
  
  if(length(unique(samples)) == 1){
    return(pca)
  }
  
  #create nested list
  pc_list    <- lapply(unique(timepoints), function(tp){
    sub_pc   <- pca[timepoints == tp, , drop = FALSE]
    sub_samp <- samples[timepoints == tp]
    list     <- lapply(unique(sub_samp), function(samp){
      sub_pc[sub_samp == samp, , drop = FALSE]
    })
    names(list) <- unique(sub_samp)
    return(list)
  })
  
  names(pc_list) <- unique(timepoints)
  
  #arrange to match timepoint order
  pc_list <- pc_list[order(match(names(pc_list), timepoint_order))]
  pc_list <- lapply(pc_list, function(x){
    x[order(match(names(x), sample_order))]
  })
  
  #perform corrections within list elements (i.e. within stages)
  correct_list <- lapply(pc_list, function(x){
    if(length(x) > 1){
      return(do.call(scran::fastMNN, c(x, "pc.input" = TRUE, BPPARAM = BPPARAM))$corrected)
    } else {
      return(x[[1]])
    }
  })
  
  #perform correction over list
  if(length(correct_list)>1){
    correct <- do.call(scran::fastMNN, c(correct_list, "pc.input" = TRUE, BPPARAM = BPPARAM))$corrected
  } else {
    correct <- correct_list[[1]]
  }
  
  correct <- correct[match(colnames(counts), rownames(correct)),]
  
  return(correct)
  
}

getHVGs <- function(sce, min.mean = 1e-3, chrY.genes.file = NULL){
  trend  <- scran::trendVar(sce, use.spikes = FALSE, loess.args = list(span = 0.05))
  decomp <- scran::decomposeVar(sce, fit = trend)
  decomp <- decomp[decomp$mean > min.mean,]
  
  #exclude sex genes
  xist <- "ENSMUSG00000086503"
  if(!is.null(chrY.genes.file)){
    #ychr <- read.table("/nfs/research1/marioni/jonny/embryos/data/ygenes.tab", stringsAsFactors = FALSE)[,1]
    ychr <- read.table(chrY.genes.file, stringsAsFactors = FALSE)[,1]
  }else{
    mouse_ensembl <- biomaRt::useMart("ensembl")
    mouse_ensembl <- biomaRt::useDataset("mmusculus_gene_ensembl", mart = mouse_ensembl)
    gene_map <- biomaRt::getBM(attributes=c("ensembl_gene_id", "chromosome_name"),
      filters = "ensembl_gene_id", values = rownames(decomp), mart = mouse_ensembl)
    ychr <- gene_map[gene_map[,2] == "Y", 1]  
  }
  
  other = c("tomato-td") #for the chimera
  decomp = decomp[!rownames(decomp) %in% c(xist, ychr, other),]
  
  decomp$FDR = p.adjust(decomp$p.value, method = "fdr")
  return(rownames(decomp)[decomp$p.value < 0.05])
}

#MAPPING FUNCTIONS

getmode <- function(v, dist) {
  tab <- table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied <- names(tab)[tab == max(tab)]
    sub  <- dist[v %in% tied]
    names(sub) <- v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
}

getcelltypes <- function(v, dist) {
  tab <- table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied <- names(tab)[tab == max(tab)]
    sub  <- dist[v %in% tied]
    names(sub) <- v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
}

getMappingScore <- function(mapping){
  celltypes_accrossK <- matrix(unlist(mapping$celltypes.mapped), 
    nrow=length(mapping$celltypes.mapped[[1]]),
    ncol=length(mapping$celltypes.mapped))
  P <- NULL
  for (i in 1:nrow(celltypes_accrossK)){
    p <- max(table(celltypes_accrossK[i,]))
    index <- which(table(celltypes_accrossK[i,]) == p)
    p <- p/length(mapping$celltypes.mapped)
    P <- c(P,p) 
  }
  return(P)  
}

mnnMap <- function(atlas_pca, atlas_meta, map_pca, map_meta, k_map = 10){
  correct <- scran::fastMNN(atlas_pca, map_pca, pc.input = TRUE)$corrected
  atlas   <- 1:nrow(atlas_pca)
  correct_atlas <- correct[atlas,]
  correct_map   <- correct[-atlas,]
  
  knns <- BiocNeighbors::queryKNN(correct_atlas, correct_map, k = k_map, get.index = TRUE,
    get.distance = FALSE)
  
  #get closest k matching cells
  k.mapped  <- t(apply(knns$index, 1, function(x) atlas_meta$cell[x]))
  celltypes <- t(apply(k.mapped, 1, function(x) atlas_meta$celltype[match(x, atlas_meta$cell)]))
  stages    <- t(apply(k.mapped, 1, function(x) atlas_meta$stage[match(x, atlas_meta$cell)]))
  celltype.mapped <- apply(celltypes, 1, function(x) getmode(x, 1:length(x)))
  stage.mapped    <- apply(stages, 1, function(x) getmode(x, 1:length(x)))
  
  out <- lapply(1:length(celltype.mapped), function(x){
    list(cells.mapped     = k.mapped[x,],
         celltype.mapped  = celltype.mapped[x],
         stage.mapped     = stage.mapped[x],
         celltypes.mapped = celltypes[x,],
         stages.mapped    = stages[x,])
  })
  
  names(out) <- map_meta$cell
  
  return(out)  
  
}

#atlas_meta MUST have a cell column, celltype column and a stage column, spelled exactly like that
#map_meta MUST have a cell column and a stage column, spelled exactly like that
mapWrap <- function(atlas_sce, atlas_meta, map_sce, map_meta, k = 30, mapstage_x=NULL, map2stage_x=NULL, return.list = FALSE){
  
  #prevent duplicate rownames
  colnames(map_sce) <- paste0("map_", colnames(map_sce))
  map_meta$cell     <- paste0("map_", map_meta$cell)

  message("Normalizing joint dataset...")
  #easier to avoid directly binding sce objects as it is a lot more likely to have issues
  sce_all <- SingleCellExperiment::SingleCellExperiment(
    list(counts=Matrix::Matrix(cbind(counts(atlas_sce),counts(map_sce)),sparse=TRUE)))
  big_sce <- scater::normalize(sce_all)
  message("Done\n")

  message("Computing highly variable genes...")
  hvgs    <- getHVGs(big_sce)
  message("Done\n")
  
  message("Performing PCA...")
  big_pca <- irlba::prcomp_irlba(t(logcounts(big_sce[hvgs,])), n = 50)$x
  rownames(big_pca) <- colnames(big_sce) 
  atlas_pca <- big_pca[1:ncol(atlas_sce),]
  map_pca   <- big_pca[-(1:ncol(atlas_sce)),]
  message("Done\n")
  
  message("Batch effect correction...")  
  #correct the atlas first
  order_df        <- atlas_meta[!duplicated(atlas_meta$sample), c("stage", "sample")]
  order_df$ncells <- sapply(order_df$sample, function(x) sum(atlas_meta$sample == x))
  order_df$stage  <- factor(order_df$stage, 
                          levels = rev(c("E8.5", 
                                         "E8.25", 
                                         "E8.0", 
                                         "E7.75", 
                                         "E7.5", 
                                         "E7.25", 
                                         "mixed_gastrulation", 
                                         "E7.0", 
                                         "E6.75", 
                                         "E6.5")))
  order_df       <- order_df[order(order_df$stage, order_df$ncells, decreasing = TRUE),]
  order_df$stage <- as.character(order_df$stage)
  
  set.seed(42)
  atlas_corrected <- doBatchCorrect(counts         = logcounts(atlas_sce[hvgs,]), 
                                   timepoints      = atlas_meta$stage, 
                                   samples         = atlas_meta$sample, 
                                   timepoint_order = order_df$stage, 
                                   sample_order    = order_df$sample, 
                                   pc_override     = atlas_pca)
  message("Done\n")
           
  message("MNN mapping...")                        
  if(!is.null(map2stage_x)){
    atlas_stage_index <- NULL
    for (i in seq(from = 1, to = length(map2stage_x))){
      index1 <- which(atlas_meta$stage==map2stage_x[i])
      atlas_stage_index <- c(atlas_stage_index,index1)
    } 
  }
  if(!is.null(mapstage_x)){
    map_stage_index   <- NULL
    for (i in seq(from = 1, to = length(mapstage_x))){    
      index2 <- which(map_meta$stage==mapstage_x[i])
      map_stage_index <- c(map_stage_index,index2)
    } 
    map_meta$stage <- map_meta$stage[map_stage_index]
    map_meta$cell  <- map_meta$cells[map_stage_index]
  }
  
  if(!is.null(map2stage_x) & !is.null(mapstage_x)){
    mapping <- mnnMap(atlas_pca = atlas_corrected[atlas_stage_index,],
                   atlas_meta = atlas_meta[atlas_stage_index,],
                   map_pca = map_pca[map_stage_index,],
                   map_meta = map_meta,k_map = k) 
  }else if(!is.null(map2stage_x) & is.null(mapstage_x)){
    mapping <- mnnMap(atlas_pca = atlas_corrected[atlas_stage_index,],
                   atlas_meta = atlas_meta[atlas_stage_index,],
                   map_pca = map_pca,
                   map_meta = map_meta,k_map = k)   
  }else if(is.null(map2stage_x) & !is.null(mapstage_x)){
    mapping <- mnnMap(atlas_pca = atlas_corrected,
                   atlas_meta = atlas_meta,
                   map_pca = map_pca[map_stage_index,],
                   map_meta = map_meta,k_map = k) 
  }else{
    mapping <- mnnMap(atlas_pca = atlas_corrected,
                   atlas_meta = atlas_meta,
                   map_pca = map_pca,
                   map_meta = map_meta,k_map = k)
  }
  message("Done\n")
             
  if(return.list){
    return(mapping)
  }
  
  message("Computing mapping scores...") 
  out <- list()
  for (i in seq(from = 1, to = k)) {
   out$closest.cells[[i]]     <- sapply(mapping, function(x) x$cells.mapped[i])
   out$celltypes.mapped[[i]]  <- sapply(mapping, function(x) x$celltypes.mapped[i])
   out$cellstages.mapped[[i]] <- sapply(mapping, function(x) x$stages.mapped[i])
  }  
  celltype.multinomial.prob <- getMappingScore(out)
  message("Done\n")
  
  message("Writing output...") 
  out$mapping <- data.frame(
  cell            = names(mapping), 
  celltype.mapped = sapply(mapping, function(x) x$celltype.mapped),
  stage.mapped    = sapply(mapping, function(x) x$stage.mapped),
  closest.cell    = sapply(mapping, function(x) x$cells.mapped[1]))

  out$mapping <- cbind(out$mapping,celltype.multinomial.prob)
  message("Done\n")
  
  return(out)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour", "Publication",
                 manual_pal(values = c(
                   "#000000", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                   "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                   "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                   "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                   "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                   "#372101", "#FFB500", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                   "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                   "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                   "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                   "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                   "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                   "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                   "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")), ...)
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill", "Publication",
                 manual_pal(values = c(
                   "#000000", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                   "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                   "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                   "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                   "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                   "#372101", "#FFB500", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                   "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                   "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                   "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                   "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                   "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                   "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                   "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")), ...)
}