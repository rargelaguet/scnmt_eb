librar
srat <- readRDS("/Users/ricard/data/gastrulation/mapping_10x/invitro_10x.rds")


srat <- NormalizeData(srat, scale.factor=1000)
srat <- ScaleData(
  srat, 
  # vars.to.regress=c("nCount_RNA"),
  model.use = "linear",
  do.scale = FALSE, do.center = TRUE, 
  do.par = TRUE, num.cores = 4
)



# Select highly variable genes
# srat <- FindVariableGenes(srat, x.low.cutoff = 0.001, y.cutoff = 0.5, do.plot = F, selection.method = "dispersion", top.genes=1000)
srat <- FindVariableFeatures(srat, nfeatures=2500)

srat <- RunPCA(srat, pcs.compute = 15)
srat <- RunUMAP(srat, n_neighbors = 15, min_dist = 0.7, dims = 2:15)

# DimPlot(srat, reduction.use="pca", group.by = "sample", cols.use = opts$colors, dim.1=13, dim.2=14)
# DimPlot(srat, reduction.use="umap", group.by = "celltype", cols.use = opts$colors)
DimPlot(srat, reduction.use="umap")



sce <- as.SingleCellExperiment(srat)

sample_metadata <- colData(sce) %>% as.data.frame %>% as.data.table(keep.rownames=T) %>%
  setnames("rn","cell") %>%
  .[,c("cell")] %>%
  .[,condition:=sapply(strsplit(cell, split = "_"), FUN = "[", 1)]
  
fwrite(sample_metadata, "/Users/ricard/data/invitro_10x/mapping/sample_metadata.txt", sep="\t", row.names=F, col.names=T)

foo <- readRDS("/Users/ricard/data/gastrulation/rna/parsed/SingleCellExperiment.rds")

dim(sce)
rownames(foo) <- rowData(foo)$symbol
rownames(foo) <- rowData(foo)$id

sce <- sce[rownames(sce) %in% rownames(foo),]
foo <- foo[rownames(sce),]
rowData(sce) <- rowData(foo)
rownames(sce) <- rowData(sce)$id
sce = calculateQCMetrics(sce)

saveRDS(sce, "/Users/ricard/data/invitro_10x/rna/SingleCellExperiment.rds")
