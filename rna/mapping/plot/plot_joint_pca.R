
source("/Users/ricard/scnmt_eb/settings.R")

#########
## I/O ##
#########

io$mapping <- paste0(io$basedir,"/rna/results/mapping/mapping.rds")
io$outdir <- paste0(io$basedir,"/rna/results/mapping/pdf")

#############
## Options ##
#############

###############
## Load data ##
###############

# Load joint PCA
pca <- readRDS(io$mapping)$pca
colnames(pca) <- paste0("PC",1:50)

# Load scNMT metadata (query)
meta_query <- sample_metadata %>%
  .[,c("id_rna","day2","lineage10x")] %>%
  setnames(c("cell","stage","celltype")) %>%
  .[cell%in%rownames(pca)] %>%
  .[,class:="scNMT"]

# Load 10x metadata (query)
io$atlas.metadata <- "/Users/ricard/data/gastrulation10x/sample_metadata.txt.gz"
meta_atlas <- fread(io$atlas.metadata) %>%
  .[,c("cell","stage","celltype")] %>%
  .[cell%in%rownames(pca)] %>%
  .[,class:="10x"]


################
## Parse data ##
################

# Concatenate metadata
meta_merged <- rbind(meta_atlas,meta_query)

# Create data.table to plot
to.plot <- pca %>% as.data.frame %>% tibble::rownames_to_column("cell") %>% 
  as.data.table %>% 
  merge(meta_merged,by="cell") 

to.plot %>% 
  .[,celltype:=stringr::str_replace_all(celltype,"_"," ")]# %>%
  # .[,celltype:=stringr::str_replace_all(celltype,"/"," ")]

unique(to.plot$celltype)

###################
## Visualise PCA ##
###################

p <- GGally::ggpairs(to.plot, 
 columns = paste0("PC",1:3),
 lower = list(continuous=GGally::wrap("points", size=0.5)), 
 diag = list(continuous='densityDiag'), 
 upper = list(continuous=GGally::wrap("points", size=0.5)), 
 mapping = aes_string(color="class"), 
 title = "",
) + 
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# for(i in 1:p$nrow) {
#   for(j in 1:p$ncol){
#     p[i,j] <- p[i,j] + scale_fill_manual(values=opts$celltype.colors)
#   }
# }

# pdf(sprintf("%s/joint_pca.pdf",io$outdir), width=7, height=7)
# png(sprintf("%s/joint_pca.png",io$outdir), width=7, height=7, units="in", res=400)
print(p)
# dev.off()

##############
## Run UMAP ##
##############

umap_embedding <- uwot::umap(pca[,2:ncol(pca)],
  n_neighbors = 20, 
  min_dist = 0.45, 
  metric = "cosine"
)


###############
## Plot UMAP ##
###############

to.plot.umap <- umap_embedding %>% as.data.table(keep.rownames = T) %>%
  .[,cell:=rownames(pca)] %>%
  merge(to.plot[,grep("PC",colnames(to.plot),invert=T),with=F], by="cell")

p <- ggplot(to.plot.umap, aes(x=V1, y=V2)) +
  # geom_point(aes(fill=class), shape=21, size=1, stroke=0.05) +
  geom_point(aes(color=class), size=1) +
  # scale_fill_manual(values=colors) +
  # labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# pdf(sprintf("%s/joint_umap.pdf",io$outdir), width=5, height=5)
print(p)
# dev.off()