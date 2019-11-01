---
title: "EB data set: preprocessing and quality control of expression data"
output:
  BiocStyle::html_document: 
    fig_width: 15
    fig_height: 8
---

```{r echo=FALSE, include=FALSE}
library(SingleCellExperiment)
library(data.table)
library(purrr)
library(scater)
library(scran)
library(ggplot2)
# source("/Users/ricard/NMT-seq_ESC/Rutils/stats_utils.R")
```

```{r define_opts, echo=FALSE, include=FALSE}

## Options ##
opts <- list()

# Stringent thresholds
# opts$coverage_threshold <- 2e5    # Minimum library size (coverage)
# opts$features_threshold <- 1000   # Minimum number of expressed features
# opts$top50_threshold <- 0.75      # Maximum fraction of reads accounting for the top 50 features
# opts$MT_threshold <- 0.25         # Maximum fraction of reads mapping to mithocondrial genes

# Lenient thresholds
opts$coverage_threshold <- 1e5    # Minimum library size (coverage)
opts$features_threshold <- 500    # Minimum number of expressed features
opts$top50_threshold <- 0.75      # Maximum fraction of reads accounting for the top 50 features
opts$MT_threshold <- 0.25         # Maximum fraction of reads mapping to mithocondrial genes

## I/O ##
io <- list()
io$in.gene_metadata <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
io$in.sample_metadata <- "/Users/ricard/data/NMT-seq_EB+ESC/sample_metadata.txt"
io$in.counts <- "/Users/ricard/data/NMT-seq_EB+ESC/rna/counts.txt.gz"
io$out.file <- "/Users/ricard/data/NMT-seq_EB+ESC/rna/SingleCellExperiment.rds"
io$out.sample_metadata <- "/Users/ricard/data/NMT-seq_EB+ESC/sample_metadata2.txt"
```

# Load counts
```{r load_data, echo=FALSE}
counts <- fread(cmd=sprintf("zcat < %s",io$in.counts)) %>% matrix.please
```

# Load cell metadata
```{r load_sample_metadata, echo=FALSE}
sample_metadata <- fread(io$in.sample_metadata) %>% setkey(id_rna)
```

# Load feature metadata
```{r load_feature_metadata, echo=FALSE}
feature_metadata <- read.csv("/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt", sep="\t", stringsAsFactors=FALSE, quote="", header=T)

# Define mithocondrial genes
mt <- feature_metadata$symbol[feature_metadata$chr == "chrMT"]

# remove duplicated genes
feature_metadata <- feature_metadata[!duplicated(feature_metadata$symbol),] %>%
  tibble::remove_rownames() %>% tibble::column_to_rownames("ens_id")
```

<!-- Parse data -->
```{r parse_data}
# Filter genes
genes <- rownames(feature_metadata[rownames(feature_metadata) %in% rownames(counts),])
feature_metadata <- feature_metadata[genes,]
counts <- counts[rownames(feature_metadata),]
```

<!-- Create SingleCellExperiment object -->
```{r echo=FALSE}

# Create rowData
rowData <- feature_metadata[rownames(counts),] %>% GRanges()

# Create colData
sample_metadata$id_rna[!sample_metadata$id_rna%in%colnames(counts)]
colData <- sample_metadata %>% as.data.frame %>% 
  tibble::remove_rownames() %>% tibble::column_to_rownames("id_rna") %>% 
  .[colnames(counts),]

# create SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), rowData=rowData, colData=colData)

# Calculate QC metrics
sce = calculateQCMetrics(sce, feature_controls=list(Mt=rownames(sce) %in% mt))
```

```{r}
sce$gene_strand <- sce$strand
sce$strand <- NULL
```

# Do QC
```{r filter_samples, echo=FALSE, include=TRUE}
# Library size
libsize.drop <- sce$total_counts < opts$coverage_threshold
libsize.drop_dt <- data.table(
  sample=colnames(sce), 
  size=sce$total_counts, 
  color=c("black","red")[as.numeric(libsize.drop)+1]
) %>% setkey(size) %>% .[,col:=size] %>% .[,sample:=factor(sample,levels=sample)]

p1 <- ggplot(libsize.drop_dt, aes(x=sample, y=size)) +
  geom_bar(stat='identity', position="dodge", fill="#3CB54E") +
  geom_hline(yintercept=opts$coverage_threshold, colour="black", linetype="dashed") +
  scale_fill_gradient(low="red", high="green") +
  labs(y="Library size") +
  barplot_theme() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=rel(1.8)),
    # axis.text.x = element_text(colour="black", color=foo$color, angle=90, size=10, vjust=0.5, hjust=1.0)
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
print(p1)

# Number of expressed genes
feature.drop <- sce$total_features_by_counts < opts$features_threshold
feature.drop_dt <- data.table(
  sample = colnames(sce), 
  features = sce$total_features_by_counts, 
  color = c("black","red")[as.numeric(feature.drop)+1]
  ) %>% setkey(features) %>% .[,col:=features] %>% .[,sample:=factor(sample,levels=sample)]

p2 <- ggplot(feature.drop_dt, aes(x=sample, y=features)) +
  geom_bar(stat='identity', position="dodge", fill="#3CB54E") +
  geom_hline(yintercept=opts$features_threshold, colour="black", linetype="dashed") +
  # scale_fill_gradient(low="red", high="green") +
  labs(y="Total number of expressed genes") +
  barplot_theme() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=rel(1.8)),
    # axis.text.x = element_text(colour="black", color=foo$color, angle=90, size=10, vjust=0.5, hjust=1.0)
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
print(p2)

# Proportion of reads accounting for the top 50 features
top50.drop <- sce$pct_counts_in_top_50_features > opts$top50_threshold*100
top50.drop_dt <- data.table(
  sample=colnames(sce), 
  pct=sce$pct_counts_in_top_50_features
  ) %>% setkey(pct) %>% .[,col:=pct] %>% .[,sample:=factor(sample,levels=sample)]

p3 <- barPlot(top50.drop_dt, ylabel="Fraction of reads accounting for the top 50 features") +
  geom_hline(yintercept=opts$top50_threshold*100, colour="blue", linetype="dashed") +
  # scale_fill_gradient(low="green", high="red") +
  theme_barplot_pub() +
  theme(
    legend.position = "none",
    # axis.text.x = element_text(colour="black", color=foo$color, angle=90, size=10, vjust=0.5, hjust=1.0)
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
print(p3)
```


```{r, echo=FALSE, include=TRUE}
# Remove cells that do not pass QC
drop.samples <- colnames(sce)[( libsize.drop | feature.drop | top50.drop )]

# Update sample metadata
# sample_metadata[,pass_rnaQC:=ifelse(id_rna%in%drop.samples,FALSE,TRUE)]

# Save updated metadata
# fwrite(sample_metadata,io$out.sample_metadata, sep="\t", col.names=T, row.names=F, na="NA", quote=F)

# Re-calculate QC statistics
sce <- sce[,!colnames(sce) %in% drop.samples]
sce <- calculateQCMetrics(sce)
```

# Normalisation and log transformation
```{r normalisation, echo=FALSE, warnings=FALSE, include=TRUE}


# Temporarily remove the lowly expressed genes
sce_filt <- sce
sce_filt <- sce_filt[!(rowMeans(counts(sce)) <= 1 | rowData(sce)$pct_dropout_by_counts > 90),]

# Compute size factors without the lowly expressed genes
sf = computeSumFactors(sce_filt, positive=TRUE, sf.out=T)

ggplot(data.frame(sf=log(sf), counts=log(sce_filt$total_counts))) +
  geom_point(aes(x=sf,y=counts)) +
  labs(y="Library size (log)", x="Size factor (log)") +
  theme_bw() +
  theme(
    axis.title = element_text(colour="black", size=15),
    axis.text = element_text(colour="black", size=12)
  )

# Normalise and log transform with the lowly expressed genes
sizeFactors(sce) <- sf; sce$sizeFactor <- sf
sizeFactors(sce_filt) <- sf; sce_filt$sizeFactor <- sf
sce <- normalize(sce, exprs_values="counts")
sce_filt <- normalize(sce_filt, exprs_values="counts")

# Update quality metrics
sce = calculateQCMetrics(sce)
```

<!-- Mean vs variance plot -->
```{r echo=TRUE, include=TRUE}
foo <- data.frame(sd=apply(exprs(sce),1,sd), mean=apply(exprs(sce),1,mean))
ggplot(foo, aes(x=mean, y=sd)) +
  geom_point() +
  stat_smooth() +
  scale_color_manual(values=c("black","red")) +
  xlab('Mean') + ylab('Standard deviation')
```

<!-- Save SingleCellExperiment object -->
```{r save, echo=FALSE, include=FALSE}
saveRDS(sce, io$out.file)
```
