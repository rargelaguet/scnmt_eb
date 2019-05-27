library(data.table)
library(purrr)
library(ggplot2)
library(scater)

matrix.please <- function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

###############
## Load data ##
###############

# Load Methylation data
met_dt.1 <- lapply(opts$met.annos, function(n) {
  data <- fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$met.dir.1,n), showProgress=F, stringsAsFactors=F, quote="") %>%
    .[V1%in%opts$met_cells.1]
}) %>% rbindlist %>% setnames(c("id_met","id","anno","Nmet","N","rate"))

met_dt.2 <- lapply(opts$met.annos, function(n) {
  data <- fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$met.dir.2,n), showProgress=F, stringsAsFactors=F, quote="") %>%
    .[V1%in%opts$met_cells.2]
}) %>% rbindlist %>% setnames(c("id_met","id","anno","Nmet","N","rate"))

met_dt <- rbind(met_dt.1,met_dt.2)

# Load Accessibility data
acc_dt.1 <- lapply(opts$acc.annos, function(n) {
  data <- fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$acc.dir.1,n), showProgress=F, stringsAsFactors=F, quote="") %>%
    .[V1%in%opts$acc_cells.1]
}) %>% rbindlist %>% setnames(c("id_acc","id","anno","Nmet","N","rate"))

acc_dt.2 <- lapply(opts$acc.annos, function(n) {
  data <- fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$acc.dir.2,n), showProgress=F, stringsAsFactors=F, quote="") %>%
    .[V1%in%opts$acc_cells.2]
}) %>% rbindlist %>% setnames(c("id_acc","id","anno","Nmet","N","rate"))

acc_dt <- rbind(acc_dt.1,acc_dt.2)

# Load RNA data
sce.1 <- readRDS(io$rna.file.1) %>% .[,opts$rna_cells.1]
sce.2 <- readRDS(io$rna.file.2) %>% .[,opts$rna_cells.2]

stopifnot(all(rownames(sce.1) == rownames(sce.2)))

sce <- SingleCellExperiment::SingleCellExperiment(
  list(counts=Matrix::Matrix(cbind(counts(sce.1),counts(sce.2)),sparse=TRUE)))
rowData(sce) <- rowData(sce.1)

sce <- scater::normalize(sce)
sce <- calculateQCMetrics(sce)
sce <- sce[rowData(sce)$pct_dropout_by_counts < 97,]

# Remove temporary variables
rm(acc_dt.1,acc_dt.2)
rm(met_dt.1,met_dt.2)
rm(sce.1,sce.2)

################
## Parse data ##
################


## Parse RNA expression data ##

# Convert to data.table
rna_dt <- as.matrix(exprs(sce)) %>% t %>% as.data.table(keep.rownames = "id_rna") %>% 
  melt(id.vars = "id_rna", value.name = "expr", variable.name = "ens_id") %>%
  merge(rowData(sce) %>% as.data.frame(row.names = rownames(sce)) %>% tibble::rownames_to_column("ens_id") %>% .[,c("symbol","ens_id")] %>% setnames("symbol","gene"))

## Parse accessibility data ##
acc_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))] # Calculate M value from Beta value

## Parse methylation data ##
met_dt[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))] # Calculate M value from Beta value


##############################
## Merge data with metadata ##
##############################

acc_dt <- merge(acc_dt, sample_metadata[,c("sample","id_acc")], by="id_acc") %>% droplevels()
met_dt <- merge(met_dt, sample_metadata[,c("sample","id_met")], by="id_met") %>% droplevels()
rna_dt <- merge(rna_dt, sample_metadata[,c("sample","id_rna")], by="id_rna") %>% droplevels()

#############################
## Filter methylation data ##
#############################

# Filter features by minimum number of CpGs
met_dt <- met_dt[N>=opts$met_min.CpGs]

# Filter features by  minimum number of cells (by stage_lineage)
met_dt[,N:=.N,by=c("id","anno")]  %>% .[N>=opts$met_min.cells] %>% .[,N:=NULL]

# Filter features by variance
met_dt <- met_dt[,var:=var(m), by=c("id","anno")] %>% .[var>0] %>% .[,var:=NULL] %>% droplevels()

###############################
## Filter accessibility data ##
###############################

# Filter features by minimum number of GpCs
acc_dt <- acc_dt[N>=opts$acc_min.GpCs]

# Filter features by  minimum number of cells (by stage_lineage)
acc_dt[,N:=.N,by=c("id","anno")]  %>% .[N>=opts$acc_min.cells] %>% .[,N:=NULL]

# Filter features by variance
acc_dt <- acc_dt[,var:=var(m), by=c("id","anno")] %>% .[var>0] %>% .[,var:=NULL] %>% droplevels()

################################
## Filter RNA expression data ##
################################

# Remove lowly expressed genes
rna_dt <- rna_dt[,mean:=mean(expr),by="ens_id"] %>% .[mean>0.1] %>% .[,mean:=NULL]

# Remove genes with constant expression levels
rna_dt <- rna_dt[,var:=var(expr),by="ens_id"] %>% .[var>0.1] %>% .[,var:=NULL]

# Filter genes with low cellular detection rate and sites with low coverage across samples
rna_dt <- rna_dt[,cdr:=sum(expr>0)/length(opts$rna_cells), by="ens_id"] %>% .[cdr>=opts$rna_min.cdr] %>% .[,cdr:=NULL]

# Extract top N highly variable genes
rna_dt <- rna_dt[,var:=var(expr), by="ens_id"] %>% .[var>0] %>% .[,var:=NULL] %>% droplevels()

############################
## Regress out covariates ##
############################

# RNA: number of expressed genes
# foo <- data.table(id_rna=colnames(sce), covariate=sce$total_features_by_counts/nrow(sce))
# rna_dt <- rna_dt %>% merge(foo, by="id_rna") %>%
#   .[,expr:=lm(formula=expr~covariate)[["residuals"]], by=c("gene")] %>%
#   .[,covariate:=NULL]


# Methylation: differences in mean methylation rate
# foo <- fread("/Users/ricard/gastrulation/met/stats/samples/out/sample_stats.txt") %>%
#   .[,mean:=log2(((mean/100)+0.01)/(1-(mean/100)+0.01))] %>%
#   .[,c("id_met","mean")]
# met_dt <- met_dt %>% merge(foo, by="id_met") %>%
#   .[,m:=lm(formula=m~mean)[["residuals"]], by=c("id","anno")]

# Accessibility: differences in global accessibility rate
# foo <- fread("/Users/ricard/gastrulation/acc/stats/samples/out/sample_stats.txt") %>%
#   .[,mean:=log2(((mean/100)+0.01)/(1-(mean/100)+0.01))] %>%
#   .[,c("id_acc","mean")]
# acc_dt <- acc_dt %>% merge(foo, by="id_acc") %>%
#   .[,m:=lm(formula=m~mean)[["residuals"]], by=c("id","anno")]

#####################################
## Select highly variable features ##
#####################################

# RNA: Extract top N highly variable genes
keep_hv_genes <- rna_dt[,.(var=var(expr)), by="ens_id"] %>% .[var>0] %>% setorder(-var) %>% head(n=opts$rna_ngenes) %>% .$ens_id
rna_dt <- rna_dt[ens_id%in%as.character(keep_hv_genes)] %>% droplevels()

# Accessibility: Extract top N most variable features
keep_hv_sites <- acc_dt %>% split(.$anno) %>% map(~ .[,.(var=var(m)), by="id"] %>% .[var>0] %>% setorder(-var) %>% head(n=opts$acc_nfeatures) %>% .$id)
acc_dt <- acc_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist %>% droplevels()

# Methylation: Extract top N most variable features
keep_hv_sites <- met_dt %>% split(.$anno) %>% map(~ .[,.(var=var(m)), by="id"] %>% .[var>0] %>% setorder(-var) %>% head(n=opts$met_nfeatures) %>% .$id)
met_dt <- met_dt %>% split(.$anno) %>% map2(.,names(.), function(x,y) x[id %in% keep_hv_sites[[y]]]) %>% rbindlist %>% droplevels()

#############################
## Create joint data.frame ##
#############################

# data1 <- rna_dt %>% .[,c("sample","gene","expr")] %>%  
#   setnames(c("sample","feature","value")) %>% .[,c("feature_group","sample_group"):=list("RNA","gastrulation")]
# data2 <- met_dt %>% .[,c("sample","id","m","anno")] %>%  
#   setnames(c("sample","feature","value","feature_group")) %>% .[,c("feature","feature_group","sample_group"):=list(paste0("met_",feature), paste0("met_",feature_group), "gastrulation")]
# data3 <- acc_dt %>% .[,c("sample","id","m","anno")] %>%  
#   setnames(c("sample","feature","value","feature_group")) %>% .[,c("feature","feature_group","sample_group"):=list(paste0("acc_",feature), paste0("acc_",feature_group), "gastrulation")]
# 
# data <- rbind(data1,data2,data3)

# outfile <- paste0(io$outdir,"/data.txt")
# fwrite(data, file=outfile, col.names=T, quote=F, sep="\t")
# system(sprintf("pigz -f %s",outfile))


###########################
## Create input matrices ##
###########################

met_cells <- as.character(unique(met_dt$sample))
rna_cells <- as.character(unique(rna_dt$sample))
acc_cells <- as.character(unique(acc_dt$sample))

rna_matrix <- rna_dt[,c("gene","expr","sample")] %>%
  .[,c("sample","gene"):=list(as.character(sample),as.character(gene))] %>%
  .[,sample:=factor(sample,levels=Reduce(union,list(rna_cells,acc_cells,met_cells)))] %>%
  dcast(sample~gene, value.var="expr", drop=F) %>% matrix.please() %>% t

met_matrix_list <- list()
for (n in unique(met_dt$anno)) {
  met_matrix_list[[paste("met",n,sep="_")]] <- met_dt[anno==n,c("id","m","sample")] %>%
    .[,c("sample","id"):=list(as.character(sample),as.character(id))] %>%
    .[,sample:=factor(sample,levels=Reduce(union,list(rna_cells,acc_cells,met_cells)))] %>%
    dcast(sample~id, value.var="m", drop=F) %>% matrix.please() %>% t

  cat(sprintf("%s methylation matrix has dim (%d,%d) with %0.02f%% missing values \n", n,
              nrow(met_matrix_list[[paste("met",n,sep="_")]]), ncol(met_matrix_list[[paste("met",n,sep="_")]]),
              100*mean(is.na(met_matrix_list[[paste("met",n,sep="_")]]))))
}

cat("\n")

acc_matrix_list <- list()
for (n in unique(acc_dt$anno)) {
  acc_matrix_list[[paste("acc",n,sep="_")]] <- acc_dt[anno==n,c("id","m","sample")] %>%
    .[,c("sample","id"):=list(as.character(sample),as.character(id))] %>%
    .[,sample:=factor(sample,levels=Reduce(union,list(rna_cells,acc_cells,met_cells)))] %>%
    dcast(sample~id, value.var="m", drop=F) %>% matrix.please() %>% t

  cat(sprintf("%s accessibility matrix has dim (%d,%d) with %0.02f%% missing values \n", n,
              nrow(acc_matrix_list[[paste("acc",n,sep="_")]]), ncol(acc_matrix_list[[paste("acc",n,sep="_")]]),
              100*mean(is.na(acc_matrix_list[[paste("acc",n,sep="_")]]))))
}

# join everything filling with missing values
all_matrix_list <- c(rna=list(rna_matrix),met_matrix_list,acc_matrix_list)
