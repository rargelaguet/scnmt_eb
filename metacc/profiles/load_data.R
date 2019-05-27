
# Load methylation data
met_list <- list()
for (cell in opts$met.cells) {
  tmp <- fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$met.dir,cell), showProgress=F, stringsAsFactors=F)
    
  if (ncol(tmp)==3) {
    tmp <- tmp %>% setnames(c("chr","pos","rate"))
  } else {
    tmp <- tmp %>% .[,c("chr","pos","rate")] %>% .[,rate:=rate*100]
  }
  tmp <- tmp  %>% .[,id_met:=cell] %>% 
    .[,c("start","end"):=list(pos,pos)] %>% setnames("pos","bp") %>% setkey("chr","start","end") %>%
    
    foverlaps(.,anno_df.met, nomatch=0) %>% .[, c("chr","i.start","i.end") := NULL] %>%
    .[,id:=as.character(id)] %>%
    .[,dist:=ifelse(strand %in% c("+","*"),bp-center,center-bp)] %>% 
    .[, dist:=opts$met.tile*round(dist/opts$met.tile)] %>%
    .[,list(rate=mean(rate), n=.N),by=.(id_met,id,dist,anno)]
  met_list[[cell]] <- tmp
}
met <- rbindlist(met_list) %>%
  .[,c("id_met","id","context"):=list(as.factor(id_met),as.factor(id),"CG")]
  
rm(met_list)

# Load accessibility data
acc_list <- list()
for (cell in opts$acc.cells) {
  tmp <- fread(cmd=sprintf("zcat < %s/%s.tsv.gz",io$acc.dir,cell), showProgress=F, stringsAsFactors=F) 
  
  if (ncol(tmp)==3) {
    tmp <- tmp %>% setnames(c("chr","pos","rate"))
  } else {
    tmp <- tmp %>% .[,c("chr","pos","rate")] %>% .[,rate:=rate*100]
  }
  
  tmp <- tmp %>% .[,id_acc:=cell] %>% 
    .[,c("start","end"):=list(pos,pos)] %>% setnames("pos","bp") %>% setkey("chr","start","end") %>%
    
    foverlaps(.,anno_df.acc, nomatch=0) %>% .[, c("chr","i.start","i.end") := NULL] %>%
    .[,id:=as.character(id)] %>%
    .[,dist:=ifelse(strand %in% c("+","*"),bp-center,center-bp)] %>% 
    .[, dist:=opts$acc.tile*round(dist/opts$acc.tile)] %>%
    .[,list(rate=mean(rate), n=.N),by=.(id_acc,id,dist,anno)]
  acc_list[[cell]] <- tmp
}
acc <- rbindlist(acc_list) %>%
  .[,c("id_acc","id","context"):=list(as.factor(id_acc),as.factor(id),"GC")]
  
rm(acc_list)

