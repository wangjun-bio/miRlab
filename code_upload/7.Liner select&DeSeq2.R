# install
#BiocManager::install('taxtree',dependencies = TRUE)
# library
library(tidyverse)
library(limma)
library(DESeq2)
library(magrittr)
library(dplyr)
library(reshape2)
library(scatterplot3d)
library(Maaslin2)
library(data.table)
library(DescTools)
#library(taxtree)
save.image('../image/2.Multi train.Rdata')
load('../image/2.Multi train.Rdata')
setwd('/mnt/data3/yiyonghao/MicroRNA')

### read the raw data ###
CA_cohort = readRDS('/mnt/data3/yiyonghao/MicroRNA/process_file/0819/CA_cohort_0820.rds')
CA_cohort = as.data.frame(CA_cohort)
rownames(CA_cohort) = CA_cohort[,1];CA_cohort = CA_cohort[,-1]
RS_cohort = readRDS('/mnt/data3/yiyonghao/MicroRNA/process_file/0819/RS_cohort_0820.rds')
RS_cohort = as.data.frame(RS_cohort)
rownames(RS_cohort) = RS_cohort[,1];RS_cohort = RS_cohort[,-1]


###### Maaslin2 Liner select ######
CA_meta = data.frame(row.names = colnames(CA_cohort))
CA_meta$Group = lapply(rownames(CA_meta),function(x) strsplit(x,'_')[[1]][[1]]) %>% as.character()
table(CA_meta$Group)
# CA cohort
dir.create('/mnt/data3/yiyonghao/MicroRNA/process_file/0819/CA_masslin2',showWarnings = F,recursive = T)
res=Maaslin2(
  input_data = CA_cohort, 
  input_metadata = CA_meta, 
  output = "/mnt/data3/yiyonghao/MicroRNA/process_file/0819/CA_masslin2/Maaslin2",
  analysis_method = "LM",
  correction = "BH",
  normalization = "TMM",
  plot_heatmap = T,
  plot_scatter = F,
  heatmap_first_n = 50,
  min_prevalence=0.1,
  max_significance=0.01,
  fixed_effects = c("Group"),
  reference = c("Group,NOR"))

res = res$results
feature_CA <- res[res$qval<0.01&res$pval<0.001,]$feature
unique(feature_CA)[!(unique(feature_CA) %in% rownames(CA_cohort))]
CA_cohort_feature = CA_cohort[unique(feature_CA),] %>% t() %>% 
  data.frame() %>% 
  mutate(group = limma::strsplit2(rownames(.),"_")[,1])
CA_cohort_feature = CA_cohort_feature[,!grepl('NA',colnames(CA_cohort_feature))]
saveRDS(CA_cohort_feature,'/mnt/data3/yiyonghao/MicroRNA/process_file/0819/CA_selected_Masslin2.rds')
write.csv(CA_cohort_feature,'/mnt/data3/yiyonghao/MicroRNA/process_file/0819/CA_selected_Masslin2.csv')
table(CA_cohort_feature$group)


# RS cohort
dir.create('/mnt/data3/yiyonghao/MicroRNA/process_file/0819/RS_masslin2',showWarnings = F,recursive = T)
RS_meta = data.frame(row.names = colnames(RS_cohort))
RS_meta$Group = lapply(rownames(RS_meta),function(x) strsplit(x,'_')[[1]][[1]]) %>% as.character()
table(RS_meta$Group)
res=Maaslin2(
  input_data = RS_cohort, 
  input_metadata = RS_meta, 
  output = "/mnt/data3/yiyonghao/MicroRNA/process_file/0819/RS_masslin2/Maaslin2",
  analysis_method = "LM",
  correction = "BH",
  normalization = "TMM",
  plot_heatmap = T,
  plot_scatter = F,
  heatmap_first_n = 50,
  min_prevalence=0.1,
  max_significance=0.01,
  fixed_effects = c("Group"),
  reference = c("Group,NOR"))
res = res$results

feature_RS <- res[res$qval<0.1,]$feature
unique(feature_RS)[!(unique(feature_RS) %in% rownames(RS_cohort))]
RS_cohort_feature = RS_cohort[unique(feature_RS),] %>% t() %>% 
  data.frame() %>% 
  mutate(group = limma::strsplit2(rownames(.),"_")[,1])
RS_cohort_feature = RS_cohort_feature[,!grepl('NA',colnames(RS_cohort_feature))]
RS_cohort_feature = RS_cohort_feature[OrderMixed(rownames(RS_cohort_feature)),]

saveRDS(RS_cohort_feature,'/mnt/data3/yiyonghao/MicroRNA/process_file/0819/RS_selected_Masslin2.rds')
write.csv(RS_cohort_feature,'/mnt/data3/yiyonghao/MicroRNA/process_file/0819/RS_selected_Masslin2.csv')


###### deseq2 ######
deseq2 = function(count,group,loopname,control,normalized) {
  #browser()
  library(DESeq2)
  library(DescTools)
  count = count[,colnames(count) %in% rownames(group)]
  group = group[rownames(group) %in% colnames(count),,drop = F]
  
  
  count = count[,OrderMixed(colnames(count))]
  sums = apply(count,1,sum)
  sums = sums[sums != 0]
  count = count[rownames(count) %in% names(sums),]
  
  
  group = group[OrderMixed(rownames(group)),,drop = F]
  dir.create('dds_result',showWarnings = F)
  
  colnames(group) <- "condition"
  group$condition <- factor(group$condition, levels = c(control, setdiff(unique(group$condition), control)))
  
  
  all_groups <- levels(group$condition)
  exp_groups <- setdiff(all_groups, control)
  
  done_pairs <- character(0)

  dds = DESeqDataSetFromMatrix(countData = count,
                               colData = group,
                               design = ~condition,
                               )
  dds = estimateSizeFactors(dds, type = "poscounts") 
  dds = DESeq(dds,quiet = T)
  if (normalized == T) {
    
    normalized_counts <- counts(dds, normalized = TRUE)
    dir.create('./dds_result/normalized',showWarnings = F)
    write.csv(normalized_counts,paste0(getwd(),'/dds_result/normalized/',loopname,'_normalized_counts.csv'))
  }
  
  for (exp in exp_groups) {
    pair_id <- paste(sort(c(exp, control)), collapse = "_vs_")
    if (!(pair_id %in% done_pairs)){
      res <- results(dds, contrast = c("condition", exp, control))
      resorder <- res[order(res$padj), ]
      dataDEG <- as.data.frame(resorder)
      
      write.csv(dataDEG,
                paste0(getwd(), '/dds_result/', exp, '_vs_', control, '_', loopname, '.csv'))
      print(paste0('The ', exp, '_vs_', control, ' is complete'))
      done_pairs <- c(done_pairs, pair_id)
    } 
    else{
      print('The combn is exist')
    }
    
    

  }
  
}

setwd('/mnt/data3/yiyonghao/MicroRNA/process_file/0819')
CA_meta_deseq = CA_meta
CA_meta_deseq$Group = ifelse(CA_meta_deseq$Group != 'NOR','PAN','NOR')
colnames(CA_meta_deseq) = 'group'
deseq2(count = CA_cohort,group = CA_meta_deseq,loopname = '',control = 'NOR',normalized = F)

RS_meta_deseq = data.frame(row.names = colnames(RS_cohort))
RS_meta_deseq$group = lapply(rownames(RS_meta_deseq),function(x) strsplit(x,'_')[[1]][[1]]) %>% as.character()
RS_meta_deseq$group = ifelse(RS_meta_deseq$group != 'NOR','Other','NOR')
deseq2(count = RS_cohort,group = RS_meta_deseq,loopname = '',control = 'NOR',normalized = F)




