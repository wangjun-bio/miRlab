library(dplyr) 
library(openxlsx)
library(data.table)
library(DescTools)
library(tibble)
library(tidyr)
library(languageserver)
## PS: this scipit is final version for generating the expression matrix of krakenuniq 0806
setwd('D:/R_project/Microbe cfRNA/MicrobeRNA_R')
## read the krakenuniq data
CA_cohort = read.csv('../process_file/1.generate the exp matrix/PAN_microbeRNA_0731.csv')
RS_cohort = read.csv('../process_file/1.generate the exp matrix/PN_microbeRNA_0731.csv',check.names = F)
colnames(RS_cohort)[1] = 'X'
NC_meta = read.csv('../process_file/1.generate the exp matrix/NC_paper_meta.csv')
PN_meta = read.xlsx('../process_file/1.generate the exp matrix/TB_meta.xlsx')
df_all = merge(CA_cohort,RS_cohort,by = 'X',all = T)

### meta  ###
NC_meta = NC_meta %>%
  select(ID,sample_type)
NC_meta$group = lapply(NC_meta$sample_type,function(x) strsplit(x,'_')[[1]][[1]])
NC_meta$group = gsub(' ','',NC_meta$group)
NC_meta$num = lapply(NC_meta$ID,function(x) strsplit(x,'_')[[1]][[3]])
NC_meta$num = gsub(' ','',NC_meta$num)
NC_meta$num = as.numeric(NC_meta$num)
# Brain
Brain_meta = read.csv('../process_file/1.generate the exp matrix/Brain_usedSamples_meta.csv',row.names = 1)
Brain_meta$num = lapply(Brain_meta$id,function(x) strsplit(x,'_')[[1]][[3]]) %>% as.numeric()
colnames(Brain_meta) = c('group','sample_type','ID','num')

# PN
PN_meta = PN_meta %>%
  select(ID,Group)
PN_meta$group = lapply(PN_meta$Group,function(x) strsplit(x,'_')[[1]][[1]])
PN_meta$group = gsub(' ','',PN_meta$group)
colnames(PN_meta)[2] = 'sample_type'
table(PN_meta$group)

### Merge RS cohort ###
set.seed(222)
NOR_selected = NC_meta[NC_meta$group == 'NOR',]
NOR_selected$num = as.numeric(NOR_selected$num)
NOR_selected = NOR_selected[sample(nrow(NOR_selected),30),]

LC_select = NC_meta[NC_meta$group == 'LC',]
LC_select = LC_select[sample(nrow(LC_select),30),]
RS_meta = bind_rows(NOR_selected,LC_select,PN_meta[PN_meta$group == 'PN',])
RS_meta = RS_meta %>%
  group_by(group) %>%
  mutate(num = ifelse(is.na(num),
                      cumsum(is.na(num)),
                      num))
RS_meta$sample_label = paste0(RS_meta$sample_type,'_',RS_meta$num)
table(RS_meta$group)


# Final output
CA_sample = intersect(colnames(df_all),PAN_meta$ID)
CA_cohort = df_all %>%
  select(X,CA_sample)
CA_meta = PAN_meta[PAN_meta$ID %in% CA_sample,]
colnames(CA_cohort)[1] = 'classification_name'
CA_cohort[is.na(CA_cohort)] = 0
table(CA_meta$group)

RS_sample = intersect(colnames(df_all),RS_meta$sample_label)
RS_meta = RS_meta[RS_meta$sample_label %in% RS_sample,]
RS_cohort = df_all %>%
  select(X,RS_sample) %>% as.data.frame()
rownames(RS_cohort) = RS_cohort[,1];RS_cohort = RS_cohort[,-1]
RS_cohort = RS_cohort[,OrderMixed(colnames(RS_cohort))]
RS_meta = RS_meta[OrderMixed(RS_meta$ID),]
colnames(RS_cohort) = RS_meta$sample_label
RS_cohort = RS_cohort %>%
  rownames_to_column(var = 'classification_name')
RS_cohort[is.na(RS_cohort)] = 0
table(RS_meta$group)

saveRDS(CA_cohort,'../process_file/1.generate the exp matrix/CA_cohort_0801.rds')
saveRDS(RS_cohort,'../process_file/1.generate the exp matrix/RS_cohort_0801.rds')
saveRDS(CA_meta,'../process_file/1.generate the exp matrix/CA_meta_0801.rds')
saveRDS(RS_meta,'../process_file/1.generate the exp matrix/RS_meta_0801.rds')
write.csv(RS_meta,'../process_file/1.generate the exp matrix/RS_meta_0801.csv')
write.csv(CA_cohort,'../process_file/1.generate the exp matrix/CA_cohort_0801.csv')
write.csv(RS_cohort,'../process_file/1.generate the exp matrix/RS_cohort_0801.csv')


# summary the exceed genus
CA_cohort = fread('../process_file/1.generate the exp matrix/CA_cohort_0801.csv')  %>% as.data.frame()
RS_cohort = fread('../process_file/1.generate the exp matrix/RS_cohort_0801.csv')  %>% as.data.frame()
rownames(CA_cohort) = CA_cohort[,2];CA_cohort = CA_cohort[,-c(1,2)]
rownames(RS_cohort) = RS_cohort[,2];RS_cohort = RS_cohort[,-c(1,2)]

filter_taxa_by_abundance <- function(abund,
                                     abundance_thr = 0.15,
                                     prop_thr      = 0.05) {
  
  n_samples      <- ncol(abund)
  passed_counts  <- rowSums(abund > abundance_thr)
  passed_prop    <- passed_counts / n_samples
  
  keep_taxa      <- passed_prop > prop_thr
  filtered_abund <- abund[keep_taxa, , drop = FALSE]
  
  return(filtered_abund)
}
CA_cohort = filter_taxa_by_abundance(CA_cohort,abundance_thr = 0.15,prop_thr = 0.05)
RS_cohort = filter_taxa_by_abundance(RS_cohort,abundance_thr = 0.15,prop_thr = 0.05)

exceed = c(rownames(CA_cohort),rownames(RS_cohort)) %>% unique()
