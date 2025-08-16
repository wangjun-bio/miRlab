## PS: this script is used to generate cmRNA abundance matrix based on krakenuniq software.
# modifation : 2025-08-11

getwd()
setwd('/mnt/data3/yiyonghao/NC_paper_rawdata/')

library(dplyr)
library(data.table)
library(pheatmap)
library(tibble)
library(here)
library(httpgd)
library(languageserver)
# save the image
save.image(file = './Lenth_hisgram.RData')
load('./Lenth_hisgram.RData')

############## Reading the contamination data ###################
## ALL_H20 samples
files = list.files('./h20_all/krakenuniq_output/',pattern = '_reportfile.tsv',full.names = T)
result = lapply(files,function(x) {
  tmpname = basename(x)
  tmpname = gsub('_reportfile.tsv','',tmpname)
  tmp = read.table(x,header = T,sep = '\t')
  tmp = tmp[tmp$rank == 'genus',]
  tmp$taxName = gsub(' ','',tmp$taxName)
  tmp = tmp[,c(9,2)]
  colnames(tmp)[2] = tmpname
  return(tmp)
})
result = Reduce(function(x,y) merge(x,y,by = 'taxName',all = T),result) %>% as.data.frame()
rownames(result) = result[,1];result = result[,-1]
result = apply(result,2,function(row){
  row[is.na(row)] = 0
  return(row)
})
colnames(result) = gsub('_cut_R1','',colnames(result))
h20_all = result


## skin-derived microbes
library(openxlsx)
skin_meta = read.xlsx('./mmc2.xlsx',check.names = F)
colnames(skin_meta) = skin_meta[1,];skin_meta = skin_meta[-1,]
skin_meta = skin_meta[skin_meta$Timepoint == 'First',]

skin = read.xlsx('./mmc3.xlsx',check.names = F)
colnames(skin) = skin[1,];skin = skin[-1,]
skin = skin[6:nrow(skin),]
skin <- skin %>%
  mutate(genus = sapply(strsplit(Taxa, ";"), function(x) x[length(x) - 1]))
skin[,1] = skin$genus
skin = skin %>%
  select(-genus)
skin = as.data.frame(skin)
skin[,-1] = lapply(skin[,-1],as.numeric)
skin = skin %>%
  group_by(Taxa) %>%
  summarise(across(everything(),sum),.groups = 'drop')
skin = skin[-1,]
skin = as.data.frame(skin)
rownames(skin) = skin[,1];skin = skin[,-1]
skin = skin[,colnames(skin) %in% skin_meta$SAMPLE]

## cutoff is existed in at least 10 % samples
min_sample_count <- ceiling(0.1 * ncol(skin))
keep_rows <- rowSums(skin > 0) >= min_sample_count
skin <- skin[keep_rows, ]
filtered_skin = rownames(skin)
filtered_skin_df <- data.frame(
  skin_derived = rep(1, length(filtered_skin)),
  row.names       = filtered_skin
)


# virus host
virus = fread('./virushostdb.daily.tsv')
virus = virus %>%
  mutate(superkingdom = sapply(strsplit(`host lineage`, ";"), function(x) x[2]))
virus = virus[virus$superkingdom == ' Eukaryota',]
virus = virus[virus$`host name` != 'Homo sapiens',]
virus = virus %>%
  mutate(genus = sapply(strsplit(`virus lineage`, ";"), function(x) x[length(x) - 1]))
virus = virus[!grepl('unclassified',virus$`virus lineage`),]
virus$genus = gsub(' ','',virus$genus)
virus = virus$genus[!duplicated(virus$genus)]
virus <- data.frame(
  virus = rep(1, length(virus)),
  row.names       = virus
)

# remove reported laboratory contaminations from "Reagent and laboratory contamination can critically impact sequence-based microbiome analyses" table 1
contaminations = c('Afipia', 'Aquabacterium', 'Asticcacaulis', 'Aurantimonas', 'Beijerinckia', 'Bosea',
                  'Bradyrhizobium', 'Brevundimonas', 'Caulobacter', 'Craurococcus', 'Devosia', 'Hoeflea',
                  'Mesorhizobium', 'Methylobacterium', 'Novosphingobium', 'Ochrobactrum', 'Paracoccus',
                  'Pedomicrobium', 'Phyllobacterium', 'Rhizobium', 'Roseomonas', 'Sphingobium', 'Sphingomonas','Sphingopyxis',
                  'Acidovorax', 'Azoarcus', 'Azospira', 'Burkholderia', 'Comamonas',
                  'Cupriavidus', 'Curvibacter', 'Delftia', 'Duganella', 'Herbaspirillum', 'Janthinobacterium', 'Kingella',
                  'Leptothrix', 'Limnobacter', 'Massilia', 'Methylophilus', 'Methyloversatilis', 'Oxalobacter', 'Pelomonas',
                  'Polaromonas', 'Ralstonia','Schlegelella', 'Sulfuritalea', 'Undibacterium', 'Variovorax',
                  'Acinetobacter','Enhydrobacter', 'Enterobacter', 'Escherichia' ,'Nevskia', 'Pseudomonas', 'Pseudoxanthomonas', 'Psychrobacter',
                  'Stenotrophomonas','Xanthomonas',
                  'Aeromicrobium', 'Arthrobacter', 'Beutenbergia', 'Brevibacterium', 'Corynebacterium', 'Curtobacterium', 'Dietzia',
                  'Geodermatophilus', 'Janibacter', 'Kocuria', 'Microbacterium', 'Micrococcus', 'Microlunatus', 'Patulibacter', 'Propionibacterium',
                  'Rhodococcus', 'Tsukamurella',
                  'Abiotrophia', 'Bacillus', 'Brevibacillus', 'Brochothrix', 'Facklamia', 'Paenibacillus', 'Streptococcus',
                  'Chryseobacterium', 'Dyadobacter', 'Flavobacterium', 'Hydrotalea', 'Niastella', 'Olivibacter', 'Pedobacter', 'Wautersiella',
                  'Deinococcus','Homo')

all_contamination = c(contaminations,filtered_skin,virus)
all_contamination = all_contamination[!duplicated(all_contamination)]

write.table(all_contamination,'Microbe_contamination.txt',col.names = F,row.names = F,quote = F,sep = '\t')


######## Reading the reportfile.tsv files ##########
files = list.files('./upload/krakenuniq_output',pattern = '_reportfile.tsv',full.names = T)
# Adding the brain samples
brain = list.files('./Brain/krakenuniq_output/',pattern = '_reportfile.tsv',full.names = T)
files = c(files,brain)
result = lapply(files,function(x) {
  tmpname = basename(x)
  tmpname = gsub('_reportfile.tsv','',tmpname)
  tmp = read.table(x,header = T,sep = '\t')
  tmp = tmp[tmp$rank == 'genus',]
  tmp$taxName = gsub(' ','',tmp$taxName)
  tmp = tmp[!(tmp$taxName %in% all_contamination), ]
  tmp = tmp[,c(9,2)]
  colnames(tmp)[2] = tmpname
  return(tmp)
})

result = Reduce(function(x,y) merge(x,y,by = 'taxName',all = T),result) %>% as.data.frame()
rownames(result) = result[,1];result = result[,-1]
result = apply(result,2,function(row){
  row[is.na(row)] = 0
  return(row)
})
colnames(result) = gsub('_cut_R1','',colnames(result))
result = as.data.frame(result)

## Merge the brain sample
Brain_meta = read.xlsx('./Brain/BrainTumor_Summary_good_samples&NewID.xlsx')
merge_two_samples <- function(df, sample_string, new_name = NULL) {
  sample_pair <- unlist(strsplit(sample_string, ","))
  if (length(sample_pair) != 2) {
    stop("样本字符串拆分后不是两个样本，请确保格式为 'sample1,sample2'")
  }
  sample1 <- sample_pair[1]
  sample2 <- sample_pair[2]
  if (!all(sample_pair %in% colnames(df))) {
    stop("某些样本ID不在数据框中")
  }
  if (is.null(new_name)) {
    new_name <- paste(sample_pair[1])
  }
  merged_col <- as.numeric(rowSums(df[, sample_pair, drop = FALSE], na.rm = TRUE))
  df_new <- df[, setdiff(colnames(df), sample_pair), drop = FALSE]
  df_new[[new_name]] <- merged_col
  rownames(df_new) <- rownames(df)
  return(as.data.frame(df_new))
}
sample = Brain_meta[sapply(strsplit(as.character(Brain_meta$LibraryID), ","), length) > 1, ]
sample = sample$LibraryID
for (i in 1:length(sample)) {
  result = merge_two_samples(df = result,sample_string = sample[i])
}

Brain_meta_filter = Brain_meta[grepl('plasma',Brain_meta$Group),] %>%
  select(Group,LibraryID)
Brain_meta_filter = Brain_meta_filter[!grepl('NOR',Brain_meta_filter$Group),]
filter_gdm_samples <- function(df, keep_list) {
  col_all <- colnames(df)
  is_gdm <- grepl("^GDM\\d{8}-[A-Z0-9]+_[ATCG]+$", col_all)
  cols_keep <- col_all[is_gdm & col_all %in% keep_list]
  cols_unaffected <- col_all[!is_gdm]
  df_filtered <- df[, c(cols_unaffected, cols_keep), drop = FALSE]
  return(df_filtered)
}
result_filter = filter_gdm_samples(df = result,keep_list = Brain_meta_filter$LibraryID)

Brain_meta_filter <- Brain_meta_filter %>%
  group_by(Group) %>%
  mutate(id = paste0(Group, "_", row_number())) %>%
  ungroup()
index = match(colnames(result_filter),Brain_meta_filter$LibraryID)
for (i in 1:length(index)) {
  tmp = index[i]
  if (!is.na(tmp)) {
    colnames(result_filter)[i] = Brain_meta_filter$id[tmp]
    
  }else{
    colnames(result_filter)[i] = colnames(result_filter)[i]
  }
  
}
PAN = result_filter
write.csv(Brain_meta_filter,'./Brain/Brain_usedSamples_meta.csv')



## reading the PN data 
files = list.files('./PN/krakenuniq_output/',pattern = '_reportfile.tsv',full.names = T)
result = lapply(files,function(x) {
  tmpname = basename(x)
  tmpname = gsub('_reportfile.tsv','',tmpname)
  tmp = read.table(x,header = T,sep = '\t')
  tmp = tmp[tmp$rank == 'genus',]
  tmp$taxName = gsub(' ','',tmp$taxName)
  tmp = tmp[!(tmp$taxName %in% all_contamination), ]
  tmp = tmp[,c(9,2)]
  colnames(tmp)[2] = tmpname
  return(tmp)
})
result = Reduce(function(x,y) merge(x,y,by = 'taxName',all = T),result) %>% as.data.frame()

rownames(result) = result[,1];result = result[,-1]
result = apply(result,2,function(row){
  row[is.na(row)] = 0
  return(row)
})
result = as.data.frame(result)
PN = result
PN <- setNames(PN, paste0("PN_SZDE_", seq_len(ncol(PN))))


## Merging all samples to decontamination 
## based on decontam packages
library(decontam)
library(phyloseq)
ALL = list(PN,PAN,h20_all)

ALL = lapply(ALL,function(x){
  x = as.data.frame(x)
  x = x  %>% rownames_to_column()
  return(x)
})  %>% Reduce(function(x,y) merge(x,y,by = 'rowname',all = T),.)
ALL[is.na(ALL)] = 0
rownames(ALL) = ALL[,1];ALL = ALL[,-1]
sample_data = data.frame(sample = colnames(ALL))
sample_data$is.neg = lapply(sample_data$sample,function(x) strsplit(x,'_')[[1]][[1]]) %>% as.character()
sample_data$is.neg[grepl('GDM',sample_data$is.neg)] = 'WaterControl'
rownames(sample_data) = sample_data[,1];sample_data = sample_data[,-1,drop = F]
sample_data$is.neg = ifelse(sample_data$is.neg == 'WaterControl',TRUE,FALSE)

#phyloseq
OTU <- otu_table(as.matrix(t(ALL)), taxa_are_rows=FALSE)
SAM <- sample_data(sample_data)
ps  <- phyloseq(OTU, SAM)

contam.df <- decontam::isContaminant(
  ps,
  method    = "prevalence",
  neg       = "is.neg"
)
table(contam.df$contaminant)
contam.df$contaminant
decontam = rownames(contam.df)[contam.df$contaminant == TRUE]

## using wilcox to decontamination
samples = rownames(sample_data)[!grepl('GDM',rownames(sample_data))]
h20_samples = rownames(sample_data)[grepl('GDM',rownames(sample_data))]

ALL_ratio = apply(ALL,2,function(x) {
  sums = sum(x)
  x = x / sums
  return(x)
})

ALL_wilcox = lapply(rownames(ALL_ratio),function(x) {
  tmp = ALL_ratio[x,]
  tmp = data.frame(sample = names(tmp),ratio = tmp)
  tmp$group = ifelse(grepl('GDM',tmp$sample),'WaterControl','samples')
  tmp$group = factor(tmp$group,levels = c('WaterControl','samples'))
  wilcox_result = wilcox.test(ratio ~ group,data = tmp)
  y = data.frame(taxa = x,p.value = wilcox_result$p.value)
  y$padj = p.adjust(y$p.value,method = 'BH')
  y$significant = ifelse(y$padj < 0.01,'Significant','Not Significant')
  y$mean_ratio = mean(tmp$ratio[tmp$group == 'samples'])
  y$mean_ratio_h20 = mean(tmp$ratio[tmp$group == 'WaterControl'])
  y$log2FoldChange = log2(y$mean_ratio / y$mean_ratio_h20)
  y$increase = ifelse(y$log2FoldChange > 0,'Increase','Decrease')
  return(y)
})
ALL_wilcox = do.call(bind_rows,ALL_wilcox)
ALL_wilcox_intersect = ALL_wilcox[ALL_wilcox$taxa %in% decontam,]
ALL_wilcox = ALL_wilcox[ALL_wilcox$significant == 'Significant' & ALL_wilcox$increase == 'Decrease',]
ALL_wilcox = ALL_wilcox[!is.na(ALL_wilcox$increase),]

#union of the decontam and wilcox
negative_contamination = c(ALL_wilcox$taxa,decontam) %>% unique()
negative_contamination_df <- data.frame(
  Negative_sample = rep(1, length(negative_contamination)),
  row.names       = negative_contamination
)


## saving the contamination matrix
contamination <- data.frame(
  lab_contamination = rep(1, length(contaminations)),
  row.names       = contaminations
)
contamination_list = list(negative_contamination_df,filtered_skin_df,virus,contamination)
for (i in 1:length(contamination_list)) {
  contamination_list[[i]] = contamination_list[[i]] %>% rownames_to_column()
  
}
contamination_list = Reduce(function(x,y) merge(x,y,by = 'rowname',all = T),contamination_list)
contamination_list[is.na(contamination_list)] = 0
write.table(contamination_list,'Microbe_contamination_df_0731.txt',col.names = T,row.names = F,quote = F,sep = '\t')

## filter microbe derived from negative controls
PN_filter = PN[!(rownames(PN) %in% negative_contamination),]
PAN_filter = PAN[!(rownames(PAN) %in% negative_contamination),]


write.csv(PN_filter,'./PN_microbeRNA_0731.csv')
write.csv(PAN_filter,'./PAN_microbeRNA_0731.csv')















