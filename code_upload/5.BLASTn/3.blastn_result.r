# reading the result of Blastn
library(languageserver)
library(dplyr)
library(data.table)
getwd()
setwd('/mnt/data3/yiyonghao/NC_paper_rawdata')

blastn = fread('/mnt/data3/yiyonghao/NC_paper_rawdata/upload/blastn/results/merged_random100_blastn.out')
head(blastn)
min(blastn$pident)
PN_blastn = fread('/mnt/data3/yiyonghao/NC_paper_rawdata/PN/blastn/results/merged_random100_blastn.out')
brain_blastn = fread('/mnt/data3/yiyonghao/NC_paper_rawdata/Brain/blastn/results/merged_random100_blastn.out')
ALL = bind_rows(blastn,PN_blastn,brain_blastn)
colnames(ALL) = c(
  "qseqid",    # Query sequence ID
  "sseqid",    # Subject (database) sequence ID
  "pident",    # Percentage of identical matches
  "length",    # Alignment length
  "mismatch",  # Number of mismatches
  "gapopen",   # Number of gap openings
  "qstart",    # Start of alignment in query
  "qend",      # End of alignment in query
  "sstart",    # Start of alignment in subject
  "send",      # End of alignment in subject
  "evalue",    # Expect value
  "bitscore"   # Bit score
)
ALL = ALL[,c('qseqid','sseqid','pident','evalue','bitscore')]
ALL_filter = ALL %>% 
    group_by(qseqid)  %>% 
    slice_min(order_by = evalue,n = 1,with_ties = F) %>% 
    ungroup()  %>% 
    as.data.frame()
species_id = unique(ALL$sseqid)
#write.table(species_id,file = '/mnt/data3/yiyonghao/NC_paper_rawdata/rRNA_region/species_id.txt',col.names = F,row.names = F,quote = F,sep = '\t')

## read species_map
species_map = fread('/mnt/data3/yiyonghao/NC_paper_rawdata/rRNA_region/accession2definition.tsv',header = F)
species_map$tax = lapply(species_map$V1,function(x) strsplit(x,'\\t',fixed = T)[[1]][[2]]) %>% as.character()
species_map$V1 = lapply(species_map$V1,function(x) strsplit(x,'\\t',fixed = T)[[1]][[1]]) %>% as.character()
head(species_map)

## mapping tax id to accessionid
index = match(ALL_filter$sseqid,species_map$V1)
ALL_filter$sseqid = species_map$tax[index]
ALL_filter = ALL_filter[,1:2]
colnames(ALL_filter)[2] = 'blastn_taxid'

## read the sequences info
paths = c('/mnt/data3/yiyonghao/NC_paper_rawdata/upload','/mnt/data3/yiyonghao/NC_paper_rawdata/PN','/mnt/data3/yiyonghao/NC_paper_rawdata/Brain')
infos = lapply(paths,function(x){
    x = list.files(paste0(x,'/blastn/info'),pattern = '.tsv',full.names = T)
    tmp = lapply(x,function(y){
        tmpname = lapply(y,function(x) strsplit(x,'\\/')[[1]][[length(strsplit(x,'\\/')[[1]])]])
        n = fread(y,header = F)
        n$sample = tmpname
        return(n)
    })
    tmp = do.call(bind_rows,tmp)
    return(tmp)
})
infos = do.call(bind_rows,infos)
colnames(infos) = c("classificied",'Sequence_id','TaxID','length','k-mer','sample')
infos = infos[,c(2:3,6)]
head(infos)
infos$sample = gsub('_random100_C.tsv','',infos$sample)


## merge the krakenuniq and blastn result
Merge_df = merge(infos,ALL_filter,by.x = 'Sequence_id',by.y = 'qseqid',all = T)

# trans taxid to genus
# write.table(Merge_df$TaxID,'/mnt/data3/yiyonghao/NC_paper_rawdata/rRNA_region/taxid.txt',col.names = F,row.names = F,quote = F,sep = '\t')
# install.packages('taxize')
library(taxize)
library(taxizedb)
library(tidyverse)
taxid = c(Merge_df$TaxID,Merge_df$blastn_taxid) %>% unique()

files = list.files('/mnt/data3/yiyonghao/NC_paper_rawdata/upload/krakenuniq_output',pattern = '_reportfile.tsv',full.names = T)
# db_download_ncbi() 
db <- db_load_ncbi()
res <- classification(taxid, db = "ncbi")
res <- res[!sapply(res, function(x) {
  length(x) == 1 && is.na(x)
})]

for(i in 1:length(res)){
    tmp = res[[i]][res[[i]]$rank == 'genus',]
    tmp = tmp$id
    res[[i]] = tmp
  }
res = unlist(res)
res = res %>% as.data.frame()
res = res %>% rownames_to_column()
colnames(res) = c('taxid','genus_taxid')

index = match(Merge_df$TaxID,res$taxid)
Merge_df$TaxID = res$genus_taxid[index]
index = match(Merge_df$blastn_taxid,res$taxid)
Merge_df$blastn_taxid = res$genus_taxid[index]

Merge_df = na.omit(Merge_df)
Merge_df$check = ifelse(Merge_df$TaxID == Merge_df$blastn_taxid,'yes','no')


Merge_df_nohomo = Merge_df[Merge_df$TaxID != '9605',]
table(Merge_df_nohomo$check)
write.csv(Merge_df_nohomo,'/mnt/data3/yiyonghao/MicroRNA/process_file/blastn_krakenuniq_compare_0827.csv',row.names = F,quote = F)

