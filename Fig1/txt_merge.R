setwd('F:\\R_script\\20230919 乳腺癌\\expr\\result_final')
library(stringr)
#a <- read.table('GDM20220313-III_AGCGTAGC_Homo-rsRNA_cal.txt',header = F,sep = '')
### 获取文件夹中的所有txt文件 ####
all_file <- list.files()
RNA <- all_file[grep('ysRNA',all_file)]############更换不同种类的cfRNA前缀

# 定义一个函数来读取txt文件
read_txt_file <- function(filename) {
  data <- read.table(filename, header = F, sep = '')  # 根据文件实际的分隔符来设置sep
  return(data)
}

# 使用lapply批量读取文件并将结果存储在一个列表中
library(tidyverse)
data_list <- lapply(RNA, read_txt_file)
###对miRNA\piRNA\rsRNA\mtt\ys表达矩阵进行处理
for (i in 1:length(data_list)) {
  colnames(data_list[[i]])[2] <- 'Count'
  colnames(data_list[[i]])[1] <- 'ID'
}
###对剩余RNA表达矩阵进行处理
for (i in 1:length(data_list)) {
     colnames(data_list[[i]])[1] <- "ID"
     colnames(data_list[[i]])[7] <- "Count"
  data_list[[i]] <- data_list[[i]][-1,]
  data_list[[i]] <- data_list[[i]] %>% select('ID','Count')
}
combine <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), data_list)
rownames(combine) <- combine$ID
combine <- combine[,-1]
library(readxl)
id <- read_xlsx('ID.xlsx')
ids <- sub('_$','',id$ID)
ids <- as.data.frame(cbind(ids,id$seq_ID))
colnames(ids)[1] <- 'ID'
colnames(ids)[2] <- 'seq-ID'
ID <- sub('_Homo-ysRNA_cal.txt','',RNA)
ids2 <- ID[match(ids$ID,ID)]
colnames(combine) <- ID
a <- combine[,match(ids$ID,ID)]
identical(ids2,ids$ID)
colnames(a) <- ids$`seq-ID`

write.csv(a,file = 'ysRNA.csv')


