################### RNA_LEN_Proportion ##################
#########rsRNA长度占比计算##############
setwd('F:\\R_script\\20230919 乳腺癌\\expr\\RS_len')
library(stringr)

all_file <- list.files()
rsRNA_len <- all_file[grepl('_sort.txt',all_file)]
read_txt_file <- function(filename) {
  data <- read.table(filename, header = FALSE, sep = '')  # 根据文件实际的分隔符来设置sep
  data$filename <- basename(filename)  # 添加文件名列
  return(data)
}
# 使用lapply批量读取文件并将结果存储在一个列表中
library(tidyverse)
data_list <- lapply(rsRNA_len, read_txt_file)
for (i in 1:length(data_list)) {
  colnames(data_list[[i]])[3] <- 'ID'
  colnames(data_list[[i]])[2] <- 'len'
  colnames(data_list[[i]])[1] <- 'reads'
}
merge_data=do.call(rbind,data_list)
merge_data$ID=gsub('_len.txt_sort.txt','',merge_data$ID)
# 使用pivot_wider将长数据形式转换为宽数据形式
merge_matrix <- merge_data %>%
  pivot_wider(names_from = ID, values_from = reads)
merge_matrix <- as.data.frame(merge_matrix)
rownames(merge_matrix) <- merge_matrix$len
merge_matrix <- merge_matrix[,-1]
merge_matrix[is.na(merge_matrix)] <- 0


rs18S_len = merge_matrix[,grepl('-18S',colnames(merge_matrix))]
colnames(rs18S_len) <- sub('_[^_]*$','',colnames(rs18S_len))
rs28S_len = merge_matrix[,grepl('-28S',colnames(merge_matrix))]
colnames(rs28S_len) <- sub('_[^_]*$','',colnames(rs28S_len))
rs5.8S_len = merge_matrix[,grepl('-5.8S',colnames(merge_matrix))]
colnames(rs5.8S_len) <- sub('_[^_]*$','',colnames(rs5.8S_len))
rs5S_len = merge_matrix[,grepl('-5S',colnames(merge_matrix))]
colnames(rs5S_len) <- sub('_[^_]*$','',colnames(rs5S_len))
rsmt12S_len = merge_matrix[,grepl('-mt12S',colnames(merge_matrix))]
colnames(rsmt12S_len) <- sub('_[^_]*$','',colnames(rsmt12S_len))
rsmt16S_len = merge_matrix[,grepl('-mt16S',colnames(merge_matrix))]
colnames(rsmt16S_len) <- sub('_[^_]*$','',colnames(rsmt16S_len))

merge_matrix <- rs18S_len + rs28S_len + rs5.8S_len + rs5S_len + rsmt12S_len + rsmt16S_len



library(readxl)
id <- read_xlsx('F:\\R_script\\20230919 乳腺癌\\result_final\\ID.xlsx')
ids <- sub('_$','',id$ID)
ids <- as.data.frame(cbind(ids,id$seq_ID))
colnames(ids)[1] <- 'ID'
colnames(ids)[2] <- 'seq-ID'
all <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\all_seq.csv')
a <- as.data.frame(all[match(ids$ID,all$X),])
rownames(a) <- ids$`seq-ID`
group <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\simple.csv')
group <- group[-c(82,83),]
a=as.data.frame(a[match(group$X,rownames(a)),])
a$group=group$condition
a=data.frame(a$X,a$ysRNA,a$group)
rownames(a) <- a$a.X
a=a[,-1]
colnames(a) <- c('count','group')

ben <- a %>% filter(a$group == 'BEN')
mal <- a %>% filter(a$group == 'MAL')


merge_matrix_ben <- merge_matrix[,match(rownames(ben),colnames(merge_matrix))]
total_reads_ben <- colSums(merge_matrix_ben)
result_ben <- data.frame(matrix(NA, nrow = nrow(merge_matrix_ben), ncol = ncol(merge_matrix_ben)))
colnames(result_ben) <- colnames(merge_matrix_ben)
rownames(result_ben) <- rownames(merge_matrix_ben)
for (i in 1:length(colnames(merge_matrix_ben))) {
  result_ben[,i] <- as.data.frame((merge_matrix_ben[,i] / total_reads_ben[i]) * 100)
}
result_ben$type <- rownames(result_ben)


merge_matrix_mal <- merge_matrix[,match(rownames(mal),colnames(merge_matrix))]
total_reads_mal <- colSums(merge_matrix_mal)
result_mal <- data.frame(matrix(NA, nrow = nrow(merge_matrix_mal), ncol = ncol(merge_matrix_mal)))
colnames(result_mal) <- colnames(merge_matrix_mal)
rownames(result_mal) <- rownames(merge_matrix_mal)
for (i in 1:length(colnames(merge_matrix_mal))) {
  result_mal[,i] <- as.data.frame((merge_matrix_mal[,i] / total_reads_mal[i]) * 100)
}
result_mal$type <- rownames(result_mal)
library(ggplot2)
library(dplyr)
library(reshape)

result_long <- melt(result_mal)
all_types <- unique(result_long$type)
# 指定要展示的ID列表
desired_IDs <- as.character(seq(15, 150))
# 将ID列转换为因子，并只包含desired_IDs
result_long$type <- factor(result_long$type, levels = desired_IDs)
#result_long <- result_long[order(result_long$type),]
# 选择15到42之间的数据子集
num_labels <- 8
step <- ceiling(length(all_types) / num_labels)
custom_order <- factor(all_types[seq(from = 1, to = length(all_types), by = step)])
custom_order <- c("15", "29",'43','57','71','86','100','114','128','142')

# 将type变量转换为因子，并且指定自定义的顺序
result_long$type <- factor(result_long$type, levels = custom_order)
result_long <- na.omit(result_long)
result_subset <- subset(result_long, type >= 15 & type <= 42)
result_subset <- subset(result_subset,type != 150)



top.mar=0.2
right.mar=0.2
bottm.mar=0.2
left.mar=0.2
mytheme <- theme_classic()+
  theme(text = element_text(family = 'sans',colour = 'black',size = 14),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.6,colour = 'black'),
        axis.ticks.length = unit(1.5,units = 'mm'),
        plot.margin = unit(x=c(top.mar,right.mar,left.mar,bottm.mar),units = 'inches'))

p <- ggplot(result_long, aes(x = type, y = value,group = variable)) +
  geom_line(color = 'grey') + geom_point(color = '#87CEEB') +
  xlab("length(nt)") +
  ylab("Proportion") +
  ggtitle("rsRNAs")+
  mytheme+scale_x_discrete(breaks = custom_order,labels = custom_order)+
  theme(plot.title = element_text(hjust = 0.5, size = 14))+ annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "transparent", color = "black")
p
ggsave(p,filename = 'MAL_rsRNA_len占比图.pdf',width = 4,height = 4)


##############ysRNA长度占比计算################
setwd('F:\\R_script\\20230919 乳腺癌\\expr\\YS_len')
library(stringr)

all_file <- list.files()
ysRNA_len <- all_file[grepl('_sort.txt',all_file)]
read_txt_file <- function(filename) {
  data <- read.table(filename, header = FALSE, sep = '')  # 根据文件实际的分隔符来设置sep
  data$filename <- basename(filename)  # 添加文件名列
  return(data)
}
# 使用lapply批量读取文件并将结果存储在一个列表中
library(tidyverse)
data_list <- lapply(ysRNA_len, read_txt_file)
for (i in 1:length(data_list)) {
  colnames(data_list[[i]])[3] <- 'ID'
  colnames(data_list[[i]])[2] <- 'len'
  colnames(data_list[[i]])[1] <- 'reads'
}
merge_data=do.call(rbind,data_list)
merge_data$ID=gsub('_len.txt_sort.txt','',merge_data$ID)
# 使用pivot_wider将长数据形式转换为宽数据形式
merge_matrix <- merge_data %>%
  pivot_wider(names_from = ID, values_from = reads)
merge_matrix <- as.data.frame(merge_matrix)
rownames(merge_matrix) <- merge_matrix$len
merge_matrix <- merge_matrix[,-1]
merge_matrix[is.na(merge_matrix)] <- 0


RNY1_len = merge_matrix[,grepl('RNY1',colnames(merge_matrix))]
colnames(RNY1_len) <- sub('_[^_]*$','',colnames(RNY1_len))
RNY3_len = merge_matrix[,grepl('RNY3',colnames(merge_matrix))]
colnames(RNY3_len) <- sub('_[^_]*$','',colnames(RNY3_len))
RNY4_len = merge_matrix[,grepl('RNY4',colnames(merge_matrix))]
colnames(RNY4_len) <- sub('_[^_]*$','',colnames(RNY4_len))
RNY5_len = merge_matrix[,grepl('RNY5',colnames(merge_matrix))]
colnames(RNY5_len) <- sub('_[^_]*$','',colnames(RNY5_len))


merge_matrix <- RNY1_len + RNY3_len + RNY4_len + RNY5_len


library(readxl)
id <- read_xlsx('F:\\R_script\\20230919 乳腺癌\\result_final\\ID.xlsx')
ids <- sub('_$','',id$ID)
ids <- as.data.frame(cbind(ids,id$seq_ID))
colnames(ids)[1] <- 'ID'
colnames(ids)[2] <- 'seq-ID'
all <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\all_seq.csv')
a <- as.data.frame(all[match(ids$ID,all$X),])
rownames(a) <- ids$`seq-ID`
group <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\simple.csv')
a=as.data.frame(a[match(group$X,rownames(a)),])
a$group=group$condition
a=data.frame(a$X,a$ysRNA,a$group)
rownames(a) <- a$a.X
a=a[,-1]
colnames(a) <- c('count','group')

ben <- a %>% filter(a$group == 'BEN')
mal <- a %>% filter(a$group == 'MAL')


merge_matrix_ben <- merge_matrix[,match(rownames(ben),colnames(merge_matrix))]
total_reads_ben <- colSums(merge_matrix_ben)
result_ben <- data.frame(matrix(NA, nrow = nrow(merge_matrix_ben), ncol = ncol(merge_matrix_ben)))
colnames(result_ben) <- colnames(merge_matrix_ben)
rownames(result_ben) <- rownames(merge_matrix_ben)
for (i in 1:length(colnames(merge_matrix_ben))) {
  result_ben[,i] <- as.data.frame((merge_matrix_ben[,i] / total_reads_ben[i]) * 100)
}
result_ben$type <- rownames(result_ben)


merge_matrix_mal <- merge_matrix[,match(rownames(mal),colnames(merge_matrix))]
total_reads_mal <- colSums(merge_matrix_mal)
result_mal <- data.frame(matrix(NA, nrow = nrow(merge_matrix_mal), ncol = ncol(merge_matrix_mal)))
colnames(result_mal) <- colnames(merge_matrix_mal)
rownames(result_mal) <- rownames(merge_matrix_mal)
for (i in 1:length(colnames(merge_matrix_mal))) {
  result_mal[,i] <- as.data.frame((merge_matrix_mal[,i] / total_reads_mal[i]) * 100)
}
result_mal$type <- rownames(result_mal)
library(ggplot2)
library(dplyr)
library(reshape)

result_long <- melt(result_ben)
all_types <- unique(result_long$type)
# 指定要展示的ID列表
desired_IDs <- as.character(seq(15, 150))
# 将ID列转换为因子，并只包含desired_IDs
result_long$type <- factor(result_long$type, levels = desired_IDs)
custom_order <- c("15", "29",'43','57','71','86','100','113')
# 指定要展示的ID列表
#desired_IDs <- as.character(seq(1, 122))
# 将ID列转换为因子，并只包含desired_IDs
#result_long$type <- factor(result_long$type, levels = desired_IDs)
#result_long <- result_long[order(result_long$type),]
# 选择15到42之间的数据子集
#result_subset <- subset(result_long, type >= 15 & type <= 42)
#result_subset <- subset(result_subset,type != 150)



top.mar=0.2
right.mar=0.2
bottm.mar=0.2
left.mar=0.2
mytheme <- theme_classic()+
  theme(text = element_text(family = 'sans',colour = 'black',size = 11),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.6,colour = 'black'),
        axis.ticks.length = unit(1.5,units = 'mm'),
        plot.margin = unit(x=c(top.mar,right.mar,left.mar,bottm.mar),units = 'inches'))

p <- ggplot(result_long, aes(x = type, y = value,group = variable)) +
  geom_line(color = 'grey') + geom_point(color = '#7B68EE') +
  xlab("length(nt)") +
  ylab("Proportion") +
  ggtitle("ysRNAs")+
  mytheme+scale_x_discrete(breaks = custom_order,labels = custom_order) +
  theme(plot.title = element_text(hjust = 0.5, size = 14))+ annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "transparent", color = "black")
p
ggsave(p,filename = 'BEN_ysRNA_len占比图.pdf',width = 4,height = 4)




#############tRNA长度占比计算##############
setwd('F:\\R_script\\20230919 乳腺癌\\expr\\tRFs_len')
library(stringr)

all_file <- list.files()
ysRNA_len <- all_file[grepl('_sort.txt',all_file)]
read_txt_file <- function(filename) {
  data <- read.table(filename, header = FALSE, sep = '')  # 根据文件实际的分隔符来设置sep
  data$filename <- basename(filename)  # 添加文件名列
  return(data)
}
# 使用lapply批量读取文件并将结果存储在一个列表中
library(tidyverse)
data_list <- lapply(ysRNA_len, read_txt_file)
for (i in 1:length(data_list)) {
  colnames(data_list[[i]])[3] <- 'ID'
  colnames(data_list[[i]])[2] <- 'len'
  colnames(data_list[[i]])[1] <- 'reads'
}
merge_data=do.call(rbind,data_list)
merge_data$ID=gsub('_tRFs_len.txt_sort.txt','',merge_data$ID)
# 使用pivot_wider将长数据形式转换为宽数据形式
merge_matrix <- merge_data %>%
  pivot_wider(names_from = ID, values_from = reads)
merge_matrix <- as.data.frame(merge_matrix)
rownames(merge_matrix) <- merge_matrix$len
merge_matrix <- merge_matrix[,-1]
merge_matrix[is.na(merge_matrix)] <- 0

library(readxl)
id <- read_xlsx('F:\\R_script\\20230919 乳腺癌\\result_final\\ID.xlsx')
ids <- sub('_$','',id$ID)
ids <- as.data.frame(cbind(ids,id$seq_ID))
colnames(ids)[1] <- 'ID'
colnames(ids)[2] <- 'seq-ID'
all <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\all_seq.csv')
a <- as.data.frame(all[match(ids$ID,all$X),])
rownames(a) <- ids$`seq-ID`
group <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\simple.csv')
a=as.data.frame(a[match(group$X,rownames(a)),])
a$group=group$condition
a=data.frame(a$X,a$ysRNA,a$group)
rownames(a) <- a$a.X
a=a[,-1]
colnames(a) <- c('count','group')

ben <- a %>% filter(a$group == 'BEN')
mal <- a %>% filter(a$group == 'MAL')


merge_matrix_ben <- merge_matrix[,match(rownames(ben),colnames(merge_matrix))]
total_reads_ben <- colSums(merge_matrix_ben)
result_ben <- data.frame(matrix(NA, nrow = nrow(merge_matrix_ben), ncol = ncol(merge_matrix_ben)))
colnames(result_ben) <- colnames(merge_matrix_ben)
rownames(result_ben) <- rownames(merge_matrix_ben)
for (i in 1:length(colnames(merge_matrix_ben))) {
  result_ben[,i] <- as.data.frame((merge_matrix_ben[,i] / total_reads_ben[i]) * 100)
}
result_ben$type <- rownames(result_ben)


merge_matrix_mal <- merge_matrix[,match(rownames(mal),colnames(merge_matrix))]
total_reads_mal <- colSums(merge_matrix_mal)
result_mal <- data.frame(matrix(NA, nrow = nrow(merge_matrix_mal), ncol = ncol(merge_matrix_mal)))
colnames(result_mal) <- colnames(merge_matrix_mal)
rownames(result_mal) <- rownames(merge_matrix_mal)
for (i in 1:length(colnames(merge_matrix_mal))) {
  result_mal[,i] <- as.data.frame((merge_matrix_mal[,i] / total_reads_mal[i]) * 100)
}
result_mal$type <- rownames(result_mal)
library(ggplot2)
library(dplyr)
library(reshape)

result_long <- melt(result_ben)
custom_order <- c("16", "19",'22','25','28','31','34','37','40','43','46','49')



top.mar=0.2
right.mar=0.2
bottm.mar=0.2
left.mar=0.2
mytheme <- theme_classic()+
  theme(text = element_text(family = 'sans',colour = 'black',size = 11),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.6,colour = 'black'),
        axis.ticks.length = unit(1.5,units = 'mm'),
        plot.margin = unit(x=c(top.mar,right.mar,left.mar,bottm.mar),units = 'inches'))

p <- ggplot(result_long, aes(x = type, y = value,group = variable)) +
  geom_line(color = 'grey') + geom_point(color = '#006400') +
  xlab("length(nt)") +
  ylab("Proportion") +
  ggtitle("tRNAs")+
  mytheme+scale_x_discrete(breaks = custom_order,labels = custom_order)+
  theme(plot.title = element_text(hjust = 0.5, size = 14))+ annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "transparent", color = "black")
p
ggsave(p,filename = 'BEN_tRNA_len占比图.pdf',width = 4,height = 4)



#############miRNA长度占比计算############
setwd('F:\\R_script\\20230919 乳腺癌\\expr\\miRNALEN')
library(stringr)

all_file <- list.files()
miRNA_len <- all_file[grepl('_sort.txt',all_file)]
read_txt_file <- function(filename) {
  data <- read.table(filename, header = FALSE, sep = '')  # 根据文件实际的分隔符来设置sep
  data$filename <- basename(filename)  # 添加文件名列
  return(data)
}
# 使用lapply批量读取文件并将结果存储在一个列表中
library(tidyverse)
data_list <- lapply(miRNA_len, read_txt_file)
for (i in 1:length(data_list)) {
  colnames(data_list[[i]])[3] <- 'ID'
  colnames(data_list[[i]])[2] <- 'len'
  colnames(data_list[[i]])[1] <- 'reads'
}
merge_data=do.call(rbind,data_list)
merge_data$ID=gsub('_mapping2maturemiRNAs.sam_sort.txt','',merge_data$ID)
# 使用pivot_wider将长数据形式转换为宽数据形式
merge_matrix <- merge_data %>%
  pivot_wider(names_from = ID, values_from = reads)
merge_matrix <- as.data.frame(merge_matrix)
rownames(merge_matrix) <- merge_matrix$len
merge_matrix <- merge_matrix[,-1]
merge_matrix[is.na(merge_matrix)] <- 0


library(readxl)
id <- read_xlsx('F:\\R_script\\20230919 乳腺癌\\result_final\\ID.xlsx')
ids <- sub('_$','',id$ID)
ids <- as.data.frame(cbind(ids,id$seq_ID))
colnames(ids)[1] <- 'ID'
colnames(ids)[2] <- 'seq-ID'
all <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\all_seq.csv')
a <- as.data.frame(all[match(ids$ID,all$X),])
rownames(a) <- ids$`seq-ID`
group <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\simple.csv')
group <- group[-c(82,83),]
a=as.data.frame(a[match(group$X,rownames(a)),])
a$group=group$condition
a=data.frame(a$X,a$ysRNA,a$group)
rownames(a) <- a$a.X
a=a[,-1]
colnames(a) <- c('count','group')
ben <- a %>% filter(a$group == 'BEN')
mal <- a %>% filter(a$group == 'MAL')
merge_matrix_ben <- merge_matrix[,match(rownames(ben),colnames(merge_matrix))]
total_reads_ben <- colSums(merge_matrix_ben)
result_ben <- data.frame(matrix(NA, nrow = nrow(merge_matrix_ben), ncol = ncol(merge_matrix_ben)))
colnames(result_ben) <- colnames(merge_matrix_ben)
rownames(result_ben) <- rownames(merge_matrix_ben)
for (i in 1:length(colnames(merge_matrix_ben))) {
  result_ben[,i] <- as.data.frame((merge_matrix_ben[,i] / total_reads_ben[i]) * 100)
}
result_ben$type <- rownames(result_ben)


merge_matrix_mal <- merge_matrix[,match(rownames(mal),colnames(merge_matrix))]
total_reads_mal <- colSums(merge_matrix_mal)
result_mal <- data.frame(matrix(NA, nrow = nrow(merge_matrix_mal), ncol = ncol(merge_matrix_mal)))
colnames(result_mal) <- colnames(merge_matrix_mal)
rownames(result_mal) <- rownames(merge_matrix_mal)
for (i in 1:length(colnames(merge_matrix_mal))) {
  result_mal[,i] <- as.data.frame((merge_matrix_mal[,i] / total_reads_mal[i]) * 100)
}
result_mal$type <- rownames(result_mal)
library(ggplot2)
library(dplyr)
library(reshape)

result_long <- melt(result_ben)
result_long <- result_long[order(result_long$type),]

result_long <- melt(result_mal)
result_long <- result_long[order(result_long$type),]

top.mar=0.2
right.mar=0.2
bottm.mar=0.2
left.mar=0.2
mytheme <- theme_classic()+
  theme(text = element_text(family = 'sans',colour = 'black',size = 11),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.6,colour = 'black'),
        axis.ticks.length = unit(1.5,units = 'mm'),
        plot.margin = unit(x=c(top.mar,right.mar,left.mar,bottm.mar),units = 'inches'))

p <- ggplot(result_long, aes(x = type, y = value,group = variable)) +
  geom_line(color = 'grey') + geom_point(color = '#FF8C00') +
  xlab("length(nt)") +
  ylab("Proportion") +
  ggtitle("miRNAs")+
  mytheme+
  theme(plot.title = element_text(hjust = 0.5, size = 14))+ annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "transparent", color = "black")
p
ggsave(p,filename = 'ben_miRNA_len占比图.pdf',width = 4,height = 4)


##############lncRNA长度占比计算################
setwd('F:\\R_script\\20230919 乳腺癌\\expr\\lncRNA_len')
library(stringr)

all_file <- list.files()
miRNA_len <- all_file[grepl('_sort.txt',all_file)]
read_txt_file <- function(filename) {
  data <- read.table(filename, header = FALSE, sep = '')  # 根据文件实际的分隔符来设置sep
  data$filename <- basename(filename)  # 添加文件名列
  return(data)
}
# 使用lapply批量读取文件并将结果存储在一个列表中
library(tidyverse)
data_list <- lapply(miRNA_len, read_txt_file)
for (i in 1:length(data_list)) {
  colnames(data_list[[i]])[3] <- 'ID'
  colnames(data_list[[i]])[2] <- 'len'
  colnames(data_list[[i]])[1] <- 'reads'
}
merge_data=do.call(rbind,data_list)
merge_data$ID=gsub('_lncRNA_S1_gene.sam_sort.txt','',merge_data$ID)
# 使用pivot_wider将长数据形式转换为宽数据形式
merge_matrix <- merge_data %>%
  pivot_wider(names_from = ID, values_from = reads)
merge_matrix <- as.data.frame(merge_matrix)
rownames(merge_matrix) <- merge_matrix$len
merge_matrix <- merge_matrix[,-1]
merge_matrix[is.na(merge_matrix)] <- 0


library(readxl)
id <- read_xlsx('F:\\R_script\\20230919 乳腺癌\\result_final\\ID.xlsx')
ids <- sub('_$','',id$ID)
ids <- as.data.frame(cbind(ids,id$seq_ID))
colnames(ids)[1] <- 'ID'
colnames(ids)[2] <- 'seq-ID'
all <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\all_seq.csv')
a <- as.data.frame(all[match(ids$ID,all$X),])
rownames(a) <- ids$`seq-ID`
group <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\simple.csv')
group <- group[-c(82,83),]
a=as.data.frame(a[match(group$X,rownames(a)),])
a$group=group$condition
a=data.frame(a$X,a$ysRNA,a$group)
rownames(a) <- a$a.X
a=a[,-1]
colnames(a) <- c('count','group')
ben <- a %>% filter(a$group == 'BEN')
mal <- a %>% filter(a$group == 'MAL')


merge_matrix_ben <- merge_matrix[,match(rownames(ben),colnames(merge_matrix))]
total_reads_ben <- colSums(merge_matrix_ben)
result_ben <- data.frame(matrix(NA, nrow = nrow(merge_matrix_ben), ncol = ncol(merge_matrix_ben)))
colnames(result_ben) <- colnames(merge_matrix_ben)
rownames(result_ben) <- rownames(merge_matrix_ben)
for (i in 1:length(colnames(merge_matrix_ben))) {
  result_ben[,i] <- as.data.frame((merge_matrix_ben[,i] / total_reads_ben[i]) * 100)
}
result_ben$type <- rownames(result_ben)


merge_matrix_mal <- merge_matrix[,match(rownames(mal),colnames(merge_matrix))]
total_reads_mal <- colSums(merge_matrix_mal)
result_mal <- data.frame(matrix(NA, nrow = nrow(merge_matrix_mal), ncol = ncol(merge_matrix_mal)))
colnames(result_mal) <- colnames(merge_matrix_mal)
rownames(result_mal) <- rownames(merge_matrix_mal)
for (i in 1:length(colnames(merge_matrix_mal))) {
  result_mal[,i] <- as.data.frame((merge_matrix_mal[,i] / total_reads_mal[i]) * 100)
}
result_mal$type <- rownames(result_mal)
library(ggplot2)
library(dplyr)
library(reshape)

result_long <- melt(result_ben)
result_long <- result_long[order(result_long$type),]
#result_long <- result_long %>% filter(value > 0.5)

result_long <- melt(result_mal)
result_long <- result_long[order(result_long$type),]
library(ggplot2)
library(dplyr)
library(reshape)

#result_long <- melt(result)
result_long <- result_long[order(result_long$type),]
all_types <- unique(result_long$type)

# 从中间抽取大约 9 个标签
num_labels <- 8
step <- ceiling(length(all_types) / num_labels)
custom_order <- factor(all_types[seq(from = 1, to = length(all_types), by = step)])
custom_order <- c("19M", "23M",'45M','56M','67M','78M','88M','104M','128M','150M')
# 将type变量转换为因子，并且指定自定义的顺序
result_long$type <- factor(result_long$type, levels = custom_order)
result_long <- na.omit(result_long)




top.mar=0.2
right.mar=0.2
bottm.mar=0.2
left.mar=0.2
mytheme <- theme_classic()+
  theme(text = element_text(family = 'sans',colour = 'black',size = 11),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.6,colour = 'black'),
        axis.ticks.length = unit(1.5,units = 'mm'),
        plot.margin = unit(x=c(top.mar,right.mar,left.mar,bottm.mar),units = 'inches'))

p <- ggplot(result_long, aes(x = type, y = value,group = variable)) +
  geom_line(color = 'grey') + geom_point(color = '#B22222') +
  xlab("length(nt)") +
  ylab("Proportion") +
  ggtitle("lncRNAs")+
  mytheme+
  theme(plot.title = element_text(hjust = 0.5, size = 14))+ annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "transparent", color = "black")
p
ggsave(p,filename = 'MAL_lncRNA_len占比图.pdf',width = 4,height = 4)



############mRNA长度占比计算##################
setwd('F:\\R_script\\20230919 乳腺癌\\expr\\mRNA_len')
library(stringr)

all_file <- list.files()
miRNA_len <- all_file[grepl('_sort.txt',all_file)]
read_txt_file <- function(filename) {
  data <- read.table(filename, header = FALSE, sep = '')  # 根据文件实际的分隔符来设置sep
  data$filename <- basename(filename)  # 添加文件名列
  return(data)
}
# 使用lapply批量读取文件并将结果存储在一个列表中
library(tidyverse)
data_list <- lapply(miRNA_len, read_txt_file)
for (i in 1:length(data_list)) {
  colnames(data_list[[i]])[3] <- 'ID'
  colnames(data_list[[i]])[2] <- 'len'
  colnames(data_list[[i]])[1] <- 'reads'
}
merge_data=do.call(rbind,data_list)
merge_data$ID=gsub('_mRNA_S1_gene.sam_sort.txt','',merge_data$ID)
# 使用pivot_wider将长数据形式转换为宽数据形式
merge_matrix <- merge_data %>%
  pivot_wider(names_from = ID, values_from = reads)
merge_matrix <- as.data.frame(merge_matrix)
rownames(merge_matrix) <- merge_matrix$len
merge_matrix <- merge_matrix[,-1]
merge_matrix[is.na(merge_matrix)] <- 0


library(readxl)
id <- read_xlsx('F:\\R_script\\20230919 乳腺癌\\result_final\\ID.xlsx')
ids <- sub('_$','',id$ID)
ids <- as.data.frame(cbind(ids,id$seq_ID))
colnames(ids)[1] <- 'ID'
colnames(ids)[2] <- 'seq-ID'
all <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\all_seq.csv')
a <- as.data.frame(all[match(ids$ID,all$X),])
rownames(a) <- ids$`seq-ID`
group <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\simple.csv')
group <- group[-c(82,83),]
a=as.data.frame(a[match(group$X,rownames(a)),])
a$group=group$condition
a=data.frame(a$X,a$ysRNA,a$group)
rownames(a) <- a$a.X
a=a[,-1]
colnames(a) <- c('count','group')
ben <- a %>% filter(a$group == 'BEN')
mal <- a %>% filter(a$group == 'MAL')


merge_matrix_ben <- merge_matrix[,match(rownames(ben),colnames(merge_matrix))]
total_reads_ben <- colSums(merge_matrix_ben)
result_ben <- data.frame(matrix(NA, nrow = nrow(merge_matrix_ben), ncol = ncol(merge_matrix_ben)))
colnames(result_ben) <- colnames(merge_matrix_ben)
rownames(result_ben) <- rownames(merge_matrix_ben)
for (i in 1:length(colnames(merge_matrix_ben))) {
  result_ben[,i] <- as.data.frame((merge_matrix_ben[,i] / total_reads_ben[i]) * 100)
}
result_ben$type <- rownames(result_ben)


merge_matrix_mal <- merge_matrix[,match(rownames(mal),colnames(merge_matrix))]
total_reads_mal <- colSums(merge_matrix_mal)
result_mal <- data.frame(matrix(NA, nrow = nrow(merge_matrix_mal), ncol = ncol(merge_matrix_mal)))
colnames(result_mal) <- colnames(merge_matrix_mal)
rownames(result_mal) <- rownames(merge_matrix_mal)
for (i in 1:length(colnames(merge_matrix_mal))) {
  result_mal[,i] <- as.data.frame((merge_matrix_mal[,i] / total_reads_mal[i]) * 100)
}
result_mal$type <- rownames(result_mal)
library(ggplot2)
library(dplyr)
library(reshape)

result_long <- melt(result_ben)
result_long <- result_long[order(result_long$type),]
#result_long <- result_long %>% filter(value > 0.5)

result_long <- melt(result_mal)
result_long <- result_long[order(result_long$type),]
library(ggplot2)
library(dplyr)
library(reshape)

all_types <- unique(result_long$type)

# 从中间抽取大约 9 个标签
num_labels <- 8
step <- ceiling(length(all_types) / num_labels)
custom_order <- factor(all_types[seq(from = 1, to = length(all_types), by = step)])
custom_order <- c("19M", "23M",'45M','56M','67M','78M','88M','104M','128M','150M')
# 将type变量转换为因子，并且指定自定义的顺序
result_long$type <- factor(result_long$type, levels = custom_order)
result_long <- na.omit(result_long)




top.mar=0.2
right.mar=0.2
bottm.mar=0.2
left.mar=0.2
mytheme <- theme_classic()+
  theme(text = element_text(family = 'sans',colour = 'black',size = 11),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.6,colour = 'black'),
        axis.ticks.length = unit(1.5,units = 'mm'),
        plot.margin = unit(x=c(top.mar,right.mar,left.mar,bottm.mar),units = 'inches'))

p <- ggplot(result_long, aes(x = type, y = value,group = variable)) +
  geom_line(color = 'grey') + geom_point(color = '#FF0000') +
  xlab("length(nt)") +
  ylab("Proportion") +
  ggtitle("mRNAs")+
  mytheme+
  theme(plot.title = element_text(hjust = 0.5, size = 14))+ annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "transparent", color = "black")
p
ggsave(p,filename = 'BEN_mRNA_len占比图.pdf',width = 4,height = 4)











