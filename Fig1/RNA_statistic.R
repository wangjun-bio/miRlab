data <- read.table('fkrm.txt',header = F,sep = '')
colnames(data)[1] <- 'Count'
colnames(data)[2] <- 'ID'
data1 <- data
####匹配定义的suffix内的名字的后缀#########
suffix <- 'ysRNA'##########替换不同种类的cfRNA
ysRNA<- data1[grepl(paste0(suffix, "$"), data1$ID), ]
####定义一个空的字符向量
new_ids <- character(length(ysRNA$ID))
library(stringr)
for (i in 1:length(ysRNA$ID)) {
  newstr <- str_split(ysRNA$ID,pattern = '_')[[i]]
  firstpart <- newstr[1]
  secondpart <- newstr[2]
  new_id <- paste(firstpart,secondpart,sep = '_')
  # 将新的ID添加到向量中
  new_ids[i] <- new_id
}
rownames(ysRNA) <- new_ids
colnames(ysRNA)[1] <- 'ysRNA'
# 删除第二列，保持第一列为一个数据框
ysRNA <- ysRNA[, "ysRNA", drop = FALSE]

######第一列的后缀与低于15reads的后缀相同，需要分步取############
suffix <- 'trim'
all <- data1[grepl(paste0(suffix, "$"), data1$ID), ]
# 隔行取130行
n_rows <- 260
step <- 2  # 隔行的步长
selected_rows <- seq(1, n_rows, by = step)
selected_rows1 <- seq(2,n_rows,by = step)
all1 <- all[selected_rows,]
low15 <- all[selected_rows1,]

new_ids <- character(length(low15$ID))
for (i in 1:length(low15$ID)) {
  newstr <- str_split(low15$ID,pattern = '_')[[i]]
  firstpart <- newstr[1]
  secondpart <- newstr[2]
  new_id <- paste(firstpart,secondpart,sep = '_')
  # 将新的ID添加到向量中
  new_ids[i] <- new_id
}

rownames(low15) <- new_ids
colnames(low15)[1] <- 'low15'
# 删除第二列，保持第一列为数据框
low15 <- low15[, "low15", drop = FALSE]



##将所有数据合并在一个矩阵
all_data <- cbind(lncRNA,low15,piRNA,miRNA,tRNA,ysRNA,rsRNA,all1)
all_data2 <- all_data
# 循环遍历每一列，从第二列开始（第一列不需要减去）
for (i in 1:7) {
  # 计算当前列的结果并赋值给下一列
  all_data2[, i] <- (all_data2[, i+1] - all_data2[, i])/4
}
all_counts <- all_data2[,8]/4
all_data2 <- all_data2[,-8]
all_data2 <- cbind(all_data2,all_counts)
unspecified <- as.data.frame(all_data2[,8]-(all_data2[,1]+all_data2[,3]+all_data2[,4]+all_data2[,5]+all_data2[,6]+all_data2[,7]))
colnames(unspecified) <- 'unspecified'
all_data2 <- cbind(all_data2,unspecified)
all_data2 <- all_data2[,-2]
colnames(all_data2) <- c('mRNA+lncRNA','piRNA','miRNA','tRNA','ysRNA','rsRNA','all_counts','unspecified')
all_data2 <- all_data2[,c('mRNA+lncRNA','piRNA','miRNA','tRNA','ysRNA','rsRNA','all_counts','unspecified')]

# 获取总counts列的名称
total_counts <- all_data2$all_counts
all_data2 <- cbind(all_data2,total_counts)
all_data2 <- all_data2[,-7]
total_counts_col <- names(all_data2)[ncol(all_data2)]
write.csv(all_data2,'RNA含量.csv')
# 循环计算每种RNA的占比
for (i in 1:(ncol(all_data2) - 1)) {  # 从第一列到倒数第二列
  col_name <- names(all_data2)[i]  # 获取当前列的名称
  all_data2[col_name] <- (all_data2[col_name] / all_data2[total_counts_col])
}
RNA <- all_data2
RNA[,8] <- 1
ID <- rownames(RNA)
a <- RNA[match(ids$ID,ID),]
a <- read.csv('all_seq.csv',row.names = 1,header = T)

library(readxl)
id <- read_xlsx('ID.xlsx')
ids <- sub('_$','',id$ID)
ids <- as.data.frame(cbind(ids,id$seq_ID))
colnames(ids)[1] <- 'ID'
colnames(ids)[2] <- 'seq-ID'
group <- read.csv('simple.csv')
group <- group[-c(82,83),]
ids <- ids[match(group$X,ids$`seq-ID`),]

a <- a[match(ids$ID,rownames(a)),]
rownames(a) <- ids$`seq-ID`
# 将数据矩阵转换为长格式

#a$all <- rowSums(a)
#for (i in 1:6) {
   # a[,i] <- (a[,i] / a$all)
 # }
a$sample <- rownames(a)
#a$group <- group$condition
a <- a[,-8]
library(dplyr)
library(tidyverse)
library(reshape2)
data_long <- melt(a, id.vars = 'sample')
data_long$log10 <- log10(data_long$value)
data_long$group <- rep(group$condition,7)

# 创建矩形图###
library(ggplot2)
p = ggplot(data_long,aes(x = sample,y = value,fill = variable)) + geom_bar(stat = 'identity',position = 'stack')
p
library(ggsci)
p <- p +scale_fill_nejm()+xlab('sample') + ylab('Proportion') + theme(axis.text.x = element_blank())+facet_wrap(~group,scales = 'free',nrow = 2)
p
ggsave(p,filename = 'RNA含量堆叠柱状图.pdf',height = 8,width = 10)
############log10转换####
top.mar=0.2
right.mar=0.2
bottm.mar=0.2
left.mar=0.2
mytheme <- theme_classic()+
  theme(text = element_text(family = 'sans',colour = 'black',size = 11),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.6,colour = 'black'),
        axis.ticks.length = unit(1.5,units = 'mm'),
        axis.text.x = element_blank(),
        plot.margin = unit(x=c(top.mar,right.mar,left.mar,bottm.mar),units = 'inches'))
pdf(file = 'RNA含量图.pdf')
# 设置变量unspecified的颜色为灰色
 ggplot(data_long, aes(x=sample, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x='Sample',y = "Proportion") +
  theme_minimal() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set3")+  # 设置颜色映射
   mytheme #隐藏行名
dev.off()
 
p <- ggplot(data_long, aes(x=value, fill=variable)) +
  geom_density(alpha=0.5) +
  theme_bw() +
  labs( x="log10 Proporition", y="Density")
p


#install.packages('PerformanceAnalytics')
#install.packages('GGally')
write.csv(df,'all_seq.csv')
#install.packages('ggplot2')
library(GGally)
library(ggplot2)
library(PerformanceAnalytics)
df <- cbind(miRNA[,1],mRNA_lncRNA[,1],piRNA,rsRNA,ysRNA,tRNA)
colnames(df) <- c('miRNA','mRNA+lncRNA','piRNA','rsRNA','ysRNA','tRNA')
chart.Correlation(df)
library(ggpubr)
ggpairs(df)
group=read.csv('simple.csv',row.names = 1,header = T)
df1 <- a[match(rownames(group),rownames(a)),]
df1 <- cbind(group$condition,df1)
colnames(df1)[1] <- 'group'
library(ggsci)
ggpairs(df1[,-1],mapping = ggplot2::aes(color=df1$group))

##############饼状图统计#################
a <- cbind(a,group$condition)
BEN <- a %>% filter(a$`group$condition` == 'BEN')
MAL <- a %>% filter(a$`group$condition` == 'MAL')
BEN <- BEN[,-7]
MAL <- MAL[,-7]

pie_percentage <- function(x){
  result <- as.data.frame(colSums(x))
  total <- sum(x)
  percentage <- (result / total) * 100
  return(percentage)
}
BEN_per <- pie_percentage(BEN)
MAL_per <- pie_percentage(MAL)
colors <- c("red", "blue", "green", "orange",'gray','yellow')
pie(BEN_per$`colSums(x)`, labels = paste(rownames(BEN_per),round(BEN_per$`colSums(x)`,3)), col = colors) + scale_fill_nejm()
colors <- c("red", "blue", "green", "orange",'gray','yellow')
pie(MAL_per$`colSums(x)`, labels = paste(rownames(MAL_per),round(MAL_per$`colSums(x)`,3)), col = colors)


#devtools::install_github('cardiomoon/moonBook')
require(moonBook)
#devtools::install_github('cardiomoon/webr')
require(webr)

a$group <- group$condition
library(tidyverse)
mal <- a %>% filter(group == 'MAL')
mal <- mal[,-9]
ben <- a %>% filter(group == 'BEN')
ben <- ben[,-9]

sums <- function(x){
  x <- as.data.frame(colSums(x))
  for (i in 1:(nrow(x)-1)) {
    x[i,] <- as.data.frame(x[i,] / x[8,])
  }
  return(x)
}
mal_sums <- sums(mal)
ben_sums <- sums(ben)
colnames(mal_sums) <- 'count'
mal_sums$MAL <- rownames(mal_sums)
mal_sums <- mal_sums[-8,]
colnames(ben_sums) <- 'count'
ben_sums$BEN <- rownames(ben_sums)
ben_sums <- ben_sums[-8,]

custom_colors <- c("#FF8C00","#B22222","#00CED1","#87CEEB","#006400", "#FFFACD","#7B68EE")
p <- PieDonut(ben_sums,aes(BEN,count = count),r0 = 0.7,labelpositionThreshold = 0.1) + scale_fill_manual(values = custom_colors)
ggsave(p,filename = 'BEN饼状图.pdf',width = 5,height = 5)






