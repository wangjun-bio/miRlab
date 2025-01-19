################ysRNA#####################
library(readxl)
group <- read.csv('simple.csv')
group <- group[-c(82,83),]
id <- read_xlsx('ID.xlsx')

library(stringr)
ysRNA_seq <- read.table('F:\\R_script\\20230919 乳腺癌\\expr\\ysrna_len.txt')
ysRNY1 = ysRNA_seq[grepl('RNY1_len.txt',ysRNA_seq$V2),]
ysRNY1$V2 <- sub('_[^_]*$','',ysRNY1$V2)
ysRNY1$V2 <- sub('_[^_]*$','',ysRNY1$V2) 
ysRNY1 <- ysRNY1[match(id$ID,ysRNY1$V2),]
colnames(ysRNY1)[1] <- 'ysRNA-RNY1'

ysRNY3 = ysRNA_seq[grepl('RNY3_len.txt',ysRNA_seq$V2),]
ysRNY3$V2 <- sub('_[^_]*$','',ysRNY3$V2)
ysRNY3$V2 <- sub('_[^_]*$','',ysRNY3$V2) 
ysRNY3 <- ysRNY3[match(id$ID,ysRNY3$V2),]
colnames(ysRNY3)[1] <- 'ysRNA-RNY3'

ysRNY4 = ysRNA_seq[grepl('RNY4_len.txt',ysRNA_seq$V2),]
ysRNY4$V2 <- sub('_[^_]*$','',ysRNY4$V2)
ysRNY4$V2 <- sub('_[^_]*$','',ysRNY4$V2) 
ysRNY4 <- ysRNY4[match(id$ID,ysRNY4$V2),]
colnames(ysRNY4)[1] <- 'ysRNA-RNY4'

ysRNY5 = ysRNA_seq[grepl('RNY5_len.txt',ysRNA_seq$V2),]
ysRNY5$V2 <- sub('_[^_]*$','',ysRNY5$V2)
ysRNY5$V2 <- sub('_[^_]*$','',ysRNY5$V2) 
ysRNY5 <- ysRNY5[match(id$ID,ysRNY5$V2),]
colnames(ysRNY5)[1] <- 'ysRNA-RNY5'


ys_propotion=as.data.frame(cbind(ysRNY1$`ysRNA-RNY1`,ysRNY3$`ysRNA-RNY3`,ysRNY4$`ysRNA-RNY4`,ysRNY5$`ysRNA-RNY5`))
colnames(ys_propotion)=c('ysRNA-RNY1','ysRNA-RNY3','ysRNA-RNY4','ysRNA-RNY5')
rownames(ys_propotion)=id$seq_ID
ys_propotion <- ys_propotion[match(group$X,rownames(ys_propotion)),]
ys_propotion$all <- rowSums(ys_propotion)
total_reads <- rowSums(ys_propotion)
head(ys_propotion)
ys_propotion <- as.data.frame(t(ys_propotion))
for (i in 1:(ncol(ys_propotion)-1)) {
  ys_propotion[,i] <- (ys_propotion[,i] / ys_propotion$all) * 100
}
for (i in 1:ncol(ys_propotion)) {
  ys_propotion[,i] <- (ys_propotion[,i] / total_reads[i]) * 1e6
}
ys_propotion <- as.data.frame(t(ys_propotion))
ys_propotion <- cbind(ys_propotion,group$condition)
colnames(ys_propotion)[5] <- 'group'


library(ggplot2)
library(reshape)
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

# 将数据从宽格式转换为长格式
ys_long <- tidyr::pivot_longer(ys_propotion, cols = starts_with("ysRNA"), names_to = "ysRNA", values_to = "value")

p <- ggplot(ys_long, aes(x=ysRNA,y=value, fill = group)) +
  geom_boxplot() +
  scale_fill_brewer(palette = 'Accent') +
  ylab('Proportion %') + xlab(NULL) +
  mytheme +
  theme(legend.position = 'right')+
 annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "transparent", color = "black")
p
ggsave(p,filename = 'RNY亚型长度表达.pdf',width = 8,height = 6)

p <- ggplot(ys_long, aes(x=ysRNA,y=value, fill = group)) +
  geom_boxplot() +
  scale_fill_brewer(palette = 'Accent') +
  ylab('RPM') + xlab(NULL) +
  mytheme +
  theme(legend.position = 'right')+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "transparent", color = "black")
p
ggsave(p,filename = 'RNY亚型长度表达.pdf',width = 8,height = 6)

#####################rsRNA######################
library(stringr)
rsRNA_seq <- read.table('F:\\R_script\\20230919 乳腺癌\\expr\\rslen.txt')
rs18S_len = rsRNA_seq[grepl('-18S_len.txt',rsRNA_seq$V2),]
rs18S_len$V2 <- sub('_[^_]*$','',rs18S_len$V2)
rs18S_len$V2 <- sub('_[^_]*$','',rs18S_len$V2) 
rs18S_len <- rs18S_len[match(id$ID,rs18S_len$V2),]
colnames(rs18S_len)[1] <- 'rsRNA-18S'

rs28S_len = rsRNA_seq[grepl('-28S_len.txt',rsRNA_seq$V2),]
rs28S_len$V2 <- gsub('_[^_]*$','',rs28S_len$V2)
rs28S_len$V2 <- gsub('_[^_]*$','',rs28S_len$V2)
rs28S_len <- rs28S_len[match(id$ID,rs28S_len$V2),]
colnames(rs28S_len)[1] <- 'rsRNA-28S'

rs5S_len = rsRNA_seq[grepl('-5S_len.txt',rsRNA_seq$V2),]
rs5S_len$V2 <- gsub('_[^_]*$','',rs5S_len$V2)
rs5S_len$V2 <- gsub('_[^_]*$','',rs5S_len$V2)
rs5S_len <- rs28S_len[match(id$ID,rs5S_len$V2),]
colnames(rs5S_len)[1] <- 'rsRNA-5S'

rs5.8S_len = rsRNA_seq[grepl('-5.8S_len.txt',rsRNA_seq$V2),]
rs5.8S_len$V2 <- gsub('_[^_]*$','',rs5.8S_len$V2)
rs5.8S_len$V2 <- gsub('_[^_]*$','',rs5.8S_len$V2)
rs5.8S_len <- rs5.8S_len[match(id$ID,rs5.8S_len$V2),]
colnames(rs5.8S_len)[1] <- 'rsRNA-5.8S'

rsmt12S_len = rsRNA_seq[grepl('-mt12S_len.txt',rsRNA_seq$V2),]
rsmt12S_len$V2 <- gsub('_[^_]*$','',rsmt12S_len$V2)
rsmt12S_len$V2 <- gsub('_[^_]*$','',rsmt12S_len$V2)
rsmt12S_len <- rsmt12S_len[match(id$ID,rsmt12S_len$V2),]
colnames(rsmt12S_len)[1] <- 'rsRNA-mt12S'

rsmt16S_len = rsRNA_seq[grepl('-mt16S_len.txt',rsRNA_seq$V2),]
rsmt16S_len$V2 <- gsub('_[^_]*$','',rsmt16S_len$V2)
rsmt16S_len$V2 <- gsub('_[^_]*$','',rsmt16S_len$V2)
rsmt16S_len <- rsmt12S_len[match(id$ID,rsmt16S_len$V2),]
colnames(rsmt16S_len)[1] <- 'rsRNA-mt16S'


rs_propotion=as.data.frame(cbind(rs18S_len$`rsRNA-18S`,rs28S_len$`rsRNA-28S`,rs5.8S_len$`rsRNA-5.8S`,rs5S_len$`rsRNA-5S`,rsmt12S_len$`rsRNA-mt12S`,rsmt16S_len$`rsRNA-mt16S`))
colnames(rs_propotion)=c('rsRNA-18S','rsRNA-28S','rsRNA-5.8S','rsRNA-5S','rsRNA-mt12S','rsRNA-mt16S')
rownames(rs_propotion)=id$seq_ID
rs_propotion[is.na(rs_propotion)] <- 0
rs_propotion <- rs_propotion[match(group$X,rownames(rs_propotion)),]
total_reads <- rowSums(rs_propotion)
rs_propotion$all <- rowSums(rs_propotion)
head(rs_propotion)
for (i in 1:(ncol(rs_propotion)-1)) {
  rs_propotion[,i] <- (rs_propotion[,i] / rs_propotion$all) * 100
}

rs_propotion <- as.data.frame(t(rs_propotion))
for (i in 1:ncol(rs_propotion)) {
  rs_propotion[,i] <- (rs_propotion[,i] / total_reads[i]) * 1e6
}
rs_propotion <- as.data.frame(t(rs_propotion))
rs_propotion <- cbind(rs_propotion,group$condition)
colnames(rs_propotion)[7] <- 'group'


top.mar=0.2
right.mar=0.2
bottm.mar=0.2
left.mar=0.2
mytheme <- theme_classic()+
  theme(text = element_text(family = 'sans',colour = 'black',size = 10),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.6,colour = 'black'),
        axis.ticks.length = unit(1.5,units = 'mm'),
        plot.margin = unit(x=c(top.mar,right.mar,left.mar,bottm.mar),units = 'inches'))

# 将数据从宽格式转换为长格式
rs_long <- tidyr::pivot_longer(rs_propotion, cols = starts_with("rsRNA"), names_to = "rsRNA", values_to = "value")

p <- ggplot(rs_long, aes(x=rsRNA,y=value, fill = group)) +
  geom_boxplot() +
  scale_fill_brewer(palette = 'Accent') +
  ylab('Proportion %') + xlab(NULL) +
  mytheme +
  theme(legend.position = 'right') +
 annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "transparent", color = "black")
p

p <- ggplot(rs_long, aes(x=rsRNA,y=value, fill = group)) +
  geom_boxplot() +
  scale_fill_brewer(palette = 'Accent') +
  ylab('RPM') + xlab(NULL) +
  mytheme +
  theme(legend.position = 'right') +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "transparent", color = "black")
p

ggsave(p,filename = 'rsRNA亚型长度表达.pdf',width = 8,height = 6)

################tsRNA_TEST#####################
#BiocManager::install('seqinr')
setwd('F:\\R_script\\20230919 乳腺癌\\expr')
library(seqinr)
library(readxl)
trfs <- as.data.frame(read_xlsx('tRFs_original_202103-1.xlsx'))
colnames(trfs) <- c('MINTbase Unique ID ','Fragment sequence','Fragment Length','5_half','5_tRF','i_tRF','3_tRF','3_half','Total number of genomic instances','10','11','12')
head(trfs)
library(tidyverse)
# 计算特定亚型出现的次数
i_tRFS <- data.frame(trfs$`MINTbase Unique ID `,trfs$i_tRF)
i_tRFS <- na.omit(i_tRFS)
i_tRFS$trfs.i_tRF <- substr(i_tRFS$trfs.i_tRF, 1, 3)
t5_tRF <- data.frame(trfs$`MINTbase Unique ID `,trfs$`5_tRF`)
t5_tRF <- na.omit(t5_tRF)
t5_tRF$trfs..5_tRF.<- substr(t5_tRF$trfs..5_tRF., 1, 3)
t5_half <- data.frame(trfs$`MINTbase Unique ID `,trfs$`5_half`)
t5_half <- na.omit(t5_half)
t5_half$trfs..5_half. <- substr(t5_half$trfs..5_half., 1, 3)
t3_tRF <- data.frame(trfs$`MINTbase Unique ID `,trfs$`3_tRF`)
t3_tRF <- na.omit(t3_tRF)
t3_tRF$trfs..3_tRF. <- substr(t3_tRF$trfs..3_tRF., 1, 3)
t3_half <- data.frame(trfs$`MINTbase Unique ID `,trfs$`3_half`)
t3_half <- na.omit(t3_half)
t3_half$trfs..3_half. <- substr(t3_half$trfs..3_half., 1, 3)
count <- read.csv('tRNA_rpm2.csv',row.names = 1,header = T)
table(i_tRFS$trfs.i_tRF)
table(t5_tRF$trfs..5_tRF.)


group = read.csv('simple.csv',header = T)
group <- group[-c(82,83),]
count <- count[,match(group$X,colnames(count))]
colnames(i_tRFS)[1] <- 'ID'
i_tRFS_count=count[match(i_tRFS$ID,rownames(count)),]
t5_tRF_count=count[match(t5_tRF$trfs..MINTbase.Unique.ID..,rownames(count)),]
t5_half_count=count[match(t5_half$trfs..MINTbase.Unique.ID..,rownames(count)),]
t3_half_count=count[match(t3_half$trfs..MINTbase.Unique.ID..,rownames(count)),]
t3_tRF_count=count[match(t3_tRF$trfs..MINTbase.Unique.ID..,rownames(count)),]

sample_mean <- function(x,group){
  sample_means=as.data.frame(apply(x, 2, mean),colnames = colnames(x))
  sample_means=cbind(sample_means,group$condition)
  colnames(sample_means)=c('values','group')
  library(reshape)
  data_long = melt(sample_means)
  return(data_long)
}
data_long <- sample_mean(i_tRFS_count,group = group)
data_long2 <- sample_mean(t5_tRF_count,group = group)
data_long3 <- sample_mean(t5_half_count,group = group)
data_long4 <- sample_mean(t3_half_count,group = group)
data_long5 <- sample_mean(t3_tRF_count,group = group)
library(dplyr)
library(ggsignif)
library(tidyr)
library(reshape)
data_long$variable = 'Internal'
data_long2$variable = '5-tRFs'
data_long3$variable = '5-half'
data_long4$variable = '3-tRFs'
data_long5$variable = '3-half'

data_all <- rbind(data_long,data_long2,data_long3,data_long4,data_long5)

top.mar=0.2
right.mar=0.2
bottm.mar=0.2
left.mar=0.2
mytheme <- theme_classic()+
  theme(text = element_text(family = 'sans',colour = 'black',size = 12),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.6,colour = 'black'),
        axis.ticks.length = unit(1.5,units = 'mm'),
        plot.margin = unit(x=c(top.mar,right.mar,left.mar,bottm.mar),units = 'inches'))
library(ggplot2)
library(gridExtra)
library(cowplot)
library(ggpubr)
ggplot(data_all,aes(variable,value,fill=group,add = 'none'))+
  geom_boxplot()+
  stat_compare_means(method = 't.test')+
  ylab('means of RPM')+
  scale_fill_brewer(palette = 'Dark2')+
  mytheme

plot_box <- function(x){
  library(ggplot2)
  library(ggpubr)
  top.mar=0.2
  right.mar=0.2
  bottm.mar=0.2
  left.mar=0.2
  mytheme <- theme_classic()+
    theme(text = element_text(family = 'sans',colour = 'black',size = 12),
          axis.line = element_blank(),
          axis.ticks = element_line(size = 0.6,colour = 'black'),
          axis.ticks.length = unit(1.5,units = 'mm'),
          axis.title.x = NULL,
          plot.margin = unit(x=c(top.mar,right.mar,left.mar,bottm.mar),units = 'inches'))
  p <- ggplot(x,aes(variable,value,fill=group,add = 'none'))+
    geom_boxplot()+
    stat_compare_means(method = 't.test',label = 'p.signif')+
    ylab('means of RPM')+
    scale_fill_brewer(palette = 'Dark2')+
    mytheme
  return(p)
}


# 创建第一个图形
plot1 <- plot_box(data_long)
# 创建第二个图形
plot2 <- plot_box(data_long2)

plot3 <- plot_box(data_long3)

plot4 <- plot_box(data_long4)
  
plot5 <- plot_box(data_long5)
# 将图形放在一个画板中
pdf('tRFs各片段表达.pdf',width = 8,height = 8)
plot_grid(plot1, plot2, plot3, plot4, plot5, align = 'v', ncol = 3)
dev.off()

colnames(i_tRFS)[2] <- 'type'
i_tRFS$condition <- 'Internal'
colnames(t5_tRF) <- c('ID','type')
t5_tRF$condition <- '5-tRF'
colnames(t5_half) <- c('ID','type')
t5_half$condition <- '5-half'
colnames(t3_tRF) <- c('ID','type')
t3_tRF$condition <- '3-tRF'
colnames(t3_half) <- c('ID','type')
t3_half$condition <- '3-half'
type <- rbind(i_tRFS,t5_tRF,t5_half,t3_tRF,t3_half)
table(type$type)
count_trf <- count[match(type$ID,rownames(count)),]
count_trf$type <- type$type
count_trf$condition <- type$condition
table(count_trf$type)

sample_means=as.data.frame(apply(count_trf[,-c(85,84)], 2, mean),colnames = colnames(count_trf))
sample_means=cbind(sample_means,group$condition)
library(tidyverse)
library(reshape)
library(ggsignif)
picDir <- './tRNA亚型箱线图1/'
if(!dir.exists(picDir)){
  dir.create(picDir)
}
for (subtype in unique(count_trf$type)) {
  print(subtype)
  subset <- filter(count_trf,type == subtype)
  sample_means=as.data.frame(apply(subset[,-c(84,85)], 2, mean),colnames = colnames(subset))
  sample_means=cbind(sample_means,group$condition)
  data_long=melt(sample_means)
  data_long$variable <- subtype
  colnames(data_long)[1] <- 'group'
  max_y <- max(data_long$value)
  p <- ggplot(data_long, aes(x = group,y = value, fill = group)) +
    geom_boxplot() +
    scale_fill_brewer(palette = 'Dark2') + 
  geom_signif(comparisons = list(c('BEN','MAL')), # 使用group数据框中的condition列
              test = t.test,
              y_position = max_y -3, # 调整为你想要的y轴位置
              tip_length = c(c(0.05, 0.05)),
              size = 0.8,
              color = "black",
              position = position_identity())+
    xlab(unique(data_long$variable)) +
    ylab('means of RPM') +
    mytheme +
    theme(legend.position = 'right')+ annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "transparent", color = "black")
  ggsave(p,filename = paste0(picDir,subtype,'.pdf'),width = 4,height = 8)####ggplot2 自带的保存图片函数，不用创pdf
}



picDir <- './tRNA亚型片段箱线图1/'
if(!dir.exists(picDir)){
  dir.create(picDir)
}
library(reshape2)
library(ggplot2)
library(ggpubr)
for (subtype in unique(count_trf$type)) {
    subset <- filter(count_trf,type == subtype & condition == 'Internal')
    sample_means=as.data.frame(apply(subset[,-c(84,85)], 2, mean),colnames = colnames(subset))
    sample_means=cbind(sample_means,group$condition)
    data_long_i = melt(sample_means)
    data_long_i$variable <- subtype
    data_long_i$condition <- 'Internal'
    colnames(data_long_i)[1] <- 'group'
    subset1 <- filter(count_trf,type == subtype & condition == '5-tRF')
    sample_means1=as.data.frame(apply(subset1[,-c(84,85)], 2, mean),colnames = colnames(subset1))
    sample_means1=cbind(sample_means1,group$condition)
    data_long_t5t = melt(sample_means1)
    data_long_t5t$variable <- subtype
    data_long_t5t$condition <- '5-tRF'
    colnames(data_long_t5t)[1] <- 'group'
    subset2 <- filter(count_trf,type == subtype & condition == '5-half')
    sample_means2=as.data.frame(apply(subset2[,-c(84,85)], 2, mean),colnames = colnames(subset2))
    sample_means2=cbind(sample_means2,group$condition)
    data_long_t5h = melt(sample_means2)
    data_long_t5h$variable <- subtype
    data_long_t5h$condition <- '5-half'
    colnames(data_long_t5h)[1] <- 'group'
    subset3 <- filter(count_trf,type == subtype & condition == '3-half')
    sample_means3=as.data.frame(apply(subset3[,-c(84,85)], 2, mean),colnames = colnames(subset3))
    sample_means3=cbind(sample_means3,group$condition)
    data_long_t3h = melt(sample_means3)
    data_long_t3h$variable <- subtype
    data_long_t3h$condition <- '3-half'
    colnames(data_long_t3h)[1] <- 'group'
    subset4 <- filter(count_trf,type == subtype & condition == '3-tRF')
    sample_means4=as.data.frame(apply(subset4[,-c(84,85)], 2, mean),colnames = colnames(subset4))
    sample_means4=cbind(sample_means4,group$condition)
    data_long_t3t = melt(sample_means4)
    data_long_t3t$variable <- subtype
    data_long_t3t$condition <- '3-tRF'
    colnames(data_long_t3t)[1] <- 'group'
    data_long <- rbind(data_long_i,data_long_t5t,data_long_t5h,data_long_t3h,data_long_t3t)
    max_y <- max(data_long$value)
    p <- ggplot(data_long, aes(x = condition,y = value, fill = group)) +
      geom_boxplot() +
      scale_fill_brewer(palette = 'Dark2') +
      ggtitle(unique(data_long$variable)) +  geom_signif(comparisons = list(c('BEN','MAL')), # 使用group数据框中的condition列
                                                         test = t.test,
                                                         y_position = max_y -3, # 调整为你想要的y轴位置
                                                         tip_length = c(c(0.05, 0.05)),
                                                         size = 0.8,
                                                         color = "black",
                                                         position = position_identity())+
      xlab(NULL) +
      ylab('means of RPM') +
      mytheme +
      theme(legend.position = 'right',
            plot.title = element_text(hjust = 0.5, size = 14))+ annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "transparent", color = "black")
    ggsave(p,filename = paste0(picDir,subtype,'.pdf'),width = 8,height = 6)####ggplot2 自带的保存图片函数，不用创pdf
    
    
}

type_index <- 84
types <- count_trf[,type_index]
# 计算每种类型的总表达值
total_expression <- sapply(unique(types), function(type) {
  subset <- count_trf[types == type, -type_index, drop = FALSE]
  total <- colSums(subset)
  return(total)
})
total_expression <- as.data.frame(total_expression)
total_expression$all <- rowSums(total_expression)
for (i in 1:(ncol(total_expression)-1)) {
  total_expression[,i] <- (total_expression[,i] / total_expression$all) * 100
}
ts_long <- melt(total_expression[,-22])

#colors <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", 
            #"#00FFFF", "#FFA500", "#800080", "#008080", "#FFC0CB", 
            #"#800000", "#008000", "#000080", "#FFD700", "#A52A2A", 
            #"#808000", "#ADD8E6", "#20B2AA", "#F08080", "#4682B4", 
            #"#808080")
ts_long$group <- rep(group$condition,21)

p <- ggplot(ts_long, aes(x = variable, y = value, fill = group)) +
  geom_boxplot() +
  #scale_fill_manual(values = c('#00FF7F','#DA70D6')) +  # 手动设置颜色
  scale_fill_brewer(palette = 'Accent')+
  #geom_hline(yintercept = mean(result_df$proportion), linetype = "dashed") +  
  xlab(NULL) +
  ylab('Proportion') +
  ggtitle('tRNAs') +
  mytheme +
  theme(legend.position = 'right',
        plot.title = element_text(hjust = 0.5, size = 14)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "transparent", color = "black")
p
ggsave(p,filename = 'tRNA各种亚型占比图.pdf',width = 8,height = 6)


#####################tsRNA相关性热图#################
count_trf <- aggregate(x=count_trf[,-c(84,85)],by=list(count_trf$type),FUN = mean)#计算亚型表达矩阵
colnames(count_trf)[1] <- 'type'
rownames(count_trf) <- count_trf$type
count_trf <- as.data.frame(t(count_trf[,-1]))
ct <- round(cor(count_trf),3)
library(psych)
library(reshape)
library(ggplot2)
library(tidyverse)
library(ggcorrplot)
data <- as.data.frame(ct) %>%
  mutate(x=rownames(ct)) %>%
  melt(id='x') %>%
  rename('y'='variable','Corr'='value')
data <- data %>% group_by(variable) %>% arrange(variable,value)

list = rownames(ct)
list = factor(list,levels = list)

list = unique(data$x)
list = factor(list,levels = list)

ggplot(data, aes(x = fct_inorder(x), y = factor(variable,levels = list), fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red',
                       limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  labs(x = NULL, y = NULL) +
  mytheme

ggplot(data,aes(factor(x,levels = list),
                factor(variable,levels = list), #定义x，y轴顺序，防止被默认改变
                fill=value))+  #根据相关性值填充颜色
  geom_tile()+  #色块函数
  scale_fill_gradient2(low = '#FF4500',mid = '#FFA500',high ='white',
                       limits=c(-1,1),breaks=c(-1,-0.5,0,0.5,1))+
  labs(x=NULL,y=NULL)+
  mytheme

ggcorrplot(ct)+mytheme








