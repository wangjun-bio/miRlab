setwd('F:\\R_script\\20230919 乳腺癌\\expr')
miRNA = read.csv('miRNA.csv',row.names = 1,header = T)
miRNA[is.na(miRNA)] <- 0
miRNA=as.data.frame(miRNA[,match(group$X,colnames(miRNA))])
total_miRNA <- colSums(miRNA)
tRNA <- read.csv('tRFs.csv',row.names = 1,header = T)
tRNA[is.na(tRNA)] <- 0
tRNA=as.data.frame(tRNA[,match(group$X,colnames(tRNA))])
total_tRNA <- colSums(tRNA)
rsRNA <- read.csv('rsRNA.csv',row.names = 1,header = T)
rsRNA[is.na(rsRNA)] <- 0
rsRNA=as.data.frame(rsRNA[,match(group$X,colnames(rsRNA))])
total_rsRNA <- colSums(rsRNA)
ysRNA <- read.csv('ysRNA.csv',row.names = 1,header = T)
ysRNA[is.na(ysRNA)] <- 0
ysRNA=as.data.frame(ysRNA[,match(group$X,colnames(ysRNA))])
total_ysRNA <- colSums(ysRNA)

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
a=a[,-1]
a$all <- rowSums(a)
proportion_miRNA <- as.data.frame((total_miRNA / a$all) * 100) 
proportion_tRNA <- as.data.frame((total_tRNA / a$all) * 100) 
proportion_rsRNA <- as.data.frame((total_rsRNA / a$all) * 100) 
proportion_ysRNA <- as.data.frame((total_ysRNA / a$all) * 100) 

td <- cbind(proportion_miRNA,proportion_rsRNA,proportion_tRNA,proportion_ysRNA)
colnames(td) <- c('miRNA','rsRNA','tRNA','ysRNA')
cor <- cor(td$miRNA,td$rsRNA)

library(ggplot2)
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
library(ggpubr)
p <- ggscatter(td, x = "miRNA", y = "rsRNA", 
          color = "blue3",fill = "lightgray",
          add = "reg.line", conf.int = TRUE, 
          add.params = list( color = "black",fill = "lightgray",fill = "lightgray"),
          cor.coef = T,
          cor.coeff.args = list(),
          cor.method = "pearson",
          cor.coef.coord = c(0.5, 0.025),
          cor.coef.size = 8)+xlab('Proportion of miRNAs')+ylab('Proportion of rsRNAs')+mytheme+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "transparent", color = "black")
ggsave(p,filename = 'miRNA-rsRNA占比相关性.pdf',width = 6,height = 10)


p <- ggscatter(td, x = "miRNA", y = "tRNA", 
               color = "orange",fill = "lightgray",
               add = "reg.line", conf.int = TRUE, 
               add.params = list( color = "black",fill = "lightgray",fill = "lightgray"),
               cor.coef = T,
               cor.coeff.args = list(),
               cor.method = "pearson",
               cor.coef.coord = c(0.5, 0.025),
               cor.coef.size = 8)+xlab('Proportion of miRNAs')+ylab('Proportion of tRNAs')+mytheme+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "transparent", color = "black")
p
ggsave(p,filename = 'miRNA-tRNA占比相关性.pdf',width = 6,height = 10)

p <- ggscatter(td, x = "miRNA", y = "ysRNA", 
               color = "#7FFFAA",fill = "lightgray",
               add = "reg.line", conf.int = TRUE, 
               add.params = list( color = "black",fill = "lightgray",fill = "lightgray"),
               cor.coef = T,
               cor.coeff.args = list(),
               cor.method = "pearson",
               cor.coef.coord = c(0.5, 0.025),
               cor.coef.size = 8)+xlab('Proportion of miRNAs')+ylab('Proportion of ysRNAs')+mytheme+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "transparent", color = "black")
p
ggsave(p,filename = 'miRNA-ysRNA占比相关性.pdf',width = 6,height = 10)








































