library(data.table)
################ TCGA ####################
ano <- fread('./TCGA/gencode.v22.annotation.gene.probeMap')
ano = ano[!duplicated(ano$gene),]

data <- as.data.frame(fread('./TCGA/TCGA-BRCA.htseq_fpkm.tsv.gz'))
rownames(data) <- data$Ensembl_ID
data <- data[,-1]
gene <- intersect(rownames(data),ano$id)
data <- data[gene,]
identical(rownames(data),ano$id)
data <- data[match(ano$id,rownames(data)),]
rownames(data) = ano$gene

sur <- as.data.frame(fread('./TCGA/TCGA-BRCA.survival.tsv'))
overlap <- intersect(colnames(data),sur$sample)
data <- data[,overlap]
rownames(sur) <- sur$sample
sur <- sur[overlap,]


group <- as.data.frame(colnames(data))
library(tidyverse)
colnames(group) <- 'id'
group <- group %>% mutate(group = ifelse(grepl('0',sur$OS),'Alive',ifelse(grepl('1',id),'Dead',NA)))
group <- na.omit(group)
table(group$group)
gene <- c('DLST','FBXO31')

data_gene <- data[gene,]
data_gene <- data_gene[,match(group$id,colnames(data_gene))]

for (i in (1:length(gene))) {
  print(i)
  data_box = data_gene[i,]
  gene_box = gene[i]
  boxplot_gene <- function(data,gene,group){
    library(reshape2)
    data <- melt(data = data)
    data$group <- group$group
    colnames(data) <- c('id','exp','Group')
    data = data[,-1]
    #data1$exp <- log2(data1$exp + 1)
    library(ggbeeswarm)
    library(ggpubr)
    p <- ggplot(data = data,aes(x = Group, y = exp, color = Group)) + ## aes用于设置x和y轴分别绘制什么,color用于给每一组分配一个颜色
      stat_boxplot(geom = "errorbar", width = 0.2, lwd = 1.0) + ## 绘制上下须子，加在geom_boxplot前，否则线会盖在箱线图上面
      geom_boxplot(width = 0.4, outlier.colour = "black", lwd =1.0,outlier.shape = NA) + ## 绘制箱线图，离群点为黑色，也可使用outlier.shape = NA来不画离群点
      #geom_beeswarm(size = 2, priority = "density") + ##绘制蜜蜂群图，这里其实是散点图
      scale_color_manual(values = c("#AF0F11", "#3372A6")) +
      theme_minimal() +
      theme(panel.grid = element_blank(),  # 去除所有网格线
            axis.line = element_line(colour = "black", size = 1),  # 设置坐标轴线的颜色和宽度
            axis.ticks = element_line(colour = "black", size = 0.5),  # 设置坐标轴刻度线的颜色和宽度
            axis.title = element_text(size = 14, face = "bold"),  # 加粗并设置坐标轴标题字号
            axis.text = element_text(size = 12, face = "bold"),  # 加粗坐标轴刻度标签
            legend.title = element_text(size = 14, face = "bold"),  # 加粗图例标题
            legend.text = element_text(size = 12, face = "bold"),   # 加粗图例项
            plot.title = element_text(size = 16, face = "bold")) + # 加粗并设置图形标题字号
      labs(x = gene, y = "log2(Exp)", title = "") +
    stat_compare_means(comparisons = list(c("Dead", "Alive")), label = "p.signif",method = 'wilcox',exact = F)
    p
    ggsave(plot = p,filename = paste('.',paste(gene,'pdf',sep = '.'),sep = '/'),width = 8,height = 4)
  }
  boxplot_gene(data = data_box,group = group,gene = gene_box)
}


############## plasma #########################
setwd('D:\\乳腺癌\\new1')
plasma_data <- read.csv('filter_mRNA_rpm.csv',row.names = 1)
group <- read.csv('./stage.csv')
plasma_data <- plasma_data[,match(group$X,colnames(plasma_data))]
identical(colnames(plasma_data),group$X)
data_gene <- plasma_data[gene,]


for (i in (1:length(gene))) {
  print(i)
  data_box = data_gene[i,]
  gene_box = gene[i]
  boxplot_gene <- function(data,gene,group){
    library(reshape2)
    data <- melt(data = data)
    data$group <- group$group
    colnames(data) <- c('id','exp','Group')
    data = data[,-1]
    data$exp <- log2(data$exp + 1)
    library(ggbeeswarm)
    p <- ggplot(data = data,aes(x = Group, y = exp, color = Group)) + ## aes用于设置x和y轴分别绘制什么,color用于给每一组分配一个颜色
      stat_boxplot(geom = "errorbar", width = 0.2, lwd = 1.0) + ## 绘制上下须子，加在geom_boxplot前，否则线会盖在箱线图上面
      geom_boxplot(width = 0.4, outlier.colour = "black", lwd =1.0,outlier.shape = NA) + ## 绘制箱线图，离群点为黑色，也可使用outlier.shape = NA来不画离群点
      #geom_beeswarm(size = 2, priority = "density") + ##绘制蜜蜂群图，这里其实是散点图
      scale_color_manual(values = c("#AF0F11", "#3372A6")) +
      theme_minimal() +
      theme(panel.grid = element_blank(),  # 去除所有网格线
            axis.line = element_line(colour = "black", size = 1),  # 设置坐标轴线的颜色和宽度
            axis.ticks = element_line(colour = "black", size = 0.5),  # 设置坐标轴刻度线的颜色和宽度
            axis.title = element_text(size = 14, face = "bold"),  # 加粗并设置坐标轴标题字号
            axis.text = element_text(size = 12, face = "bold"),  # 加粗坐标轴刻度标签
            legend.title = element_text(size = 14, face = "bold"),  # 加粗图例标题
            legend.text = element_text(size = 12, face = "bold"),   # 加粗图例项
            plot.title = element_text(size = 16, face = "bold")) + # 加粗并设置图形标题字号
      labs(x = gene, y = "log2(Exp)", title = "") +
    stat_compare_means(comparisons = list(c("MAL", "BEN")), label = "p.signif",method = 'wilcox',exact = F)
    ggsave(plot = p,filename = paste('.',paste(gene,'pdf',sep = '.'),sep = '/'),width = 8,height = 4)
  }
  boxplot_gene(data = data_box,group = group,gene = gene_box)
}

############## tissue ##################
setwd('D:\\乳腺癌\\tissue')
library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(readxl)
library(ggsci)
group <- read_xlsx('tissue_ID.xlsx')
group <- group %>% filter(sample_type == 'MAL' | sample_type == 'BEN')
data <- read.csv('./tissue_mRNA.csv',row.names = 1,check.names = F)
gene <- c('DLST','FBXO31')
data <- data[,match(group$seq_ID,colnames(data))]
filter_genes <- function(expression_matrix,threshold = 0.5) {
  # Replace NA values with 0
  expression_matrix[is.na(expression_matrix)] <- 0
  # Calculate the row sums of zeros
  row_sum <- rowSums(expression_matrix == 0)
  # Calculate the threshold number of columns
  n <- threshold * ncol(expression_matrix)
  # Filter the rows based on the threshold
  filtered_matrix <- expression_matrix[row_sum <= (ncol(expression_matrix) - n), ]
  return(filtered_matrix)
}
#data <- filter_genes(data,threshold = 0.5)
data_gene <- data[gene,]
identical(group$seq_ID,colnames(data_gene))

library(reshape2)
for (i in (1:length(gene))) {
  print(i)
  data_box = data_gene[i,]
  gene_box = gene[i]
  boxplot_gene <- function(data,gene,group){
    library(reshape2)
    data <- melt(data = data)
    data$group <- group$sample_type
    colnames(data) <- c('id','exp','Group')
    data = data[,-1]
    data$exp <- log2(data$exp + 1)
    library(ggbeeswarm)
    p <- ggplot(data = data,aes(x = Group, y = exp, color = Group)) + ## aes用于设置x和y轴分别绘制什么,color用于给每一组分配一个颜色
      stat_boxplot(geom = "errorbar", width = 0.2, lwd = 1.0) + ## 绘制上下须子，加在geom_boxplot前，否则线会盖在箱线图上面
      geom_boxplot(width = 0.4, outlier.colour = "black", lwd =1.0,outlier.shape = NA) + ## 绘制箱线图，离群点为黑色，也可使用outlier.shape = NA来不画离群点
      #geom_beeswarm(size = 2, priority = "density") + ##绘制蜜蜂群图，这里其实是散点图
      scale_color_manual(values = c("#AF0F11", "#3372A6")) +
      theme_minimal() +
      theme(panel.grid = element_blank(),  # 去除所有网格线
            axis.line = element_line(colour = "black", size = 1),  # 设置坐标轴线的颜色和宽度
            axis.ticks = element_line(colour = "black", size = 0.5),  # 设置坐标轴刻度线的颜色和宽度
            axis.title = element_text(size = 14, face = "bold"),  # 加粗并设置坐标轴标题字号
            axis.text = element_text(size = 12, face = "bold"),  # 加粗坐标轴刻度标签
            legend.title = element_text(size = 14, face = "bold"),  # 加粗图例标题
            legend.text = element_text(size = 12, face = "bold"),   # 加粗图例项
            plot.title = element_text(size = 16, face = "bold")) + # 加粗并设置图形标题字号
      labs(x = gene, y = "log2(Exp)", title = "") +
    stat_compare_means(comparisons = list(c("MAL", "BEN")), label = "p.signif",method = 'wilcox',exact = F)
    ggsave(plot = p,filename = paste('.',paste(paste(gene,'new',sep = '_'),'pdf',sep = '.'),sep = '/'),width = 8,height = 4)
  }
  boxplot_gene(data = data_box,group = group,gene = gene_box)
}


