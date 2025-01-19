setwd('F:\\R_script\\20230919 乳腺癌\\expr\\new1')
tRNA <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\tRFs.csv',row.names = 1)
rsRNA <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\rsRNA.csv',row.names = 1)
ysRNA <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\ysRNA.csv',row.names = 1)
mRNA <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\mRNA.csv',row.names = 1)
miRNA <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\miRNA.csv',row.names = 1)
lncRNA <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\lncRNA.csv',row.names = 1)
piRNA <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\piRNA.csv',row.names = 1)

library(readxl)
group <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\simple.csv')
group <- group[-c(82,83),]
#ids <- group[match(group$X,colnames(data)),]
ben <- group[group$condition == 'BEN',]
mal <- group[group$condition == 'MAL',]
ids <- rbind(mal,ben)
#write.table(ids,file = 'malvsben.csv')
##############过滤掉低表达基因
filter_genes <- function(expression_matrix, ids, threshold = 0.75) {
  # Reorder the columns of the expression matrix to match the order in ids
  expression_matrix <- expression_matrix[, match(ids$X, colnames(expression_matrix))]
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
rsRNA <- filter_genes(rsRNA,ids = ids)
tRNA <- filter_genes(tRNA,ids)
ysRNA <- filter_genes(ysRNA,ids)
mRNA <- filter_genes(mRNA,ids)
miRNA <- filter_genes(miRNA,ids)
lncRNA <- filter_genes(lncRNA,ids)
piRNA <- filter_genes(piRNA,ids)

write.table(rsRNA,file = 'filter_rsRNA.csv',sep = ',',col.names = NA)
write.table(ysRNA,file = 'filter_ysRNA.csv',sep = ',',col.names = NA)
write.table(tRNA,file = 'filter_tRNA.csv',sep = ',',col.names = NA)
write.table(mRNA,file = 'filter_mRNA.csv',sep = ',',col.names = NA)
write.table(miRNA,file = 'filter_miRNA.csv',sep = ',',col.names = NA)
write.table(lncRNA,file = 'filter_lncRNA.csv',sep = ',',col.names = NA)
write.table(piRNA,file = 'filter_piRNA.csv',sep = ',',col.names = NA)

##################DESeq2差异分析
library(DESeq2)
library(ggplot2)
ids2 <- ids
ids2 <- as.data.frame(ids2[,-1])
rownames(ids2) <- ids$X
colnames(ids2) <- 'condition'
ids2$condition <- as.factor(ids2$condition)
library(DESeq2)
######数据打包成DESeq格式的数据集 基因名在行名，没有数字列
dds <- DESeqDataSetFromMatrix(countData = ysRNA, colData = ids2, design= ~condition)
######过滤掉低表达的counts值，count函数过滤
#dds <- dds[rowSums(counts(dds))>=0,]###可以定义为1
########## 差异分析
dds <- DESeq(dds)
######### 构建contrast对象，用于后续的差异结果的提取，谁是肿瘤/普通
table(ids$condition)
contrast <- c('condition','MAL','BEN')
res <- results(dds,contrast = contrast)
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
library(tidyverse)
res1 <- res1 %>% arrange(res1$pvalue)


deg1 <- res1
colnames(deg1)[2] <- 'logFC'
deg1$logP <- -log10(deg1$pvalue)
deg1$group <- 'NS'
deg1$group[which(deg1$logFC >= log2(2) & deg1$pvalue <0.05)] <- 'up'
deg1$group[which(deg1$logFC <= -log2(2) & deg1$pvalue < 0.05)] <- 'down'
table(deg1$group)

#####新增一列找出前10基因
library(tidyverse)
deg1$Lable <- ''
deg1 <- arrange(deg1,desc(abs(logFC)))##排序
deg1$id <- rownames(deg1)
up <- head(deg1$id[which(deg1$group == 'up')],25)#提取p值最小的前10个
down <- head(deg1$id[which(deg1$group == 'down')],25)
top10_gene <- c(as.character(up),as.character(down))##合并
deg1$Lable[match(top10_gene,deg1$id)] <- top10_gene
write.table(deg1, 'MALvsBEN_ysRNA.DESeq2.csv',sep = ',',col.names = NA)

lnc_deg <- deg1
mi_deg <- deg1
pi_deg <- deg1
m_deg <- deg1
t_deg <- deg1
rs_deg <- deg1
ys_deg <- deg1

lnc_deg$group2 <- 'lncRNA'
mi_deg$group2 <- 'miRNA'
m_deg$group2 <- 'mRNA'
pi_deg$group2 <- 'piRNA'
rs_deg$group2 <- 'rsRNA'
ys_deg$group2 <- 'ysRNA'
t_deg$group2 <- 'tRNA'
all2 <- rbind(lnc_deg,mi_deg,m_deg,pi_deg,rs_deg,ys_deg,t_deg)

library(ggplot2)
library(RColorBrewer)
library(grid)
library(scales)
library(dplyr)

all2$group2 <- as.factor(all2$group2)
all2$group <- as.factor(all2$group)
all2 <- na.omit(all2)
all_bg <- all2 %>% group_by(group2) %>% summarise(max_log2FC = max(all2$logFC),min_log2FC = min(all2$logFC))
#setwd('F:\\R_script\\20230919 乳腺癌\\expr\\分期new')
pdf('F:\\R_script\\20230919 乳腺癌\\expr\\new1\\差异火山图.pdf',width = 8,height = 6) 
p <- ggplot() +
  geom_col(data = all_bg,mapping = aes(group2,max_log2FC),
           fill = 'grey85',width = 0.8,alpha = 0.5)+
  geom_col(data = all_bg,mapping = aes(group2,min_log2FC),
           fill = 'grey85',width = 0.8,alpha = 0.5)
p1 <- p + geom_jitter(data = all2,mapping = aes(x = group2,y = logFC,color = group),
                      size = 2,width = 0.4,alpha = 0.7)

p1 <- p1 +geom_col(data = all_bg,mapping = aes(x = group2,y = 0.3,fill = group2))+
  geom_col(data = all_bg,mapping = aes(x = group2,y = -0.3,fill = group2))

p1 <- p1 + geom_text(data = all_bg,mapping = aes(x = group2,y = 0, label = group2),
                     size = 4,color = '#dbebfa') 
p1 + theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.8),
        axis.text.y = element_text(size = 12,color = 'black'),
        axis.title = element_text(size = 14,color = 'black'),
        axis.ticks.y = element_line(linewidth = 0.8)) +
  labs(x = 'type',y='log2FoldChange',fill = NULL,color = NULL)+
  guides(color = guide_legend(override.aes = list(size = 6,alpha = 1)))
dev.off()

###############rpm标准化
normalize_to_RPM_optimized <- function(expression_matrix) {
  # Calculate the sum of reads for each column (sample)
  col_sums <- colSums(expression_matrix)
  # Normalize each value to RPM
  for (i in 1:ncol(expression_matrix)) {
    expression_matrix[,i] <- (expression_matrix[,i] / col_sums[i]) * 1e6
  }
  return(expression_matrix)
}

lnc_rpm <- normalize_to_RPM_optimized(lncRNA)
m_rpm <- normalize_to_RPM_optimized(mRNA)
mi_rpm <- normalize_to_RPM_optimized(miRNA)
t_rpm <- normalize_to_RPM_optimized(tRNA)
rs_rpm <- normalize_to_RPM_optimized(rsRNA)
ys_rpm <- normalize_to_RPM_optimized(ysRNA)
pi_rpm <- normalize_to_RPM_optimized(piRNA)

a = t_deg[t_deg$group =='up'|t_deg$group == 'down',]
b = t_rpm[a$id,]
b = na.omit(b)
b = log2(b+1)

##############分期test###############
stage <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\group2.csv',row.names = 1)
s1 <- stage[stage$stage == 'I',]
s2 <- stage[stage$stage == 'II',]
s3 <- stage[stage$stage == 'III',]
s4 <- stage[stage$stage == 'IV',]
missing <- stage[stage$stage == 'Missing',]
ben <- stage[stage$stage == 'BEN',]
stage <- rbind(s1,s2,s3,s4,missing,ben)
#write.table(stage,file = 'stage.csv',sep = ',',col.names = NA)
b <- b[,match(rownames(stage),colnames(b))]
g <- as.data.frame(stage$stage)
rownames(g) <- rownames(stage)
colnames(g) <- 'stage'



##############热图
library(pheatmap)
p <- pheatmap(b, 
              #annotation_row=dfGene, # （可选）指定行分组文件
              annotation_col=g, # （可选）指定列分组文件
              #annotation_row = annotation_row,
              show_colnames = F, # 是否显示列名
              show_rownames = F,
              # 是否显示行名
              fontsize=9, # 字体大小
              color = colorRampPalette(c('navy','#ffffff','#d53e4f'))(50), # 指定热图的颜色
              annotation_legend=TRUE, # 是否显示图例
              breaks = seq(-2,2,length.out =50),
              #border_color='black',  # 边框颜色 NA表示没有
              scale = 'row', # 指定归一化的方式。"row"按行归一化，"column"按列归一化，"none"不处理
              cluster_rows = T, # 是否对行聚类
              cluster_cols = F ,# 是否对列聚类
              main = 'heatmap of differential tRNA',
              border = F,
              #clustering_distance_rows = "euclidean", # 聚类距离方法
              #clustering_method = "ward.D2" # 聚类方法
)
p
ggsave(p,filename = 'F:\\R_script\\20230919 乳腺癌\\expr\\new1\\tRNA-heatmap.pdf',width = 8,height = 6)


ys_lable <- ys_rpm[match(ys_deg$Lable,rownames(ys_rpm)),]
ys_lable <- na.omit(ys_lable)
rs_lable <- rs_rpm[match(rs_deg$Lable,rownames(rs_rpm)),]
rs_lable <- na.omit(rs_lable)
t_lable <- t_rpm[match(t_deg$Lable,rownames(t_rpm)),]
t_lable <- na.omit(t_lable)
m_lable <- m_rpm[match(m_deg$Lable,rownames(m_rpm)),]
m_lable <- na.omit(m_lable)
mi_lable <- mi_rpm[match(mi_deg$Lable,rownames(mi_rpm)),]
mi_lable <- na.omit(mi_lable)
pi_lable <- pi_rpm[match(pi_deg$Lable,rownames(pi_rpm)),]
pi_lable <- na.omit(pi_lable)
lnc_lable <- lnc_rpm[match(lnc_deg$Lable,rownames(lnc_rpm)),]
lnc_lable <- na.omit(lnc_lable)

lnc_lable$type <- 'lncRNA'
mi_lable$type <- 'miRNA'
m_lable$type <- 'mRNA'
pi_lable$type <- 'piRNA'
rs_lable$type <- 'rsRNA'
ys_lable$type <- 'ysRNA'
t_lable$type <- 'tRNA'

all <- rbind(lnc_lable[,-84],m_lable[,-84],mi_lable[,-84],pi_lable[,-84],ys_lable[,-84],rs_lable[,-84],t_lable[,-84])
type <- data.frame(Gene=rownames(lnc_lable),type = lnc_lable$type)
type1 <- data.frame(Gene=rownames(m_lable),type = m_lable$type)
type2 <- data.frame(Gene=rownames(mi_lable),type = mi_lable$type)
type3 <- data.frame(Gene=rownames(pi_lable),type = pi_lable$type)
type4 <- data.frame(Gene=rownames(ys_lable),type = ys_lable$type)
type5 <- data.frame(Gene=rownames(rs_lable),type = rs_lable$type)
type6 <- data.frame(Gene=rownames(t_lable),type = t_lable$type)
annotation_row=rbind(type,type1,type2,type3,type4,type5,type6)
annotation_row=as.data.frame(annotation_row[,-1])
colnames(annotation_row) <- 'type'
rownames(annotation_row) <- rownames(all)
all1 <- cbind(ids2,as.data.frame(t(all)))
all1 <- t(all1[,-c(1,2)])

color_gradient <- colorRampPalette(c('navy', '#ffffff', '#d53e4f'))(50)
col_annotation_colors <- list(
  ids2 = c("MAL" = "blue", "BEN" = "red")
)
all1 <- all1[,match(rownames(stage),colnames(all1))]
all1 <- log2(all1+1)


P <- pheatmap(all1, 
              #annotation_row=dfGene, # （可选）指定行分组文件
              annotation_col=g, # （可选）指定列分组文件
              annotation_row = annotation_row,
              show_colnames = F, # 是否显示列名
              show_rownames = F,
              color = color_gradient,
              # 是否显示行名
              fontsize=9, # 字体大小
              # 指定热图的颜色
              annotation_legend=TRUE, # 是否显示图例
              border_color='black',  # 边框颜色 NA表示没有
              scale = 'row', # 指定归一化的方式。"row"按行归一化，"column"按列归一化，"none"不处理
              cluster_rows = F, # 是否对行聚类
              cluster_cols = F ,# 是否对列聚类
              main = 'heatmap of RNA',method = 'average',gaps_row = c(21,71,121,171,221,271),
              breaks = seq(-2,2,length.out =50),
              annotation_colors = col_annotation_colors
              
)
P
library(ggplot2)
ggsave(P,filename = 'F:\\R_script\\20230919 乳腺癌\\expr\\new1\\RNA-heatmap.pdf',width = 8,height = 6)

a = lnc_deg[lnc_deg$group =='up'|lnc_deg$group == 'down',]
b = lnc_rpm[a$id,]
write.table(b,file = 'deg_lncRNA.csv',sep = ',',col.names = NA)





