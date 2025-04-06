getwd()
setwd('/Users/willow/Desktop/2025/cfRNA乳腺癌/DATA-cfRNA')
tRNA <- read.csv('./表达矩阵/tRFs.csv',row.names = 1)
rsRNA <- read.csv('./表达矩阵/rsRNA.csv',row.names = 1)
ysRNA <- read.csv('./表达矩阵/ysRNA.csv',row.names = 1)
mRNA <- read.csv('./表达矩阵/mRNA.csv',row.names = 1)
miRNA <- read.csv('./表达矩阵/miRNA.csv',row.names = 1)
lncRNA <- read.csv('./表达矩阵/lncRNA.csv',row.names = 1)
piRNA <- read.csv('./表达矩阵/piRNA.csv',row.names = 1)

library(readxl)
group <- read.csv('./simple.csv')
group <- group[-c(82,83),]
#ids <- group[match(group$X,colnames(data)),]
ben <- group[group$condition == 'BEN',]
mal <- group[group$condition == 'MAL',]
ids <- rbind(mal,ben)
#write.table(ids,file = 'malvsben.csv')
##############过滤掉低表达基因############
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
rsRNA <- filter_genes(rsRNA,ids)
tRNA <- filter_genes(tRNA,ids)
ysRNA <- filter_genes(ysRNA,ids)
mRNA <- filter_genes(mRNA,ids)
miRNA <- filter_genes(miRNA,ids)
lncRNA <- filter_genes(lncRNA,ids)
piRNA <- filter_genes(piRNA,ids)

# write.table(rsRNA,file = 'filter_rsRNA.csv',sep = ',',col.names = NA)
# write.table(ysRNA,file = 'filter_ysRNA.csv',sep = ',',col.names = NA)
# write.table(tRNA,file = 'filter_tRNA.csv',sep = ',',col.names = NA)
# write.table(mRNA,file = 'filter_mRNA.csv',sep = ',',col.names = NA)
# write.table(miRNA,file = 'filter_miRNA.csv',sep = ',',col.names = NA)
# write.table(lncRNA,file = 'filter_lncRNA.csv',sep = ',',col.names = NA)
# write.table(piRNA,file = 'filter_piRNA.csv',sep = ',',col.names = NA)

##################DESeq2差异分析############
library(DESeq2)
library(ggplot2)
ids2 <- ids
ids2 <- as.data.frame(ids2[,-1])
rownames(ids2) <- ids$X
colnames(ids2) <- 'condition'
ids2$condition <- as.factor(ids2$condition)
######数据打包成DESeq格式的数据集 基因名在行名，没有数字列############
dds_zip <- function(expression_matrix){
  dds <- DESeqDataSetFromMatrix(countData = expression_matrix, colData = ids2, design= ~condition)
  ######过滤掉低表达的counts值，count函数过滤############
  #dds <- dds[rowSums(counts(dds))>=0,]###可以定义为1
  ########## 差异分析############
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
  
  #####新增一列找出前25基因############
  library(tidyverse)
  deg1$Lable <- ''
  deg1 <- arrange(deg1,desc(abs(logFC)))##排序
  deg1$id <- rownames(deg1)
  up <- head(deg1$id[which(deg1$group == 'up')],25)#提取p值最小的前25个
  down <- head(deg1$id[which(deg1$group == 'down')],25)
  top10_gene <- c(as.character(up),as.character(down))##合并
  deg1$Lable[match(top10_gene,deg1$id)] <- top10_gene
  return(deg1)
  #write.table(deg1, 'MALvsBEN_.DESeq2.csv',sep = ',',col.names = NA)
}

lnc_deg <- dds_zip(lncRNA)
mi_deg <- dds_zip(miRNA)
pi_deg <- dds_zip(piRNA)
m_deg <- dds_zip(mRNA)
t_deg <- dds_zip(tRNA)
rs_deg <- dds_zip(rsRNA)
ys_deg <- dds_zip(ysRNA)

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

##########差异火山图###########
# pdf('./figure/fig2/差异火山图.pdf',width = 8,height = 6) 
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

p
# dev.off()

####----plot火山图-2----####
library(ggplot2)
library(dplyr)

# 假设 all2 数据框已经存在，其中包含以下主要变量：
# group: 数据的分组（例如 "NS", "up", "down"）
# group2: RNA类型或其它分组，用作 x 轴及底部矩形
# logFC: 对数 fold change
# Lable: 文本标签（如果需要显示，可使用，但这里不再针对 top 数据做额外处理）

# 1. 对 all2 数据按照 group2 分组处理：如果某个 RNA 类型数据点超过 5000，则随机保留50%
plot_data <- all2 %>%
  group_by(group2) %>%
  do(if(nrow(.) > 5000) sample_frac(., 0.5) else .) %>%
  ungroup()

# 2. 构建颜色梯度（如散点用的颜色，此处与之前不同，可根据需要调整）
# 如果需要使用连续调色板可如下构造（这里作为例子）
# color_gradient <- colorRampPalette(c("blue", "white", "red"))(50)
# 但下文我们主要使用手动设置散点颜色 scale_color_manual() 与 scale_fill_manual()

# 3. 创建 ggplot 对象，只绘制两组抖动点（不再额外处理每组中 logFC 最大的数据）
p <- ggplot() +
  
  # 针对 group 为 "NS" 的数据
  geom_jitter(data = plot_data %>% filter(group == "NS"),
              aes(x = group2, y = logFC, 
                  color = group, size = abs(logFC),
                  alpha = abs(logFC)),
              width = 0.4) +
  
  # 针对非 "NS" 的数据
  geom_jitter(data = plot_data %>% filter(group != "NS"),
              aes(x = group2, y = logFC, 
                  color = group, size = abs(logFC),
                  alpha = abs(logFC)),
              width = 0.4) +
  
  # 底部绘制矩形，颜色根据 group2 填充
  geom_tile(data = plot_data,
            mapping = aes(x = group2, y = 0, fill = group2),
            height = 0.4) +
  
  # 在矩形上添加分组标签（仅显示每个 group2 一次）
  geom_text(data = plot_data %>% distinct(group2, .keep_all = TRUE),
            aes(x = group2, y = 0, label = group2),
            size = 6) +
  
  # 设置 y 轴范围
  scale_y_continuous(limits = c(-2, 3)) +
  
  # 设置散点大小与透明度范围
  scale_size(range = c(0.5, 4)) +
  scale_alpha(range = c(0.1, 1)) +
  
  # 手动设置散点颜色
  scale_color_manual(values = c("up" = "#f46d43",
                                "NS" = "#bdbdbd",
                                "down" = "#3288bd")) +
  
  # 手动设置底部矩形的填充颜色（根据 group2 分组，颜色可自定义）
  scale_fill_manual(values = c('#d4c953','#93c665','#aea7c6','#7eace0','#3673b1','#849989','#d8744f')) +

  # 使用黑白主题，并自定义其它主题设置
  theme_bw() +
  theme(
    axis.text = element_text(color = "#000000", size = 12),
    axis.title = element_text(color = "#000000", size = 15),
    panel.grid = element_blank(),
    legend.background = element_rect(color = "#969696"),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank()
  ) +
  labs(title = "heatmap of RNA")

p

# 保存绘制好的图形为 PDF 文件
ggsave(filename = "./figure/fig2/差异火山图.pdf", plot = p, height = 9, width = 12)


###############rpm标准化############
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

##############分期############
stage <- read.csv('./group2.csv',row.names = 1)
s1 <- stage[stage$stage == 'I',]
s2 <- stage[stage$stage == 'II',]
s3 <- stage[stage$stage == 'III',]
s4 <- stage[stage$stage == 'IV',]
missing <- stage[stage$stage == 'Missing',]
ben <- stage[stage$stage == 'BEN',]
stage <- rbind(s1,s2,s3,s4,missing,ben)
g <- as.data.frame(stage$stage)
rownames(g) <- rownames(stage)
colnames(g) <- 'stage'
#############RNA单独绘图###########
Group_draw <- function(deg_matrix, rpm_matrix, type){
  a = deg_matrix[deg_matrix$group =='up'|deg_matrix$group == 'down',]
  b = rpm_matrix[a$id,]
  b = na.omit(b)
  b = log2(b+1)

  #write.table(stage,file = 'stage.csv',sep = ',',col.names = NA)
  b <- b[,match(rownames(stage),colnames(b))]
  col_groups <- as.character(stage[,1])
  gaps_col <- which(col_groups[-length(col_groups)] != col_groups[-1])
  ##############热图##########
  library(pheatmap)
  p <- pheatmap(b, 
                #annotation_row=dfGene, # （可选）指定行分组文件
                annotation_col=g, # （可选）指定列分组文件
                #annotation_row = annotation_row,
                show_colnames = F, # 是否显示列名
                show_rownames = F,
                # 是否显示行名
                fontsize=9, # 字体大小
                color = colorRampPalette(c('#4c71d0','#ffffff','#d50000'))(50), # 指定热图的颜色
                annotation_legend=TRUE, # 是否显示图例
                breaks = seq(-3,3,length.out =50),
                border_color = 'gray',      # 单元格边框颜色
                scale = 'row', # 指定归一化的方式。"row"按行归一化，"column"按列归一化，"none"不处理
                cluster_rows = T, # 是否对行聚类
                cluster_cols = F ,# 是否对列聚类
                
                main = paste0('heatmap of differential ', type),
                border = F,
                gaps_col = gaps_col,        # 根据 g 进行列分割
                #clustering_distance_rows = "euclidean", # 聚类距离方法
                #clustering_method = "ward.D2" # 聚类方法
  )
  p
  filename = paste0('.././Figure/Fig2/', type, '-heatmap.pdf')
  
  ggsave(p,filename = filename, width = 8,height = 6)
}
getwd()
Group_draw(t_deg, t_rpm, 'tRNA')
Group_draw(ys_deg, ys_rpm, 'ysRNA')
Group_draw(pi_deg, pi_rpm, 'piRNA')
Group_draw(rs_deg, rs_rpm, 'rsRNA')
Group_draw(mi_deg, mi_rpm, 'miRNA')
Group_draw(m_deg, m_rpm, 'mRNA')
Group_draw(lnc_deg, lnc_rpm, 'lncRNA')


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


#all <- rbind(lnc_lable,m_lable,mi_lable,pi_lable,ys_lable,rs_lable,t_lable)
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
# write.csv(all1, './差异分析/deg_allRNA.csv', row.names= TRUE)
all1 <- t(all1[,-1])

annotation_colors <- list(
  group = c("MAL" = "#dd641a", "BEN" = "#2d555f"),
  stage = c("BEN"= "#2d555f", "I"="#e48ea0", "II"= "#79acaf", "III"= "#9b84e4", "IV"= "#d362a6", "Missing" = "#addbee"),
  type =  c("tRNA" = '#849989', "rsRNA" = '#3673b1', "ysRNA" = '#d8744f', 
            "mRNA" = '#aea7c6', "lncRNA" = '#d4c953', "miRNA" = '#93c665', 
            "piRNA" = '#7eace0')  # 新分组Treatment
)

all1 <- all1[,match(rownames(stage),colnames(all1))]
all1 <- log2(all1+1)


# 构建颜色梯度和断点（这里颜色从 navy 到 white 到 #d53e4f，断点固定为 -2 到 2）
color_gradient <- colorRampPalette(c("#4c71d0","#ffffff",  "#d50000"))(50)
breaks_custom <- seq(-3, 3, length.out = 51)
# 根据 annotation_row 计算行分割位置：
# 假设 annotation_row 第一列存放分组信息
row_groups <- as.character(annotation_row[,1])
gaps_row <- which(row_groups[-length(row_groups)] != row_groups[-1])
# gaps_row 保存的是组别转换时的行号（注意 pheatmap 在每个gap位置后插入空白行）

# 根据 g 计算列分割位置：
# 假设 g 第一列存放列分组信息
col_groups <- as.character(stage[,1])
gaps_col <- which(col_groups[-length(col_groups)] != col_groups[-1])

# 绘制热图
P <- pheatmap(all1, 
              annotation_col = stage,         # 列注释（例如存放分组信息）
              annotation_row = annotation_row,  # 行注释（例如存放 RNA 类型）
              show_colnames = FALSE,      # 是否显示列名
              show_rownames = FALSE,      # 是否显示行名
              color = color_gradient,     # 自定义颜色梯度
              fontsize = 9,               # 字体大小
              annotation_legend = TRUE,   # 显示注释图例

              border_color = 'gray',      # 单元格边框颜色
              scale = 'row',              # 行归一化
              cluster_rows = FALSE,       # 不对行聚类
              cluster_cols = FALSE,       # 不对列聚类
              main = 'heatmap of RNA',
              # method = 'average',         # 聚类方法（未启用聚类时可忽略）
              gaps_row = gaps_row,        # 根据 annotation_row 进行行分割
              gaps_col = gaps_col,        # 根据 g 进行列分割
              breaks = breaks_custom,     # 使用自定义断点
              annotation_colors = annotation_colors
)
P

ggsave(P,filename = './figure/fig2/RNA-heatmap.pdf',width = 8,height = 6)



