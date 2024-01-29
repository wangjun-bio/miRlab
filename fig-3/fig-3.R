#####################################################################################################################################################
#DEseq2
library(DESeq2)
setwd("C:/Users/树金/Desktop/代码/fig-3/DEseq2")
getwd()

if (!require("circlize")) {
  BiocManager::install("circlize")
  library(circlize)
}


file_list <- list.files(pattern = "_full_join.csv")
data <- list()
for (file in file_list) {
  data[[file]] <- read.csv(file = file, header = T, sep = ",")
}

NP <- data[[1]]
NT <- data[[2]]
PT <- data[[3]]

mycounts <- NP
rownames(mycounts)<-mycounts$ID
mycounts <- mycounts[,-1]

#生成新的grop的矩阵
a <- mycounts[1:60,59:60]
rownames(a) <- NULL
a$KK <- colnames(mycounts)
b <- a[,-1]
rownames(b) <- NULL
b$condition  = factor(rep(c("Normal" ,"tumor"), c(30,  30)))
colnames(b)[1] <- "ID"
grop <- b
grop <- grop[,-1]
row.names(grop) <- grop$KK
grop$condition <- as.factor(grop$condition)
class(grop$condition)

#检查count文件和colData文件中的样本顺序是否一致
all(rownames(grop) == colnames(mycounts))
condition <- factor(grop$condition)
condition = relevel( condition, "Normal")
#用于检查对象 mycounts 中是否存在任何缺失（NA）值,s输出FALSE即为没有
any(is.na(mycounts))
mycounts[is.na(mycounts)] <- 0
dds <- DESeqDataSetFromMatrix(countData=mycounts, colData=grop, design=~condition)
dds <- DESeq(dds)
#查看基因的差异倍数，和显著性p值。T，N 在后，意为 T 相较于 N中哪些基因上调/下调
res <- results(dds, contrast = c('condition', 'tumor', 'Normal'))
res = res[order(res$pvalue),]
head(res)
write.csv(res,file="circRNA-NP-DEseq2.csv",quote = FALSE)

dds <- estimateSizeFactors(dds)
tpm <- t(apply(counts(dds, normalized = TRUE), 1, function(x) (x / sum(x)) * 1e6))
write.csv(tpm,file="circRNA-NT-tpm.csv",quote = FALSE)

#####################################################################################################################################################
#火山图
library(ggplot2)
library(scales)
library(ggrepel)
library(svglite)
library(dplyr)
setwd("C:/Users/树金/Desktop/代码/fig-3/A")
getwd()

file_list <- list.files(pattern = "DEseq2.csv")
data <- list()
for (file in file_list) {
  data[[file]] <- read.csv(file = file, header = T, sep = ",")
}
NP <- data[[1]]
NT <- data[[2]]
PT <- data[[3]]

m <- NT
rownames(m) <- m$X
m <- m[,-1]
#去除NA值
m <- m[complete.cases(m),]
cut_off_pvalue = 0.05  
cut_off_padj = 0.05
cut_off_basemean = 3.51
cut_off_logFC = 1

m$change = ifelse(m$pvalue <= cut_off_pvalue & 
                  m$padj <= cut_off_padj & 
                  m$baseMean > cut_off_basemean & 
                  abs(m$log2FoldChange) >= cut_off_logFC, 
                  ifelse(m$log2FoldChange > cut_off_logFC ,'Up','Down'),'Not sig')

p <- ggplot(m, aes(x = log2FoldChange, y = -log10(pvalue), color = change)) +
  geom_point(alpha = 1, size = 3.5) +   # 添加散点并设置透明度和大小
  scale_color_manual(values = c("blue", "black", "red")) +   # 自定义颜色
  theme_bw() +   # 使用白底黑字主题
  theme(panel.grid = element_blank()) +   # 去除网格线
  geom_vline(xintercept = c(-1, 1), lty = 2, col = "black", lwd = 0.4) +   # 添加垂直虚线
  geom_hline(yintercept = -log10(cut_off_pvalue), lty = 4, col = "black", lwd = 0.4) +   # 添加水平虚线
  ylab('-log10 (p-value)') +   # 设置y轴标签
  xlab('log2 (FoldChange)') +   # 设置x轴标签
  labs(col = '') +   # 设置图例标签为空
  coord_cartesian(xlim = c(-12, 12), ylim = c(0, 30)) +   # 设置坐标轴范围
  theme(legend.position = "bottom", axis.text = element_text(colour = 'black')) +   # 设置图例位置和坐标轴文本颜色
  ggtitle("Normal vs. Paracancerous") +   # 设置图表标题
  theme(plot.title = element_text(hjust = 0.5)) +   # 标题居中
  scale_x_continuous(expand = c(0, 0), limits = c(-12, 12), breaks = seq(-30, 30, by = 10)) +   # 自定义x轴
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 30, by = 5)) +   # 自定义y轴
  theme(text = element_text(size = 30))   # 自定义文本大小

count_summary <- m %>% group_by(change) %>% summarize(count = n())
p1 <- p + geom_text(aes(x = -9, y = 25, label = paste("Up: ", count_summary$count[1])), size = 5, color = "red") +
          geom_text(aes(x = -9, y = 23, label = paste("Down: ", count_summary$count[2])), size = 5, color = "blue") +
          geom_text(aes(x = -9, y = 21, label = paste("Not sig: ", count_summary$count[3])), size = 5, color = "black")
p1
ggsave("circRNA-Volcano-NP.pdf",width =10, height =9)


#####################################################################################################################################################
#热图
library(ggplot2)
library(pheatmap)
setwd("C:/Users/树金/Desktop/代码/fig-3/B")
getwd()

file_list <- list.files(pattern = ".csv")
data <- list()
for (file in file_list) {
  data[[file]] <- read.csv(file = file, header = T, sep = ",")
}
circRNA_55_TPM <- data[[1]]
group_55 <- data[[2]]

#设置函数，赋值文件的名称(NP\NT\PT)
set_value <- function(var_name, value) {
  assign(var_name, value, envir = .GlobalEnv)
}
set_value("file_name", circRNA_55_TPM)

#处理表达矩阵
m <- file_name
row.names(m) <- m$ID
m <- m[,-1]
top_circRNA <- m
#修改输入文件的列明
colnames_list <- paste0("Test_", 1:90)
colnames(top_circRNA) <- colnames_list
#处理分组信息
pheatmap_grop <- read.csv(file = "grop-pheatmap.csv",header = T, sep = ",")
rownames(pheatmap_grop) <- pheatmap_grop$ID
pheatmap_grop <- pheatmap_grop[,-1]
# 绘制热图
bk <- c(seq(-4,-0.001,by=0.001),seq(0,4,by=0.001))
# 使用 pheatmap 函数时，设置 annotation_colors 参数
p <- pheatmap(
  top_circRNA,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  scale = "row",
  treeheight_col = 2,
  annotation_col = pheatmap_grop,
  color = c(
    colorRampPalette(colors = c("blue", "white"))(length(bk) / 2),
    colorRampPalette(colors = c("white", "red"))(length(bk) / 2)
  ),
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "DEcircRNAs-NT",
  border = "white",
  cellwidth = 6, cellheight = 6,
  fontsize = 10,
  fontsize_col = 5,
  fontsize_row = 5,
  legend_breaks=seq(-4,4,2),
  breaks=bk)
p

pdf_file <- "A.pdf"
ggsave(pdf_file, plot = p, device = "pdf", width = 18, height = 16, units = "in")

#########################################################################
library(ggpubr)
library(ggplot2)
setwd("C:/Users/树金/Desktop/代码/fig-3/C")

file_list <- list.files(pattern = ".csv")
data <- list()
for (file in file_list) {
  data[[file]] <- read.csv(file = file, header = T, sep = ",")
}
NT <- data[[1]]
PT <- data[[2]]

data <- PT[,1:4]
data[is.na(data)] <- 0
colnames(data)[3] <- "circRNA"
colnames(data)[4] <- "mRNA"

p <- ggscatter(data = data, x = "circRNA", y = "mRNA", 
               color = "red",  # 将红色改为黑色，也可以设置其他颜色
               fill = "red",     # 填充颜色为红色
               shape = 21, size = 5,        # 16代表圆圈的形状
               add = "reg.line", conf.int = TRUE, 
               add.params = list(color = "black", fill = "blue"),
               cor.coef = TRUE,
               cor.method = "pearson",
               cor.coef.coord = c(-3, 1)) + 
  ggtitle("PT") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(-3, 2), ylim = c(-3, 2)) +
  annotate("text", x = -2.2, y = 1.5, 
           label = paste("Pearson Correlation:", round(cor(m$circRNA, m$mRNA), 2)))

p

ggsave("PT.pdf", plot = p, width = 7, height = 7)




N_bed <- read.csv(file = "55_N_bed.csv", row.names = 1)
T_bed <- read.csv(file = "55_T_bed.csv", row.names = 1)

# 初始化 Circos 图，设置染色体 ideogram
circos.initializeWithIdeogram(species = "hg38")

circos.genomicDensity(T_bed, col = "#d9942e", track.height = 0.3) 
cols <-  c("red", "blue")
bed_list <- list(N_bed, T_bed)
circos.genomicTrack(bed_list, 
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicPoints(region, value, pch = 16, cex = 1, col = cols[i], ...)
                    })

