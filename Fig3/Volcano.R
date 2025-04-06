# 之前的代码保持不变
setwd('/Users/willow/Desktop/2025/cfRNA乳腺癌/DATA-public/')
p <- read.table('./TCGA/TCGA-BRCA.GDC_phenotype.tsv.gz',header = T,
                sep = '\t',quote = '')

colnames(p)[grep("receptor_status", colnames(p))]

table(p$metastatic_breast_carcinoma_estrogen_receptor_status == 'Negative' &
        p$metastatic_breast_carcinoma_progesterone_receptor_status == 'Negative' &
        p$metastatic_breast_carcinoma_lab_proc_her2_neu_immunohistochemistry_receptor_status == 'Negative')

library(stringr)

tab1 <- p[1:2] # includes 'sampleID' & "AJCC_Stage_nature2012"  
tumor <- tab1[substr(tab1$submitter_id.samples,14,15) < 10,]         
tumor$TCGAID <- str_sub(tumor$submitter_id.samples,1,12)
normal <- tab1[!substr(tab1$submitter_id.samples,14,15) <10,]
normal$TCGAID <- str_sub(normal$submitter_id.samples,1,12)
dim(tumor)
dim(normal)

library(data.table)
ano <- fread('./TCGA/gencode.v22.annotation.gene.probeMap')
ano = ano[!duplicated(ano$gene),]

data <- as.data.frame(fread('./TCGA/TCGA-BRCA.htseq_fpkm.tsv.gz'))
rownames(data) <- data$Ensembl_ID
data <- data[,-1]
gene <- intersect(rownames(data),ano$id)
data <- data[gene,]
identical(rownames(data),ano$id)

data <- data[match(ano$id,rownames(data)),]
identical(rownames(data),ano$id)

rownames(data) = ano$gene

sur <- as.data.frame(fread('./TCGA/TCGA-BRCA.survival.tsv'))
overlap <- intersect(colnames(data),sur$sample)
data <- data[,overlap]
rownames(sur) <- sur$sample
sur <- sur[overlap,]


gene <- read.csv('./ROC-top10_mRNA.csv')
gene <- gene$id
gene_exp <- as.data.frame(t(data[gene,]))
gene_exp <- gene_exp[,-62]
gene_sur <- cbind(sur$OS,sur$OS.time,gene_exp)
colnames(gene_sur)[1] <- 'OS'
colnames(gene_sur)[2] <- 'OS.time'
gene_sur <- na.omit(gene_sur)

rt <- gene_sur


var1 <- 'DLST'
var2 <- 'EGFL7'
var3 <- 'DOCK4'

# 根据基因表达中位数分组
a <- ifelse(rt[, var1] < median(rt[, var1]), paste0(var1, " low"), paste0(var1, " high"))
b <- ifelse(rt[, var2] < median(rt[, var2]), paste0(var2, " low"), paste0(var2, " high"))
c <- ifelse(rt[, var3] < median(rt[, var3]), paste0(var3, " low"), paste0(var3, " high"))

# 组合分组
Type <- paste(a, "+", b, "+", c)
rt <- cbind(rt, Type)

# 计算组与组之间的 P 值
fit <- survdiff(Surv(OS.time, OS) ~ Type, data = rt)
p_value <- 1 - pchisq(fit$chisq, length(fit$n) - 1)

library(survival)
library(survminer)

# 1. 整体生存差异检验
fit <- survfit(Surv(OS.time, OS) ~ Type, data = rt)
log_rank <- survdiff(Surv(OS.time, OS) ~ Type, data = rt)
p_overall <- 1 - pchisq(log_rank$chisq, df = length(log_rank$n)-1)

# 2. 两两比较（使用BH法校正p值）
pairwise_res <- pairwise_survdiff(
  Surv(OS.time, OS) ~ Type,
  data = rt,
  p.adjust.method = "BH"
)

# 3. 整理输出结果
p_matrix <- round(pairwise_res$p.value, 4)
p_matrix_formatted <- format(p_matrix, scientific = FALSE, na.encode = FALSE)

cat("Overall log-rank test p-value:", p_overall, "\n\n")
cat("Pairwise comparison adjusted p-values:\n")
print(p_matrix_formatted, quote = FALSE)



library(survival)
library(survminer)
library(stringr)

# 创建新分组
rt$high_count <- str_count(rt$Type, "high")
rt$low_count <- str_count(rt$Type, "low")

# 定义分组规则：
# High_group = 包含2个或3个high表达
# Low_group = 包含2个或3个low表达
rt$Group <- ifelse(rt$high_count >= 2, "High_group",
                   ifelse(rt$low_count >= 2, "Low_group", NA))

library(DESeq2)
library(dplyr)
getwd()

##############火山图###########
data_gene <- rt[rt$Group == 'High_group' | rt$Group == 'Low_group',]

data_gene <- data_gene %>% mutate(group = case_when(Group == 'High_group' ~ 'High',
                                                    Group == 'Low_group' ~ 'Low'))
data_exp <- data[,rownames(data_gene)]


data_exp[is.na(data_exp)] <- 0
row_sum <- rowSums(data_exp == 0)
# Calculate the threshold number of columns
n <- 0.75 * ncol(data_exp)
# Filter the rows based on the threshold
filtered_mRNA <- data_exp[row_sum <= (ncol(data_exp) - n), ]

group_list <- as.data.frame(data_gene$group)
rownames(group_list) <- rownames(data_gene)
colnames(group_list) <- 'condition'
group_list$condition <- as.factor(group_list$condition)

# str(filtered_mRNA)
filtered_mRNA <- round(filtered_mRNA)
# write.csv(filtered_mRNA,file = 'fitermRNA.csv')
# write.csv(group_list,file = 'group_list.csv')

# filtered_mRNA <- read.csv('./fitermRNA.csv',row.names = 1,check.names = F)
# group_list <- read.csv('./group_list.csv',row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = filtered_mRNA, colData = group_list, design= ~condition)
######过滤掉低表达的counts值，count函数过滤
dds <- dds[rowSums(counts(dds))>=1,]###可以定义为1
########## 差异分析
dds <- DESeq(dds)
######### 构建contrast对象，用于后续的差异结果的提取，谁是肿瘤/普通
table(group_list$condition)
contrast <- c('condition','High','Low')
res <- results(dds,contrast = contrast)
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
library(tidyverse)
res1 <- res1 %>% arrange(res1$pvalue)


deg1 <- res1
colnames(deg1)[2] <- 'logFC'
deg1$logP <- -log10(deg1$padj)
deg1$group <- 'NS'
deg1$group[which(deg1$logFC >= 1 & deg1$padj <0.05)] <- 'up'
deg1$group[which(deg1$logFC <= -1 & deg1$padj < 0.05)] <- 'down'
write.csv(deg1, './TCGA-deg-mRNA.csv')
table(deg1$group)

bullcol <- c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25")

logfc_low <- -1
logfc_high <- 1
P <- 0.05
table(deg1$group)
up_num = 499
down_num = 18

library(ggplot2)
library(ggrepel)
p <- ggplot(deg1, aes(logFC, logP)) +
  geom_hline(yintercept = -log10(P), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(logfc_low, logfc_high), linetype = "dashed", color = "black") +
  geom_point(aes(size = logP, color = logP)) +
  scale_color_gradientn(colors = bullcol) +
  scale_size_continuous(range = c(1, 3)) +
  theme_bw(base_size = 12) +
  annotate("text", x = Inf, y = Inf, label = paste("Up = ", up_num), hjust = 1.85, vjust = 2.5, cex = 5, 
           arrow = arrow(type = "closed", length = unit(0.2, "inches"))) +
  annotate("text", x = -Inf, y = Inf, label = paste("Down = ", down_num), hjust = -0.1, vjust = 2.5, cex = 5, 
           arrow = arrow(type = "closed", length = unit(0.2, "inches"))) +
  theme(panel.grid = element_blank(), legend.position = 'right', legend.justification = c(0, 1)) +
  guides(color = guide_colorbar(title = "-Log10(pvalue)", barheight = 8, barwidth = 1), size = "none") +
  xlab("Log2FC") + ylab("-Log10(FDR)") 

# 提取排名前 5 的基因
# 这里假设按照 logP 从大到小排序，取前 5 个基因
top_5_genes <- rownames(deg1)[order(deg1$logP, decreasing = TRUE)[1:5]]
deg1$label <- ifelse(rownames(deg1) %in% top_5_genes, rownames(deg1), NA)

# 添加基因名称标签
p <- p + geom_text_repel(data = deg1, aes(label = label), size = 3.5, max.overlaps = 15, force = 10, 
                         box.padding = unit(0.5, "lines"), nudge_x = 0.5, nudge_y = 1)

# 显示图形
print(p)

# 保存图形
ggsave(plot = p, filename = 'volcano.pdf', height = 8, width = 10, dpi = 600)




p = ggplot(deg1, aes(logFC, -log10(padj))) +
  # 横向水平参考线：
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color="#999999") +
  # 纵向垂直参考线：
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="#999999") +
  # 散点图:
  geom_point(aes(size=-log10(padj), color=-log10(padj))) +
  # 指定颜色渐变模式：
  scale_color_gradientn(values=seq(0,1,0.2),
                        colors=c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25")) +
  # 指定散点大小渐变模式：
  scale_size_continuous(range=c(0.5,4)) +
  # 主题调整：
  theme_bw() +
  # 调整主题和图例位置：
  theme(panel.grid=element_blank(),
        legend.position=c(0.01,0.7),
        legend.justification=c(0,1)) +
  # 设置部分图例不显示：
  guides(col=guide_colourbar(title="-Log10(FDR)"),
         size="none") +
  # 修改坐标轴：
  xlab("Log2FC") +
  ylab("-Log10(FDR)") +
  scale_x_continuous(limits = c(-3,3)) +
  # 添加对称的注释
  annotate("text", x = Inf, y = Inf, label = paste("Up = ", up_num), hjust = 2, vjust = 2.5, cex = 5.5,
           arrow = arrow(type = "closed", length = unit(0.2, "inches"))) +
  annotate("text", x = -Inf, y = Inf, label = paste("Down = ", down_num), hjust = -1, vjust = 2.5, cex = 5.5,
           arrow = arrow(type = "closed", length = unit(0.2, "inches")))
p

# ggsave(plot = p,filename = 'volcano.pdf',height = 8,width = 11,dpi = 600)

