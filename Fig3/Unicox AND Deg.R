rm(list = ls())
setwd('D:/r_project/乳腺癌/')
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
rownames(data) = ano$gene

sur <- as.data.frame(fread('./TCGA/TCGA-BRCA.survival.tsv'))
overlap <- intersect(colnames(data),sur$sample)
data <- data[,overlap]
rownames(sur) <- sur$sample
sur <- sur[overlap,]
#saveRDS(data,file = './TCGA/rawdata.rds')

gene <- read.csv('./ROC-top10mRNA.csv')
gene <- gene$X
gene_exp <- as.data.frame(t(data[gene,]))
gene_exp <- gene_exp[,-62]
gene_sur <- cbind(sur$OS,sur$OS.time,gene_exp)
colnames(gene_sur)[1] <- 'OS'
colnames(gene_sur)[2] <- 'OS.time'
gene_sur <- na.omit(gene_sur)

######Univariate-COX
library(survival)
library(survminer)

set.seed(1234)
fit <- coxph(Surv(OS.time,OS) ~ .,data = gene_sur)
summary(fit)
single_result=data.frame()
for(i in colnames(gene_sur[,3:ncol(gene_sur)])){
  cox <- coxph(Surv(OS.time,OS) ~ get(i), data = gene_sur)
  coxSummary = summary(cox)
  single_result=rbind(single_result,
                      cbind(id=i,
                            HR=coxSummary$conf.int[,"exp(coef)"],
                            HR.95L=coxSummary$conf.int[,"lower .95"],
                            HR.95H=coxSummary$conf.int[,"upper .95"],
                            pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
}
single_result[,2:5]<-apply(single_result[,2:5],2,as.numeric)
single_result
single_result <- single_result[single_result$pvalue < 0.1,]
gene <- single_result$id
hr <- sprintf("%.3f", single_result$HR) # 将surv.data数据框中的HR列取出，保留3位小数，并用sprintf函数进行格式化，赋值给hr变量
hrLow  <- sprintf("%.3f", single_result$HR.95L) # 将surv.data数据框中的lower95列取出，保留3位小数，并用sprintf函数进行格式化，赋值给hrLow变量
hrHigh <- sprintf("%.3f", single_result$HR.95H) # 将surv.data数据框中的upper95列取出，保留3位小数，并用sprintf函数进行格式化，赋值给hrHigh变量
Hazard.ratio <- paste0(hr, "(", hrLow, "-", hrHigh, ")") # 将hr、hrLow和hrHigh变量中的内容合并为一个字符，使用paste0函数，括号前面是hr的值，括号内是hrLow和hrHigh值的范围，赋值给Hazard.ratio变量
pValue <- ifelse(single_result$pvalue < 0.001, "<0.001",round(single_result$pvalue,3))#,sprintf("%.3f", single_result$pValue)) # 如果surv.data数据框中的pValue列的值小于0.001，则将pValue变量赋值为"<0.001"，否则将pValue变量赋值为pValue的值，保留3位小数，并用sprintf函数进行格式化
write.csv(single_result,'top10mRNA_unicox_result.csv')

pdf('top10_UniCox.pdf',width = 13.5,height = 9)
png('top10_UniCox.png',width = 13.5*300,height = 9*300,res = 600)
n <- nrow(single_result) # 计算surv.data数据框的行数，即数据的观测数量
nRow <- n+1 # 将行数加一，用于设置y轴的范围，以确保能够容纳全部数据
ylim <- c(1,nRow) # 设置y轴的范围，范围从1到nRow
layout(matrix(c(1,2),nc=2),width=c(2.5,2)) # 使用layout函数设置图形的布局，将画布分为两列，并设置宽度比例为2.5和2

xlim = c(0, 2.0) # 设置x轴的范围，范围从0到2.5
par(mar = c(4, 2.5, 2, 1)) # 使用par函数设置图形参数，mar参数用于设置边距，上下左右依次为4, 2.5, 2, 1
plot(0, xlim = xlim, ylim = ylim, type = "n", axes = T, xlab = "", ylab = "") # 创建一个空白的绘图区域，使用plot函数，设置x轴范围、y轴范围、绘图类型为“n”，不绘制坐标轴，并设置x轴和y轴标签为空字符串
text.cex = 0.8 # 设置文本的缩放比例
text(0, n:1, gene, adj = 0, cex = text.cex) # 在坐标轴上添加gene变量的文本，位置为(0, n:1)，水平方向对齐方式为0，文本缩放比例为
text(0, n + 1, 'Gene', cex = text.cex, font = 2, adj = 0)
text(0.9, n:1, pValue, adj = 1, cex = text.cex) # 在坐标轴上添加pValue变量的文本，位置为(1.2 - 0.5 * 0.2, n:1)，水平方向对齐方式为1，文本缩放比例为
text(0.9, n + 1, 'pValue', cex = text.cex, font = 2, adj = 1) # 在坐标轴上添加文本'pValue'，位置为(1.2 - 0.5 * 0.2, n + 1)，水平方向对齐方式为1，文本缩放比例为text.cex，设置字体为粗体
text(1.8, n:1, Hazard.ratio, adj = 1, cex = text.cex) # 在坐标轴上添加Hazard.ratio变量的文本，位置为(2.5, n:1)，水平方向对齐方式为1，文本缩放比例为
text(1.7, n + 1, 'Hazard ratio', cex = text.cex, font = 2, adj = 1) # 在坐标轴上添加文本'Hazard ratio'，位置为(2.5, n + 1)，水平方向对齐方式为1，文本缩放比例为text.cex，设置字体为粗体

par(mar = c(4, 1, 2, 1), mgp = c(2, 0.5, 0)) # 使用par函数设置图形参数，其中mar参数用于设置边距，mgp参数用于设置刻度标签与轴线的距离
xlim = c(0, max(as.numeric(hrLow), as.numeric(hrHigh))) # 设置x轴的范围，范围从0到hrLow和hrHigh中的最大值
plot(0, xlim = xlim, ylim = ylim, type = "n", axes = T, ylab = "", xaxs = "i", xlab = "Hazard ratio") # 创建一个空白的绘图区域，使用plot函数，设置x轴范围、y轴范围、绘图类型为“n”，绘制坐标轴，y轴标签为空，x轴使用整数刻度类型，x轴标签为"Hazard ratio"
arrows(as.numeric(hrLow), n:1, as.numeric(hrHigh), n:1, angle = 90, code = 3, length = 0.05, col = "darkblue", lwd = 2.5) # 使用arrows函数绘制箭头，从hrLow到hrHigh，位置为(as.numeric(hrLow), n:1)到(as.numeric(hrHigh), n:1)，箭头角度为90度，箭头类型为3，箭头长度为0.05，箭头颜色为深蓝色，线宽为2.5
abline(v = 1, col = "black", lty = 2, lwd = 2) # 使用abline函数绘制垂直于x轴的虚线，位置为x = 1，线颜色为黑色，线型为2，线宽为2
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green') # 根据hr的数值是否大于1来设置盒子颜色，当hr大于1时，盒子颜色为红色，否则为绿色
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex = 1.3) # 使用points函数绘制散点图，位置为(as.numeric(hr), n:1)，点的形状为方块，颜色为盒子颜色，大小为1.3倍
axis(1) # 绘制x轴刻度线

dev.off()



###########生存分析 DLST
# 根据基因表达，对数据分组
rt <- gene_sur
var1 <- 'DLST'
var2 <- 'FBXO31'
a=ifelse(rt[,var1]<median(rt[,var1]),paste0(var1," low"),paste0(var1," high"))
b=ifelse(rt[,var2]<median(rt[,var2]),paste0(var2," low"),paste0(var2," high"))
Type = paste(a,"+",b)
rt = cbind(rt, Type)  # 自定义调色板

# 生存差异统计
length=length(levels(factor(Type)))# 计算分组的数量
diff=survdiff(Surv(OS.time,OS)~Type,data=rt)# 进行生存差异分析
pValue=1 - pchisq(diff$chisq,df = length - 1)# 计算p值
if(pValue<0.001){
  pValue="p<0.001"# 如果p值小于0.001，则设定为"p<0.001"
}else{
  pValue = paste0("p=", sprintf("%.03f",pValue))# 否则格式化p值
}
fit <- survfit(Surv(OS.time,OS) ~ Type, data = rt)# 拟合生存曲线


ggsurvplot(fit,
           data = rt,
           conf.int = F,# 不包含置信区间
           pval = pValue,
           pval.size = 5,
           legend.labs = levels(factor(rt[,"Type"])),# 图例标签
           legend.title = "Type", # 图例标题
           legend = "right",# 图例位置
           xlab = "Time(days)",# x轴标签
           break.time.by = 730, # x轴刻度间隔
           risk.table.title = "", # 风险表格标题
           #surv.median.line = "hv", # 添加中位生存时间线
           risk.table = F)# 风险表高度
pdf(file = 'surplot.pdf',height = 8,width = 12)
png(filename = 'surplot.png',height = 12*300,width = 16.5*300,res = 600)
#print(surPlot)
dev.off()
##############
data_gene <- rt[rt$Type == 'DLST high + FBXO31 high' | rt$Type == 'DLST low + FBXO31 low',]
library(dplyr)
data_gene <- data_gene %>% mutate(group = case_when(Type == 'DLST high + FBXO31 high' ~ 'High',
                                                    Type == 'DLST low + FBXO31 low' ~ 'Low'))
data_exp <- data[,rownames(data_gene)]

library(DESeq2)
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

str(filtered_mRNA)
filtered_mRNA <- round(filtered_mRNA)
write.csv(filtered_mRNA,file = 'fitermRNA.csv')
write.csv(group_list,file = 'group_list.csv')

filtered_mRNA <- read.csv('./fitermRNA.csv',row.names = 1,check.names = F)
group_list <- read.csv('./group_list.csv',row.names = 1)

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
  geom_hline(yintercept = -log10(P),linetype= "dashed", color = "black") +
  geom_vline(xintercept = c(logfc_low,logfc_high), linetype = "dashed", color = "black") +
  geom_point(aes(size = logP,color = logP)) +
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
#geom_text_repel(data = core_gene_data, aes(label = gene, color = pvalue_log), size = 3.5,max.overlaps = 15,force = 10,box.padding = unit(0.5, "lines"), nudge_x = 0.5, nudge_y = 1)
p
ggsave(plot = p,filename = 'volcano.pdf',height = 8,width = 10,dpi = 600)

deg1$label<-c(rownames(deg1)[1:10],rep(NA,(nrow(deg1)-10)))




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

ggsave(plot = p,filename = 'volcano.pdf',height = 8.5,width = 8,dpi = 600)
ggsave(plot = p,filename = 'volcano.png',height = 8.5,width = 8,dpi = 600)



































