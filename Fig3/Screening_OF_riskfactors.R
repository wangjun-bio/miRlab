rm(list = ls())
getwd()
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
#saveRDS(data,file = './TCGA/rawdata.rds')

gene <- read.csv('./All_Model_mRNA_Top_10_Features.csv')
gene <- gene$Feature
gene_exp <- as.data.frame(t(data[gene,]))
# 找出列名包含 "NA." 的列的索引
cols_to_remove <- grep("NA", colnames(gene_exp))

# 如果存在这样的列，则删除它们
if (length(cols_to_remove) > 0) {
  gene_exp <- gene_exp[, -cols_to_remove]
}
gene_exp <- subset(gene_exp, select = -ANKRD30BL)

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
# single_result <- single_result[single_result$pvalue < 0.1,]
gene <- single_result$id
hr <- sprintf("%.3f", single_result$HR) # 将surv.data数据框中的HR列取出，保留3位小数，并用sprintf函数进行格式化，赋值给hr变量
hrLow  <- sprintf("%.3f", single_result$HR.95L) # 将surv.data数据框中的lower95列取出，保留3位小数，并用sprintf函数进行格式化，赋值给hrLow变量
hrHigh <- sprintf("%.3f", single_result$HR.95H) # 将surv.data数据框中的upper95列取出，保留3位小数，并用sprintf函数进行格式化，赋值给hrHigh变量
Hazard.ratio <- paste0(hr, "(", hrLow, "-", hrHigh, ")") # 将hr、hrLow和hrHigh变量中的内容合并为一个字符，使用paste0函数，括号前面是hr的值，括号内是hrLow和hrHigh值的范围，赋值给Hazard.ratio变量
pValue <- ifelse(single_result$pvalue < 0.001, "<0.001",round(single_result$pvalue,3))#,sprintf("%.3f", single_result$pValue)) # 如果surv.data数据框中的pValue列的值小于0.001，则将pValue变量赋值为"<0.001"，否则将pValue变量赋值为pValue的值，保留3位小数，并用sprintf函数进行格式化
#write.csv(single_result,'top10mRNA_unicox_result.csv')

pdf('./figure/top10_UniCox.pdf',width = 13.5,height = 9)
# png('./figure/top10_UniCox.png',width = 13.5*300,height = 9*300,res = 600)
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


library(dplyr)
library(caret)
df <- gene_sur
inTrain<-createDataPartition(y=c(1:nrow(df)), p=0.7, list=F)
train<-df[inTrain,]      
test<-df[-inTrain,]      
trainOut=cbind(id=row.names(train),train)#结果输出的格式      
testOut=cbind(id=row.names(test),test)#结果输出的格式
table(train$OS)
table(test$OS)
library(glmnet)
x <- as.matrix(train[,-c(1,2,25)])
y <- data.matrix(Surv(time = train$OS.time,event = train$OS))

set.seed(1005)
fit <- glmnet(x,y,family = 'cox',alpha = 1)
plot(fit,xvar = 'lambda')
cv_fit <- cv.glmnet(x,y,family = 'cox',type.measure = 'deviance')
cv_fit$lambda.min
plot(cv_fit)

coefficients <- coef(cv_fit,s=cv_fit$lambda.min)
Active.index <- which(as.numeric(coefficients) != 0)
active.coefficients <- as.numeric(coefficients)[Active.index]
sig_gene_lasso <- rownames(coefficients)[Active.index]
sig_gene_lasso
sig_gene_exp <- train[,sig_gene_lasso]
sig_gene_exp$OS <- train$OS
sig_gene_exp$OS.time <- train$OS.time

########multiCox
sum(is.na(sig_gene_exp))
str(sig_gene_exp)
multi_cox <- coxph(Surv(OS.time,OS)~., data = sig_gene_exp)
summary(multi_cox)
multi_cox <- step(multi_cox,direction = 'both')
multi_coxsum = summary(multi_cox)
nrow(multi_coxsum$coefficients)

outMultiTab = data.frame()
outMultiTab = cbind(
  coef = multi_coxsum$coefficients[,"coef"],
  HR = multi_coxsum$conf.int[,"exp(coef)"],
  HR.95L = multi_coxsum$conf.int[,'lower .95'],
  HR.95H = multi_coxsum$conf.int[,'upper .95'],
  pvalue = multi_coxsum$coefficients[,'Pr(>|z|)']
)
outMultiTab <- as.data.frame(outMultiTab)
gene <- rownames(outMultiTab)
hr <- sprintf("%.3f", outMultiTab$HR) # 将surv.data数据框中的HR列取出，保留3位小数，并用sprintf函数进行格式化，赋值给hr变量
hrLow  <- sprintf("%.3f", outMultiTab$HR.95L) # 将surv.data数据框中的lower95列取出，保留3位小数，并用sprintf函数进行格式化，赋值给hrLow变量
hrHigh <- sprintf("%.3f", outMultiTab$HR.95H) # 将surv.data数据框中的upper95列取出，保留3位小数，并用sprintf函数进行格式化，赋值给hrHigh变量
Hazard.ratio <- paste0(hr, "(", hrLow, "-", hrHigh, ")") # 将hr、hrLow和hrHigh变量中的内容合并为一个字符，使用paste0函数，括号前面是hr的值，括号内是hrLow和hrHigh值的范围，赋值给Hazard.ratio变量
pValue <- ifelse(outMultiTab$pvalue < 0.001, "<0.001",round(outMultiTab$pvalue,3))#,sprintf("%.3f", single_result$pValue)) # 如果surv.data数据框中的pValue列的值小于0.001，则将pValue变量赋值为"<0.001"，否则将pValue变量赋值为pValue的值，保留3位小数，并用sprintf函数进行格式化
write.csv(outMultiTab,'Multicox_result.csv')


pdf('./figure/MultiCox.pdf',width = 12,height = 9)
# png('./figure/MultiCox.png',width = 12*300,height = 9*300,res = 300)
n <- nrow(outMultiTab) # 计算surv.data数据框的行数，即数据的观测数量
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