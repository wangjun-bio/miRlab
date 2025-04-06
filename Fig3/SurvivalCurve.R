# 之前的代码保持不变
setwd('/Users/willow/Desktop/2025/cfRNA乳腺癌/DATA-public/')
p <- read.table('./TCGA/TCGA-BRCA.GDC_phenotype.tsv.gz',header = T, sep = '\t',quote = '')

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


gene <- read.csv('./All_Model_mRNA_Top_10_Features.csv')
gene <- gene$Feature
gene_exp <- as.data.frame(t(data[gene,]))
gene_exp <- gene_exp[, colSums(is.na(gene_exp)) != nrow(gene_exp)]

gene_sur <- cbind(sur$OS,sur$OS.time,gene_exp)
colnames(gene_sur)[1] <- 'OS'
colnames(gene_sur)[2] <- 'OS.time'
gene_sur <- na.omit(gene_sur)


# 生存分析部分
library(survival)
library(survminer)

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

write.csv(rt, 'TCGA-生存分析分组.csv')
#绘制生存分析图-----------------------------------------------------------------
rt <- read.csv('./TCGA-生存分析分组.csv')
# 检查分组结果
table(rt$Group, useNA = "always")

# 生存分析
fit <- survfit(Surv(OS.time, OS) ~ Group, data = rt)
surv_diff <- survdiff(Surv(OS.time, OS) ~ Group, data = rt)

# 计算p值
p_value <- 1 - pchisq(surv_diff$chisq, df = length(surv_diff$n)-1)

# 可视化
surPlot <- ggsurvplot(fit,
                      data = rt,
                      pval = paste0("Log-Rank p = ", signif(p_value,3)),
                      risk.table = TRUE,
                      title = "Survival by Combined Gene Expression",
                      legend.labs = c("≥2 High", "≥2 Low"),
                      palette = c("#c5b4d8", "#377EB890"))

surPlot <- ggsurvplot(fit,
                      data = rt,
                      pval = paste0("Log-Rank p = ", signif(p_value,3)),
                      conf.int = TRUE,                # 显示置信区间
                      conf.int.style = "ribbon",      # 使用带状显示（默认）
                      conf.int.alpha = 0.2,           # 设置透明度（0-1）
                      risk.table = TRUE,
                      title = "Survival by Combined Gene Expression",
                      legend.labs = c("≥2 High", "≥2 Low"),
                      palette = c("#b7a4d0", "#377EB890"))

print(surPlot)  # 显示带置信区间的生存曲线

# 输出统计结果
cat("Comparison between ≥2 High vs ≥2 Low groups:\n")
cat("Log-Rank p-value:", p_value, "\n")
cat("Number in High_group:", sum(rt$Group == "High_group", na.rm = T), "\n")
cat("Number in Low_group:", sum(rt$Group == "Low_group", na.rm = T), "\n")

# 可选：计算风险比(HR)
cox_model <- coxph(Surv(OS.time, OS) ~ Group, data = rt)
summary_cox <- summary(cox_model)
hr <- summary_cox$coefficients[,"exp(coef)"]
hr_ci <- summary_cox$conf.int[, "lower .95"]
hr_ci_up <- summary_cox$conf.int[, "upper .95"]

# 保存图形
# pdf(file = '生存曲线.pdf', height = 8, width = 10)
# print(surPlot)
# dev.off() 