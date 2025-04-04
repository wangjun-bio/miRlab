################
data <- read.csv('../fitermRNA.csv',row.names = 1,check.names = F)
group <- read.csv('./new_var/TCGA-生存分析分组DLST+ EGFL7 + DOCK4(1).csv',row.names = 1)
library(tidyverse)
group$Tumor_Sample_Barcode <- rownames(group)
group <- group %>% select(Tumor_Sample_Barcode,Group)

library(tidyverse)
library(maftools)
library(data.table)
mafFilePath <- dir(path = "./SNV/",
                   pattern="masked.maf.gz$",
                   full.names = T,
                   recursive=T)
head(mafFilePath)
maf.list <- lapply(mafFilePath, fread, 
                   skip = 7, 
                   sep = "\t", 
                   header = T)
maf.merge <- do.call(rbind,maf.list)
maf.merge$Tumor_Sample_Barcode <- gsub('01A.*','01A',maf.merge$Tumor_Sample_Barcode)
maf <- read.maf(maf = maf.merge,clinicalData = group)
head(maf@clinical.data)
output <- file.path('D:/乳腺癌/var/new_var')

png(paste(output,'summary.png',sep = '/'),width = 10*300,height = 10*300,res = 300)
pdf(paste(output,'summary.pdf',sep = '/'),width = 10,height = 10)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

#oncoplot for top ten mutated genes.
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
png(paste(output,'oncoplot1.png',sep = '/'),width = 10*300,height = 10*300,res = 300)
pdf(paste(output,'oncoplot1.pdf',sep = '/'),width = 10,height = 10)
oncoplot(maf = maf,
         clinicalFeatures = c("Group"),
         top = 10,
         colors = vc_cols,
         sortByAnnotation=T,borderCol=NULL
)
dev.off()

group1_samples <- group[group$Group == "High_group",]
group1_samples <- group1_samples[maf@clinical.data$Tumor_Sample_Barcode %in% group1_samples$Tumor_Sample_Barcode, ]$Tumor_Sample_Barcode
group2_samples <- group[group$Group == "Low_group",]
group2_samples <- group2_samples[maf@clinical.data$Tumor_Sample_Barcode %in% group2_samples$Tumor_Sample_Barcode, ]$Tumor_Sample_Barcode
getFields(maf)
luad.high <- subsetMaf(maf = maf, tsb = group1_samples, fields = c('Missense_Mutation','Frame_Shift_Del',
                                                                   'In_Frame_Del','Frame_Shift_Ins','Nonsense_Mutation',
                                                                   'Nonstop_Mutation','Splice_Site','Translation_Start_Site','HGVSp_Short',
                                                                   'Protein_position','t_ref_count','t_alt_count'))
luad.low <- subsetMaf(maf = maf, tsb = group2_samples, fields = c('Missense_Mutation','Frame_Shift_Del',
                                                                  'In_Frame_Del','Frame_Shift_Ins','Nonsense_Mutation',
                                                                  'Nonstop_Mutation','Splice_Site','Translation_Start_Site','HGVSp_Short',
                                                                  'Protein_position','t_ref_count','t_alt_count'))
maf_subset <- subset(maf@data, Tumor_Sample_Barcode %in% group1_samples)
table(maf_subset$Variant_Classification)
plotVaf(maf = maf_group1,top = 10)
plotVaf(maf = maf_group2,top = 10)

###手动计算VAF VAF = t_alt_count / (t_ref_count + t_alt_count)
if ("t_ref_count" %in% colnames(luad.high@data) & "t_alt_count" %in% colnames(luad.high@data)) {
  luad.high@data$tumor_vaf <- with(luad.high@data, t_alt_count / (t_ref_count + t_alt_count))
}

if ("t_ref_count" %in% colnames(luad.low@data) & "t_alt_count" %in% colnames(luad.low@data)) {
  luad.low@data$tumor_vaf <- with(luad.low@data, t_alt_count / (t_ref_count + t_alt_count))
}

# 提取高 VAF 组的 VAF 数据
vaf_data_high <- data.frame(
  Gene = luad.high@data$Hugo_Symbol,
  VAF = luad.high@data$tumor_vaf,
  Group = rep("High", nrow(luad.high@data))
)

# 提取低 VAF 组的 VAF 数据
vaf_data_low <- data.frame(
  Gene = luad.low@data$Hugo_Symbol,
  VAF = luad.low@data$tumor_vaf,
  Group = rep("Low", nrow(luad.low@data))
)

# 合并数据
vaf_data_combined <- rbind(vaf_data_high, vaf_data_low)

top10_gene <- c('TP53','PIK3CA','CDH1','MAP3K1','GATA3','KMT2C','MUC16','TTN','HMCN1','FLG')

# 根据选定的基因筛选数据
filtered_vaf_data <- vaf_data_combined[vaf_data_combined$Gene %in% top10_gene,]

# 绘制 VAF 分布图，展示前十的基因
library(ggpubr)
library(ggplot2)
png(paste(output,'VAF.png',sep = '/'),width = 8*300,height = 6*300,res = 300)
pdf(paste(output,'VAF.pdf',sep = '/'),width = 8,height = 6)
ggplot(filtered_vaf_data, aes(x = Gene, y = VAF, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),  # 去除主要网格线
    panel.grid.minor = element_blank(),  # 去除次要网格线
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1,size = 20), # 旋转x轴标签45度
    axis.text = element_text(size = 20),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5) # 调整图表边距
  ) +
  labs(title = "", y = "Tumor VAF", x = "") +
  scale_fill_manual(values = c("High" = "#BF1D2D","Low" = "#293890")) +
  stat_compare_means(aes(label = ..p.signif..), label = "p.signif", label.x = 1) +
  theme_bw() 
dev.off()

library(rstatix)

stat_results <- filtered_vaf_data %>%
  group_by(Gene) %>%
  t_test(VAF ~ Group) %>%
  adjust_pvalue(method = "BH") %>%  # 可选：调整 p 值
  add_significance()

# 查看结果
print(stat_results)
write.csv(stat_results,file = './new_var/tumor_vaf_stat.csv')

# 计算 Ti/Tv 比率
group1_titv <- titv(maf = luad.high, plot = FALSE)
group2_titv <- titv(maf = luad.low, plot = FALSE)

# 假设有多个样本的 Ti/Tv 比率数据
group1_fraction <- group1_titv$TiTv.fractions %>% mutate(Group = 'High')
group2_fraction <- group2_titv$TiTv.fractions %>% mutate(Group = 'Low')


#### C/T
group1_ct <- group1_titv$fraction.contribution %>% mutate(Group = 'High')
group2_ct <- group2_titv$fraction.contribution %>% mutate(Group = 'Low')


# 合并数据
combined_titv_data <- rbind(group1_fraction, group2_fraction)
combined_ct_data <- rbind(group1_ct,group2_ct)

# 将数据转换为长格式
titv_long <- tidyr::pivot_longer(combined_ct_data, cols = c('C>A', 'C>G','C>T','T>C','T>A','T>G'), names_to = "Mutation_Type", values_to = "Value")

# 绘制箱线图
png(paste(output,'mutation.png',sep = '/'),width = 8*300,height = 6*300,res = 300)
pdf(paste(output,'mutation.pdf',sep = '/'),width = 8,height = 6)
p <- ggplot(titv_long, aes(x = Mutation_Type, y = Value, fill = Group)) +
  geom_boxplot(outlier.shape = NA) + # 不显示异常值
  scale_fill_manual(values = c("High" = "#BF1D2D", "Low" = "#293890")) + # 根据Group设置颜色
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(y = "% Mutations", x = "") +
  theme(axis.text.x = element_text(size = 14)) + # 旋转x轴标签45度
  stat_compare_means(aes(label = ..p.signif..), label = "p.signif", method = "t.test") + # 添加显著性标记
theme_bw() 
print(p)
dev.off()

library(rstatix)

stat_results <- titv_long %>%
  group_by(Mutation_Type) %>%
  t_test(Value ~ Group) %>%
  adjust_pvalue(method = "BH") %>%  # 可选：调整 p 值
  add_significance()

# 查看结果
print(stat_results)

# 保存为 CSV
write.csv(stat_results, "./new_var/stat_test_TITV_by_Mutation_Type.csv", row.names = FALSE)


pdf(paste(output,"genecast.vs.know.maf.diff.coBarplot.pdf",sep = '/'),width=5.5,height=5.5)
coBarplot(m1 = luad.high, m2 = luad.low, m1Name = 'High', m2Name = 'Low',colors = vc_cols)
dev.off()

luad.sig <- oncodrive(maf=maf, minMut=5, AACol="HGVSp_Short", pvalMethod="zscore")

# 使用plotOncodrive函数绘制散点图
plotOncodrive(res = luad.sig, fdrCutOff = 0.1, useFraction = TRUE)
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
genes <- c('TP53','PIK3CA','CDH1','MAP3K1','GATA3','KMT2C','MUC16','TTN','HMCN1','FLG')
pdf(paste(output,"genecast.vs.know.maf.diff.coOncoplot.pdf",sep = '/'),width=7,height=5.5)
coOncoplot(m1=luad.high, m2=luad.low, m1Name="High", m2Name="Low", genes=genes,colors = vc_cols)
dev.off()

pdf(paste(output,"genecast.vs.know.maf.diff.collipoplot.pdf",sep = '/'),width=7,height=5.5)
lollipopPlot2(m1=luad.high, m2=luad.low,m1_name="High", m2_name="Low", gene="PIK3CA",AACol1 = "HGVSp_Short", AACol2 = "HGVSp_Short")
dev.off()

##########TMB###############

get_TMB <- function(file) {
  # 需要用到的列
  use_cols <- c(
    "Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode", 
    "HGVSc", "t_depth", "t_alt_count"
  )
  # 删除这些突变类型
  mut_type <- c(
    "5'UTR", "3'UTR", "3'Flank", "5'Flank", "Intron", "IGR",
    "Silent", "RNA", "Splice_Region"
  )
  # 读取文件
  df <- file
  data <- df %>% subset(!Variant_Classification %in% mut_type) %>%
    # 计算 VAF
    mutate(vaf = t_alt_count / t_depth) %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(mut_num = n(), TMB = mut_num / 30, MaxVAF = max(vaf))
  return(data)
}
#计算
tmb_table_wt_log = tmb(maf = maf1_list)
quantile(tmb_table_wt_log$total_perMB)



paad.filtered <- get_TMB(maf1)
max(paad.filtered$TMB)
#计算tmb值
paad.filtered <- as.data.frame(paad.filtered)
rownames(paad.filtered) <- paad.filtered$Tumor_Sample_Barcode
clin <- maf1_list@clinical.data
id = intersect(clin$Tumor_Sample_Barcode,paad.filtered$Tumor_Sample_Barcode)
clin <- as.data.frame(clin)
rownames(clin) <- clin$Tumor_Sample_Barcode
clin <- clin[id,]
paad.filtered <- paad.filtered[id,]
paad.filtered <- inner_join(paad.filtered,clin,by = 'Tumor_Sample_Barcode')

library(ggplot2)
library(ggpubr)

paad.filtered <- paad.filtered %>% 
  filter(TMB != max(TMB))
max(paad.filtered$TMB)

ggboxplot(paad.filtered, x = "Group", y = "TMB",
          fill = "Group") +
  stat_compare_means(comparisons = list(
    c("High", "Low")
  )) +
  stat_compare_means(label.y = 4.5) 








