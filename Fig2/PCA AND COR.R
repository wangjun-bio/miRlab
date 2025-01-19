tRNA <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\new1\\filter_tRNA.csv',row.names = 1)
rsRNA <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\new1\\filter_rsRNA.csv',row.names = 1)
ysRNA <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\new1\\filter_ysRNA.csv',row.names = 1)
mRNA <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\new1\\filter_mRNA.csv',row.names = 1)
miRNA <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\new1\\filter_miRNA.csv',row.names = 1)
lncRNA <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\new1\\filter_lncRNA.csv',row.names = 1)
piRNA <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\new1\\filter_piRNA.csv',row.names = 1)

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
all_rpm <- rbind(lnc_rpm,m_rpm,mi_rpm,t_rpm,rs_rpm,ys_rpm,pi_rpm)

write.table(rs_rpm,file = 'filter_rsRNA_rpm.csv',sep = ',',col.names = NA)
write.table(ys_rpm,file = 'filter_ysRNA_rpm.csv',sep = ',',col.names = NA)
write.table(t_rpm,file = 'filter_tRNA_rpm.csv',sep = ',',col.names = NA)
write.table(m_rpm,file = 'filter_mRNA_rpm.csv',sep = ',',col.names = NA)
write.table(mi_rpm,file = 'filter_miRNA_rpm.csv',sep = ',',col.names = NA)
write.table(lnc_rpm,file = 'filter_lncRNA_rpm.csv',sep = ',',col.names = NA)
write.table(pi_rpm,file = 'filter_piRNA_rpm.csv',sep = ',',col.names = NA)

group <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\group2.csv')
ids <- as.data.frame(group$group)
rownames(ids) <- group$X
colnames(ids) <- 'group'


data <- as.data.frame(t(rbind(ids$group,all_rpm)))
colnames(data)[1] <- 'group'

library(tidyverse)
BRE_BEN <- data %>% filter(group == 'BEN')
BRE_MAL <- data %>% filter(group == 'MAL')

BRE_BEN <- BRE_BEN %>% select(-group)
BRE_BEN <- as.data.frame(t(BRE_BEN))
BEN <- rownames(BRE_BEN)
BRE_BEN <- apply(BRE_BEN, 2, as.numeric)
rownames(BRE_BEN) <- BEN
rowmeans_BEN <- as.data.frame(rowMeans(BRE_BEN))
lben <- log2(rowmeans_BEN+1)

BRE_MAL <- BRE_MAL %>% select(-group)
BRE_MAL <- as.data.frame(t(BRE_MAL))
MAL <- rownames(BRE_MAL)
BRE_MAL <- apply(BRE_MAL, 2, as.numeric)
rownames(BRE_MAL) <- MAL
rowmeans_MAL <- as.data.frame(rowMeans(BRE_MAL))
lmal <- log2(rowmeans_MAL+1)

data1 <- cbind(lben,lmal)
colnames(data1)[1] <- 'BEN'
colnames(data1)[2] <- 'MAL'
###############correlation##########
library(ggpubr)
scatter <- ggscatter(data1, x = "MAL", y = "BEN",
                     xlab = "log2(RPM+1) of BRE-MAL(n=41)", ylab = "log2(RPM+1) of BEN(n=42)",
                     cor.coef = TRUE, # 显示相关性系数
                     cor.coeff.args = list(method = "pearson"), # 选择相关性方法
                     conf.int = TRUE, # 显示置信区间
                     add = "reg.line",# 添加回归线
                     color = 'red'
) + ggtitle("allRNA") +
  theme(plot.title = element_text(hjust = 0.5))  # 将标题水平居中
scatter

ggsave(scatter,filename = 'cor of allRNA.pdf',width = 8,height = 6)
library(ggplot2)

library(FactoMineR)
library(factoextra)
######################PCA##########
data_pca <- PCA(t(all_rpm),graph = F)
p <- fviz_pca_ind(data_pca,geom.ind = 'point',col.ind = ids$group,palette = 'Dark2',addEllipses = T,legend.title='Group',title="PCA of allRNA",
                  mean.point =F,ellipse.level = 0.7,pointsize =3)+
  theme_bw()+theme(plot.title = element_text(size = 16,hjust = 0.5))+
theme(text=element_text(size=14,face="plain",color="black"),
      axis.title=element_text(size=16,face="plain",color="black"),
      axis.text = element_text(size=14,face="plain",color="black"),
      legend.title = element_text(size=16,face="plain",color="black"),
      legend.text = element_text(size=14,face="plain",color="black"),
      legend.background = element_blank()
)
p
ggsave(p,filename = 'PCA of allRNA.pdf',width = 8,height = 6)




















