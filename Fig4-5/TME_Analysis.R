library(CIBERSORT)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggpubr)


sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
expr <- read.csv('./var/new_var/免疫浸润_new/fitermRNA.csv',row.names = 1,check.names = F)
group <- read.csv('./var/new_var//TCGA-生存分析分组DLST+ EGFL7 + DOCK4(1).csv',row.names = 1)
library(tidyverse)
group$Tumor_Sample_Barcode <- rownames(group)
group <- group %>% select(Tumor_Sample_Barcode,Group)

GeneSymbol <- rownames(expr)
expr <- cbind(GeneSymbol,expr)
exp2 <- sweep(expr[,-1],1,apply(expr[,-1],1,mean,na.rm=T))
exp2 <- exp2[!duplicated(exp2),]
exp2 <- exp2[match(rownames(exp2),GeneSymbol),]
GeneSymbol <- GeneSymbol[match(rownames(exp2),GeneSymbol)]
exp2 <- cbind(GeneSymbol,exp2)
write.table(expr,"./var/new_var/免疫浸润_new/exp.txt",row.names=FALSE,sep="\t",quote=FALSE)

res <- cibersort(sig_matrix = sig_matrix,'./var/new_var/免疫浸润_new/exp.txt',
                 perm = 50,QN = T
)
res1 <- as.data.frame(res[,1:22])
#write.csv(res1,'./var/new_var/免疫浸润_new/cibersort.csv')
#Group <- read.csv('./group_list.csv',row.names = 1)
#Group <- Group$condition
#col_sum <- apply(res1, 2, sum)
#col_sum_index <- which(col_sum > 10)
#res1 <- res1[,col_sum_index]  
#res1 <- read.csv('./var/new_var/免疫浸润_new/cibersort.csv',row.names = 1)
group <- group %>% mutate(Group = case_when(Group == 'High_group' ~ 'High',
                                            Group == 'Low_group' ~ 'Low'))

Group <- group$Group
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
dat <- res1 %>% 
  as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  mutate(group = Group) %>% 
  gather(key = Cell_type,value = Proportion,-Sample,-group) %>% 
  arrange(group)

dat$Sample = factor(dat$Sample,ordered = T,levels = unique(dat$Sample)) #定横坐标顺序
# 先把group排序，然后将sample设为了因子，确定排序后的顺序为水平，所以两图的顺序是对应的。
dat2 = data.frame(a = 1:ncol(expr[,-1]),
                  b = 1, 
                  group = sort(Group))
output <- file.path('./var/new_var/免疫浸润_new/')

####### 堆叠柱状图
library(ggh4x)
library(reshape2)
library(ggalluvial)
library(ggplot2)

p1 = ggplot(dat2,aes(x = a, y = b)) + 
  geom_tile(aes(fill = group)) + 
  scale_fill_manual(values = c("#BF1D2D","#293890")) +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        axis.title = element_blank()) + 
  scale_x_continuous(expand = c(0, 0)) +
  labs(fill = "Group")
p1
p2 = ggplot(dat,aes(Sample, Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") + 
  labs(fill = "Cell Type",x = "",y = "Cell Proportion") + xlab(NULL)+
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 4)
  ) + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(22))
p2
library(patchwork)
p3 <- p1 / p2 + plot_layout(heights = c(1,10),guides = "collect" ) &
  theme(legend.position = "bottom")
p3
ggsave(plot = p3,filename = paste(output,'堆叠图.png',sep = '/'),width = 12,height = 7,dpi = 600)
ggsave(plot = p3,filename = paste(output,'堆叠图.pdf',sep = '/'),width = 12,height = 7)

p <- ggplot(dat, aes(x = Sample, y = Proportion, fill = Cell_type,
                     stratum = Cell_type, alluvium = Cell_type)) +
  geom_col(position = 'stack', width = 0.6) +
  geom_stratum(width = 0.6, color = 'white') +
  geom_alluvium(alpha = 0.4, width = 0.6, color = 'white', linewidth = 1, curve_type = "linear") +
  scale_fill_manual(values = mypalette) +
  xlab('') + 
  ylab('') +
  scale_y_continuous(expand = c(0, 0))+
  theme_bw(base_size = 12) + 
  theme(
    axis.text = element_text(color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )
p

ggsave(plot = p,filename = paste(output,'堆叠图.png',sep = '/'),width = 12,height = 10,dpi = 300)
ggsave(plot = p,filename = paste(output,'堆叠图.pdf',sep = '/'),width = 12,height = 10)




#############箱线图
res2 <- data.frame(res[,1:22])%>%
  mutate(group = Group)%>%
  rownames_to_column("sample")%>%
  pivot_longer(cols = colnames(.)[2:23],
               names_to = "cell.type",
               values_to = 'value')


library(ggplot2)
library(ggpubr)
library(RColorBrewer)
theme<-theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=13,face = 'bold'),
        axis.text.y=element_text(size=12,face = 'bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=16,face = 'bold'),
        axis.line=element_line(size=1),
        plot.title=element_blank(),
        legend.text=element_text(size=16,face = 'bold'),
        legend.key=element_rect(fill='transparent'),
        legend.background=element_rect(fill='transparent'),
        legend.position="top",
        legend.title=element_blank())
p<-ggplot(dat,aes(x=Cell_type,y=Proportion,fill=factor(group)))+
  geom_violin(position=position_dodge(width=1),scale='width')+
  geom_boxplot(position=position_dodge(width=1),outlier.shape=NA,
               width=0.25,alpha=0.2,show.legend=FALSE)+
  #scale_fill_manual(values=alpha(brewer.pal(8,"Set1")[1:2],0.6))+
  scale_fill_manual(values = c("High" = "#BF1D2D", "Low" = "#293890"))+
  scale_y_continuous(expand=c(0,0),limits=c(0,0.8),breaks=seq(0,0.8,0.1),
                     labels=seq(0,0.8,0.1))+
  theme+
  labs(y='Immune Cell Relative Proportion',fill = NULL)
p <- ggplot(dat, aes(x = Cell_type, y = Proportion, fill = factor(group))) +
  geom_boxplot(position = position_dodge(width = 1), outlier.shape = NA,
               width = 0.7, alpha = 0.8) +
  scale_fill_manual(values = c("High" = "#BF1D2D", "Low" = "#293890")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.8), breaks = seq(0, 0.8, 0.1),
                     labels = seq(0, 0.8, 0.1)) +
  theme+
  labs(y = 'Immune Cell Relative Proportion', fill = NULL)

###添加显著性检验
Data_summary<-as.data.frame(compare_means(Proportion~group,dat,method="wilcox.test",
                                          paired=FALSE,group.by="Cell_type"))

stat.test<-Data_summary[,-c(2,9)]
stat.test$xmin=c(1:22)-0.2
stat.test$xmax=c(1:22)+0.2
stat.test$p.signif[stat.test$p.signif=="ns"]<-NA
write.csv(stat.test,'./var/new_var/免疫浸润_new/cibersort_box_stat.csv')
p2<-p+geom_signif(xmin=stat.test$xmin,xmax=stat.test$xmax,
                  annotations=stat.test$p.signif,margin_top=0.00,
                  y_position=0.72,size=0.5,textsize=7.5,
                  tip_length=0)#横线两侧折线长度
p2

ggsave(plot = p2,filename = paste(output,'箱线图.png',sep = '/'),width = 12,height = 8,dpi = 600)
ggsave(plot = p2,filename = paste(output,'箱线图.pdf',sep = '/'),width = 12,height = 8)



genelist <- c("DLST","EGFL7",'DOCK4')
#提取基因集的表达矩阵
goal_exp<-filter(expr,rownames(expr) %in% genelist)
goal_exp <- goal_exp[,-1]
#goal_exp <- as.data.frame(expr)
combine<-rbind(goal_exp,as.data.frame(t(res1)))
combine <- as.data.frame(t(apply(combine, 1, as.numeric)))
colnames(combine)  <- colnames(goal_exp)
comcor<-cor(t(combine))

imm <- res1
exp <- expr[,-1]
#相关性计算
immuscore <- function(gene){  y <- as.numeric(exp[gene,])  
colnames <- colnames(imm)  
do.call(rbind,lapply(colnames, function(x){    
  dd  <- cor.test(as.numeric(imm[,x]),y,type="spearman")    
  data.frame(gene=gene,immune_cells=x,
             cor=dd$estimate,
             p.value=dd$p.value )}))}

gene <- genelist
data1 <- do.call(rbind,lapply(gene,immuscore))
head(data1)
write.csv(data1, "./var/new_var/免疫浸润_new/cibersort_correlation.csv", quote = F, row.names = F)


data1$pstar <- ifelse(data1$p.value < 0.05,ifelse(data1$p.value < 0.01,
                                                  ifelse(data1$p.value < 0.001,"***", "**"), "*"), "")

#相关性热图的绘制

cor_plot <- ggplot(data1, aes(immune_cells, gene)) + 
  geom_tile(aes(fill = cor), colour = "grey",size=1) + 
  scale_fill_gradient2(low = "#293890",mid = "white",high = "#BF1D2D") + 
  geom_text(aes(label=pstar),col ="black",size = 5) + 
  theme_minimal() + 
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1,size = 20,family = 'SimHei'),
        axis.text.y = element_text(size = 20))+
  labs(fill =paste0('Pvalue','\n\n',"* p < 0.05","\n\n","** p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation"))+
  coord_flip()+theme_bw()
cor_plot
ggsave(plot = cor_plot,filename = paste(output,'相关性图.png',sep = '/'),width = 9,height = 9,dpi = 600)
ggsave(plot = cor_plot,filename = paste(output,'相关性图.pdf',sep = '/'),width = 12,height = 10)

library(utils)
rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
exp <- read.csv('./var/new_var/免疫浸润_new/fitermRNA.csv',row.names = 1,check.names = F)
#Group <- read.csv('./group_list.csv',row.names = 1)
Group <- read.csv('./var/new_var/TCGA-生存分析分组DLST+ EGFL7 + DOCK4(1).csv',row.names = 1)
library(tidyverse)
Group$Tumor_Sample_Barcode <- rownames(Group)
Group <- group %>% select(Tumor_Sample_Barcode,Group)

estimate <- function(dat,pro){
  input.f=paste0(pro,'_estimate_input.txt')  
  output.f=paste0(pro,'_estimate_gene.gct')  
  output.ds=paste0(pro,'_estimate_score.gct')  
  write.table(dat,file = input.f,sep = '\t',quote = F)  
  library(estimate)  
  filterCommonGenes(input.f=input.f,
                    output.f=output.f,
                    id="GeneSymbol")  
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="illumina")  ## platform  
  scores=read.table(output.ds,
                    skip = 2,
                    header = T,
                    check.names = F)
  rownames(scores)=scores[,1]  
  scores=t(scores[,3:ncol(scores)]) 
  library(stringr)  
  rownames(scores)=str_replace_all(rownames(scores),'[.]','-') # 这里TCGA样本名里面的-变成.了，进行恢复  
  write.csv(scores,file="Stromal_Immune_ESTIMATE.Score.csv") #  
  return(scores)
}

pro = 'BRCA'
scores = estimate(exp,pro)
scores <- as.data.frame(scores)
scores$Group = Group$Group
scores$TumorPurity = cos(0.6049872018+0.0001467884 * scores[,3])

#boxplot
library(ggpubr)
library(ggsci)


p1 <- ggplot(scores, aes(x = Group, y = TumorPurity, fill = Group)) +
  geom_boxplot(position = position_dodge(0.8),
               outlier.shape = NA,  # 确保异常值不显示
               outlier.size = 0) +  # 进一步确保异常值不显示
  scale_fill_nejm() +
  labs(x = "", y = 'Tumor Purity') +
  scale_fill_manual(values = c("High" = "#BF1D2D", "Low" = "#293890"))+
  stat_compare_means(
    comparisons = combn(unique(scores$Group), 2, simplify = FALSE),
    label = "p.signif",  # 使用星号表示显著性水平
    method = "t.test"  # 选择合适的统计方法
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.title.x = element_blank(),  # 删除 x 轴标题
    axis.text.x = element_text(angle = 45, hjust = 1)  # 如果需要调整 x 轴标签角度
  )
p1

ttest_result <- t.test(TumorPurity ~ Group, data = scores)

# 创建结果表格
results_table <- data.frame(
  Test = "Student's t-test",
  Group1 = "High",
  Group2 = "Low",
  t_statistic = ttest_result$statistic,
  df = ttest_result$parameter,
  p_value = ttest_result$p.value,
  Alternative = ttest_result$alternative,
  group = 'tumtorPurity'
)

# 显示表格
print(results_table)


p2 <- ggplot(scores, aes(x = Group, y = StromalScore, fill = Group)) +
  geom_boxplot(position = position_dodge(0.8),
               outlier.shape = NA,  # 确保异常值不显示
               outlier.size = 0) +  # 进一步确保异常值不显示
  scale_fill_nejm() +
  labs(x = "", y = 'StromalScore') +
  scale_fill_manual(values = c("High" = "#BF1D2D", "Low" = "#293890"))+
  stat_compare_means(
    comparisons = combn(unique(scores$Group), 2, simplify = FALSE),
    label = "p.signif",  # 使用星号表示显著性水平
    method = "t.test"  # 选择合适的统计方法
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.title.x = element_blank(),  # 删除 x 轴标题
    axis.text.x = element_text(angle = 45, hjust = 1)  # 如果需要调整 x 轴标签角度
  )
p2

# 进行独立样本 t 检验
ttest_result <- t.test(StromalScore ~ Group, data = scores)

# 提取关键统计量并整理成表格
result_table1 <- data.frame(
  Test = "Student's t-test",
  Group1 = "High",
  Group2 = "Low",
  t_statistic = ttest_result$statistic,
  df = ttest_result$parameter,
  p_value = ttest_result$p.value,
  Alternative = ttest_result$alternative,
  group = 'StromalScore'
)

# 查看结果
print(result_table)


p3 <- ggplot(scores, aes(x = Group, y = ImmuneScore, fill = Group)) +
  geom_boxplot(position = position_dodge(0.8),
               outlier.shape = NA,  # 确保异常值不显示
               outlier.size = 0) +  # 进一步确保异常值不显示
  scale_fill_nejm() +
  labs(x = "", y = 'ImmuneScore') +
  scale_fill_manual(values = c("High" = "#BF1D2D", "Low" = "#293890"))+
  stat_compare_means(
    comparisons = combn(unique(scores$Group), 2, simplify = FALSE),
    label = "p.signif",  # 使用星号表示显著性水平
    method = "t.test"  # 选择合适的统计方法
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.title.x = element_blank(),  # 删除 x 轴标题
    axis.text.x = element_text(angle = 45, hjust = 1)  # 如果需要调整 x 轴标签角度
  )
p3
# 进行独立样本 t 检验
ttest_result <- t.test(ImmuneScore ~ Group, data = scores)

# 提取关键统计量并整理成表格
result_table2 <- data.frame(
  Test = "Student's t-test",
  Group1 = "High",
  Group2 = "Low",
  t_statistic = ttest_result$statistic,
  df = ttest_result$parameter,
  p_value = ttest_result$p.value,
  Alternative = ttest_result$alternative,
  group = 'ImmuneScore'
)

result_estimate <- rbind(result_table1,result_table2,results_table)
write.csv(result_estimate,file = './var/new_var/免疫浸润_new/estimate_stat.csv')

# 查看结果
print(result_table)

output <- file.path('./新建文件夹')
ggsave(plot = p3,filename = paste(output,'免疫细胞评分.pdf',sep = '/'),height = 8,width = 10)
ggsave(plot = p3,filename = paste(output,'免疫细胞评分.png',sep = '/'),height = 8,width = 10,dpi = 300)

ggsave(plot = p2,filename = paste(output,'基质细胞评分.pdf',sep = '/'),height = 8,width = 10)
ggsave(plot = p2,filename = paste(output,'基质细胞评分.png',sep = '/'),height = 8,width = 10,dpi = 300)

ggsave(plot = p1,filename = paste(output,'肿瘤纯度评分.pdf',sep = '/'),height = 8,width = 10)
ggsave(plot = p1,filename = paste(output,'肿瘤纯度评分.png',sep = '/'),height = 8,width = 10,dpi = 300)

########TIDE#######################
############TIDE网站得到计算结果#######################
exp2 <- sweep(expr,1,apply(expr,1,mean,na.rm=T))
#write.table(exp2,"TIDE.txt",sep="\t",quote=F)
  
#TIDE <- as.data.frame(read.table('./var/new_var/免疫浸润_new/TIDE.csv',header = T,sep = '\t',check.names = F,row.names = 1))
TIDE <- read.csv('./var/new_var/免疫浸润_new/TIDE.csv',row.names = 1)
TIDE$id <- rownames(TIDE)
Group <- group_list
TIDE <- TIDE[match(rownames(Group),rownames(TIDE)),]
TIDE$Group <- Group$condition
group = levels(factor(TIDE$Group))
comp = combn(group,2)
my_comparisons = list()
for (i in 1:ncol(comp)) {
  my_comparisons[[i]] <- comp[,i]
}


library(ggplot2)
library(ggpubr)
p <- ggviolin(TIDE,x = 'Group',y = 'TIDE',fill = 'Group',
              xlab = '',ylab = 'TIDE_Scores',
              palette = c('Firebrick2','DodgerBlue1'),
              legend.title = 'Group',
              add = 'boxplot',add.params = list(fill = 'white')) +
  stat_compare_means(comparisons = my_comparisons)
p3 <- ggplot(TIDE, aes(x = Group, y = TIDE, fill = Group)) +
  geom_boxplot(position = position_dodge(0.8),
               outlier.shape = NA,  # 确保异常值不显示
               outlier.size = 0) +  # 进一步确保异常值不显示
  scale_fill_nejm() +
  labs(x = "", y = 'TIDE Scores') +
  scale_fill_manual(values = c("High_group" = "#BF1D2D", "Low_group" = "#293890"))+
  stat_compare_means(
    comparisons = combn(unique(TIDE$Group), 2, simplify = FALSE),
    label = "p.signif",  # 使用星号表示显著性水平
    method = "t.test"  # 选择合适的统计方法
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.title.x = element_blank(),  # 删除 x 轴标题
    axis.text.x = element_text(angle = 45, hjust = 1)  # 如果需要调整 x 轴标签角度
  )
p3


# 进行独立样本 t 检验
ttest_result <- t.test(TIDE ~ Group, data = TIDE)

# 提取关键统计量并整理成表格
result_table <- data.frame(
  Test = "Student's t-test",
  Group1 = "High_group",
  Group2 = "Low_group",
  t_statistic = ttest_result$statistic,
  df = ttest_result$parameter,
  p_value = ttest_result$p.value,
  Alternative = ttest_result$alternative
)

# 查看结果
print(result_table)

write.csv(result_table,file = './var/new_var/免疫浸润_new/TIDE_stat.csv')


ggsave(plot = p,filename = './新建文件夹/TIDE评分.pdf',height = 8,width = 10)
ggsave(plot = p,filename = './新建文件夹/TIDE评分.png',height = 8,width = 10,dpi = 300)

res <- as.data.frame(read.csv('./TIDE.csv',row.names = 1))
res <- res[match(rownames(Group),rownames(res)),]
res$Group <- Group$Group
table(res$Responder,res$Group)
f=fisher.test(table(res$Responder,res$Group))
label=paste("fisher.test pvalue=",round(f$p.value,3))
label

library(ggplot2)
library(dplyr)
res=arrange(res,desc(TIDE))
p1=ggplot(res,aes(x=1:nrow(res),
                  y=TIDE,
                  fill=Responder))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#e04030","#6cb8d2"))+
  xlab("patient")+
  ylab("TIDE value")+
  annotate("text",x=300,y=-1.1,label=label,size=5)+
  theme_bw()+
  theme(legend.position="none") +#把P1的图注去掉了
  theme(axis.title.x = element_text(size = 14,face = 'bold'),
        axis.title.y = element_text(size = 14,face = 'bold'),
        axis.text.x = element_text(size = 10,face = 'bold'),
        axis.text.y = element_text(size = 10,face = 'bold'))
p1

########### ICB Scores #############
library(IMvigor210CoreBiologies)
#devtools::install_local("C:/Users/Administrator.DESKTOP-9U8QK3S/Downloads/IOBR-master.zip")
library(IOBR)
library(easier)
data <- read.csv('./var/new_var/免疫浸润_new/fitermRNA.csv',row.names = 1)
group_list$SampleID <- rownames(group_list)
colnames(data) <- rownames(group_list)
ips <- deconvo_tme(eset = data, method = "ips", plot= FALSE)
ips$group <- group_list$condition

hallmarks_of_immune_response <- c("CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS")

immune_response_scores <- compute_scores_immune_response(RNA_tpm = data, 
                                                         selected_scores = hallmarks_of_immune_response)
write.csv(immune_response_scores,file = './var/new_var/免疫浸润_new/immuneresponsescores_tcga.csv')


identical(rownames(immune_response_scores),group_list$SampleID)
library(tidyverse)
immune_response_scores = immune_response_scores %>% mutate(SampleID = group_list$SampleID,Group = group_list$condition)
df_long <- immune_response_scores %>%
  pivot_longer(
    cols = -c(SampleID, Group),  # 排除ID和分组列
    names_to = "ScoreType",
    values_to = "ScoreValue"
  )
df_long <- df_long %>% mutate(Group = case_when(Group == 'High_group' ~ 'High',
                                                Group == 'Low_group' ~ 'Low',
                                                TRUE ~ NA))
immune_boxplot <- ggplot(df_long, aes(x = ScoreType, y = ScoreValue, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +  # 隐藏离群点
  #geom_jitter(width = 0.2, alpha = 0.5, aes(color = Group)) +  # 添加数据点
  stat_compare_means(aes(group = Group), 
                     method = "wilcox.test", 
                     label = "p.signif") +  # 添加统计检验
  scale_fill_brewer(palette = "Set2") +
  labs(title = "") + ylab('Scores') + xlab('') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
immune_boxplot
ggsave(plot = immune_boxplot,filename = './var/new_var/免疫浸润_new/immune_type_tcga.pdf',height = 6,width = 8)
library(tidyr)
library(dplyr)

calculate_PCA_ICB <- function(data) {
  icb_features <- c("CYT", "TLS", "IFNy", "Ayers_expIS", "Tcell_inflamed",
                    "Roh_IS", "Davoli_IS", "chemokines")
  
  available_features <- intersect(icb_features, colnames(data))
  
  # 执行PCA
  pca_result <- prcomp(data[, available_features], scale. = TRUE)
  
  # 使用第一主成分作为ICB评分
  data$ICB_Score_PC1 <- pca_result$x[,1]
  
  # 解释方差
  var_explained <- summary(pca_result)$importance[2,1] * 100
  message(sprintf("PC1 explains %.1f%% of variance", var_explained))
  
  return(data)
}

df_with_PCA_ICB <- calculate_PCA_ICB(immune_response_scores)
df_with_PCA_ICB <- df_with_PCA_ICB[,c(12,13,14)]
# 将宽格式转换为长格式
df_with_PCA_ICB <- df_with_PCA_ICB %>% mutate(Group = case_when(Group == 'High_group' ~ 'High',
                                                Group == 'Low_group' ~ 'Low',
                                                TRUE ~ NA))
write.csv(df_with_PCA_ICB,file = './var/new_var/免疫浸润_new/ICBScores_tcga.csv')
advanced_boxplot <- ggplot(df_with_PCA_ICB, aes(x = Group, y = ICB_Score_PC1, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +  # 隐藏离群点
  #geom_jitter(width = 0.2, alpha = 0.5, aes(color = Group)) +  # 添加数据点
  stat_compare_means(aes(group = Group), 
                     method = "wilcox.test", 
                     label = "p.signif") +  # 添加统计检验
  scale_fill_brewer(palette = "Set2") +
  labs(title = "") + ylab('ICB Scores') + xlab('') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
advanced_boxplot
ggsave(plot = advanced_boxplot,filename = './var/new_var/免疫浸润_new/ICBScores_tcga.pdf',height = 6,width = 6.5)

library(tidyverse)
library(data.table)
exp_seq <- read_delim("./ICGC/exp_array.BRCA-FR.tsv.gz", "\t", escape_double = FALSE, trim_ws = TRUE)

exp <- unique(exp_seq[,c("icgc_donor_id","gene_id","normalized_expression_value")])
colnames(exp) <- c("icgc_donor_id","gene_id","normalized_expression_value")
dat <- as.data.frame(dcast(data.table(exp),gene_id~icgc_donor_id,
             value.var="normalized_expression_value",fun.aggregate = max))
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("./ICGC/gencode.v38lift37.annotation.gtf.gz",format="gtf")
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
gene <- as.matrix(exons_gene_lens) %>% as.data.frame() %>% rownames_to_column()
colnames(gene)<-c('gene_id','length')
gene$gene_id <- str_split(gene$gene,"[.]",simplify = T)[,1] # 删除版本的“.”和后的数字
gene$length <- as.numeric(gene$length)
gene <- unique(gene) # 去重
library(org.Hs.eg.db)
gene$symbol <- mapIds(org.Hs.eg.db,
                      keys = gene$gene_id,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")
colnames(gene) <- c('ENSID','length','gene_id')
tmp <- inner_join(dat,gene,by = 'gene_id') # 匹配gene_count与gene_length
tmp<- tmp[!duplicated(tmp$gene_id),]

gene <- gene[match(tmp$gene_id,gene$gene_id),]
gene <- gene[,-1]
tmp <- inner_join(dat,gene,by = 'gene_id') # 匹配gene_count与gene_length
# count2tpm
tpm <- data.frame(row.names = tmp$gene_id)
for (i in 2:(dim(tmp)[2]-1)){
  col <- tmp[[i]]
  len <- tmp[[dim(tmp)[2]]]
  rate <- col/len
  N <- sum(col) # 计算每个样本的mapped reads数
  TPMi <- (rate*1e6)/(sum(rate)) # 计算tpm值
  print(sum(rate))
  TPMi <- pmax(TPMi,0) %>% as.data.frame() # 去掉矫正带来的负值
  colnames(TPMi) <- colnames(tmp)[i]
  tpm <- cbind(tpm,TPMi)
}
dat <- tpm
write.csv(dat, "./ICGC/row_count.csv",row.names = F)

phe <- read_delim("./ICGC/donor.BRCA-FR.tsv.gz", "\t", escape_double = FALSE, trim_ws = TRUE)
phe <- phe[match(colnames(dat),phe$icgc_donor_id),]
write.csv(phe, "./ICGC/phe.csv",row.names = F)

dat2 <- as.data.frame(t(dat))
target_genes <- c("DLST", "EGFL7", "DOCK4")
gene_medians <- apply(dat2[, target_genes], 2, median)
is_high <- t(apply(dat2[, target_genes], 1, function(x) x > gene_medians))
high_count <- rowSums(is_high)
dat2$Group <- ifelse(high_count >= 2, "High", "Low")
table(dat2$Group)
phe <- phe[match(rownames(dat2),phe$icgc_donor_id),]
phe$Group <- dat2$Group
age_means <- phe %>%
  group_by(Group) %>%
  summarise(
    Mean_Age = mean(donor_age_at_diagnosis, na.rm = TRUE),
    SD_Age = sd(donor_age_at_diagnosis, na.rm = TRUE),
    N = n()
  )
print(age_means)
#write.csv(dat2, "./ICGC/count.csv",row.names = F)
#write.csv(phe, "./ICGC/metadata.csv",row.names = F)
tnm_data <- data.frame(group = phe$Group,tnm = phe$donor_tumour_stage_at_diagnosis)

immune_response_scores <- compute_scores_immune_response(RNA_tpm = dat, 
                                                         selected_scores = hallmarks_of_immune_response)
write.csv(immune_response_scores,file = './var/new_var/免疫浸润_new/immuneresponsescores_icgc.csv')

identical(rownames(immune_response_scores),phe$icgc_donor_id)
library(tidyverse)
immune_response_scores = immune_response_scores %>% mutate(SampleID = phe$icgc_donor_id,Group = dat2$Group)
df_long <- immune_response_scores %>%
  pivot_longer(
    cols = -c(SampleID, Group),  # 排除ID和分组列
    names_to = "ScoreType",
    values_to = "ScoreValue"
  )

immune_boxplot <- ggplot(df_long, aes(x = ScoreType, y = ScoreValue, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +  # 隐藏离群点
  #geom_jitter(width = 0.2, alpha = 0.5, aes(color = Group)) +  # 添加数据点
  stat_compare_means(aes(group = Group), 
                     method = "wilcox.test", 
                     label = "p.signif") +  # 添加统计检验
  scale_fill_brewer(palette = "Set2") +
  labs(title = "") + ylab('Scores') + xlab('') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
immune_boxplot
ggsave(plot = immune_boxplot,filename = './var/new_var/免疫浸润_new/immune_type_icgc.pdf',height = 6,width = 8)
library(tidyr)
library(dplyr)

calculate_PCA_ICB <- function(data) {
  icb_features <- c("CYT", 'TLS',"IFNy", "Ayers_expIS", "Tcell_inflamed",
                    "Roh_IS", "Davoli_IS", "chemokines")
  
  available_features <- intersect(icb_features, colnames(data))
  
  # 执行PCA
  pca_result <- prcomp(data[, available_features], scale. = TRUE)
  
  # 使用第一主成分作为ICB评分
  data$ICB_Score_PC1 <- pca_result$x[,1]
  
  # 解释方差
  var_explained <- summary(pca_result)$importance[2,1] * 100
  message(sprintf("PC1 explains %.1f%% of variance", var_explained))
  
  return(data)
}

df_with_PCA_ICB <- calculate_PCA_ICB(immune_response_scores)
df_with_PCA_ICB <- df_with_PCA_ICB[,c(12,13,14)]
write.csv(df_with_PCA_ICB,file = './var/new_var/免疫浸润_new/ICBScores_icgc.csv')
advanced_boxplot <- ggplot(df_with_PCA_ICB, aes(x = Group, y = ICB_Score_PC1, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +  # 隐藏离群点
  #geom_jitter(width = 0.2, alpha = 0.5, aes(color = Group)) +  # 添加数据点
  #stat_compare_means(aes(group = Group), 
   #                  method = "wilcox.test", 
    #                 label = "p.signif") +  # 添加统计检验
  scale_fill_brewer(palette = "Set2") +
  labs(title = "") + ylab('ICB Scores') + xlab('') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
advanced_boxplot
ggsave(plot = advanced_boxplot,filename = './var/new_var/免疫浸润_new/ICBScores_icgc.pdf',height = 6,width = 6.5)


