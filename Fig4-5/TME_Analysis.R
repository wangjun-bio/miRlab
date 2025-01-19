library(CIBERSORT)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggpubr)


sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
expr <- read.csv('./fitermRNA.csv',row.names = 1,check.names = F)
GeneSymbol <- rownames(expr)
expr <- cbind(GeneSymbol,expr)
exp2 <- sweep(expr[,-1],1,apply(expr[,-1],1,mean,na.rm=T))
exp2 <- exp2[!duplicated(exp2),]
exp2 <- exp2[match(rownames(exp2),GeneSymbol),]
GeneSymbol <- GeneSymbol[match(rownames(exp2),GeneSymbol)]
exp2 <- cbind(GeneSymbol,exp2)
write.table(expr,"./新建文件夹/exp.txt",row.names=FALSE,sep="\t",quote=FALSE)

res <- cibersort(sig_matrix = sig_matrix,'./新建文件夹/exp.txt',
                 perm = 50,QN = T
)
res1 <- as.data.frame(res[,1:22])
#write.csv(res1,'cibersort.csv')
Group <- read.csv('./group_list.csv',row.names = 1)
Group <- Group$condition
#col_sum <- apply(res1, 2, sum)
#col_sum_index <- which(col_sum > 10)
#res1 <- res1[,col_sum_index]  
  
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
output <- file.path('./新建文件夹/')

####### 堆叠柱状图
library(ggh4x)
library(reshape2)
library(ggalluvial)
library(ggplot2)

p1 = ggplot(dat2,aes(x = a, y = b)) + 
  geom_tile(aes(fill = group)) + 
  scale_fill_manual(values = c("Red","Blue")) +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        axis.title = element_blank()) + 
  scale_x_continuous(expand = c(0, 0)) +
  labs(fill = "Group")

p2 = ggplot(dat,aes(Sample, Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") + 
  labs(fill = "Cell Type",x = "",y = "Cell Proportion") + xlab(NULL)+
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 11)
  ) + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(22))

library(patchwork)
p3 <- p1 / p2 + plot_layout(heights = c(1,10),guides = "collect" ) &
  theme(legend.position = "bottom")
p3
ggsave(plot = p3,filename = paste(output,'堆叠图.png',sep = '/'),width = 12,height = 10,dpi = 600)
ggsave(plot = p3,filename = paste(output,'堆叠图.pdf',sep = '/'),width = 12,height = 10)

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
  scale_fill_manual(values=alpha(brewer.pal(8,"Set1")[1:2],0.6))+
  scale_y_continuous(expand=c(0,0),limits=c(0,0.8),breaks=seq(0,0.8,0.1),
                     labels=seq(0,0.8,0.1))+
  theme+
  labs(y='Immune Cell Relative Proportion',fill = NULL)

###添加显著性检验
Data_summary<-as.data.frame(compare_means(Proportion~group,dat,method="wilcox.test",
                                          paired=FALSE,group.by="Cell_type"))

stat.test<-Data_summary[,-c(2,9)]
stat.test$xmin=c(1:22)-0.2
stat.test$xmax=c(1:22)+0.2
stat.test$p.signif[stat.test$p.signif=="ns"]<-NA
p2<-p+geom_signif(xmin=stat.test$xmin,xmax=stat.test$xmax,
                  annotations=stat.test$p.signif,margin_top=0.00,
                  y_position=0.72,size=0.5,textsize=7.5,
                  tip_length=0)#横线两侧折线长度
p2

ggsave(plot = p2,filename = paste(output,'箱线图.png',sep = '/'),width = 12,height = 10,dpi = 600)
ggsave(plot = p2,filename = paste(output,'箱线图.pdf',sep = '/'),width = 12,height = 10)



genelist <- c("DLST","FBXO31")
#提取基因集的表达矩阵
goal_exp<-filter(expr,rownames(expr) %in% genelist)
#goal_exp <- as.data.frame(expr)
combine<-rbind(goal_exp,as.data.frame(t(res1)))
combine <- as.data.frame(t(apply(combine, 1, as.numeric)))
colnames(combine)  <- colnames(goal_exp)
comcor<-cor(t(combine))

imm <- res1
exp <- expr
#相关性计算
immuscore <- function(gene){  y <- as.numeric(exp[gene,])  
colnames <- colnames(imm)  
do.call(rbind,lapply(colnames, function(x){    
  dd  <- cor.test(as.numeric(imm[,x]),y,type="spearman")    
  data.frame(gene=gene,immune_cells=x,
             cor=dd$estimate,
             p.value=dd$p.value )}))}

gene <- genelist
data1 <- do.call(rbind,lapply(genelist,immuscore))
head(data1)
write.csv(data1, "correlation.csv", quote = F, row.names = F)


data1$pstar <- ifelse(data1$p.value < 0.05,ifelse(data1$p.value < 0.01,
                                                  ifelse(data1$p.value < 0.001,"***", "**"), "*"), "")

#相关性热图的绘制

cor_plot <- ggplot(data1, aes(immune_cells, gene)) + 
  geom_tile(aes(fill = cor), colour = "grey",size=1) + 
  scale_fill_gradient2(low = "blue",mid = "white",high = "red") + 
  geom_text(aes(label=pstar),col ="black",size = 5) + 
  theme_minimal() + 
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1,size = 20),
        axis.text.y = element_text(size = 20))+
  labs(fill =paste0("* p < 0.05","\n\n","** p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation"))+
  coord_flip()+theme_bw()
cor_plot
ggsave(plot = cor_plot,filename = paste(output,'相关性图.png',sep = '/'),width = 12,height = 10,dpi = 300)
ggsave(plot = cor_plot,filename = paste(output,'相关性图.pdf',sep = '/'),width = 12,height = 10)

######相关性

genelist <- c("DLST","FBXO31")
data <- as.data.frame(t(exp2[genelist,-1,drop = F]))
rt <- cbind(res1,data)

outtab = data.frame()
for (i in colnames(rt)[1:(ncol(rt)-2)]) {
  for (j in 1:length(genelist)) {
    x = as.numeric(rt[,genelist[j]])
    y = as.numeric(rt[,i])
    cor = cor.test(x,y,method = 'spearman')
    outVector = cbind(Gene = genelist[j],Cell = i,cor = cor$estimate,pvalue = cor$p.value)
    outtab = rbind(outtab,outVector)
  }
}

outtab$cor = as.numeric(outtab$cor)
outtab$pvalue = as.numeric(outtab$pvalue)

library(ggplot2)
library(patchwork)
library(dplyr)
library(reshape2)
library(grid)
library(cowplot)
library(gridExtra)
library(ggcharts)

color=c(
  "#FF5733", # 红色
  "#33FF57", # 绿色
  "#3357FF", # 蓝色
  "#FFFF33", # 黄色
  "#FF33FF", # 紫色
  "#33FFFF", # 青色
  "#FF8C00", # 暗橙色
  "#C70039", # 暗红色
  "#9CB570", # 黄绿色
  "#698B69", # 暗青色
  "#7B7C9C", # 石板灰
  "#FFC300", # 金色
  "#C700FF", # 亮紫色
  "#00E676", # 亮绿色
  "#536DFE", # 亮蓝色
  "#FF6F00", # 橙色
  "#8A2BE2", # 蓝紫色
  "#DA70D6", # 兰花紫
  "#FF4500", # 珊瑚色
  "#00CED1", # 浅青色
  "#FF6347", # 番茄色
  "#4682B4"  # 钢蓝色
)
unique_names <- unique(outtab$Cell)
# 创建颜色映射字典
color_dict <- setNames(color, unique_names)

a <- outtab %>% filter(Gene == 'DLST')
a$size = -log10(a$pvalue)
b <- outtab %>% filter(Gene == 'FBXO31')
b$size = -log10(a$pvalue)


p1 <- diverging_lollipop_chart(a,Cell,cor,                                         
                              lollipop_colors = c("#be8738", "#5569a5"),                                          
                              line_size = 1,                                          
                              point_size = 2,                                          
                              text_size = 13,                                          
                              text_color = c("#be8738", "#5569a5"))+  
                              labs(x = 'CellType', y = "Cor of DLST")+ 
  theme_light()+  
  theme(axis.title = element_text(size = 14),axis.text.x = element_text(size = 13),
        panel.grid = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),        
        panel.border = element_blank(),panel.background = element_blank())+ 
  ylim(-0.3, 0.3)

p1

p2 <- diverging_lollipop_chart(b,Cell,cor,                                         
                               lollipop_colors = c("#be8738", "#5569a5"),                                          
                               line_size = 1,                                          
                               point_size = 2,                                          
                               text_size = 13,                                          
                               text_color = c("#be8738", "#5569a5"))+  
  labs(x = '', y = "Cor of FBXO31")+ 
  theme_light()+  
  theme(axis.title = element_text(size = 14),axis.text.x = element_text(size = 13),
        panel.grid = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),        
        panel.border = element_blank(),panel.background = element_blank())+ 
  ylim(-0.3, 0.3)

p2

p3 <- p1 + p2
p3
ggsave(plot = p3,filename = '免疫细胞与基因相关性.pdf',width = 12,height = 8)
ggsave(plot = p3,filename = '免疫细胞与基因相关性.png',width = 12,height = 8,dpi = 600)


#################Estimate###################
library(utils)
rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
exp <- read.csv('./fitermRNA.csv',row.names = 1,check.names = F)
Group <- read.csv('./group_list.csv',row.names = 1)

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
scores = estimate(exp2,pro)
scores <- as.data.frame(scores)
scores$Group = Group$condition
scores$TumorPurity = cos(0.6049872018+0.0001467884 * scores[,3])

#boxplot
library(ggpubr)
library(ggsci)
p1 = ggplot(scores,aes(x=Group,y=TumorPurity,fill=Group))+
  geom_boxplot(position=position_dodge(0.8))+
  #geom_point(aes = (color = Group),position=position_jitterdodge(),size = 3)+
  scale_fill_nejm()+
  labs(x="",y='Tumor Purity')+
  stat_compare_means()+
  stat_compare_means(comparisons=combn(unique(scores$Group),2,simplify=FALSE))+
  theme_bw(base_size=16)
p1

p2 = ggplot(scores,aes(x=Group,y=StromalScore,fill=Group))+
  geom_boxplot(position=position_dodge(0.8))+
  #geom_point(aes = (color = Group),position=position_jitterdodge(),size = 3)+
  scale_fill_nejm()+
  labs(x="",y='Stromal Score')+
  stat_compare_means()+
  stat_compare_means(comparisons=combn(unique(scores$Group),2,simplify=FALSE))+
  theme_bw(base_size=16)
p2

p3 = ggplot(scores,aes(x=Group,y=ImmuneScore,fill=Group))+
  geom_boxplot(position=position_dodge(0.8))+
  #geom_point(aes = (color = Group),position=position_jitterdodge(),size = 3)+
  scale_fill_nejm()+
  labs(x="",y='Immune Score')+
  stat_compare_means()+
  stat_compare_means(comparisons=combn(unique(scores$Group),2,simplify=FALSE))+
  theme_bw(base_size=16)
p3

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

TIDE <- as.data.frame(read.table('./TIDE.txt',header = T,sep = '\t',check.names = F,row.names = 1))
TIDE$id <- rownames(TIDE)
TIDE <- TIDE[match(rownames(Group),rownames(TIDE)),]
TIDE$Group <- Group$condition
group = levels(factor(TIDE$Group))
comp = combn(group,2)
my_comparisons = list()
for (i in 1:ncol(comp)) {
  my_comparisons[[i]] <- comp[,i]
}


library(ggplot2)
p <- ggviolin(TIDE,x = 'Group',y = 'TIDE',fill = 'Group',
         xlab = '',ylab = 'TIDE_Scores',
         palette = c('Firebrick2','DodgerBlue1'),
         legend.title = 'Group',
         add = 'boxplot',add.params = list(fill = 'white')) +
  stat_compare_means(comparisons = my_comparisons)
ggsave(plot = p,filename = './新建文件夹/TIDE评分.pdf',height = 8,width = 10)
ggsave(plot = p,filename = './新建文件夹/TIDE评分.png',height = 8,width = 10,dpi = 300)

res <- as.data.frame(read.csv('./TIDE.csv',row.names = 1))
res <- res[match(rownames(Group),rownames(res)),]
res$Group <- Group$condition
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
########免疫反应与亚型
library(dplyr)
dat=count(res,Group,Responder)
dat=dat%>%group_by(Group)%>%
  summarise(Responder=Responder,n=n/sum(n))
dat$Responder=factor(dat$Responder,levels=c("False","True"))
dat

library(ggplot2)
p2=ggplot(data=dat)+
  geom_bar(aes(x=Group,y=n,
               fill=Responder),
           stat="identity")+
  scale_fill_manual(values=c("#e04030","#6cb8d2"))+
  geom_label(aes(x=Group,y=n,
                 label=scales::percent(n),
                 fill=Responder),
             color="white",
             size=6,label.size=0,
             show.legend=FALSE,
             position=position_fill(vjust=0.5))+
  ylab("Percentage")+
  theme_bw()+
  guides(fill=guide_legend(title="Responder")) +#仅保留一个图例
  theme(axis.text.x = element_text(size = 12,face = 'bold'),
        axis.title.x = element_text(size = 14,face = 'bold'),
        axis.text.y = element_text(size = 10,face = 'bold'),
        axis.title.y = element_text(size = 14,face = 'bold'))
  
  
  
  
library(patchwork)
p3 <- p1+p2+plot_layout(widths=c(3,2),guides="collect")
p3
ggsave(plot = p3,filename = './新建文件夹/total.png',width=14,height=8,dpi = 600)
ggsave(plot = p3,filename = './新建文件夹/total.pdf',width=14,height=8)

###############ICI免疫治疗反应评分##########################
library(reshape2)
library(ggpubr)
tcia = read.table('./TCIA-ClinicalData.tsv',header = T,sep = '\t',check.names = F,row.names = 1)
tcia = tcia[,c('ips_ctla4_neg_pd1_neg','ips_ctla4_neg_pd1_pos','ips_ctla4_pos_pd1_neg','ips_ctla4_pos_pd1_pos')]

#rownames(Group) = gsub("(.*?)\\-(.*?)-(.*?)\\-(.*?)","\\1\\-\\2\\-\\3",rownames(Group))
Group <- read.csv('./group_list.csv',row.names = 1)
rownames(Group) = gsub('-01A','',rownames(Group))
overlap <- intersect(rownames(Group),rownames(tcia))
tcia <- tcia[overlap,]
Group$id = rownames(Group)
Group <- Group[overlap,]
tcia$Group = Group$condition
tcia$Group = factor(tcia$Group,levels = c('Low','High'))
tcia = na.omit(tcia)



for (i in 1:4) {
  a <- tcia
  colnames(a)[i] <- 'nn'
  col <- colnames(tcia)[i]
  p <- ggboxplot(a,           
                 x="Group",           
                 y="nn",          
                 xlab = "",         
                 ylab = col,          
                 color = "Group",          
                 order=c("High","Low"),         
                 bxp.errorbar=T,          
                 bxp.errorbar.width = 0.1,          
                 palette = "lancet",           
                 add = "jitter",           
                 lwd=0.7,           
                 fatten=2,)+  
    scale_color_manual(values = c("#e60033", "#38a1db"))+  
    geom_signif(comparisons = list(c("High","Low")),              
                map_signif_level = T,               
                test = "wilcox.test",               
                y_position = c(11),              
                tip_length = c(c(0.05,0.05)),              
                size=0.8,color="black")+  
    theme_bw()
  theme(legend.position = "top",        
        axis.title.x =element_text(size=16,face = 'bold.italic'),        
        axis.title.y=element_text(size=16,face = 'bold.italic'),        
        axis.text.x=element_text(size=16,face = 'bold'),        
        axis.text.y=element_text(size=16,face = 'bold'),
        plot.margin=unit(c(2,2,2,2),'cm'))
  ggsave(plot = p,filename = paste(col,'.pdf'),height = 9,width = 7)
  ggsave(plot = p,filename = paste(col,'.png'),height = 7,width = 7,dpi = 600)
}






