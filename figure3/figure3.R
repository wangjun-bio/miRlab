# 将DCC、CICR2数据预处理
if (!require("ggvenn")) {
  install.packages("ggvenn")
  library(ggvenn)
}
# 加载 ggplot2 包
if (!require("ggplot2")) {
  BiocManager::install("ggplot2")
  library("ggplot2")
}
# 加载 ggpubr 包
if (!require("ggpubr")) {
  BiocManager::install("ggpubr")
  library("ggpubr")
}
# 加载 ggpmisc 包
if (!require("ggpmisc")) {
  BiocManager::install("ggpmisc")
  library("ggpmisc")
}
if (!require("VennDiagram")) {
  BiocManager::install("VennDiagram")
  library("VennDiagram")
}
if (!require("circlize")) {
  BiocManager::install("circlize")
  library(circlize)
}




# a
ciri2_count_distribution <- read.csv("ciri2_count_distribution.csv", row.names = 1)
dcc_count_distribution <- read.csv("dcc_count_distribution.csv", row.names = 1)
ciri2_filtered_count_distribution <- read.csv("ciri2_filtered_count_distribution.csv", row.names = 1)
dcc_filtered_count_distribution <- read.csv("dcc_filtered_count_distribution.csv", row.names = 1)
# b
NC <- read.csv("dcc_vs_ciri2_NC_venn_plot.csv", header = TRUE, row.names = 1)
Hyp <- read.csv("dcc_vs_ciri2_Hyp_venn_plot.csv", header = TRUE, row.names = 1)
Su <- read.csv("dcc_vs_ciri2_Su_venn_plot.csv", header = TRUE, row.names = 1)
MCT <- read.csv("dcc_vs_ciri2_MCT_venn_plot.csv", header = TRUE, row.names = 1)
dcc_vs_ciri2_NC_correlation_plot <- read.csv(file = "dcc_vs_ciri2_NC_correlation_plot.csv", row.names = 1)
dcc_vs_ciri2_Hyp_correlation_plot <- read.csv(file = "dcc_vs_ciri2_Hyp_correlation_plot.csv", row.names = 1)
dcc_vs_ciri2_Su_correlation_plot <- read.csv(file = "dcc_vs_ciri2_Su_correlation_plot.csv", row.names = 1)
dcc_vs_ciri2_MCT_correlation_plot <- read.csv(file = "dcc_vs_ciri2_MCT_correlation_plot.csv", row.names = 1)
# c
circRNA_type_barplot <- read.csv(file = "circRNA_type_barplot.csv", row.names = 1)
circRNA_hostGene_number <- read.csv("circRNA_hostGene_number.csv", row.names = 1)
# d
NC_bed<- read.csv("NC_bed.csv", row.names = 1)
HA_bed<- read.csv("HA_bed.csv", row.names = 1)
SA_bed<- read.csv("SA_bed.csv", row.names = 1)
MA_bed<- read.csv("MA_bed.csv", row.names = 1)
# e
five_set_venn_plot <- read.csv(file = "5_set_venn_plot.csv", row.names = 1)

# 堆叠柱状图
barplot(as.matrix(ciri2_count_distribution),                                        
        col = c("#2878B5", "#9AC9DB", "#F8AC8C", "#C82423"), # 柱子颜色
        border = "NA",
        ylab = "Number of CircRNAs",
        xlab = "Samples", # 柱子边框为无色
        main = "ciri2",
        cex.axis = 1,
        cex.names = 0.8,
        axisnames = T,
        xlim = c(0,15),
        ylim = c(0,3500),
        xpd = T,
        width = 0.8,
        names.arg = colnames(ciri2_count_distribution),  # 柱子名称
        #xlab = 'class',  # X轴标题
        # ylab = 'value',  # Y轴标题
        #horiz = TRUE # 水平柱状图
)    
#添加图注
legend(
  "topright", # 指定图例的位置为右上角
  legend = c("1 reads", "2 to 5 readds", "6 to 10 readds", ">10 readds"),
  fill = c("#2878B5", "#9AC9DB", "#F8AC8C", "#C82423"),
  y.intersp = 0.9
)

# 堆叠柱状图
barplot(as.matrix(dcc_count_distribution),                                        
        col = c("#2878B5", "#9AC9DB", "#F8AC8C", "#C82423"), # 柱子颜色
        border = "NA",
        ylab = "Number of CircRNAs",
        xlab = "Samples", # 柱子边框为无色
        main = "dcc",
        cex.axis = 1,
        cex.names = 0.8,
        axisnames = T,
        xlim = c(0,15),
        ylim = c(0,55000),
        xpd = T,
        width = 0.8,
        names.arg = colnames(dcc_count_distribution),  # 柱子名称
        #xlab = 'class',  # X轴标题
        # ylab = 'value',  # Y轴标题
        #horiz = TRUE # 水平柱状图
)    
#添加图注
legend(
  "topright", # 指定图例的位置为右上角
  legend = c("1 reads", "2 to 5 readds", "6 to 10 readds", ">10 readds"),
  fill = c("#2878B5", "#9AC9DB", "#F8AC8C", "#C82423"),
  y.intersp = 0.9
)



# 

# 堆叠柱状图
barplot(as.matrix(ciri2_filtered_count_distribution),                                        
        col = c("#2878B5", "#9AC9DB", "#F8AC8C", "#C82423"), # 柱子颜色
        border = "NA",
        ylab = "Number of CircRNAs",
        xlab = "Samples", # 柱子边框为无色
        main = "ciri2_filtered",
        cex.axis = 1,
        cex.names = 0.8,
        axisnames = T,
        xlim = c(0,15),
        ylim = c(0,2000),
        xpd = T,
        width = 0.8,
        names.arg = colnames(ciri2_filtered_count_distribution),  # 柱子名称
        #xlab = 'class',  # X轴标题
        # ylab = 'value',  # Y轴标题
        #horiz = TRUE # 水平柱状图
)    
#添加图注
legend(
  "topright", # 指定图例的位置为右上角
  legend = c("1 reads", "2 to 5 readds", "6 to 10 readds", ">10 readds"),
  fill = c("#2878B5", "#9AC9DB", "#F8AC8C", "#C82423"),
  y.intersp = 0.9
)

barplot(as.matrix(dcc_filtered_count_distribution),                                        
        col = c("#2878B5", "#9AC9DB", "#F8AC8C", "#C82423"), # 柱子颜色
        border = "NA",
        ylab = "Number of CircRNAs",
        xlab = "Samples", # 柱子边框为无色
        main = "dcc_filtered",
        cex.axis = 1,
        cex.names = 0.8,
        axisnames = T,
        xlim = c(0,15),
        ylim = c(0,3000),
        xpd = T,
        width = 0.8,
        names.arg = colnames(dcc_filtered_count_distribution),  # 柱子名称
        #xlab = 'class',  # X轴标题
        # ylab = 'value',  # Y轴标题
        #horiz = TRUE # 水平柱状图
)    
#添加图注
legend(
  "topright", # 指定图例的位置为右上角
  legend = c("1 reads", "2 to 5 readds", "6 to 10 readds", ">10 readds"),
  fill = c("#2878B5", "#9AC9DB", "#F8AC8C", "#C82423"),
  y.intersp = 0.9
)


# NC

#输出分组取交集venn
x <- list(ciri2 = na.omit(NC$ciri2_ids),
          dcc = na.omit(NC$dcc_ids)
)
ggvenn(
  x, 
  fill_color = c("#CDE9EF", "#FAF8CB"),
  stroke_size = 0.5, set_name_size = 5
)


# Hyp

#输出分组取交集venn
x <- list(ciri2 = na.omit(Hyp$ciri2_ids),
          dcc = na.omit(Hyp$dcc_ids)
)
ggvenn(
  x, 
  fill_color = c("#CDE9EF", "#FAF8CB"),
  stroke_size = 0.5, set_name_size = 5
)


# Su

#输出分组取交集venn
x <- list(ciri2 = Su$ciri2_ids,
          dcc = Su$dcc_ids
)
ggvenn(
  x, 
  fill_color = c("#CDE9EF", "#FAF8CB"),
  stroke_size = 0.5, set_name_size = 5
)


# MCT

#输出分组取交集venn
x <- list(ciri2 = MCT$ciri2_ids,
          dcc = MCT$dcc_ids
)
ggvenn(
  x, 
  fill_color = c("#CDE9EF", "#FAF8CB"),
  stroke_size = 0.5, set_name_size = 5
)



#分组筛选后circRNA相关性图输出



xlimit <- 20
ylimit <- 15

ggplot(dcc_vs_ciri2_NC_correlation_plot,aes(x = ciri2_mean, y = dcc_mean)) +
  geom_point(color = "black") +
  geom_smooth( method = "lm", se = TRUE, color = "blue") +
  stat_cor(method = "pearson",label.x = 0, label.y = ylimit, size = 4.7) +
  coord_cartesian(xlim = c(0, xlimit), ylim = c(0, ylimit))+#限制显示的视野，不会舍弃视野外的数值
  scale_x_continuous(limits = c(0,max(dcc_vs_ciri2_NC_correlation_plot$ciri2_mean))) +
  scale_y_continuous(limits = c(0, max(dcc_vs_ciri2_NC_correlation_plot$dcc_mean))) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  geom_text(  
    aes(x = xlimit/2, y = ylimit-1, label = "NC"),  
    hjust = -0.5, vjust = -2, size = 4.7  
  )+
  xlab("mean_count of circ2") +
  ylab("mean_count of dcc")  






ggplot(dcc_vs_ciri2_Hyp_correlation_plot,aes(x = ciri2_mean, y = dcc_mean)) +
  geom_point(color = "black") +
  geom_smooth( method = "lm", se = TRUE, color = "blue") +
  stat_cor(method = "pearson",label.x = 0, label.y = ylimit, size = 4.7) +
  coord_cartesian(xlim = c(0, xlimit), ylim = c(0, ylimit))+#限制显示的视野，不会舍弃视野外的数值
  scale_x_continuous(limits = c(0,max(dcc_vs_ciri2_Hyp_correlation_plot$ciri2_mean))) +
  scale_y_continuous(limits = c(0, max(dcc_vs_ciri2_Hyp_correlation_plot$dcc_mean))) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  geom_text(  
    aes(x = xlimit/2, y = ylimit-1, label = "Hyp"),  
    hjust = -0.5, vjust = -2, size = 4.7  
  )+
  xlab("mean_count of circ2") +
  ylab("mean_count of dcc")  






xlimit <- 20
ylimit <- 15

ggplot(dcc_vs_ciri2_Su_correlation_plot,aes(x = ciri2_mean, y = dcc_mean)) +
  geom_point(color = "black") +
  geom_smooth( method = "lm", se = TRUE, color = "blue") +
  stat_cor(method = "pearson",label.x = 0, label.y = ylimit, size = 4.7) +
  coord_cartesian(xlim = c(0, xlimit), ylim = c(0, ylimit))+#限制显示的视野，不会舍弃视野外的数值
  scale_x_continuous(limits = c(0,max(dcc_vs_ciri2_Su_correlation_plot$ciri2_mean))) +
  scale_y_continuous(limits = c(0, max(dcc_vs_ciri2_Su_correlation_plot$dcc_mean))) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  geom_text(  
    aes(x = xlimit/2, y = ylimit-1, label = "Su"),  
    hjust = -0.5, vjust = -2, size = 4.7  
  )+
  xlab("mean_count of circ2") +
  ylab("mean_count of dcc")   



xlimit <- 20
ylimit <- 15

ggplot(dcc_vs_ciri2_MCT_correlation_plot,aes(x = ciri2_mean, y = dcc_mean)) +
  geom_point(color = "black") +
  geom_smooth( method = "lm", se = TRUE, color = "blue") +
  stat_cor(method = "pearson",label.x = 0, label.y = ylimit, size = 4.7) +
  coord_cartesian(xlim = c(0, xlimit), ylim = c(0, ylimit))+#限制显示的视野，不会舍弃视野外的数值
  scale_x_continuous(limits = c(0,max(dcc_vs_ciri2_MCT_correlation_plot$ciri2_mean))) +
  scale_y_continuous(limits = c(0, max(dcc_vs_ciri2_MCT_correlation_plot$dcc_mean))) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  geom_text(  
    aes(x = xlimit/2, y = ylimit-1, label = "MCT"),  
    hjust = -0.5, vjust = -2, size = 4.7  
  )+
  xlab("mean_count of circ2") +
  ylab("mean_count of dcc")  


# circRNA种类统计

data_ouput <- as.matrix(circRNA_type_barplot)
# 堆叠柱状图
barplot(data_ouput,
        col = c("#9AC9DB", "#F8AC8C", "#C82423"), # 柱子颜色
        border = "NA",
        xlab = "Number of CircRNAs",
        ylab = "Samples", # 柱子边框为无色
        main = "circRNA种类统计",
        cex.axis = 1,
        cex.names = 0.8,
        axisnames = TRUE,
        ylim = c(0, 5),
        xlim = c(0, 2500),
        xpd = TRUE,
        width = 0.8,
        names.arg = colnames(data_ouput),  # 柱子名称
        horiz = TRUE
)   
#添加图注
legend(
  "topright", # 指定图例的位置为右上角
  legend = rownames(result),
  fill = c("#9AC9DB", "#F8AC8C", "#C82423"),
  y.intersp = 0.9
)

# circRNA与母源基因数量统计



# 堆叠柱状图
barplot(as.matrix(circRNA_hostGene_number),
        col = c("#9AC9DB", "#F8AC8C", "#C82423","red"), # 柱子颜色
        border = "NA",
        xlab = "Number of CircRNAs",
        ylab = "Samples", # 柱子边框为无色
        main = "circRNA_hostGene_number",
        cex.axis = 1,
        cex.names = 0.8,
        axisnames = TRUE,
        ylim = c(0, 5),
        xlim = c(0, 1000),
        xpd = TRUE,
        width = 0.8,
        names.arg = colnames(circRNA_hostGene_number),  # 柱子名称
        horiz = TRUE
)   
#添加图注
legend(
  "topright", # 指定图例的位置为右上角
  legend = rownames(circRNA_hostGene_number),
  fill = c("#9AC9DB", "#F8AC8C", "#C82423","red"),
  y.intersp = 0.9
)


# 分布圈图



# 初始化 Circos 图，设置染色体 ideogram
circos.initializeWithIdeogram(species = "rn6")

# 绘制染色体密度图
circos.genomicDensity(MA_bed, col = "#407434",track.height = 0.15)
circos.genomicDensity(SA_bed, col = "#019ED5",track.height = 0.15)
circos.genomicDensity(HA_bed, col = "#D9742B",track.height = 0.15)
circos.genomicDensity(NC_bed, col = "#BA2835",track.height = 0.2)





# 输出与外部数据库cicrAltas的交集


x <- list(NC = na.omit(five_set_venn_plot$NC),
          Hyp = na.omit(five_set_venn_plot$Hyp),
          Su = na.omit(five_set_venn_plot$Su),
          MCT = na.omit(five_set_venn_plot$MCT),
          circAltas = na.omit(five_set_venn_plot$circAltas)
)

data <- venn.diagram(x=x,
                     
                     scaled = F, # 根据比例显示大小
                     
                     alpha= 0.5, #透明度
                     
                     lwd=1,lty=1,col=c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC"), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
                     
                     label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
                     
                     cex = 2, # 数字大小
                     
                     fontface = "bold",  # 字体粗细；加粗bold
                     
                     fill=c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC"), # 填充色 配色https://www.58pic.com/
                     
                     # category.names = c("Set1", "Set2","Set3","Set4","Set5") , #标签名
                     
                     cat.dist = c(0.2, 0.2, 0.2, 0.2, 0.2), # 标签距离圆圈的远近
                     
                     cat.pos = c(0, -10, 240, 120, 20), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
                     
                     cat.cex = 2, #标签字体大小
                     
                     cat.fontface = "bold",  # 标签字体加粗
                     
                     cat.col=c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC"),   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
                     
                     cat.default.pos = "outer",  # 标签位置, outer内;text 外
                     
                     output=FALSE,
                     
                     filename= NULL,# 文件保存
                     
                     # imagetype="png",  # 类型（tiff png svg）
                     
                     resolution = 400,  # 分辨率
                     
                     compression = "lzw"# 压缩算法
                     
)

grid.draw(data)










