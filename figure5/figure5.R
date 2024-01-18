if (!require("ggplot2")) {
  BiocManager::install("ggplot2")
  library("ggplot2")
}
GO_KEGG <- read.csv(file = "TargetGenesGO_KEGG.csv", row.names = 1)

p1 <- ggplot(data = GO_KEGG, mapping = aes(x = GeneRatio ,y = Description)) + # 锁定数据框，X轴Y轴数据
  geom_point(aes(size=Count,color=p.adjust))+ #绘制气泡图
  # coord_flip() + # 将X轴与Y轴调换
  scale_color_gradient(low="#55668E",high ="#862657")+ #颜色范围
  theme(legend.title = element_text(size = 15, face = 2))+ #图例标题大小
  theme(legend.key.size=unit(1,'cm'))+ #图例标题大小
  theme(legend.text = element_text( size = 15,face = 'bold'))+ #图例图形里文字大小
  labs(x="GeneRatio",y = "",title = "target genes GO & KEGG") +# X轴、 Y轴 、图的biao
  theme_bw()+
  theme(plot.background = element_rect(fill = "White"))+ #图片背景颜色设置
  theme(panel.grid.major=element_line(colour="grey")) +
  scale_fill_gradient(low = "pink", high = "red")+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 12, color = "black", face = "bold"))+
  # theme(axis.text.y = element_text(size = 12, color =c("red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue"),face = "bold"))+
  theme(plot.title = element_text(size = 15,face = 4, hjust =0.5))+
  facet_grid(group ~ ., scales = "free_y") 

p1