if (!require("ggplot2")) {
  BiocManager::install("ggplot2")
  library("ggplot2")
}
if(!require("ggforce")){
  BiocManager::install("ggforce")
  library("ggforce")
}


HA_VS_NC <- read.csv("HA_vs_NC_output_matrix.csv")
SA_vs_NC <- read.csv("SA_vs_NC_output_matrix.csv")
MA_VS_NC <- read.csv("MA_vs_NC_output_matrix.csv")

DE_GO_result <- read.csv("DE_GO_GSEA.csv", header = TRUE, check.names = FALSE, row.names = 1)
DE_KEGG_result <- read.csv("DE_KEGG_GSEA.csv", header = TRUE, check.names = FALSE, row.names = 1)
intersect_GO_result <- read.csv("DE_GO_GSEA.csv", header = TRUE, check.names = FALSE, row.names = 1)
intersect_KEGG_result <- read.csv("DE_KEGG_GSEA.csv", header = TRUE, check.names = FALSE, row.names = 1)

pic <- read.csv("DE_PCA.csv", row.names = 1)

result_list <-  list(HA_vs_NC = HA_VS_NC,
                   SA_vs_NC = SA_vs_NC,
                   MA_vs_NC = MA_VS_NC)

# 创建结果列表
contrast_names <- names(result_list)
# 创建pvalue火山图并保存为PDF文件
for (i in 1:length(result_list)) {
  x_scale <- max(na.omit(abs(result_list[[i]]$log2FoldChange))) + 1
  y_scale <- max(na.omit(-log10(result_list[[i]]$pvalue))) + 1.5
  cut_off_pvalue <-  0.05
  cut_off_logFC <- 1
  volcano_plot <- ggplot(result_list[[i]], aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = ifelse(pvalue < cut_off_pvalue & abs(log2FoldChange) > cut_off_logFC,
                                  ifelse(log2FoldChange > 0, "Upregulated", "Downregulated"),
                                  "Non-significant")),
               alpha = 0.4, size = 2) +
    scale_color_manual(values = c("Upregulated" = "#ff4757", "Downregulated" = "#546de5", "Non-significant" = "#d2dae2"),
                       labels = c("Upregulated" = "Upregulated Genes",
                                  "Downregulated" = "Downregulated Genes",
                                  "Non-significant" = "Non-significant Genes")) +
    theme_bw()+    #去除背景色
    theme(panel.grid=element_blank())+    #去除网格线
    geom_vline(xintercept=c(-1,1),lty=2,col="black",lwd=0.4)+    #添加横线|FoldChange|>2
    geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.4)  + #添加竖线padj<0.05
    labs(title = paste("Volcano Plot -", contrast_names[i]),
         x = "log2(Fold Change)", y = "-log10(Adjusted p-value)",
         color = "Gene Type") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text = element_text(colour = 'black'))+    #图例位置
    guides(color = guide_legend(title = "Gene Type"))+
    coord_cartesian(xlim = c(-17, 17), ylim = c(0, 6))+
    scale_x_continuous(limits = c(-x_scale, x_scale))  # 设置x轴刻度范围为-2到2
  print(volcano_plot)
  
  
  # 保存为PDF文件
  pdf_file_name <- paste("Volcano_Plot_pvalue_", contrast_names[i], ".pdf", sep = "")
  ggsave(file = pdf_file_name, plot = volcano_plot, width = 9.03, height = 7.16)
}



GO_plot <- ggplot(data = DE_GO_result, mapping = aes(x = group ,y = Description)) + # 锁定数据框，X轴Y轴数据
  geom_point(aes(size=-log10(pvalue),color=NES))+ #绘制气泡图
  # coord_flip()+ # 将X轴与Y轴调换
  scale_color_gradient(low="#55668E",high ="#862657")+ #颜色范围
  theme(legend.title = element_text(size = 15, face = 2))+ #图例标题大小
  theme(legend.key.size=unit(1,'cm'))+ #图例标题大小
  theme(legend.text = element_text( size = 15,face = 'bold'))+ #图例图形里文字大小
  labs(x="PAH_RAT_MODEL",y = "Term",title = "GO") +# X轴、 Y轴 、图的biao
  theme_bw()+
  theme(plot.background = element_rect(fill = "White"))+ #图片背景颜色设置
  theme(panel.grid.major=element_line(colour="grey")) +
  scale_fill_gradient(low = "pink", high = "red")+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 12, color = "black", face = "bold"))+
  # theme(axis.text.y = element_text(size = 12, color =c("red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue"),face = "bold"))+
  theme(plot.title = element_text(size = 15,face = 4, hjust =0.5))
GO_plot


KEGG_plot <- ggplot(data = DE_KEGG_result, mapping = aes(x = group ,y = Description)) + # 锁定数据框，X轴Y轴数据
  geom_point(aes(size=-log10(pvalue),color=NES))+ #绘制气泡图
  # coord_flip()+ # 将X轴与Y轴调换
  scale_color_gradient(low="#55668E",high ="#862657")+ #颜色范围
  theme(legend.title = element_text(size = 15, face = 2))+ #图例标题大小
  theme(legend.key.size=unit(1,'cm'))+ #图例标题大小
  theme(legend.text = element_text( size = 15,face = 'bold'))+ #图例图形里文字大小
  labs(x="PAH_RAT_MODEL",y = "Term",title = "KEGG") +# X轴、 Y轴 、图的biao
  theme_bw()+
  theme(plot.background = element_rect(fill = "White"))+ #图片背景颜色设置
  theme(panel.grid.major=element_line(colour="grey")) +
  scale_fill_gradient(low = "pink", high = "red")+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 12, color = "black", face = "bold"))+
  # theme(axis.text.y = element_text(size = 12, color =c("red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue"),face = "bold"))+
  theme(plot.title = element_text(size = 15,face = 4, hjust =0.5))

KEGG_plot

GO_plot <- ggplot(data = intersect_GO_result, mapping = aes(x = group ,y = Description)) + # 锁定数据框，X轴Y轴数据
  geom_point(aes(size=-log10(pvalue),color=NES))+ #绘制气泡图
  # coord_flip()+ # 将X轴与Y轴调换
  scale_color_gradient(low="#55668E",high ="#862657")+ #颜色范围
  theme(legend.title = element_text(size = 15, face = 2))+ #图例标题大小
  theme(legend.key.size=unit(1,'cm'))+ #图例标题大小
  theme(legend.text = element_text( size = 15,face = 'bold'))+ #图例图形里文字大小
  labs(x="PAH_RAT_MODEL",y = "Term",title = "GO") +# X轴、 Y轴 、图的biao
  theme_bw()+
  theme(plot.background = element_rect(fill = "White"))+ #图片背景颜色设置
  theme(panel.grid.major=element_line(colour="grey")) +
  scale_fill_gradient(low = "pink", high = "red")+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 12, color = "black", face = "bold"))+
  # theme(axis.text.y = element_text(size = 12, color =c("red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue"),face = "bold"))+
  theme(plot.title = element_text(size = 15,face = 4, hjust =0.5))
GO_plot


KEGG_plot <- ggplot(data = intersect_KEGG_result, mapping = aes(x = group ,y = Description)) + # 锁定数据框，X轴Y轴数据
  geom_point(aes(size=-log10(pvalue),color=NES))+ #绘制气泡图
  # coord_flip()+ # 将X轴与Y轴调换
  scale_color_gradient(low="#55668E",high ="#862657")+ #颜色范围
  theme(legend.title = element_text(size = 15, face = 2))+ #图例标题大小
  theme(legend.key.size=unit(1,'cm'))+ #图例标题大小
  theme(legend.text = element_text( size = 15,face = 'bold'))+ #图例图形里文字大小
  labs(x="PAH_RAT_MODEL",y = "Term",title = "KEGG") +# X轴、 Y轴 、图的biao
  theme_bw()+
  theme(plot.background = element_rect(fill = "White"))+ #图片背景颜色设置
  theme(panel.grid.major=element_line(colour="grey")) +
  scale_fill_gradient(low = "pink", high = "red")+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 12, color = "black", face = "bold"))+
  # theme(axis.text.y = element_text(size = 12, color =c("red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","orange","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue"),face = "bold"))+
  theme(plot.title = element_text(size = 15,face = 4, hjust =0.5))

KEGG_plot




pic$BATCH <- factor(pic$BATCH, levels = c("NC", "HYPOXIA", "HYPOXIA+sugen5416", "MCT"))

ggplot(data = pic, aes(x = PC1, y = PC2, color = BATCH)) +
  geom_mark_ellipse(aes(fill = BATCH), expand = unit(0.5, "mm"), color = "transparent") +
  geom_point(alpha = 0.5, size = 3) +
  labs(x = 'PC1: 0.3381', y = 'PC2: 0.1326', color = "BATCH", title = 'PCA-PLOT') +
  guides(fill = FALSE) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("NC" = "#BA2835", "HYPOXIA" = "#D9742B", "HYPOXIA+sugen5416" = "#019ED5", "MCT" = "#407434")) +
  scale_fill_manual(values = c("NC" = "#BA2835", "HYPOXIA" = "#D9742B", "HYPOXIA+sugen5416" = "#019ED5", "MCT" = "#407434"))





