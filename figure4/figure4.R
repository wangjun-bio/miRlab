if (!require("ggplot2")) {
  BiocManager::install("ggplot2")
  library("ggplot2")
}
if (!require("ggVennDiagram")) {
  devtools::install_github("gaospecial/ggVennDiagram")
  library("ggVennDiagram")
}
if (!require("clusterProfiler")) {
  BiocManager::install("clusterProfiler")
  library("clusterProfiler")
}
if (!require("circlize")) {
  BiocManager::install("circlize")
  library(circlize)
}
if (!require("pheatmap")) {
  BiocManager::install("pheatmap")
  library("pheatmap")
}



# 火山图
HA_VS_NC <- read.csv("HA_vs_NC_output_matrix.csv")
SA_vs_NC <- read.csv("SA_vs_NC_output_matrix.csv")
MA_VS_NC <- read.csv("MA_vs_NC_output_matrix.csv")
DE_venn <- read.csv("DE_intersect_venn.csv", row.names = 1)
DE_HA_bed<- read.csv("DE_HA_bed.csv", row.names = 1)
DE_SA_bed<- read.csv("DE_SA_bed.csv", row.names = 1)
DE_MA_bed<- read.csv("DE_MA_bed.csv", row.names = 1)
all_heatmap <- read.csv(file = "all_heatmap.csv", row.names = 1)
all_heatmap_annotation_col <- read.csv(file = "all_heatmap_annotation_col.csv", row.names = 1)
ann_colors <- read.csv(file = "all_heatmap_ann_colors.csv", row.names = 1)
anno_summary <- read.csv(file = "anno_summary.csv", row.names = 1)
host_and_circ_FC <- read.csv(file = "host_circ_FC_correlation.csv")
filtered_data <- read.csv("ciri2_filtered_annotion_matrix.csv", header = TRUE, check.names = FALSE)

result_list <- list(HA_VS_NC = HA_VS_NC, SA_vs_NC = SA_vs_NC, MA_VS_NC = MA_VS_NC)
# 创建结果列表
contrast_names <- names(result_list)
cut_off_pvalue <-  0.05
cut_off_logFC <- 1
for (i in 1:length(result_list)) {
  x_scale <- max(na.omit(abs(result_list[[i]]$log2FoldChange))) + 1
  y_scale <- max(na.omit(-log10(result_list[[i]]$pvalue))) + 1.5
  ggplot(result_list[[i]], aes(x = log2FoldChange, y = -log10(pvalue))) +
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
}

# 差异基因venn

vennlist <- list(Hyp_vs_Nor = DE_venn$HA_vs_NC,
                 HySu_vs_Nor = DE_venn$SA_vs_NC,
                 MCT_vs_Nor = DE_venn$MA_vs_NC)

venn <- Venn(vennlist)
data <- process_data(venn)
ggplot() +
  geom_sf(aes(fill = count), 
          data = venn_region(data)) +
  geom_sf(color="grey", 
          size = 1, 
          data = venn_setedge(data), 
          show.legend = FALSE) +
  scale_fill_gradient(low="blue",high = "yellow",name = "gene")+
  geom_sf_text(aes(label = name), 
               data = venn_setlabel(data),
               size = 3) +
  geom_sf_label(aes(label = count), 
                data = venn_region(data),
                size = 4) +
  theme_void()



anno_summary$group <- factor(group, levels = c( "unannotation","has_hostGene","in_circAtlas"))
ggplot() + 
  geom_bar(data =anno_summary, 
           aes(x = group, y = nums),
           stat = "identity",
           fill = c("purple","blue","gray"))+ 
  coord_flip()






# 每种模型的差异基因分布



# 初始化 Circos 图，设置染色体 ideogram
circos.initializeWithIdeogram(species = "rn6")


# 绘制染色体密度图
circos.genomicDensity(DE_MA_bed, col = "#407434",track.height = 0.2)
circos.genomicDensity(DE_SA_bed, col = "#019ED5",track.height = 0.2)
circos.genomicDensity(DE_HA_bed, col = "#D9742B",track.height = 0.3)



# heatmap


all_heatmap_ann_colors <- list(treatment = unlist(ann_colors$treatment))
names(all_heatmap_ann_colors[[1]]) <- rownames(ann_colors)

# 使用 rpm 数据制作热图，只包括 "condition" 列的注释
heatmap_output <-pheatmap(all_heatmap,
                          cluster_rows = T,
                          show_rownames = FALSE,
                          show_colnames = F,
                          cluster_cols = F,
                          scale = "row",
                          treeheight_col = 5,
                          fontsize =11,
                          fontsize_col = 10,
                          annotation_col = all_heatmap_annotation_col,
                          annotation_colors = all_heatmap_ann_colors,
                          main = "heatmap_pvalue_ciri2",
                          filename = "heatmap_pvalue_ciri2.pdf",width = 15, height = 6)


# circRNA_hostgene_FC


# HA_vs_NC
xlimit <- 6
ylimit <- 2
ggplot(host_and_circ_FC,aes(x = HA_vs_NC_log2FoldChange, y = mRNA_HA_Vs_NC_log2FoldChange)) +
  geom_point(color = "black") +
  geom_smooth( method = "lm", se = TRUE, color = "blue") +
  stat_cor(method = "pearson",label.x = -xlimit, label.y = ylimit, size = 4.7) +
  coord_cartesian(xlim = c(-xlimit, xlimit), ylim = c(-ylimit, ylimit))+#限制显示的视野，不会舍弃视野外的数值
  # scale_x_continuous(limits = c(-10, 10)) +
  # scale_y_continuous(limits = c(-10, 10)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  geom_text(  
    aes(x = 0, y = ylimit-0.3, label = "Hyp"),  
    hjust = -0.5, vjust = -2, size = 4.7  
  )+
  xlab("log2foldchange of circRNA") +
  ylab("log2foldchange of host_gene")  

# SA_vs_NC
ggplot(host_and_circ_FC,aes(x = SA_vs_NC_log2FoldChange, y = mRNA_SA_Vs_NC_log2FoldChange)) +
  geom_point(color = "black") +
  geom_smooth( method = "lm", se = TRUE, color = "blue") +
  stat_cor(method = "pearson",label.x = -xlimit, label.y = ylimit, size = 4.7) +
  coord_cartesian(xlim = c(-xlimit, xlimit), ylim = c(-ylimit, ylimit))+#限制显示的视野，不会舍弃视野外的数值
  # scale_x_continuous(limits = c(-10, 10)) +
  # scale_y_continuous(limits = c(-10, 10)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  geom_text(  
    aes(x = 0, y = ylimit-0.3, label = "Su"),  
    hjust = -0.5, vjust = -2, size = 4.7  
  )+
  xlab("log2foldchange of circRNA") +
  ylab("log2foldchange of host_gene")  

# MA_vs_NC
xlimit <- 7
ylimit <- 3
ggplot(host_and_circ_FC,aes(x = MA_vs_NC_log2FoldChange, y = mRNA_MA_Vs_NC_log2FoldChange)) +
  geom_point(color = "black") +
  geom_smooth( method = "lm", se = TRUE, color = "blue") +
  stat_cor(method = "pearson",label.x = -xlimit, label.y = ylimit, size = 4.7) +
  coord_cartesian(xlim = c(-xlimit, xlimit), ylim = c(-ylimit, ylimit))+#限制显示的视野，不会舍弃视野外的数值
  # scale_x_continuous(limits = c(-10, 10)) +
  # scale_y_continuous(limits = c(-10, 10)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  geom_text(  
    aes(x = 0, y = ylimit-0.3, label = "MCT"),  
    hjust = -0.5, vjust = -2, size = 4.7  
  )+
  xlab("log2foldchange of circRNA") +
  ylab("log2foldchange of host_gene")  

# GO



GO_database <- "org.Rn.eg.db"
KEGG_database <- "rno"
GO<-enrichGO( filtered_data$host_gene_entrezgene_id,#GO富集分析
              OrgDb = GO_database ,
              keyType = "ENTREZID",#设定读取的gene ID类型
              ont = "ALL",#(ont为ALL因此包括 Bioal Process,Cellular Component,Mollecular Function三部分）
              pvalueCutoff = 1,#设定p值阈值
              qvalueCutoff = 1,#设定q值阈值
              readable = T)


KEGG<-enrichKEGG(filtered_data$host_gene_entrezgene_id,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 1,
                 qvalueCutoff = 1)
pathway <- c("GO:0048738","GO:0003012","GO:0034329","GO:0003712",	
             "GO:0055007","GO:0006338","GO:0022804", "GO:0098984","GO:0019838","GO:0015631")
dev.new()
pathway <- GO@result[GO@result$ID %in% pathway, "Description"]
enrichplot::cnetplot(GO,circular=TRUE,colorEdge = TRUE,showCategory = pathway)
barplot(KEGG,showCategory = 10,title = 'KEGG Pathway')

