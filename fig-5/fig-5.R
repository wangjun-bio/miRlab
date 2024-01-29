library(openxlsx)#读取.xlsx文件
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(ggridges)
#library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
#library(ComplexHeatmap)#绘制图例
#install.packages("enrichplot")
#BiocManager::install("HPO.db")
setwd("C:/Users/树金/Desktop/代码/fig-5")
getwd()

file_list <- list.files(pattern = "GSEA.csv")
data <- list()
for (file in file_list) {
  data[[file]] <- read.csv(file = file, header = T, sep = ",")
}
info <- data[[1]]
colnames(info)[1] <- "gene_symbol"

GO_database <- "org.Hs.eg.db"  #加载准换ID的数据库
KEGG_database <- "hsa"         #使用人的KEGG数据库
gene <- bitr(info$gene_symbol,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)    #将母源基因转换为GO数据库的代号

GO<-enrichGO(gene$ENTREZID,#GO富集分析
             OrgDb = GO_database,
             keyType = "ENTREZID",#设定读取的gene ID类型
             ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
             pvalueCutoff = 1,#设定p值阈值
             qvalueCutoff = 1,#设定q值阈值
             readable = T)

KEGG<-enrichKEGG(gene$ENTREZID,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 1,
                 qvalueCutoff = 1) 

names(info) <- c('SYMBOL','a','baseMean','Log2FoldChange','b','c','pvalue','padj','d')
info_merge <- merge(info,gene,by='SYMBOL')
GSEA_input <- info_merge$Log2FoldChange
names(GSEA_input) = info_merge$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE)
GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, pvalueCutoff = 1)#GSEA富集分析

barplot(GO, showCategory = 25,split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#柱状图
barplot(KEGG,showCategory = 50,title = 'KEGG Pathway')

dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#点状图
dotplot(KEGG)

ridgeplot(GSEA_KEGG,showCategory = 21) 
gseaplot2(GSEA_KEGG,4)
gseaplot2(GSEA_KEGG,1:9)#30是根据ridgeplot中有30个富集通路得到的
ggsave("GSEA.pdf", plot = last_plot(), device = "pdf", width = 12, height = 15)

a <- GO[,,]
b <- KEGG[,,]
c <- GSEA_KEGG[,,]
