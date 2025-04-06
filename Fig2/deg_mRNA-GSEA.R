# 安装并加载必要的包
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!require("clusterProfiler")) {
  BiocManager::install("clusterProfiler")
  library("clusterProfiler")
}
if (!require("org.Hs.eg.db")) {
  BiocManager::install("org.Hs.eg.db")
  library("org.Hs.eg.db")
}
if (!require("ggplot2")) {
  install.packages("ggplot2")
  library("ggplot2")
}
setwd('/Users/willow/Desktop/2025/cfRNA乳腺癌/DATA-cfRNA/差异分析/')
# 读取 CSV 文件
# 请将下面的文件路径替换为你实际的 CSV 文件路径
data <- read.csv("./MALvsBEN_mRNA.DESeq2.csv")

# 提取基因列表和对应的表达值
# 假设 CSV 文件中有 "gene_name" 和 "log2FC" 两列，可按需修改
gene_list <- data$logFC
names(gene_list) <- data$id

# 去除重复基因并按表达值排序
gene_list <- gene_list[!duplicated(names(gene_list))]
gene_list <- sort(gene_list, decreasing = TRUE)

# 将基因名称转换为 Entrez ID
gene_entrez <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_list <- gene_list[match(gene_entrez$SYMBOL, names(gene_list))]
names(gene_list) <- gene_entrez$ENTREZID

# 3. **基于 GO (ALL) 进行 GSEA**
gsea_go_official <- gseGO(
  geneList = gene_list, 
  OrgDb = org.Hs.eg.db, 
  ont = "ALL",  
  keyType = "ENTREZID", 
  minGSSize = 1,
  maxGSSize = 500,
  pvalueCutoff = 1,
  verbose = TRUE,
)


gsea_results <- as.data.frame(gsea_go_official)
name = 'GO'

# 过滤显著性结果
gsea_res_sig <- subset(gsea_results, p.adjust < 0.05)
title <- paste0(name, " (p.adjust < 0.05)")

# 绘制 dotplot
p <- dotplot(gsea_go_official, showCategory = gsea_res_sig$Description, title = title)
p

# 保存 dotplot
pdf_name <- paste0("./", name, "_dot.pdf")
#ggsave(p, filename = pdf_name, width = 6, height = 8)

# 取显著通路中 NES 最高的 3 个通路绘制 GSEA 曲线
top3_ids <- head(gsea_res_sig[order(gsea_res_sig$p.adjust, decreasing = TRUE), "ID"], 3)

# 检查 top3_ids 是否为空
if (length(top3_ids) > 0) {
  p_gsea <- gseaplot2(
    gsea_res, 
    geneSetID = top3_ids, 
    title = title,
    pvalue_table = TRUE, 
    subplots = 1:3,
    rel_heights = c(1.5, 0.3, 1)
  )
  
  print(p_gsea)
  
  # 保存 GSEA 曲线图
  pdf_name <- paste0("./", name, ".pdf")
  ggsave(p_gsea, filename = pdf_name, width = 10, height = 8)
} else {
  message("未找到足够的显著通路来绘制 GSEA 曲线。")
}



# 进行 KEGG GSEA 分析
kegg_gsea <- gseKEGG(geneList = gene_list,
                     organism = "hsa",
                     nPerm = 1000,
                     minGSSize = 1,
                     maxGSSize = 500,
                     pvalueCutoff = 1,
                     verbose = FALSE)

# 进行 GO GSEA 分析（这里以生物过程 "BP" 为例，也可选择 "MF" 或 "CC"）
go_gsea <- gseGO(geneList = gene_list,
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 nPerm = 1000,
                 minGSSize = 1,
                 maxGSSize = 500,
                 pvalueCutoff = 1,
                 verbose = FALSE)

# 绘制 KEGG 通路富集图
kegg_plot <- dotplot(kegg_gsea, showCategory = 15) + ggtitle("KEGG Pathway GSEA")

# 绘制 GO 通路富集图
go_plot <- dotplot(go_gsea, showCategory = 15) + ggtitle("GO Pathway GSEA (Biological Process)")

# 显示图形
print(kegg_plot)
print(go_plot)
