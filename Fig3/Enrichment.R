library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(GOplot)



gene <- as.data.frame(read.csv('./de_tcga_filter.csv'))
gene <- gene %>% filter(group != 'NS')
gene = gene$X


#转换为ENTREZID
entrezIDs = bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb= "org.Hs.eg.db", drop = TRUE)
gene = entrezIDs$ENTREZID

##GO####
go <- enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, qvalueCutoff = 0.05,ont="all",readable =T,minGSSize = 10,
               maxGSSize = 500)
GO_result <- go@result
#write.csv(GO_result,'GO_all.csv')
GO_result <- go@result %>% filter(Count > 15)
GO_result$ONTOLOGY <- as.factor(GO_result$ONTOLOGY)
GOCircle(GO_result)

#############kegg###########
kk <- enrichKEGG(gene = gene,keyType = "kegg",organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05, pAdjustMethod = "fdr")
kegg_result <- kk@result
write.csv(kegg_result,'kegg_all.csv')
kegg_result <- kegg_result %>% filter(Count > 7)
kegg_result <- na.omit(kegg_result)

###############输出表格，可视化由CNS KNOWN网站实现##################





