library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
###############################################################################################################################################
# 获取数据文件夹下的所有样本文件列表
samples_list <- list.files("./")
head(samples_list ,n=50)
# 创建一个空的列表来存储Seurat对象
seurat_list <- list()
# 读取每个样本的10x数据并创建Seurat对象
for (sample in samples_list) {
  data.path <- paste0("./", sample)  # 拼接文件路径
  seurat_data <- Read10X(data.dir = data.path) # 读取10x数据，data.dir参数指定存放文件的路径
  seurat_obj <- CreateSeuratObject(counts = seurat_data, min.cells = 3, min.features = 250, project = sample)  # 创建Seurat对象，并指定项目名称为样本文件名
  #过滤检测少于200个基因的细胞（min.features = 200）和少于3个细胞检测出的基因（min.cells = 3）
  seurat_list <- append(seurat_list, seurat_obj)  # 将Seurat对象添加到列表中
}
# 输出所有的Seurat对象列表
seurat_list
# 合并Seurat对象，将所有Seurat对象合并到一个对象中
seurat_combined <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list[-1]))
# 打印合并后的Seurat对象
print(seurat_combined)
table(seurat_combined$orig.ident)
#############################################################################################################################################
#在每一个分组中添加线粒体信息
seurat_combined[["percent.mt"]] <- PercentageFeatureSet(seurat_combined, pattern = "^MT-")
#展示前五列
head(seurat_combined@meta.data, 5)
#使用小提琴图可视化QC指标
#nFeature_RNA代表每个细胞测到的基因数目。
#nCount_RNA代表每个细胞测到所有基因的表达量之和。
#percent.mt代表测到的线粒体基因的比例。
p <- VlnPlot(seurat_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", ncol = 1,pt.size = 0)
ggsave("VlnPlot_mt.pdf", plot = p, width = 11, height = 11)
#FeatureScatter通常用于可视化 feature-feature 相关性，
#nCount_RNA 与percent.mt的相关性
plot1 <- FeatureScatter(seurat_combined, feature1 = "nFeature_RNA", feature2 = "nCount_RNA")
plot2 <- FeatureScatter(seurat_combined, feature1 = "percent.mt", feature2 = "nCount_RNA")
plot3 <- FeatureScatter(seurat_combined, feature1 = "percent.mt", feature2 = "nFeature_RNA")
p <- plot1 + plot2 + plot3 
ggsave("A.pdf", plot = p, width = 11, height = 9)
#过滤线粒体基因表达比例过高的细胞，和一些极值细胞（可以根据小提琴图判断，查看两端离群值）。
seurat_combined <- subset(seurat_combined, subset = nFeature_RNA > 100 & nFeature_RNA < 7500 & percent.mt < 35 & nCount_RNA > 1000)
print(seurat_combined)
p <- VlnPlot(seurat_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", ncol = 1,pt.size = 0)
ggsave("VlnPlot_mt-2.pdf", plot = p, width = 11, height = 11)
plot1 <- FeatureScatter(seurat_combined, feature1 = "nFeature_RNA", feature2 = "nCount_RNA")
plot2 <- FeatureScatter(seurat_combined, feature1 = "percent.mt", feature2 = "nCount_RNA")
plot3 <- FeatureScatter(seurat_combined, feature1 = "percent.mt", feature2 = "nFeature_RNA")
p <- plot1 + plot2 + plot3 
ggsave("B.pdf", plot = p, width = 20, height = 9)
##############################################################################################################################################
#鉴定高变基因
#计算数据集中表现出高细胞间变异的特征基因(即，它们在某些细胞中高表达，而在其他细胞中低表达)。这些基因有助于突出单细胞数据集中的生物信号
#每个数据集返回2000个features 。这些将用于下游分析，如PCA。
seurat_combined <- NormalizeData(seurat_combined, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_combined <- FindVariableFeatures(seurat_combined, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seurat_combined), 10)
head(top10)
p <- VariableFeaturePlot(seurat_combined)
ggsave("K.pdf", plot = p, width = 20, height = 9)
all.genes <- rownames(seurat_combined)
seurat_combined <- ScaleData(seurat_combined, features = all.genes)
#线性降维,接下来，对缩放的数据执行PCA。
seurat_combined <- RunPCA(seurat_combined, features = VariableFeatures(object = seurat_combined))
#查看PCA结果
print(seurat_combined[["pca"]], dims = 1:10, nfeatures = 5)
p <- VizDimLoadings(seurat_combined, dims = 1:2, reduction = "pca")
ggsave("M.pdf", plot = p, width = 7, height = 7)
DimPlot(seurat_combined, reduction = "pca", group.by = "orig.ident")
#DimHeatmap()可以方便地探索数据集中异质性的主要来源，并且可以确定哪些PC维度可以用于下一步的下游分析。细胞和基因根据PCA分数来排序。
DimHeatmap(seurat_combined, dims = 1, cells = 500, balanced = TRUE) #1个PC 500个细胞
DimHeatmap(seurat_combined, dims = 1:15, cells = 500, balanced = TRUE) #15个PC
##############################################################################################################################################
#筛子最佳dim的数值，
seurat_combined <- JackStraw(seurat_combined, num.replicate = 200)
seurat_combined <- ScoreJackStraw(seurat_combined, dims = 1:20)
p <- JackStrawPlot(seurat_combined, dims = 1:20)
ggsave("D.pdf", plot = p, width = 20, height = 9)
p <- ElbowPlot(seurat_combined,ndims = 40 )
ggsave("C.pdf", plot = p, width = 20, height = 9)
##############################################################################################################################################
seurat_combined <- FindNeighbors(object = seurat_combined, dims = 1:30)
seurat_combined <- FindClusters(seurat_combined, resolution = 0.4)
seurat_combined <- RunUMAP(seurat_combined, dims = 1:30)
seurat_combined <- RunTSNE(seurat_combined, dims = 1:30)
#展示绘图
p1 <- DimPlot(object = seurat_combined, dims = c(1,2), reduction = "umap", group.by = c("seurat_clusters"), label = TRUE)
ggsave("umap_Plot.pdf", plot = p1, width = 10, height = 9)
p2 <- DimPlot(object = seurat_combined, dims = c(1,2), reduction = "tsne", group.by = c("seurat_clusters"), label = TRUE)
ggsave("tsne_Plot.pdf", plot = p2, width = 10, height = 9)
##############################################################################################################################################
seurat_combined <- JoinLayers(seurat_combined, features = c("RNA"))
pbmc.markers <- FindAllMarkers(seurat_combined, only.pos = TRUE, min.pct = 0.15, logfc.threshold = 0.35)
colnames(pbmc.markers)[colnames(pbmc.markers) == "p_val"] <- "P_value"
colnames(pbmc.markers)[colnames(pbmc.markers) == "p_val_adj"] <- "padj"
significant_markers <- pbmc.markers[abs(pbmc.markers$P_value) <= 0.05 & pbmc.markers$padj <= 0.05, ]
write.table(significant_markers,file="markers.xls",sep="\t",row.names=F,quote=F)
cluster_0_markers <- significant_markers[significant_markers$cluster == 0, ]
cluster_1_markers <- significant_markers[significant_markers$cluster == 1, ]
cluster_2_markers <- significant_markers[significant_markers$cluster == 2, ]
cluster_3_markers <- significant_markers[significant_markers$cluster == 3, ]
cluster_4_markers <- significant_markers[significant_markers$cluster == 4, ]
cluster_5_markers <- significant_markers[significant_markers$cluster == 5, ]
cluster_6_markers <- significant_markers[significant_markers$cluster == 6, ]
cluster_7_markers <- significant_markers[significant_markers$cluster == 7, ]
cluster_8_markers <- significant_markers[significant_markers$cluster == 8, ]
cluster_9_markers <- significant_markers[significant_markers$cluster == 9, ]
cluster_10_markers <- significant_markers[significant_markers$cluster == 10, ]
cluster_11_markers <- significant_markers[significant_markers$cluster == 11, ]
cluster_12_markers <- significant_markers[significant_markers$cluster == 12, ]
cluster_13_markers <- significant_markers[significant_markers$cluster == 13, ]
cluster_14_markers <- significant_markers[significant_markers$cluster == 14, ]
cluster_15_markers <- significant_markers[significant_markers$cluster == 15, ]
cluster_16_markers <- significant_markers[significant_markers$cluster == 16, ]
cluster_17_markers <- significant_markers[significant_markers$cluster == 17, ]
cluster_18_markers <- significant_markers[significant_markers$cluster == 18, ]
cluster_19_markers <- significant_markers[significant_markers$cluster == 19, ]
cluster_20_markers <- significant_markers[significant_markers$cluster == 20, ]
cluster_21_markers <- significant_markers[significant_markers$cluster == 21, ]
#######################################################################################################################################
seurat_combined$cell_type <- "Red blood cell (erythrocyte)" # 默认为Red blood cell (erythrocyte)
#细胞注释
seurat_combined$cell_type <- "NA" # 默认为NA
seurat_combined$cell_type[seurat_combined$seurat_clusters == 2] <- "Tcm+Tfh+Treg"
seurat_combined$cell_type[seurat_combined$seurat_clusters %in% c(0, 3, 5, 16)] <- "CD8+T"
seurat_combined$cell_type[seurat_combined$seurat_clusters == 9] <- "plasma B"
seurat_combined$cell_type[seurat_combined$seurat_clusters %in% c(4, 19)] <- "NKT"
seurat_combined$cell_type[seurat_combined$seurat_clusters %in% c(17, 11)] <- "goblet"
seurat_combined$cell_type[seurat_combined$seurat_clusters %in% c(13, 1)] <- "follicular B"
seurat_combined$cell_type[seurat_combined$seurat_clusters == 18] <- "fibroblasts"
seurat_combined$cell_type[seurat_combined$seurat_clusters %in% c(12, 10, 7, 6)] <- "epithelial_carcinoma"
seurat_combined$cell_type[seurat_combined$seurat_clusters == 8] <- "epithelial_adenoma"
seurat_combined$cell_type[seurat_combined$seurat_clusters %in% c(14, 15)] <- "enterocyte"
seurat_combined$cell_type[seurat_combined$seurat_clusters == 20] <- "endothelial"
p3 <- DimPlot(object = seurat_combined, dims = c(1,2), reduction = "tsne", group.by = c("cell_type"), label = F)
ggsave("tsne_Plot-2.pdf", plot = p3, width = 11, height = 9)
########################################################################################################################################
#组织类型的注释
seurat_combined$Tissue_Type <- "normal" 
seurat_combined$Tissue_Type[seurat_combined$orig.ident %in% c("Patient_0", "Patient_1_carcinoma", "Patient_2_carcinoma")] <- "CRC"
seurat_combined$Tissue_Type[seurat_combined$orig.ident %in% c("Patient_1_adenoma", "Patient_2_adenoma", "Patient_3_adenoma_1", "Patient_3_adenoma_2")] <- "adenoma"
seurat_combined$Tissue_Type[seurat_combined$orig.ident == "Patient_3_blood"] <- "blood"
tissue_type_counts <- table(seurat_combined@meta.data$Tissue_Type)
print(tissue_type_counts)
########################################################################################################################################
#患者类型的注释
seurat_combined$Pation_Type <- "Patient_0" # 默认为Red blood cell (erythrocyte)
seurat_combined$Pation_Type[seurat_combined$orig.ident %in% c("Patient_1_adenoma", "Patient_1_carcinoma", "Patient_1_normal")] <- "Patient_1"
seurat_combined$Pation_Type[seurat_combined$orig.ident %in% c("Patient_2_adenoma", "Patient_2_carcinoma", "Patient_2_normal")] <- "Patient_2"
seurat_combined$Pation_Type[seurat_combined$orig.ident %in% c("Patient_3_adenoma_1","Patient_3_adenoma_2", "Patient_3_carcinoma", "Patient_3_normal", "Patient_3_blood")] <- "Patient_3"
Pation_type_counts <- table(seurat_combined@meta.data$Pation_Type)
print(Pation_type_counts)
#######################################################################################################################################
# 提取Tissue_Type列中内容为CRC的行
crc_rows <- subset(seurat_combined, subset = Tissue_Type == "CRC")
colnames(crc_rows@meta.data)[colnames(crc_rows@meta.data) == "cell_type"] <- "CRC_cell_type"
NC_rows <- subset(seurat_combined, subset = Tissue_Type == "normal")
colnames(NC_rows@meta.data)[colnames(NC_rows@meta.data) == "cell_type"] <- "normal_cell_type"
adenoma_rows <- subset(seurat_combined, subset = Tissue_Type == "adenoma")
colnames(adenoma_rows@meta.data)[colnames(adenoma_rows@meta.data) == "cell_type"] <- "adenoma_cell_type"
p3 <- DimPlot(object = adenoma_rows, dims = c(1,2), reduction = "tsne", group.by = c("adenoma_cell_type"), label = F)
blood_rows <- subset(seurat_combined, subset = Tissue_Type == "blood")
colnames(blood_rows@meta.data)[colnames(blood_rows@meta.data) == "cell_type"] <- "blood_cell_type"
#################################################################################################################################################
#绘制三合一
p1 <- DimPlot(object = seurat_combined, dims = c(1,2), reduction = "tsne", group.by = c("Tissue_Type"), label = F)
p2 <- DimPlot(object = seurat_combined, dims = c(1,2), reduction = "tsne", group.by = c("cell_type"), label = F)
p3 <- DimPlot(object = seurat_combined, dims = c(1,2), reduction = "tsne", group.by = c("Pation_Type"), label = F)
p4 <- p1+p2+p3
ggsave("tsne-3.pdf", plot = p4, width = 27, height = 7)
p4 <- DimPlot(object = seurat_combined, dims = c(1,2), reduction = "tsne", group.by = c("seurat_clusters"), label = TRUE)
ggsave("tsne.pdf", plot = p, width = 15, height = 8)
###########################################################################################################################
genes_of_interest <- c("CAMSAP1","LPAR1","FAT1","EDIL3","PLCE1","SKA3","PTPN22","TMEM181","SPATA13","EPB41L2","PRKAR1B")
p5 <- DotPlot(NC_rows, features = genes_of_interest,dot.scale = 13,dot.min = 0,
              assay = "RNA" , group.by = "normal_cell_type") +
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 90, hjust = 0.5,vjust=0.5))+ #轴标签
  coord_flip() + #翻转
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("N Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('blue','red')) #颜色
df_N<- p5$data
write.csv(df_N,file = "./scRNA_N.csv")

p6 <- DotPlot(adenoma_rows, features = genes_of_interest,dot.scale = 13,dot.min = 0,
              assay = "RNA" , group.by = "adenoma_cell_type") +
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 90, hjust = 0.5,vjust=0.5))+ #轴标签
  coord_flip() + #翻转
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("P Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('blue','red')) #颜色
df_P<- p6$data
write.csv(df_P,file = "./scRNA_P.csv")
p7 <- DotPlot(crc_rows, features = genes_of_interest,dot.scale = 13,dot.min = 0,
              assay = "RNA" , group.by = "CRC_cell_type") +
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 90, hjust = 0.5,vjust=0.5))+ #轴标签
  coord_flip() + #翻转
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("T Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('blue','red')) #颜色
df_T<- p7$data
write.csv(df_T,file = "./scRNA_T.csv")
p <- ggplot(scRNA_N, aes(x = Group, y = id, color = avg.exp.scaled, size = pct.exp)) +
  geom_point() +
  scale_size_continuous(range = c(1, 11)) +  # 调整气泡大小范围
  scale_color_gradient2(low = "gray", mid = "white", high = "red", midpoint = 0.5) +  # 调整颜色范围
  labs(title = "Bubble Plot", x = "Group", y = "ID", color = "Avg.Exp.Scaled", size = "Pct.Exp") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0))
p
getwd()
ggsave("circRNA_11.pdf", plot = p, width = 15, height = 8)
###########################################################################################################################################
#绘制每一个聚类模块的top5的气泡图
DotPlot(seurat_combined, features = cluster_11_markers$gene, assay = "RNA", group.by = "seurat_clusters")
unique(seurat_combined@meta.data$Tissue_Type)

cell_type_levels <- c("CD8+T","Tcm+Tfh+Treg","endothelial", "epithelial_adenoma","epithelial_carcinoma", 
                      "enterocyte","goblet","fibroblasts","plasma B", "NKT","follicular B")

seurat_combined$cell_type <- factor(seurat_combined$cell_type, levels = cell_type_levels)
seurat_combined$cell_type_level <- as.integer(seurat_combined$cell_type)

genes_of_interest <- c("CD8A", "CD8B", "CD3G", "CD3D", "CD3E", "KLRD1", "TRBC1", "TRAC",
                       "VWF","IL7R","CD2","PTPRC","LTB","SLC26A3","TFF3","ICA1","GUCA2A","GUCA2B","EPCAM",
                       "DCN","MUC2", "KLRF1", "FCGR3A","MZB1","MS4A1", "CD79A", "CD79B")

p10 <- DotPlot(seurat_combined, features = genes_of_interest,
               assay = "RNA" , group.by = "cell_type") + 
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  coord_flip() + #翻转
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色

ggsave("MARKER.pdf", plot = p10, width = 11, height = 9)
#############################################################################################################################################
p11 <- FeaturePlot(seurat_combined, features = c("CAMSAP1","LPAR1","FAT1","EDIL3","PLCE1","SKA3","PTPN22","TMEM181","SPATA13","EPB41L2","PRKAR1B"))
ggsave("FeaturePlot.pdf", plot = p11, width = 24, height = 11)

p11 <- FeaturePlot(object = seurat_combined, cols = c("#336699", "red"), features = c("CAMSAP1","LPAR1","FAT1","EDIL3","PLCE1","SKA3","PTPN22",
                                                                                      "TMEM181","SPATA13","EPB41L2","PRKAR1B"))
#############################################################################################################################################

genes_of_interest <- c("CAMSAP1","LPAR1","FAT1","EDIL3","PLCE1","SKA3","PTPN22","TMEM181","SPATA13","EPB41L2","PRKAR1B")

p <- VlnPlot(seurat_combined, features = genes_of_interest, group.by = "cell_type", split.by ="Tissue_Type", sort = T,adjust = 1,log = T,combine = TRUE, ncol = 4,pt.size = 0.01)
ggsave("w4.pdf", plot = p, width = 20, height = 15)




