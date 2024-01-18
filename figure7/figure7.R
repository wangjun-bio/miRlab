if(!require("multtest")){
  BiocManager::install("multtest")
  library("multtest")
}
if(!require("Seurat")){
  BiocManager::install("Seurat")
  library("Seurat")
}
if(!require("dplyr")){
  BiocManager::install("dplyr")
  library("dplyr")
}
if(!require("patchwork")){
  BiocManager::install("patchwork")
  library("patchwork")
}
if(!require("R.utils")){
  BiocManager::install("R.utils")
  library("R.utils")
}
if(!require("ggplot2")){
  BiocManager::install("ggplot2")
  library("ggplot2")
}
if(!require("sctransform")){
  BiocManager::install("sctransform")
  library("sctransform")
}
if (!require("biomaRt")) {
  BiocManager::install("biomaRt")
  library("biomaRt")
}
if (!require("metap")) {
  BiocManager::install("metap")
  library("metap")
}
if (!require("readxl")) {
  BiocManager::install("readxl")
  library(readxl)
}

devtools::install_github('immunogenomics/presto')

output_cell_distribution <- function(umap, name){
  allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
              "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
              "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
              "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  p <- ggplot(umap,aes(x= umap_1 , y = umap_2 ,color = cell_type)) +  
    geom_point(size = 1 , alpha =1 )  +  
    scale_color_manual(values = allcolour)
  output_distribution(p = p,name = name)
}

output_MA_distribution <- function(umap, name){
  allcolour=c("#0000FF","grey","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
              "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
              "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
              "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  p <- ggplot(umap,aes(x= umap_1 , y = umap_2 ,color = MA)) +  
    geom_point(size = 1 , alpha =1 )  +  
    scale_color_manual(values = allcolour)
  output_distribution(p = p,name = name)
}
output_SA_distribution <- function(umap, name){
  allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
              "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
              "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
              "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  p <- ggplot(umap,aes(x= umap_1 , y = umap_2 ,color = SA)) +  
    geom_point(size = 1 , alpha =1 )  +  
    scale_color_manual(values = allcolour)
  output_distribution(p = p,name = name)
}
output_HA_distribution <- function(umap, name){
  allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
              "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
              "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
              "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  p <- ggplot(umap,aes(x= umap_1 , y = umap_2 ,color = HA)) +  
    geom_point(size = 1 , alpha =1 )  +  
    scale_color_manual(values = allcolour)
  output_distribution(p = p,name = name)
}
output_distribution <- function(p, name){
  p2 <- p  +
    theme(panel.grid.major = element_blank(), #主网格线
          panel.grid.minor = element_blank(), #次网格线
          panel.border = element_blank(), #边框
          axis.title = element_blank(),  #轴标题
          axis.text = element_blank(), # 文本
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = 'white'), #背景色
          plot.background=element_rect(fill="white"))
  
  p3 <- p2 +         
    theme(
      legend.title = element_blank(), #去掉legend.title 
      legend.key=element_rect(fill='white'), #
      legend.text = element_text(size=20), #设置legend标签的大小
      legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
    guides(color = guide_legend(override.aes = list(size=5))) #设置legend中 点的大小 
  
  p4 <- p3 + 
    geom_segment(aes(x = min(umap$umap_1) , y = min(umap$umap_2) ,
                     xend = min(umap$umap_1) +3, yend = min(umap$umap_2) ),
                 colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
    geom_segment(aes(x = min(umap$umap_1)  , y = min(umap$umap_2)  ,
                     xend = min(umap$umap_1) , yend = min(umap$umap_2) + 3),
                 colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
    annotate("text", x = min(umap$umap_1) +1.5, y = min(umap$umap_2) -1, label = "UMAP_1",
             color="black",size = 3, fontface="bold" ) + 
    annotate("text", x = min(umap$umap_1) -1, y = min(umap$umap_2) + 1.5, label = "UMAP_2",
             color="black",size = 3, fontface="bold" ,angle=90) 
  
  ggsave(p4,file=paste(name,"_cell_distribution.pdf",sep="")) 
  return(p4)
  
}


MA <- read.csv("MA_geneInfo.csv", header = TRUE, check.names = FALSE)
SA <- read.csv("SA_geneInfo.csv", header = TRUE, check.names = FALSE)
HA <- read.csv("HA_geneInfo.csv", header = TRUE, check.names = FALSE)


# 读取注释
GSE210248_Human_PA_Metadata <- read.csv2("GSE210248_Human_PA_Metadata.csv", row.names = 1)

# markers
cell_markers <- read.csv("cell_markers.csv")
# 数据读取
PA.data <- Read10X(data.dir = "Human_PA_Integrated")

PA <- CreateSeuratObject(counts = PA.data, project = "human_PA", min.cells = 3, min.features = 200)




# 添加注释
PA <- AddMetaData(object = PA, metadata = GSE210248_Human_PA_Metadata, col.name = colnames(GSE210248_Human_PA_Metadata))

# 提取单个细胞检测基因在200-2500个的细胞，太少测序深度低，太多可能为多个细胞建库(mt占比反应细胞死亡)
PA_Set <- subset(PA, subset = nFeature_RNA >= 200 & percent.mito <= 5)

# 对单个细胞转录组表达进行标准化
PA_Set <- NormalizeData(PA, normalization.method = "LogNormalize", scale.factor = 10000)

#  寻找高变基因
PA_Set <- FindVariableFeatures(PA_Set, selection.method = "vst", nfeatures = 2000)

g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(PA_Set))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(PA_Set))
PA_Set <- CellCycleScoring(object=PA_Set,  g2m.features=g2m_genes,  s.features=s_genes)

# 单个基因在所有细胞中进行标准化
PA_Set <- ScaleData(PA_Set, features = rownames(PA_Set))

# 计算PCA，进行降维
PA_Set <- RunPCA(PA_Set, features = VariableFeatures(object = PA_Set))
DimHeatmap(PA_Set, dims = 1:15, cells = 2000, balanced = TRUE)

plot1 <- DimPlot(PA_Set, reduction = "pca", group.by="orig.ident")
#（左图）根据主成分1和2的值将细胞在平面上展示出来
plot2 <- ElbowPlot(PA_Set, ndims=20, reduction="pca") 
#（右图）展示前20个主成分的解释量
plot1+plot2

# 选取主成分数量，上图斜率变小几乎不变的点开始截取
pc.num=1:8
# dims参数，需要指定哪些pc轴用于分析；这里利用上面的分析，选择18
PA_Set <- FindNeighbors(PA_Set, dims = pc.num) 
PA_Set <- FindClusters(PA_Set, resolution = 0.2)

# 2、UMAP降维（tSNE、UMAP二选一）
PA_Set <- RunUMAP(PA_Set, dims = pc.num)
embed_umap <- Embeddings(PA_Set, 'umap')
PA_Set$Cell_annotation <- gsub("2","",PA_Set$Cell_annotation )
PA_Set$Cell_annotation <- gsub("1","",PA_Set$Cell_annotation )
Idents(PA_Set) <- PA_Set$Cell_annotation

DimPlot(PA_Set, reduction = "umap")


umap <-  PA_Set@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = PA_Set@meta.data$Cell_annotation) # 注释后的label信息 ，改为cell_type

a <- output_cell_distribution(umap = umap, name = "all")

Object_list <- SplitObject(PA_Set, split.by = "disease")

PAH_object <- Object_list$PAH
Donor_object <- Object_list$Donor


PA_umap <-  PAH_object@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = PAH_object@meta.data$Cell_annotation) # 注释后的label信息 ，改为cell_type

control_umap <-  Donor_object@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = Donor_object@meta.data$Cell_annotation) # 注释后的label信息 ，改为cell_type

# 总的
output_cell_distribution(umap = PA_umap, name = "PA")
output_cell_distribution(umap = control_umap, name = "control")

# 计算基因在每种细胞平均表达量
PA_cell_average <- AverageExpression(PAH_object, group.by = "Cell_annotation")

# 提取矩阵
PA_cell_average <-  as.data.frame(PA_cell_average$RNA)

markers <- FindAllMarkers(object = PAH_object, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25)  


avg_expr_cut_off <- 0.5
# 高表达基因
high_expr_gene <- PA_cell_average %>%
  filter(rowSums(. > avg_expr_cut_off) > 0) %>% #任意细胞高表达
  cbind(geneName = rownames(.), .) # 添加列


# 高表达的marker基因
high_expr_markers <- markers %>%
  filter(p_val < 0.05, pct.1 > 0.5, pct.2 < 0.5) %>% # 在单种细胞中与其他细胞类型差异显著，变化倍数在2倍markers
  semi_join(high_expr_gene, by = c("gene" = "geneName")) #属于高表达基因


GS_cut_off <- 0
MM_cut_off <- 0

modules <- list(MA = c("brown"), SA = "yellow", HA = c("turquoise"))
PA_list <- list(MA = MA, SA = SA, HA = HA)

name_list <- names(modules)
for (i in 1:length(modules)) {
  mod <- modules[[i]]
  # 筛选对应颜色基因
  dat <- PA_list[[i]]
  dat <- dat[dat$moduleColor %in% mod,]
  
  # 根据筛选条件进行筛选
  dat <- dat[abs(dat[[paste0("GS.", name_list[i])]]) >= GS_cut_off ,]
  filter_1 <- dat[[paste0("MM.", mod[1])]] >= MM_cut_off
  if(length(mod) == 2){
    filter_2 <- dat[[paste0("MM.", mod[2])]] >= MM_cut_off
    dat <- dat[filter_1 | filter_2, ]
    
  }else{
    dat <- dat[filter_1, ]
  }
  
  PA_list[[i]] <- dat
}



# 连接Ensembl数据库，跨物种转换，构建数据集需要指定网站
ensembl <- useMart("ensembl",host = "https://dec2021.archive.ensembl.org/")
human_data <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
rat_data <- useDataset("rnorvegicus_gene_ensembl", mart = ensembl)


MA_human_genes <- getLDS(attributes = c("entrezgene_id"),
                         filters = "ensembl_gene_id",
                         values = unique(PA_list[["MA"]]$host_gene_id),
                         mart = rat_data, 
                         attributesL = c("external_gene_name"),
                         martL = human_data)

MA_genes <- base::merge(MA_human_genes,PA_list[["MA"]],by.x = "NCBI.gene..formerly.Entrezgene..ID", by.y = "host_gene_entrezgene_id", all = FALSE)


SA_human_genes <- getLDS(attributes = c("entrezgene_id"),
                         filters = "ensembl_gene_id",
                         values = unique(PA_list[["SA"]]$host_gene_id),
                         mart = rat_data, 
                         attributesL = c("external_gene_name"),
                         martL = human_data)

SA_genes <- base::merge(SA_human_genes,PA_list[["SA"]],by.x = "NCBI.gene..formerly.Entrezgene..ID", by.y = "host_gene_entrezgene_id", all = FALSE)


HA_human_genes <- getLDS(attributes = c("entrezgene_id"),
                         filters = "ensembl_gene_id",
                         values = unique(PA_list[["HA"]]$host_gene_id),
                         mart = rat_data, 
                         attributesL = c("external_gene_name"),
                         martL = human_data)
HA_genes <- base::merge(HA_human_genes,PA_list[["HA"]],by.x = "NCBI.gene..formerly.Entrezgene..ID", by.y = "host_gene_entrezgene_id", all = FALSE)


PA_list$MA <- MA_genes %>%
  inner_join(high_expr_markers, by = c("Gene.name" = "gene"))
PA_list$SA <- SA_genes %>%
  inner_join(high_expr_markers, by = c("Gene.name" = "gene"))
PA_list$HA <- HA_genes %>%
  inner_join(high_expr_markers, by = c("Gene.name" = "gene"))


# 设置排序顺序
PAH_object$Cell_annotation <- factor(PAH_object$Cell_annotation, levels = c("T / NK cells",
                                                                            "T cells",
                                                                            "B cells",
                                                                            "NK cells",
                                                                            "DC",
                                                                            "Granulocytes",
                                                                            "Mast cells",
                                                                            "Mono / Macs",
                                                                            "Epithelial",
                                                                            "Endo ",
                                                                            "Fibro",
                                                                            "SMC "))


names <- names(PA_list)
cut_off <- 10
vlnList <- list()
DotPlot_list <- list()

genes <- c(sort(unique(PA_list$MA$Gene.name)), sort(unique(PA_list$SA$Gene.name)), sort(unique(PA_list$HA$Gene.name)))


MA_num <- length(sort(unique(PA_list[[1]]$Gene.name)))
SA_num <- length(sort(unique(PA_list[[2]]$Gene.name)))
HA_num <- length(sort(unique(PA_list[[3]]$Gene.name)))

group <- c(rep("MA", MA_num), rep("SA", SA_num), rep("HA", HA_num))
group <- factor(group, levels = c("HA","SA","MA"))


p1 <- DotPlot(PAH_object, features = split(genes, group),
              group.by = "Cell_annotation",
              dot.scale = 9)+RotatedAxis()


p1+ theme(
  # 面板
  panel.border = element_rect(color="black"), #面板边框
  panel.spacing = unit(1, "mm"), #面板间距
  
  # 分面标题
  #strip.background = element_rect(color="red"),
  strip.text = element_text(margin=margin(b=3, unit="mm")),
  strip.placement = 'outlet', #
  
  # 坐标轴线
  axis.line = element_blank(),
)+labs(x="", y="")+
  scale_color_gradient(low="#EEC03A",high ="#691E4D")


vln_plots <- VlnPlot(PA_Set, features = genes,
                     group.by = "Cell_annotation", split.by = "disease", ncol = 1)


featureplots <- list()
for (i in 1:length(genes)) {
  gene <- genes[i]
  featureplots[[i]] <- FeaturePlot(PA_Set, features = gene)
}


for (i in 1:length(PA_list)) {
  m <- PA_list[[i]] %>% 
    arrange(desc(sum_weight)) 
  vlnList[[names(PA_list)[i]]] <- VlnPlot(PA_Set, features = m$Gene.name,group.by = "Cell_annotation", split.by = "disease")
}


cell_markers <- cell_markers %>% 
  group_by(Annotation) %>% 
  slice_head(n = 3) %>% 
  ungroup()


markers_plot <- DotPlot(PAH_object, features = cell_markers$gene,
                        group.by = "Cell_annotation",
                        dot.scale = 9)

markers_plot2 <- markers_plot + 
  scale_color_gradient(low="#EEC03A",high ="#691E4D")+
  theme(axis.text.x = element_text(angle = 60,vjust = 0.85,hjust = 0.75))
markers_plot2
