rm(list=ls())
options(stringsAsFactors=F)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(data.table)
library(stringr)
library(patchwork)
library(harmony)

dir = list.dirs("D:/乳腺癌/singlecell_data2/GSE161529_RAW/")[-1]#设置数据存放路径
dir
names(dir) <- list.files('./singlecell_data2/GSE161529_RAW/',recursive = F)
dir
sce.all=CreateSeuratObject(counts = Read10X(dir),                            
                           min.cells = 5,                           
                           min.features = 300,)
dim(sce.all)
table(sce.all@meta.data$orig.ident)

library(dplyr)
sce.all@meta.data$group <- case_when(
  sce.all@meta.data$orig.ident == "GSM4909315" ~ "ER",
  sce.all@meta.data$orig.ident == "GSM4909290" ~ "HER2",
  sce.all@meta.data$orig.ident == "GSM4909282" ~ "TN"
)
table(sce.all$group)


####计算线粒体
sce.all[["mt_percent"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
####计算红细胞
HB_genes <- c('HBA1','HBA2','HBB','HBD','HBE1','HBG1','HBG2','HBM','HBQ1','HBZ')
HB_m <- match(HB_genes,rownames(sce.all@assays$RNA))
HB_genes <- rownames(sce.all@assays$RNA)[HB_m]
HB_genes <- HB_genes[!is.na(HB_genes)]
sce.all[['HB_percent']] <- PercentageFeatureSet(sce.all,features = HB_genes)

sce_filter <- subset(sce.all,subset = nFeature_RNA > 300 & 
                       nFeature_RNA < 5000 & 
                       #nCount_RNA < quantile(nCount_RNA,0.97) &
                       nCount_RNA < 7000 &
                       mt_percent < 20 &
                       HB_percent < 5)
table(sce_filter$group)
dim(sce_filter)

all.gene = rownames(sce_filter)

violin = VlnPlot(sce_filter,
                 features = c('nFeature_RNA','nCount_RNA','mt_percent','HB_percent'),
                 pt.size = 0.01,
                 ncol = 4,
                 assay = "RNA", 
                 slot = "counts")
violin
ggsave(plot = violin,filename = './singlecell2/violin.png',width = 15,height = 7,dpi = 300)
ggsave(plot = violin,filename = './singlecell2/violin.pdf',width = 15,height = 7)

#########标准化
sce_filter = NormalizeData(sce_filter,normalization.method="LogNormalize",scale.factor=10000) %>%
  FindVariableFeatures(selection.method = "vst",nfeatures=2000) %>% ScaleData(features = all.gene) %>%
  RunPCA(verbose = T,npcs = 30)
scRNA_harmony = RunHarmony(sce_filter,group.by.vars = 'orig.ident')

#######聚类
ElbowPlot(scRNA_harmony,ndims = 50,reduction = 'harmony')
scRNA_harmony <- FindNeighbors(scRNA_harmony,reduction = 'harmony',dims = 1:7) %>% FindClusters(resolution = 0.2)
scRNA_harmony <- RunUMAP(scRNA_harmony,reduction = 'harmony',dims = 1:7)

############注释
markers <- c('CD3D','CD3E','CD8A' ,#T
             'CD68','CD14','CD163',###TAM
             'PDGFRA',#fibroblasts
             'MCAM',#pericytes
             'NKG7',#NK
             "CD19","CD79A","MS4A1",'IGHG3', # B
             'IGLC7','IGHA1',##plasma 
             'VWF','CD34',### endothelial
            'EPCAM','EGFR','KRT8',###epithelial
            'COL1A1','COL3A1','DCN',##CAFs
            "CLEC9A","XCR1",#DC
            "KIT","TPSAB1","TPSB2",'CPA3'#Mast
             )

library(reshape2)
library(dplyr)
library(tibble)

scRNA_harmony <- readRDS('./singlecell2/scRNA_harmony.rds')
Idents(scRNA_harmony) <- 'seurat_clusters'
p <- DotPlot(scRNA_harmony, features = markers, cols = "RdYlBu") +  RotatedAxis() +  
  xlab(NULL) + ylab(NULL) +
    theme(strip.text = element_text(angle = 50,size = 10,face = 'bold.italic'),
        axis.text.y = element_text(size = 10,face = 'bold.italic'),        
        axis.text.x = element_text(angle = 45,size = 10,face = 'bold.italic'),        
        legend.text = element_text(size = 9,face = 'bold.italic'),        
        legend.title = element_text(size = 9,face = 'bold.italic'))
p
ggsave(plot = p,filename = './singlecell2/markers.pdf',height = 6,width = 8)
ggsave(plot = p,filename = './singlecell2/markers.png',height = 6,width = 8,dpi = 600)

source("./custom_seurat_functions.R")
library(data.table)
library(stringr)
writeLines(paste0(0:7,','))
celltype = read.table('./singlecell2/celltype.txt',sep = ',')
celltype

#重新规定
scRNA_harmony$celltype <- plyr::mapvalues(scRNA_harmony$seurat_clusters,
                                          from = 0:7,
                                          to = celltype$V2
)
umap2 <- DimPlot(scRNA_harmony, reduction = "umap",group.by = 'group',
                 label = T,label.box = T,label.size = 4,repel = T) + theme_bw() +
  labs( x= "Umap 1",y= "Umap 2",title = NULL) +
  theme(panel.grid=element_blank(), # 去网格线
        plot.title = element_text(size = 15,color="black",hjust = 0.5),
        axis.text.x = element_text(size = 13, color = 'black',face = 'bold.italic'),
        axis.text.y = element_text(size = 13, color = 'black',face = 'bold.italic'),
        axis.title.x = element_text(size = 14, color = 'black',face = 'bold.italic'),
        axis.title.y = element_text(size = 14, color = 'black',face = 'bold.italic'),
        axis.ticks = element_line(color = 'black', lineend = 'round'),
        legend.position = 'bottom',
        legend.text = element_text(size = 13, color = 'black',face = 'bold.italic'),
        legend.title = element_text(size = 13, color = 'black',face = 'bold.italic'),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  scale_fill_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.65)) +
  scale_color_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.65))


p2 = DimPlot(scRNA_harmony, reduction = "umap",group.by = "celltype",
             label = T,label.box = T,label.size = 4,repel = T) + theme_bw() +
  labs( x= "Umap 1",y= "Umap 2",title = NULL) +
  theme(panel.grid=element_blank(), # 去网格线
        plot.title = element_text(size = 15,color="black",hjust = 0.5),
        axis.text.x = element_text(size = 13, color = 'black',face = 'bold.italic'),
        axis.text.y = element_text(size = 13, color = 'black',face = 'bold.italic'),
        axis.title.x = element_text(size = 14, color = 'black',face = 'bold.italic'),
        axis.title.y = element_text(size = 14, color = 'black',face = 'bold.italic'),
        axis.ticks = element_line(color = 'black', lineend = 'round'),
        legend.position = 'bottom',
        legend.text = element_text(size = 13, color = 'black',,face = 'bold.italic'),
        legend.title = element_text(size = 13, color = 'black',face = 'bold.italic'),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  scale_fill_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.65)) +
  scale_color_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.65))
p3 <- p2 + umap2
p3
ggsave(plot = p3,filename = './singlecell2/umap.pdf',height = 6,width = 8)
ggsave(plot = p3,filename = './singlecell2/umap.png',height = 6,width = 8,dpi = 600)

##############
Idents(scRNA_harmony) <- 'celltype'
scRNA_harmony<-JoinLayers(scRNA_harmony)
brac.markers<-FindAllMarkers(scRNA_harmony,only.pos=T,min.pct=0.25,logfc.threshold=0.25,verbose=FALSE)
top5=brac.markers%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)
g = unique(top5$gene)
p<-DotPlot(scRNA_harmony,
           features=split(top5$gene,top5$cluster),
           cols=c("#ffffff","firebrick3"))
p

#将数据存储在dotplot_data中
library(tidyr)
library(dplyr)
library(pheatmap)
dotplot_data<-p$data
heatmap_data<-dotplot_data%>%
  select(features.plot,id,avg.exp.scaled)%>%
  pivot_wider(names_from=features.plot,values_from=avg.exp.scaled)
# 将细胞分群ID列设置为行名
library(tidyverse)
heatmap_data=column_to_rownames(heatmap_data,var="id")
pheatmap(heatmap_data,
         cluster_rows=F,
         cluster_cols=F,
         color=colorRampPalette(c("#AECDD7","#c7e0ed","#6a0624"))(100),
         border_color = 'white')
pdf(file = './var/new_var/new_singlecell/top5gene.pdf',height = 6,width = 8)
pheatmap(heatmap_data,
         cellwidth=10,
         cellheight=20,
         cluster_rows=F,
         cluster_cols=F,
         color=colorRampPalette(c("#AECDD7","#c7e0ed","#C85D4D"))(100),
         border_color = 'white')
dev.off()

#重新整理细胞亚群的排列，倒序排列
p$data$feature.groups2<-factor(p$data$feature.groups,
                               levels = c('T cells','CAFs','B cells','TAM','Epithelial','T/NK cells'))
library(ggh4x)
library(RColorBrewer)
strip<-strip_themed(
  background_x=elem_list_rect(fill=brewer.pal(7,"Paired")))
p1<-p$data %>%
  ggplot(aes(x=features.plot,
             y=id))+
  geom_point(aes(size=pct.exp,
                 color=avg.exp.scaled))+
  facet_wrap2(~feature.groups2,
              scales="free_x",
              strip=strip,
              nrow=1)+
  theme_classic()+
  RotatedAxis()+
  theme(strip.text.x=element_text(size=8),
        axis.text.x=element_text(color="black",size=10,face = 'bold.italic',hjust = 1),
        axis.title=element_blank(),
        strip.background=element_rect(color="white"),
        axis.text.y=element_blank())+
  scale_color_gradient(low="#ffffff",
                       high="firebrick3",
                       name="Avg.exp")
p1

df<-data.frame(x=0,y=levels(scRNA_harmony),stringsAsFactors=F)
df$y<-factor(df$y,levels=df$y)

#通过shape选择自己想要的图例形状
p2<-ggplot(df,aes(x,y,color=factor(y)))+
  geom_point(size=7,shape=18,show.legend=F)+
  scale_color_manual(values=rev(brewer.pal(7,"Paired")))+
  theme_classic()+
  scale_x_continuous(expand=c(0,0))+
  theme(
    plot.margin=margin(r=0),
    axis.title=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=14,color = 'black',face = 'bold.italic'),
    axis.ticks=element_blank(),
    axis.line=element_blank()
  )
p2
library(cowplot)
p3 = plot_grid(p2,p1,align="h",axis="bt",rel_widths=c(1.5,9))
p3
ggsave(plot = p3,filename = './singlecell2/top5marker.pdf',height = 6,width = 9)
ggsave(plot = p3,filename = './singlecell2/top5marker.png',height = 8,width = 11,dpi = 600)
#saveRDS(brac.markers,'./singlecell2/brac.markers.rds')
#########统计不同分组细胞含量
Idents(scRNA_harmony) <- 'celltype'
table(Idents(scRNA_harmony),scRNA_harmony$group)
library(reshape2)
library(ggplot2)
library(ggalluvial)
library(ggh4x)
Cellratio <- prop.table(table(Idents(scRNA_harmony), scRNA_harmony$group), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c('Cell','Group','value')
colors=c('gold',"#DC143C","#2E8B57","#00008B","#00BFFF",'purple')###"#FFB6C1","#7FFFAA",
p <- ggplot(Cellratio, aes(x = Group, y = value, fill = Cell,
                           stratum = Cell, alluvium = Cell)) +
  geom_col(position = 'stack', width = 0.6) +
  geom_stratum(width = 0.6, color = 'white') +
  geom_alluvium(alpha = 0.4, width = 0.6, color = 'white', linewidth = 1, curve_type = "linear") +
  scale_fill_manual(values = colors) +
  xlab('') + 
  ylab('') +
  scale_y_continuous(expand = c(0, 0))+
  theme_bw(base_size = 12) + 
  theme(
    axis.text = element_text(color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    legend.text = element_text(size = 13.4,face = 'bold.italic'),
    #legend.title = element_text(size = 14,face = 'bold.italic'),
    legend.title = element_blank(),
    axis.text.x = element_text(size = 14,face = 'bold.italic'),
    axis.text.y = element_text(size = 12,face = 'bold.italic')
  )
p

p2 <- DimPlot(scRNA_harmony, reduction = 'umap', label = F, group.by = 'group') +
  labs(title = "") +
  theme(plot.title = element_blank()) +
  theme(
    axis.text = element_text(color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    legend.text = element_text(size = 13,face = 'bold.italic'),
    legend.title = element_text(size = 14,face = 'bold.italic'),
    axis.text.x = element_text(size = 12,face = 'bold.italic'),
    axis.text.y = element_text(size = 12,face = 'bold.italic'),
    #  axis.text.x = element_blank(),
    #  axis.text.y = element_blank(),
    axis.title.x = element_text(size = 15,face = 'bold.italic'),
    axis.title.y = element_text(size = 15,face = 'bold.italic')
  )
p2

p3 <- p + p2
p3

ggsave(plot = p3,filename = './singlecell2/细胞含量.pdf',height = 6,width = 15)
ggsave(plot = p3,filename = './singlecell2/细胞含量.png',height = 6,width = 15,dpi = 600)


##########
library(Nebulosa)
gene <- c('DLST','EGFL7','DOCK4')
gene <- c('DLST','FBXO31')
p = plot_density(scRNA_harmony,gene)
p
p = plot_density(scRNA_harmony,'DLST')
ggsave(plot = p,filename = './var/new_var/new_singlecell/DLST.pdf',height = 6,width = 8)
p = plot_density(scRNA_harmony,'EGFL7')
ggsave(plot = p,filename = './var/new_var/new_singlecell/EGFL7.pdf',height = 6,width = 8)
p = plot_density(scRNA_harmony,'DOCK4')
ggsave(plot = p,filename = './var/new_var/new_singlecell/DOCK4.pdf',height = 6,width = 8)
library(patchwork)
p/umap2

ggsave(plot = p,filename = './singlecell/gene.pdf',height = 10,width = 12)
ggsave(plot = p,filename = './singlecell/gene.png',height = 10,width = 12,dpi = 300)


library(reshape2)
library(dplyr)
library(tibble)
Idents(scRNA_harmony) <- 'celltype'
AverageExpression <- AverageExpression(scRNA_harmony,assays = 'RNA',features = gene)
Ave_df <- melt(AverageExpression)
Ave_df <- Ave_df[,-7]
colnames(Ave_df) <- AverageExpression$RNA@Dimnames[[2]]
Ave_df_melt <- Ave_df %>% rownames_to_column(var = 'Gene') %>% melt(id.vars = 'Gene')
colnames(Ave_df_melt) <- c('Gene','Celltype','Expression')
Ave_df_melt$Group <- paste(Ave_df_melt$Celltype,Ave_df_melt$Gene,sep = '_')

######
library(ggplot2)
pb1 <- ggplot(Ave_df_melt,aes(x = Gene,y = Expression)) +
  geom_hline(yintercept = seq(0,2,2.5),linetype = 2,color = 'lightgrey',size = 1) + ######添加水平线
  geom_line() + #添加折线
  geom_segment(aes(x = Gene,xend = Gene,y = 0,yend = Expression),color = 'lightgrey',size = 1.5) + ##添加线段
  geom_point(size = 3,aes(color = Gene)) + 
  scale_color_manual(values = c('#BB8ABE','#EA9C9D','#B9E5FA')) +
  theme_bw() + 
  theme(panel.grid = element_blank()) + ###移除面板网络
  labs(x = '',y = 'Ave_Exp')
pb1
pb1 <- pb1 + facet_wrap(~Celltype, ncol = length(unique(Ave_df_melt$Celltype)))
pb1
# 设置分面标签的字体大小和背景颜色
#pb1 <- pb1 + scale_y_continuous(position = 'right')+ #更改y轴位置
#  theme(axis.text.y = element_text(size = 12,colour = 'black')) +
#  theme(axis.text.x = element_blank()) +
 # theme(axis.title.y = element_text(size = 12,colour = 'black')) +
#  theme(axis.title.x = element_blank()) +
#  theme(legend.position = 'right',      #设置图例
#        panel.border = element_blank(), #######移除面板
 #       axis.ticks.x = element_line(colour = NA)) ###移除x轴刻度线
pb1
ggsave(plot = pb1,filename = './var/new_var/new_singlecell/Ave_exp.pdf',height = 4,width = 8)
ggsave(plot = pb1,filename = './singlecell2/DLST_FBXO31.png',height = 9,width = 11,dpi = 600)

levels(Idents(scRNA_harmony))
scRNA <- subset(scRNA_harmony,idents = c('T/NK cells','T cells','Epithelial','TAM','B cells','CAFs'))
obj_markers <- FetchData(scRNA,vars = gene,layer = 'data')
#colmean <- colMeans(obj_markers)
obj_markers <- obj_markers %>%
  rowwise() %>%
  mutate(Group = if_else(sum(c_across(everything()) > 0) > 0, 'High', 'Low'))
table(obj_markers$Group)
scRNA <- AddMetaData(scRNA,metadata = obj_markers)
Idents(scRNA) <- scRNA@meta.data$Group
table(Idents(scRNA))
DEG_sc <- FindMarkers(scRNA,ident.1 = 'High',ident.2 = 'Low',
                      verbose = F,
                      min.pct = 0.1,
                      test.use = 'wilcox')
write.csv(DEG_sc,file = './singlecell2/DEG_sc.csv',row.names = T)
DEG_sc <- read.csv('./var/new_var/new_singlecell/DEG_sc.csv',row.names = 1)
log2FC <- 1
padj <- 0.05
DEG_sc1 <- DEG_sc %>% 
  mutate(threshold = case_when(
    avg_log2FC > log2FC & p_val_adj < padj ~ 'up',
    avg_log2FC < -log2FC & p_val_adj < padj ~ 'down',
    TRUE ~ 'ns'
  ))
table(DEG_sc1$threshold)
DEG_sc1$symbol <- rownames(DEG_sc1)
DEG_sc1 <- DEG_sc1[-c(1,2,3),]
write.csv(DEG_sc1,file = './var/new_var/new_singlecell/deg_sc.csv')
#DEG_sc1 <- DEG_sc1 %>% filter(DEG_sc1$symbol != 'IGLC3')
DEG_sc1 <- DEG_sc1 %>% filter(p_val != 1)
library(ggplot2)
library(ggrepel)
pdf(file = './singlecell2/scvolcano.pdf',height = 6,width = 9)
ggplot(DEG_sc1, aes(x = avg_log2FC, y = -log10(p_val), colour=threshold)) +
  geom_point(alpha=1, size=4.5) +
  scale_color_manual(values=c("#B9E5FA", "#d2dae2","#EA9C9D"))+
  # 辅助线
  geom_vline(xintercept=c(-0.585,0.585),lty=4,col="grey",lwd=0.8) +
  geom_hline(yintercept = -log10(padj),
             lty=4,col="grey",lwd=0.8) +
  # 坐标轴
  labs(x="",
       y="")+
  theme_bw()+
  ggtitle("")+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )
dev.off()
Symbol <- DEG_sc %>% filter(threshold != 'ns') %>% rownames()

library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
#library(ggstatsplot)
organism = 'hsa' 
OrgDb = 'org.Hs.eg.db'
need_DEG <- DEG_sc1
need_DEG$SYMBOL <- rownames(need_DEG)
need_DEG <- need_DEG[,c(2,8)] 
df <- bitr(rownames(need_DEG), 
           fromType = "SYMBOL",
           toType =  "ENTREZID",
           OrgDb = OrgDb)
need_DEG <- merge(need_DEG, df, by='SYMBOL')
geneList <- need_DEG$avg_log2FC
names(geneList) <- need_DEG$ENTREZID
geneList <- sort(geneList, decreasing = T)   #从大到小排序
KEGG_kk_entrez <- gseKEGG(geneList = geneList,
                          organism = organism,
                          pvalueCutoff = 1,
                          minGSSize = 5,
                          maxGSSize = 500)
KEGG_kk <- DOSE::setReadable(KEGG_kk_entrez, 
                             OrgDb=OrgDb,
                             keyType='ENTREZID')
gsea_result <- KEGG_kk@result
gsea_result = gsea_result %>% arrange(desc(setSize))
write.csv(gsea_result,file = './var/new_var/new_singlecell/gsea_result_singlecell.csv')
gsea_result <- gsea_result %>% filter(p.adjust < 0.05)
dotplot(KEGG_kk)+scale_color_continuous_c4a_seq('rd_pu') + mytheme
mytheme<- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11),
                 plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5))

pdf(file = './var/new_var/new_singlecell/scgseakegg.pdf',height = 10,width = 6)
dotplot(KEGG_kk)+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
dev.off()

###########T细胞亚群注释
levels(Idents(scRNA_harmony)) #打出来细胞类型供复制
scRNA = subset(scRNA_harmony,idents = c('CAFs')) #要修改的
table(Idents(scRNA))
all.gene <- rownames(scRNA)
scRNA = NormalizeData(scRNA,normalization.method="LogNormalize",scale.factor=10000) %>%
  FindVariableFeatures(selection.method = "vst",nfeatures=2000) %>% ScaleData(features = all.gene) %>%
  RunPCA(verbose = T,npcs = 30)
ElbowPlot(scRNA)
scRNA <- FindNeighbors(scRNA,dims = 1:10) %>% FindClusters(resolution = 0.1)
scRNA <- RunUMAP(scRNA,dims = 1:15)
mycolors<-c('#E64A35','#4DBBD4','#01A187','#3C5588','#F29F80',
             '#8491B6','#91D0C1','#7F5F48','#AF9E85','#4F4FFF','#CE3D33',
             '#739B57','#EFE685','#446983','#BB6239','#5DB1DC','#7F2268','#6BD66B','#800202','#D8D8CD','pink'
)
DimPlot(scRNA,reduction="umap",label=T,cols=mycolors)
library(SingleR)
library(celldex)
library(dplyr)
library(stringr)
library(pheatmap)
library(ReactomeGSA)
library(ggplot2)
#library(singleseqgset)
library(devtools)
levels(Idents(scRNA))
Idents(scRNA) <- 'seurat_clusters'
scRNA<-JoinLayers(scRNA)
brac.markers<-FindAllMarkers(scRNA,only.pos=T,min.pct=0.25,logfc.threshold=0.25,verbose=FALSE)
top5=brac.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
g = unique(top5$gene)
p<-DotPlot(scRNA,
           features=split(top5$gene,top5$cluster),
           cols=c("#ffffff","firebrick3"))
p
write.csv(top5,file = './var/new_var/new_singlecell/top5CAFs.csv')
#将数据存储在dotplot_data中
library(tidyr)
library(dplyr)
library(pheatmap)
dotplot_data<-p$data
heatmap_data<-dotplot_data%>%
  select(features.plot,id,avg.exp.scaled)%>%
  pivot_wider(names_from=features.plot,values_from=avg.exp.scaled)
# 将细胞分群ID列设置为行名
library(tidyverse)
heatmap_data=column_to_rownames(heatmap_data,var="id")
pheatmap(heatmap_data,
         cluster_rows=F,
         cluster_cols=F,
         color=colorRampPalette(c("#AECDD7","#c7e0ed","#6a0624"))(100),
         border_color = 'white')
pdf(file = './var/new_var/new_singlecell/CAFsMarkers.pdf',height = 6,width = 8)
pheatmap(heatmap_data,
         cellwidth=10,
         cellheight=20,
         cluster_rows=F,
         cluster_cols=F,
         color=colorRampPalette(c("#AECDD7","#c7e0ed","#C85D4D"))(100),
         border_color = 'white')
dev.off()


##############
markers <- c(
  'CD27','CD28','CCR7','SELL','TCF7', #### Naive T
  'GZMK','CXCR4','CXCR3','CD44',### CD8 Tem 
  'PRF1','GZMA', 'CCL5', 'GPR183', ######## CD8 Tcm
  'PTGER2', 'ICAM2', 'ANXA1','GNLY',############## CD4 Tcm
  'IFNG','RUNX3','EOMES',############## CD4 Tem
  'FOXP3','IL2RA','IL10RA','IKZF2','RTKN2', 'COC25B','S1PR4' ,### Normal Treg
  'IL23R','IL17A',###TH17
  'CXCL13',###########TH1
  'IL10','CTLA4',
  'LAG3' #Exhausted T
)

scRNA <- readRDS('./singlecell2/scRNA_亚群.rds')
Idents(scRNA) <- 'seurat_clusters'
p <- DotPlot(scRNA, features = markers, cols = "RdYlBu") +  RotatedAxis() +  
  xlab(NULL) + ylab(NULL) +
  theme(strip.text = element_text(angle = 50,size = 10,face = 'bold.italic'),
        axis.text.y = element_text(size = 10,face = 'bold.italic'),        
        axis.text.x = element_text(angle = 45,size = 10,face = 'bold.italic'),        
        legend.text = element_text(size = 9,face = 'bold.italic'),        
        legend.title = element_text(size = 9,face = 'bold.italic'))
p
ggsave(plot = p,filename = './singlecell2/markers2.pdf',height = 6,width = 8)
ggsave(plot = p,filename = './singlecell2/markers2.png',height = 8,width = 10,dpi = 600)

source("./custom_seurat_functions.R")
library(data.table)
library(stringr)
writeLines(paste0(0:5,','))
celltype = read.table('./var/new_var/new_singlecell/celltype2.txt',sep = ',')
celltype

#重新规定
scRNA$celltype <- plyr::mapvalues(scRNA$seurat_clusters,
                                          from = 0:4,
                                          to = celltype$V2
)
umap = scRNA@reductions$umap@cell.embeddings %>%    #坐标信息
  as.data.frame() %>% 
  cbind(celltype = scRNA@meta.data$celltype) # 注释后的cell_type信息
head(umap)
#计算每个细胞类型的median 坐标位置，留作在图上加标签
cell_type_med <- umap %>%
  group_by(celltype) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2)
  )
p <- ggplot(umap,aes(x= umap_1 , y = umap_2 ,color = celltype)) +  
  geom_point(size = 1 , alpha =1 )  + 
  scale_color_manual(values = allcolour)
p1 <- p  +
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=20), #设置legend标签的大小
    legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5)))  + #设置legend中 点的大小 
  #缩放坐标轴
  #添加圈圈
  #ggunchull::stat_unchull(alpha = 0.1, size = 0.5,lty = 2)  +
  #添加注释 repel防止标签重叠，把注释加到图片的点上去
  ggrepel:: geom_label_repel(aes(label=celltype), fontface="bold",data = cell_type_med,
                             point.padding=unit(0.3, "lines")) +
  # #去掉legend
  ggrepel::geom_label_repel(aes(label=celltype), fontface="bold",data = cell_type_med,
                            point.padding=unit(0.5, "lines")) +
  theme(legend.position = "none") +
  #去掉网格
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03))
p1
ggsave(plot = p1,filename = './var/new_var/new_singlecell/CAFs.pdf',height = 6,width = 8)



umap2 <- DimPlot(scRNA, reduction = "umap",group.by = 'group',
                 label = T,label.box = T,label.size = 4,repel = T) + theme_bw() +
  labs( x= "Umap 1",y= "Umap 2",title = NULL) +
  theme(panel.grid=element_blank(), # 去网格线
        plot.title = element_text(size = 15,color="black",hjust = 0.5),
        axis.text.x = element_text(size = 13, color = 'black',face = 'bold.italic'),
        axis.text.y = element_text(size = 13, color = 'black',face = 'bold.italic'),
        axis.title.x = element_text(size = 14, color = 'black',face = 'bold.italic'),
        axis.title.y = element_text(size = 14, color = 'black',face = 'bold.italic'),
        axis.ticks = element_line(color = 'black', lineend = 'round'),
        legend.position = 'bottom',
        legend.text = element_text(size = 13, color = 'black',face = 'bold.italic'),
        legend.title = element_text(size = 13, color = 'black',face = 'bold.italic'),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  scale_fill_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.65)) +
  scale_color_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.65))


p2 = DimPlot(scRNA, reduction = "umap",group.by = "celltype",
             label = T,label.box = T,label.size = 4,repel = T) + theme_bw() +
  labs( x= "Umap 1",y= "Umap 2",title = NULL) +
  theme(panel.grid=element_blank(), # 去网格线
        plot.title = element_text(size = 15,color="black",hjust = 0.5),
        axis.text.x = element_text(size = 13, color = 'black',face = 'bold.italic'),
        axis.text.y = element_text(size = 13, color = 'black',face = 'bold.italic'),
        axis.title.x = element_text(size = 14, color = 'black',face = 'bold.italic'),
        axis.title.y = element_text(size = 14, color = 'black',face = 'bold.italic'),
        axis.ticks = element_line(color = 'black', lineend = 'round'),
        legend.position = 'bottom',
        legend.text = element_text(size = 13, color = 'black',,face = 'bold.italic'),
        legend.title = element_text(size = 13, color = 'black',face = 'bold.italic'),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  scale_fill_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.65)) +
  scale_color_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.65))
p3 <- p2 + umap2
p3
ggsave(plot = p3,filename = './singlecell2/umap2.pdf',height = 6,width = 8)
ggsave(plot = p3,filename = './singlecell2/umap2.png',height = 8,width = 12,dpi = 600)

##############拟时序
library(Biobase)
library(monocle)
library(data.table)
library(ggsci)
library(reshape2)
library(magrittr)

Idents(scRNA) <- scRNA$celltype
expr <- GetAssayData(object = scRNA,layer = 'counts',assay = 'RNA')
sample_sheet <- scRNA@meta.data
sample_sheet <- cbind(rownames(sample_sheet),sample_sheet)
colnames(sample_sheet)[1] = 'cell'
gene_annotation = data.frame(gene_short_name = rownames(scRNA))
rownames(gene_annotation) = rownames(scRNA)
pd<-new("AnnotatedDataFrame",data=sample_sheet)
fd<-new("AnnotatedDataFrame",data=gene_annotation)
cds<-newCellDataSet(expr,phenoData=pd,featureData=fd,expressionFamily=negbinomial.size())
cds

###########################
#估计尺度因子
cds<-estimateSizeFactors(cds)
#估计离散度
cds<-estimateDispersions(cds)
####过滤
cds <- detectGenes(cds,min_expr = 0.1)
expressed_genes<-row.names(subset(fData(cds),num_cells_expressed>=10))

#计算Cluster间的marker基因
diff_test_res<-differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~celltype",cores = 6)
head(diff_test_res)
write.csv(diff_test_res,file = './singlecell2/degForcellOrdering.csv')

##############
diff_test_res <- diff_test_res[order(diff_test_res$qval),]
ordering_gene = row.names(diff_test_res[1:1000,])
#diff_test_res <- row.names(subset(diff_test_res$qval < -.-1))
cds = setOrderingFilter(cds = cds,ordering_genes = ordering_gene)
p = plot_ordering_genes(cds)
p
ggsave(plot = p,filename = './singlecell2/ordering_gene.pdf',height = 10,width = 8)
ggsave(plot = p,filename = './singlecell2/ordering_gene.png',height = 10,width = 8,dpi = 300)

#############降维、绘制细胞轨迹
cds <- reduceDimension(cds,method = 'DDRTree')
cds <- orderCells(cds = cds)

save(cds,file = './var/new_var/new_singlecell/OrderingCell_cds_done.Rdata')

p <- plot_cell_trajectory(cds,color_by = 'celltype') + theme(text = element_text(size = 18))
p

p1 <- plot_cell_trajectory(cds,color_by="Pseudotime")+theme(legend.position="right")
p1


GM_state <- function(x){
  if (length(unique(pData(cds)$State)) > 1) {
    T0_counts <- table(pData(cds)$State,pData(cds)$celltype)[,'eCAFs']
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  } else {return(1)}
}
cds <- orderCells(cds,root_state = GM_state(cds))
#plot_cell_trajectory(cds = cds,color_by = 'celltype')

mycolors<-c('#E64A35','#4DBBD4','#01A187','#3C5588','#F29F80',
            '#8491B6','#91D0C1','#7F5F48','#AF9E85','#4F4FFF','#CE3D33',
            '#739B57','#EFE685','#446983','#BB6239','#5DB1DC','#7F2268','#6BD66B','#800202','#D8D8CD','pink'
)
colors <- c('#6BD66B','#EFE685','#8491B6','#800202','#E64A35','#91D0C1')

p1<-plot_cell_trajectory(cds = cds,color_by="celltype") +
    #theme_bw() +
    theme(legend.position="right")+scale_color_manual(values=colors,name="celltype") +
    theme(axis.title.x = element_text(size = 14,face = 'bold.italic'),
          axis.title.y = element_text(size = 14,face = 'bold.italic'),
          axis.text.x = element_text(size = 13,face = 'bold.italic'),
          axis.text.y = element_text(size = 13,face = 'bold.italic'),
          legend.title = element_text(size = 14,face = 'bold.italic'),
          legend.text = element_text(size = 13,face = 'bold.italic'))
p1
ggsave(filename="./newCytoTRACE/monocle.pdf",width=8,height=6,plot=p1)

##########CytoTRACE
#devtools::install_local("./CytoTRACE_0.3.3.tar.gz")
library(CytoTRACE)

expr <- as.matrix(GetAssayData(scRNA,assay = 'RNA',layer = 'counts'))
####过滤掉在少于5个细胞表达的基因
expr <- expr[apply(expr>0,1,sum) >= 5,]
results <- CytoTRACE(expr,ncores = 1)
phenot <- scRNA$celltype
phenot <- as.character(phenot)
names(phenot) <- rownames(scRNA@meta.data)
####提取umap降维信息
emb <- scRNA@reductions[['umap']]@cell.embeddings
####
gene <- c('DLST','FBXO31')
gene <- c('DLST','EGFL7','DOCK4')
plotCytoTRACE(results,phenotype=phenot,emb = emb,gene = ,outputDir = 'newCytoTRACE')

############## CellChat
#devtools::install_github('sqjin/CellChat')
library(CellChat)
library(Seurat)
library(patchwork)
library(tidyverse)
library(dplyr)
library(magrittr)
scRNA_harmonY <- readRDS('./singlecell2/scRNA_harmony.rds')

expr <- GetAssayData(scRNA_harmonY,assay = 'RNA',layer = 'counts')
meta <- scRNA_harmonY@meta.data[,c('celltype','group')]
meta$celltype %<>% as.vector(.)

sc.cellchat <- createCellChat(object = expr)
sc.cellchat <- addMeta(sc.cellchat,meta = meta)
sc.cellchat <- setIdent(sc.cellchat,ident.use = 'celltype')

levels(sc.cellchat@idents)
groupSize <- as.numeric(table(sc.cellchat@idents))
                                  
####配受体数据库
dplyr::glimpse(CellChatDB.human$interaction)
sc.cellchat@DB <- CellChatDB.human
sc.cellchat <- subsetData(sc.cellchat,features = NULL)
future::plan('multisession',workers = 4)
options(future.globals.maxSize = 2000 * 1024^2)
sc.cellchat <- identifyOverExpressedGenes(sc.cellchat)####寻找高表达基因
sc.cellchat <- identifyOverExpressedInteractions(sc.cellchat)###########寻找其相互作用
sc.cellchat <- projectData(sc.cellchat,PPI.human)###将基因投射到PPI中
###########计算信号通路的配受体通讯概率
sc.cellchat <- computeCommunProb(sc.cellchat,raw.use = F,type = "triMean")
##########过滤小于10个细胞的胞间通讯网络
sc.cellchat <- filterCommunication(sc.cellchat,min.cells = 10)
#########汇总所有相关配体和受体，计算通路水平上的通信概率
sc.cellchat <- computeCommunProbPathway(sc.cellchat)
sc.cellchat <- aggregateNet(sc.cellchat) ####计算聚合网络
##############'netP'表示推断的信号通路的细胞间通讯网络
sc.cellchat <- netAnalysis_computeCentrality(sc.cellchat,slot.name = 'netP')
####################细胞通讯结果
group1.net <- subsetCommunication(sc.cellchat)
write.csv(group1.net,file = './singlecell2/cellchat.csv')


#########互作网络
par(mfrow = c(1,2))
netVisual_circle(sc.cellchat@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge = F,
                 title.name = '')
pdf(file = './cellchat/细胞互作网络.pdf',height = 6,width = 6)
netVisual_circle(
  sc.cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = '',
  vertex.label.cex = 1.5  # 增大节点标签的字体大小
)
dev.off()


############互作强度

mat <- sc.cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), label.edge= F, title.name = rownames(mat)[i])
}

mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
mat2[4,] <- mat[4,]
netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.size = .5,arrow.width = 1,
                 edge.weight.max = max(mat), label.edge= F, title.name = rownames(mat)[4])

#############
sc.cellchat@netP$pathways
levels(sc.cellchat@idents)
vertex.receiver = c(1,4,5) 
pathways.show <- c("CD45") 
par(mfrow=c(1,1))
pdf('./cellchat/CXCL and CD45.pdf',height = 6,width = 8)
netVisual_aggregate(sc.cellchat, signaling = pathways.show,layout = 'hierarchy',vertex.receiver = vertex.receiver)

dev.off()
par(mfrow=c(1,1))
netVisual_heatmap(sc.cellchat, signaling = pathways.show, color.heatmap = "Reds")



#############计算配受体对信号通路的贡献
pdf('./cellchat/配受体信号通路贡献.pdf',height = 6,width = 8)
netAnalysis_contribution(sc.cellchat,signaling = pathways.show)
dev.off()
#############查看细胞配体和其他细胞的受体
#netVisual_bubble(sc.cellchat,sources.use = 'T cells',remove.isolate = F,font.size = 14)
pathways.show <- c("CXCL")####CXCL 
sc.cellchat <- netAnalysis_computeCentrality(sc.cellchat, slot.name = "netP")
################
pdf('./cellchat/net热图.pdf',height = 6,width = 8)
netAnalysis_signalingRole_network(sc.cellchat, signaling = pathways.show, 
                                  width = 15, height = 6, font.size = 10)
dev.off()

############cellchat聚类
library(NMF)
library(ggalluvial)
selectK(sc.cellchat,pattern = 'incoming')
nPatterns = 2
sc.cellchat <- identifyCommunicationPatterns(sc.cellchat, pattern = "incoming", k = nPatterns,width = 6,height = 12,font.size = 10)
pdf('./cellchat/incoming聚类.pdf',height = 9,width = 8)
netAnalysis_river(sc.cellchat, pattern = "incoming",font.size = 3,cutoff = 0.5,do.order = F,slot.name = 'netP')
dev.off()

library(uwot)
sc.cellchat <- computeNetSimilarity(sc.cellchat, type = "functional")
sc.cellchat <- netEmbedding(sc.cellchat, type = "functional",umap.method = 'uwot')
sc.cellchat <- netClustering(sc.cellchat, type = "functional",do.parallel = F)
pdf('./cellchat/功能相似性识别信号基团.pdf',height = 8,width = 8)
netVisual_embedding(sc.cellchat, type = "functional", label.size = 3.5)
dev.off()


sc.cellchat <- computeNetSimilarity(sc.cellchat, type = "structural")
sc.cellchat <- netEmbedding(sc.cellchat, type = "structural",umap.method = 'uwot')
sc.cellchat <- netClustering(sc.cellchat, type = "structural",do.parallel = F)
pdf('./cellchat/结构相似性识别信号基团.pdf',height = 8,width = 8)
netVisual_embedding(sc.cellchat, type = "structural", label.size = 3.5)
dev.off()

saveRDS(sc.cellchat,'./singlecell2/cellchat.rds')





#bpe.se = BlueprintEncodeData()
#hpca.se = HumanPrimaryCellAtlasData()

#unique(hpca.se$label.main)
#unique(hpca.se$label.fine)
#unique(bpe.se$label.main)
#unique(bpe.se$label.fine)

#anno<-SingleR(scRNA@assays$RNA$data,
#                ref=list(BP=bpe.se,HPCA=hpca.se),
#                labels=list(bpe.se$label.fine,hpca.se$label.main),
#                clusters=scRNA@meta.data$seurat_clusters
#)

#plotScoreHeatmap(anno,clusters=anno@rownames,show_colnames=T)
#table(anno$labels)

#Celltype=data.frame(ClusterID=rownames(anno),
#                      celltype=anno$labels,
#                      stringsAsFactors=F)








