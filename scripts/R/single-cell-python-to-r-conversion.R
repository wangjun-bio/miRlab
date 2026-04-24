library(scplotter)
library(DoubletFinder)
if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(patchwork))install.packages("patchwork")
if(!require(R.utils))install.packages("R.utils")
if(!require(tidyverse))install.packages("tidyverse")
if(!require(hdf5r))install.packages("hdf5r")
if(!require(ggplot2))install.packages("ggplot2")
library(ggpubr)
library(viridis)
library(patchwork)
library(harmony)
library(scales)
library(devtools)
library(presto)
library(ggsci)
library(monocle3)
library(ggrepel)
library(tidydr)
library(reshape2)
library(ComplexHeatmap)
library(rcartocolor)
library(ggh4x)
library(cowplot)
library(ggsignif)
library(DoubletFinder)
library("FactoMineR")
library(hexbin)
library(RColorBrewer)
library(GSVA)
library(AUCell)
library(SeuratDisk)
library("factoextra")

# Set project paths
setwd("./data")
RESULT_DIR <- "../result"
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

# Convert h5ad to h5Seurat
Convert(
  file.path(".", "filtered_subset_for_convert.h5ad"),
  dest       = "h5seurat",
  overwrite  = TRUE,
  metadata   = TRUE,
  assays     = list(RNA = list(counts = "raw_counts", data = "X")),
  reductions = list(pca = "X_pca", scvi = "X_scvi", umap = "X_umap")
)

# Load converted Seurat object and metadata
seu <- LoadH5Seurat(file.path(".", "filtered_subset_for_convert.h5seurat"), meta.data = FALSE, misc = FALSE)
META=read.csv(file.path(".", "adata_obs_metadata.csv"), header = T, row.names = 1)

all(rownames(seu@meta.data) %in% rownames(META))

META2 <- META[rownames(seu@meta.data), , drop = FALSE]

seu <- AddMetaData(
  object   = seu,
  metadata = META2
)

# Subset key tissues
Idents(seu)=seu$organ_tissue

target_tissues <- c("Vasculature", "Blood", "Heart", "Skin")
cells_to_keep <- seu@meta.data$organ_tissue %in% target_tissues
sc_data <- subset(seu, cells = colnames(seu)[cells_to_keep])
print(paste("原始数据细胞数:", ncol(seu)))
print(paste("子集数据细胞数:", ncol(sc_data)))
print("子集中包含的组织类型:")
print(table(sc_data@meta.data$organ_tissue))

# Merge original annotations into simplified cell-type labels
maps <- tibble::tribble(
  ~original,                                                      ~merged,
  "Blood_cd4-positive, alpha-beta memory t cell",               "Blood_Lymphoid_cell",
  "Blood_cd8-positive, alpha-beta cytokine secreting effector t cell", "Blood_Lymphoid_cell",
  "Blood_cd8-positive, alpha-beta t cell",                      "Blood_Lymphoid_cell",
  "Blood_naive thymus-derived cd4-positive, alpha-beta t cell", "Blood_Lymphoid_cell",
  "Blood_cd4-positive, alpha-beta t cell",                      "Blood_Lymphoid_cell",
  "Blood_naive b cell",                                         "Blood_Lymphoid_cell",
  "Blood_memory b cell",                                        "Blood_Lymphoid_cell",
  "Blood_plasma cell",                                          "Blood_Lymphoid_cell",
  "Blood_nk cell",                                              "Blood_Lymphoid_cell",
  "Blood_type i nk t cell",                                     "Blood_Lymphoid_cell",
  "Blood_classical monocyte",                                   "Blood_Myeloid_cell",
  "Blood_monocyte",                                             "Blood_Myeloid_cell",
  "Blood_macrophage",                                           "Blood_Myeloid_cell",
  "Blood_neutrophil",                                           "Blood_Myeloid_cell",
  "Blood_erythrocyte",                                          "Blood_Erythrocyte",
  "Blood_platelet",                                             "Blood_Platelet",

  "Skin_t cell",                                                "Skin_Lymphoid_cell",
  "Skin_cd8-positive, alpha-beta memory t cell",                "Skin_Lymphoid_cell",
  "Skin_cd8-positive, alpha-beta cytotoxic t cell",             "Skin_Lymphoid_cell",
  "Skin_cd4-positive, alpha-beta memory t cell",                "Skin_Lymphoid_cell",
  "Skin_regulatory t cell",                                     "Skin_Lymphoid_cell",
  "Skin_nk cell",                                               "Skin_Lymphoid_cell",
  "Skin_nkt cell",                                              "Skin_Lymphoid_cell",
  "Skin_macrophage",                                            "Skin_Myeloid_cell",
  "Skin_cd1c-positive myeloid dendritic cell",                  "Skin_Myeloid_cell",
  "Skin_mast cell",                                             "Skin_Myeloid_cell",
  "Skin_stromal cell",                                          "Skin_Stromal_cell",
  "Skin_muscle cell",                                           "Skin_Muscle_cell",
  "Skin_endothelial cell",                                      "Skin_Endothelial_cell",
  "Skin_epithelial cell",                                       "Skin_Epithelial_cell",

  "Heart_cardiac endothelial cell",                             "Heart_Endothelial_cell",
  "Heart_cardiac muscle cell",                                 "Heart_Muscle_cell",
  "Heart_smooth muscle cell",                                   "Heart_Muscle_cell",
  "Heart_fibroblast of cardiac tissue",                         "Heart_Fibroblast",
  "Heart_hepatocyte",                                           "Heart_Hepatocyte",

  "Vasculature_t cell",                                       "Vasculature_Lymphoid_cell",
  "Vasculature_nk cell",                                      "Vasculature_Lymphoid_cell",
  "Vasculature_macrophage",                                   "Vasculature_Myeloid_cell",
  "Vasculature_mast cell",                                    "Vasculature_Myeloid_cell",
  "Vasculature_fibroblast",                                   "Vasculature_Fibroblast",
  "Vasculature_smooth muscle cell",                           "Vasculature_Muscle_cell",
  "Vasculature_pericyte cell",                                "Vasculature_Pericyte_cell",
  "Vasculature_artery endothelial cell",                       "Vasculature_Endothelial_cell",
  "Vasculature_endothelial cell",                              "Vasculature_Endothelial_cell"
)

# Add integrated cell-type labels
meta <- sc_data@meta.data
meta$celltype_integration <- maps$merged[match(meta$tissue_celltype, maps$original)]

sc_data@meta.data <- meta
Idents(sc_data)=sc_data$celltype_integration

# Re-cluster and visualize the integrated single-cell subset
sc_data <- FindVariableFeatures(sc_data, selection.method = "vst", nfeatures = 4000)
gc()
sc_data <- ScaleData(sc_data,verbose = T)
gc()
sc_data <- RunPCA(sc_data,features = VariableFeatures(object = sc_data),verbose = T)
sc_data=RunHarmony(sc_data,"donor", plot_convergence = F)
sc_data <- FindNeighbors(sc_data, dims = 1:30,reduction = "harmony")
sc_data <- FindClusters(sc_data, resolution = 1)
sc_data=RunUMAP(sc_data,reduction = "harmony", dims = 1:30)
Idents(sc_data)=sc_data$celltype_integration

# Save re-clustered object
dir.create(file.path(RESULT_DIR, "GSE201333_RAW/GSM6058681_TabulaSapiens.h5ad"), recursive = TRUE, showWarnings = FALSE)
saveRDS(sc_data, file.path(RESULT_DIR, "GSE201333_RAW/GSM6058681_TabulaSapiens.h5ad/recluster_key_celltype_tissue.rds"))
