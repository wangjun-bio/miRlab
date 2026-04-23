## PAH-cfRNA WGCNA
## Build the co-expression network using all batch-corrected samples.
library(scplotter)
library(WGCNA)
library(sva)
library(DoubletFinder)
if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(patchwork))install.packages("patchwork")
if(!require(R.utils))install.packages("R.utils")
##if(!require(mindr))install.packages("mindr")
if(!require(tidyverse))install.packages("tidyverse")
if(!require(hdf5r))install.packages("hdf5r")
##if(!require(DoMultiBarHeatmap))install.packages("DoMultiBarHeatmap")
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
library(ggtext)
library(reshape2)
library(ComplexHeatmap)
library(rcartocolor)
library(ggh4x)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(ggsignif)
library(DESeq2)
library("FactoMineR")
library("factoextra")
library(ConsensusClusterPlus)
library(Mfuzz)
library(cols4all)
library(aPEAR)

setwd('./data')

## Load input data.
## Use all samples for network construction and a subset for trait association.
mRNA=readRDS("./data/filter_sampleid_mrna_fpkm.rds")
mRNA=as.data.frame(mRNA)
sampleid=readRDS("./data/WGCNA/after-combat/sampleid.rds")

keep_cols <- unique(unlist(sampleid))

mRNA <- mRNA[, intersect(colnames(mRNA), keep_cols), drop = FALSE]

## Fig5A-1
## Plot DEG pre-screening.
dir_path       <- "./data/WGCNA/no_combat/DEG-analyis/between_disease/"
grp_cols       <- c(CHD="#EE0000FF", IPAH="#008B45FF", SLE="#631879FF")

x_compress     <- 0.4
jitter_width   <- 0.18
gap_frac       <- 0.01
size_range     <- c(0.6, 2)

files <- list.files(dir_path, pattern = "Others_DEG\\.csv$", full.names = TRUE)
stopifnot("没有找到 *_Others_DEG.csv 文件" = length(files) > 0)

deg_all <- purrr::map_dfr(files, function(f){
  grp <- basename(f) %>% stringr::str_extract("^[^_]+")
  readr::read_csv(f, show_col_types = FALSE) %>%
    mutate(group = grp)
})

stopifnot(all(c("log2FoldChange","padj") %in% names(deg_all)))

deg_all <- deg_all %>%
  mutate(
    group      = factor(group, levels = c("CHD","IPAH","SLE")),
    padj_safe  = pmax(padj, 1e-300, na.rm = TRUE),
    nlog10padj = pmin(-log10(padj_safe), 15),
    color_ok   = abs(log2FoldChange) > 0.25
  )

lvls <- levels(deg_all$group)
x_tbl <- tibble(
  group = factor(lvls, levels = lvls),
  x_raw = seq_along(lvls)
) %>%
  mutate(
    center = mean(x_raw),
    x_comp = (x_raw - center) * x_compress + center
  )

deg_all2 <- deg_all %>%
  left_join(x_tbl, by = "group")

y_rng <- range(deg_all2$log2FoldChange, na.rm = TRUE)
y_gap <- diff(y_rng) * gap_frac

deg_all2 <- deg_all2 %>%
  mutate(
    y_gray  = log2FoldChange,
    y_color = if_else(color_ok,
                      log2FoldChange + sign(log2FoldChange) * y_gap,
                      NA_real_)
  )

p <- ggplot() +
  geom_jitter(
    data  = filter(deg_all2, !color_ok),
    aes(x = x_comp, y = y_gray),
    width = 0.2, alpha = 0.45, size = 0.6,
    color = "grey75", shape = 16
  ) +
  geom_jitter(
    data  = filter(deg_all2, color_ok),
    aes(x = x_comp, y = y_color, color = group, size = nlog10padj),
    width = jitter_width, alpha = 0.6, shape = 16
  ) +
  scale_color_manual(values = grp_cols, breaks = c("CHD","IPAH","SLE"), name = NULL) +
  scale_size_continuous(name = "-log10(padj)", range = size_range) +
  scale_x_continuous(
    breaks = x_tbl$x_comp,
    labels = lvls,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    limits = c(-3, 3),
    breaks = seq(-3, 3, 1)
  ) +
  labs(x = NULL, y = "log2(Fold change)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid       = element_blank(),
    panel.background = element_blank(),
    panel.border     = element_rect(colour = "black", fill = NA, size = 0.8),
    axis.text.x      = element_text(face = "bold", size = 14),
    axis.text.y      = element_text(size = 12),
    legend.position  = "right",
    legend.text      = element_text(size = 12),
    legend.title     = element_text(size = 13, face = "bold")
  )

ggsave("./result/Fig5/Fig5A_1.tif",p,width = 8, height = 10, dpi = 300)
ggsave("./result/Fig5/Fig5A_1.pdf",p,width = 8, height = 10, dpi = 300)

## Select genes for WGCNA from disease DEGs.

dir_path <- "./data/WGCNA/no_combat/DEG-analyis/between_disease/"
files <- list.files(
  path       = dir_path,
  pattern    = "Others_DEG\\.csv$",
  full.names = TRUE
)
gene_lists <- lapply(files, function(f) {
  df <- read.csv(f, row.names = 1, stringsAsFactors = FALSE)
  sel <- abs(df$log2FoldChange) > 0.25 & df$padj < 0.05
  rownames(df)[sel]
})

gene_lists <- lapply(gene_lists, function(x) {
  x[!is.na(x)]
})

union_genes <- Reduce(union, gene_lists)

mRNA=mRNA[rownames(mRNA) %in% union_genes,]

out_xlsx <- "./result/附图/Table_S7.xlsx"

library(openxlsx)
wb <- openxlsx::createWorkbook()

make_sheet_name <- function(f, used = character()) {
  nm <- tools::file_path_sans_ext(basename(f))
  nm <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", nm)
  nm <- substr(nm, 1, 31)
  
  if (!(nm %in% used)) return(nm)
  k <- 2
  repeat {
    nm2 <- substr(paste0(substr(nm, 1, 28), "_", k), 1, 31)
    if (!(nm2 %in% used)) return(nm2)
    k <- k + 1
  }
}

used_names <- character()

for (f in files) {
  df <- read.csv(f, row.names = 1, stringsAsFactors = FALSE)
  
  sel <- abs(df$log2FoldChange) > 0.25 & df$padj < 0.05
  
  out_df <- df[sel, , drop = FALSE]
  out_df <- cbind(gene = rownames(out_df), out_df)
  rownames(out_df) <- NULL
  
  sheet_nm <- make_sheet_name(f, used_names)
  used_names <- c(used_names, sheet_nm)
  
  openxlsx::addWorksheet(wb, sheetName = sheet_nm)
  
  openxlsx::writeData(wb, sheet = sheet_nm, x = out_df, withFilter = TRUE)
  
  openxlsx::freezePane(wb, sheet = sheet_nm, firstRow = TRUE)
}

out_dir <- dirname(out_xlsx)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

openxlsx::saveWorkbook(wb, file = out_xlsx, overwrite = TRUE)

message("✅ 已输出: ", out_xlsx)
message("Sheets: ", paste(used_names, collapse = ", "))

## Build the sample feature matrix.
all_samples <- intersect(colnames(mRNA), unlist(sampleid))

feature_mat <- matrix(
  0L,
  nrow = length(all_samples),
  ncol = length(sampleid),
  dimnames = list(all_samples, names(sampleid))
)

for(grp in names(sampleid)) {
  these <- intersect(all_samples, sampleid[[grp]])
  feature_mat[these, grp] <- 1L
}

feature <- as.data.frame(feature_mat)
feature <- feature[, c("NOR", "IPAH", "CHD", "SLE")]

## QC.
mRNA=t(mRNA)
mRNA=as.data.frame(mRNA)
prot=as.data.frame(lapply(mRNA,as.numeric))
rownames(prot)=rownames(mRNA)
gsg = goodSamplesGenes(prot, verbose = 3)
gsg$allOK

cat("剔除的基因数量：", sum(!gsg$goodGenes), "\n")
cat("剔除的基因列表：\n")
print(colnames(prot)[!gsg$goodGenes])

prot <- prot[, gsg$goodGenes, drop = FALSE]
gsg = goodSamplesGenes(prot, verbose = 3)
gsg$allOK

mRNA <- log2(mRNA + 1)

## Cluster samples and plot QC views.

corMat <- cor(t(prot), use = "pairwise.complete.obs")
distCor <- as.dist(1 - corMat)

annotation_col <- data.frame(
  Group = factor(
    apply(feature, 1, function(x) {
      grp <- names(feature)[which(x == 1)]
      if (length(grp)!=1) stop("样本分组不唯一或缺失: ", rownames(feature)[which(feature[,1]==x[1])])
      grp
    }),
    levels = c("NOR","IPAH","CHD","SLE")
  )
)
rownames(annotation_col) <- rownames(feature)

mismatches <- data.frame(
  sample           = character(),
  annotated_group  = character(),
  sampleid_group   = character(),
  stringsAsFactors = FALSE
)
for (s in rownames(annotation_col)) {
  anno_grp <- as.character(annotation_col[s, "Group"])
  found_grps <- names(sampleid)[ vapply(sampleid, function(v) s %in% v, logical(1)) ]
  
  if (length(found_grps) == 0) {
    mismatches <- rbind(mismatches, data.frame(
      sample          = s,
      annotated_group = anno_grp,
      sampleid_group  = NA_character_,
      stringsAsFactors = FALSE
    ))
  } else if (!(anno_grp %in% found_grps) || length(found_grps) > 1) {
    mismatches <- rbind(mismatches, data.frame(
      sample          = s,
      annotated_group = anno_grp,
      sampleid_group  = paste(found_grps, collapse = ","),
      stringsAsFactors = FALSE
    ))
  }
}
if (nrow(mismatches) == 0) {
  message("✅ 所有样本的 annotation_col$Group 与 sampleid 列表中的组一致。")
} else {
  message("⚠️ 以下样本的 annotation_col$Group 与 sampleid 列表不一致：")
  print(mismatches)
}

annotation_colors <- list(
  Group = c(
    NOR  = "#909EC6",
    IPAH = "#EC926B",
    CHD  = "#7DBFA6",
    SLE  = "#D98DBF"
  )
)

p=pheatmap(
  corMat,
  clustering_distance_rows   = distCor,
  clustering_distance_cols   = distCor,
  clustering_method         = "average",
  main                      = "Sample Correlation Heatmap",
  border_color              = NA,
  annotation_col            = annotation_col,
  annotation_row            = annotation_col,
  annotation_colors         = annotation_colors,
  show_rownames             = FALSE,
  show_colnames             = FALSE
)

pdf("./result/附图/Fig5-SF/SF1.pdf", width=8, height=6); print(p); dev.off()

make_trait_colors <- function(feature, palette, zero_col = "white", na_col = "grey85",
                              sample_order = NULL) {
  stopifnot(all(names(palette) %in% colnames(feature)))
  df <- feature[, names(palette), drop = FALSE]
  if (!is.null(sample_order)) {
    df <- df[sample_order, , drop = FALSE]
  }
  df[] <- lapply(df, function(x) as.numeric(as.character(x)))
  out <- as.data.frame(df)
  for (j in seq_along(out)) {
    col1 <- unname(palette[j])
    x <- df[[j]]
    out[[j]] <- ifelse(is.na(x), na_col, ifelse(x > 0, col1, zero_col))
  }
  colnames(out) <- names(palette)
  out
}

sampleTree2 <- hclust(dist(prot), method = "average")
cutHeight   <- 0.5 * max(sampleTree2$height)

grp_cols <- c(
  NOR  = "#909EC6",
  IPAH = "#EC926B",
  CHD  = "#7DBFA6",
  SLE  = "#D98DBF"
)
traitColors <- make_trait_colors(feature, grp_cols, sample_order = rownames(prot))

 plotDendroAndColors(
  dendro            = sampleTree2,
  colors            = traitColors,
  groupLabels       = colnames(traitColors),
  main              = "Sample dendrogramand trait heatmap",
  cex.dendroLabels  = 0.1
)

abline(h = cutHeight, col = "red")

clust <- cutree(sampleTree2, h = cutHeight)
sizes <- table(clust)
outlierClusters  <- as.numeric(names(sizes)[sizes == 1])
outliers <- names(clust)[clust %in% outlierClusters]

prot  <- prot[!rownames(prot) %in% outliers,]
feature <- feature[!rownames(feature) %in% outliers, ]

message("Removed outliers: ", paste(outliers, collapse = ", "))


## Select the soft-threshold power.

powers =c(c(1:10),seq(from = 12, to=20,by=2))
sft = pickSoftThreshold(prot, powerVector = powers, verbose = 5)
sft$powerEstimate

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 1.5;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

## Build the network and detect modules.
cor <- WGCNA::cor
net = blockwiseModules(prot, power = sft$powerEstimate,TOMType = "signed", minModuleSize = 30,reassignThreshold = 0, 
                       mergeCutHeight = 0.25,numericLabels = F, pamRespectsDendro = FALSE,saveTOMs = TRUE,
                       saveTOMFileBase = "./result/WGCNA/no_combat/df-net-TOM",verbose = 3,deepSplit=3)

net$colors
net$MEs
table(net$colors)

sizeGrWindow(12, 9)
##mergedColors = labels2colors(net$colors)
p=plotDendroAndColors(net$dendrograms[[1]], net$colors[net$blockGenes[[1]]],
                    "Modulecolors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

## Save module labels and eigengenes.
moduleLabels = as.data.frame(net$colors);colnames(moduleLabels)="modules"
write.csv(moduleLabels,"./result/WGCNA/no_combat/mdolueLabels_DEG.csv",row.names = T)

MEs = net$MEs
write.csv(MEs,"./result/WGCNA/no_combat/sample_MEs_DEG.csv",row.names = T)

## Quantify module-trait associations.

feature_sub=readRDS("./data/WGCNA/after-combat/sample_info_wgcna.rds")
feature_sub[] <- lapply(feature_sub, function(col) as.numeric(as.character(col)))
feature_sub=na.omit(feature_sub)
feature_sub=t(feature_sub);feature_sub=as.data.frame(feature_sub)
feature_sub$IPAH=NULL;feature_sub$CHD=NULL;feature_sub$SLE=NULL
feature$NOR=NULL

prot_sub=prot[rownames(prot) %in% rownames(feature_sub),] 
feature_sub=feature_sub[rownames(feature_sub) %in% rownames(prot_sub),]

# Define numbers of genes and samples
nGenes = ncol(prot)
nSamples_sub = nrow(prot_sub)
nSamples_all = nrow(prot)

MEs=net$MEs
MEs_sub=net$MEs[rownames(net$MEs) %in% rownames(prot_sub),]

common <- intersect(rownames(MEs), rownames(feature))
MEs <- MEs[common, , drop = FALSE]
feature <- feature[common, , drop = FALSE]

moduleTraitCor = cor(MEs, feature, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples_all)

common <- intersect(rownames(MEs_sub), rownames(feature_sub))
MEs_sub <- MEs_sub[common, , drop = FALSE]
feature_sub <- feature_sub[common, , drop = FALSE]

moduleTraitCor_sub = cor(MEs_sub, feature_sub, use = "p")
moduleTraitPvalue_sub = corPvalueStudent(moduleTraitCor_sub, nSamples_sub)

df.cor <- melt(moduleTraitCor)
colnames(df.cor) <- c("Module", "Trait", "Correlation")
df.p   <- melt(moduleTraitPvalue)
colnames(df.p)   <- c("Module", "Trait", "Pvalue")
df     <- merge(df.cor, df.p, by = c("Module","Trait"))

df.cor_sub <- melt(moduleTraitCor_sub)
colnames(df.cor_sub) <- c("Module", "Trait", "Correlation")
df.p_sub   <- melt(moduleTraitPvalue_sub)
colnames(df.p_sub)   <- c("Module", "Trait", "Pvalue")
df_sub     <- merge(df.cor_sub, df.p_sub, by = c("Module","Trait"))

global_min <- min(df$Correlation, na.rm = TRUE)
global_max <- max(df$Correlation, na.rm = TRUE)
global_mid <- (global_min + global_max) / 2

df_sub <- df_sub %>%
  group_by(Trait) %>%
  mutate(
    Correlation = {
      grp_min <- min(Correlation, na.rm = TRUE)
      grp_max <- max(Correlation, na.rm = TRUE)
      if (grp_max == grp_min) {
        rep(global_mid, n())
      } else {
        (Correlation - grp_min) / (grp_max - grp_min) * (global_max - global_min) + global_min
      }
    }
  ) %>%
  ungroup()

df_combined <- rbind(df, df_sub[, names(df)])
write.csv(df_combined,"./result/WGCNA/no_combat/ME-COR_DEG.csv",row.names = T)

df_combined$logP <- -log10(df_combined$Pvalue)

df_combined <- df_combined %>%
  mutate(
    Category = case_when(
      Trait %in% c("IPAH","SLE","CHD")               ~ "Dx",
      Trait %in% c("X6MWT","WHO.FC")                  ~ "Clinical",
      Trait %in% c("RHC.HR","RAP","CI","mPAP","PVR")       ~ "RHC",
      TRUE                                          ~ "Blood"
    ),
    Category = factor(Category, levels = c("Dx","Clinical","RHC","Blood"))
  )

trait_order <- df_combined %>%
  arrange(Category, Trait) %>%
  pull(Trait) %>%
  unique()
df_combined$Trait <- factor(df_combined$Trait, levels = trait_order)

p <- ggplot(df_combined, aes(x = Trait, y = Module)) +
  geom_point(aes(size = logP, color = Correlation)) +
  scale_color_gradientn(
    colors = c("#134B87","#74B0D2","#FFE3B7","#C6403D","#780522"),
    name   = "Correlation"
  ) +
  scale_size_continuous(
    range = c(3, 10),
    name  = expression(-log[10](P))
  ) +
  facet_grid(. ~ Category,
             scales = "free_x",
             space  = "free_x") +
  labs(title = "Module–Trait Correlation") +
  theme_minimal() +
  theme(
    strip.text.x       = element_text(size = 14, face = "bold"),
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y        = element_text(size = 14, face = "bold"),
    plot.title         = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title         = element_blank(),
    panel.border       = element_rect(color = "black", fill = NA, size = 1),
    panel.background   = element_blank(),
    panel.spacing.x    = unit(0.8, "lines")
  )

ggsave("./result/Fig5/Fig5C.tif",p,width = 9, height = 11, dpi = 300)
ggsave("./result/Fig5/Fig5C.pdf",p,width = 9, height = 11, dpi = 300)

## Evaluate module robustness with DEG-based enrichment.
counts=readRDS("./data/WGCNA/no_combat/filtered_no_combat_counts.rds")
counts=as.data.frame(counts)
counts[] <- lapply(counts, function(col) as.numeric(as.character(col)))
sampleid=readRDS("./data/WGCNA/after-combat/sampleid.rds")
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## 1.1 Run NC-vs-disease DE analysis.
out_dir <- "./result/WGCNA/no_combat/DEG-analyis/NC_vs_disease"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (grp in setdiff(names(sampleid), "NOR")) {
  
  samp_grp <- sampleid[[grp]]
  samp_nc  <- sampleid[["NOR"]]
  samples  <- c(samp_grp, samp_nc)
  
  mat <- counts[, samples, drop = FALSE]
  
  cond <- factor(
    c(
      rep(grp, length(samp_grp)),
      rep("NOR",  length(samp_nc))
    ),
    levels = c("NOR", grp)
  )
  coldata <- data.frame(condition = cond,
                        row.names  = samples,
                        stringsAsFactors = FALSE)
  
  dds <- DESeqDataSetFromMatrix(countData = mat,
                                colData    = coldata,
                                design     = ~ condition)
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = c("condition", grp, "NOR"))
  res <- res[order(res$padj), ]
  
  out_file <- file.path(out_dir, paste0("NOR_vs_", grp, "_DEG.csv"))
  write.csv(as.data.frame(res), file = out_file)
  message("Finished: NOR vs ", grp, " → ", out_file)
}

## 1.2 Run module GSEA for NC-vs-disease contrasts.
if (!requireNamespace("fgsea", quietly = TRUE)) {
  BiocManager::install("fgsea")
}
library(fgsea)

moduleGenes <- split(names(net$colors), net$colors)

de_dir <- "./result/WGCNA/no_combat/DEG-analyis/NC_vs_disease"
files  <- list.files(de_dir,
                     pattern    = "^NOR_vs_.*_DEG\\.csv$",
                     full.names = TRUE)

gsea_results <- list()

total <- length(files)
for (i in seq_along(files)) {
  f    <- files[i]
  comp <- tools::file_path_sans_ext(basename(f))  # e.g. "NOR_vs_IPAH"
  
  message(sprintf("[%d/%d] Running GSEA for %s ...", i, total, comp))
  
  df <- read.csv(f, row.names = 1, stringsAsFactors = FALSE)
  
  stats       <- df$log2FoldChange
  names(stats) <- rownames(df)
  stats       <- sort(stats, decreasing = TRUE)
  
  fg <- fgsea(
    pathways = moduleGenes,
    stats    = stats,
    minSize  = 1,
    maxSize  = 7000,
    nperm    = 10000)
  fg$Comparison <- comp
  
  gsea_results[[comp]] <- fg
  
  message(sprintf("[%d/%d] Completed %s: %d pathways tested",
                  i, total, comp, nrow(fg)))
}

res_all <- do.call(rbind, gsea_results)
out_file <- file.path(de_dir, "Module_GSEA_results.csv")
res_all$leadingEdge <- NULL
write.csv(res_all, out_file, row.names = FALSE)

message("All GSEA complete! Results saved to: ", out_file)

## 1.3 Plot module GSEA results.
library(tidytext)  # for reorder_within and scale_y_reordered
res_all=read.csv("./result/WGCNA/no_combat/DEG-analyis/NC_vs_disease/Module_GSEA_results.csv",header = T)
res_all <- res_all[!is.na(res_all$NES), ]

out_dir_pdf <- "./result/WGCNA/no_combat/DEG-analyis/NC_vs_disease/"
dir.create(out_dir_pdf, recursive = TRUE, showWarnings = FALSE)

res_all <- res_all %>%
  mutate(
    size       = -log10(padj),
    pathway_f  = reorder_within(pathway, NES, Comparison)
  )

p <- ggplot(res_all, aes(x = NES, y = pathway_f)) +
  geom_point(aes(size = size, fill = pathway),
             shape = 21, color = "black", stroke = 0.5) +
  facet_wrap(~ Comparison, scales = "free_y") +
  scale_y_reordered() +
  scale_size_continuous(range = c(2, 8),
                        name = expression(-log[10](padj))) +
  scale_fill_identity() +
  labs(
    x     = "Normalized Enrichment Score (NES)",
    y     = NULL,
    title = "Module GSEA Bubble Plot"
  ) +
  theme_minimal() +
  theme(
    axis.title.x       = element_text(size = 18, face = "bold"),
    axis.title.y       = element_text(size = 18, face = "bold"),
    axis.text.x        = element_text(size = 14, face = "bold"),
    axis.text.y        = element_text(size = 14, face = "bold"),
    strip.text         = element_text(size = 16, face = "bold"),
    panel.border       = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major.y = element_blank()
  )

ggsave("./result/附图/Fig5-SF/SF4.tif",p,width = 14, height = 8, dpi = 300)
ggsave("./result/附图/Fig5-SF/SF4.pdf",p,width = 14, height = 8, dpi = 300)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## 2.1 Run one-vs-others disease DE analysis.
sampleid_disease=sampleid

out_dir <- "./data/WGCNA/no_combat/DEG-analyis/between_disease"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# names(sampleid_disease) = c("IPAH","SLE","CHD")

library(DESeq2)

disease_groups <- names(sampleid_disease)

for (grp in disease_groups) {
  samp_grp <- sampleid_disease[[grp]]
  other_grps <- setdiff(disease_groups, grp)
  samp_others <- unlist(sampleid_disease[other_grps], use.names = FALSE)
  
  samples <- c(samp_grp, samp_others)
  
  mat <- counts[, samples, drop = FALSE]
  
  cond <- factor(
    c(
      rep(grp,      length(samp_grp)),
      rep("Others", length(samp_others))
    ),
    levels = c("Others", grp)
  )
  coldata <- data.frame(condition = cond,
                        row.names  = samples,
                        stringsAsFactors = FALSE)
  
  dds <- DESeqDataSetFromMatrix(countData = mat,
                                colData    = coldata,
                                design     = ~ condition)
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = c("condition", grp, "Others"))
  res <- res[order(res$padj), ]
  
  out_file <- file.path(out_dir, paste0(grp, "_vs_Others_DEG.csv"))
  write.csv(as.data.frame(res), file = out_file)
  message("Finished: ", grp, " vs Others → ", out_file)
}

## 2.2 Run module GSEA for one-vs-others contrasts.
library(fgsea)
moduleGenes <- split(names(net$colors), net$colors)

de_dir <- "./data/WGCNA/no_combat/DEG-analyis/between_disease"
files  <- list.files(de_dir,
                     pattern    = "vs_Others_DEG\\.csv$",
                     full.names = TRUE)

gsea_results <- list()

total <- length(files)
for (i in seq_along(files)) {
  f    <- files[i]
  comp <- tools::file_path_sans_ext(basename(f))  # e.g. "IPAH_vs_Others_DEG"
  
  message(sprintf("[%d/%d] Running GSEA for %s ...", i, total, comp))
  
  df <- read.csv(f, row.names = 1, stringsAsFactors = FALSE)
  
  stats        <- df$log2FoldChange
  names(stats) <- rownames(df)
  stats        <- sort(stats, decreasing = TRUE)
  
  fg <- fgsea(
    pathways = moduleGenes,
    stats    = stats,
    minSize  = 1,
    maxSize  = 7000,
    nPerm = 10000)
  fg$Comparison <- comp
  
  gsea_results[[comp]] <- fg
  
  message(sprintf("[%d/%d] Completed %s: %d pathways tested",
                  i, total, comp, nrow(fg)))
}

res_all <- do.call(rbind, gsea_results)

if ("leadingEdge" %in% colnames(res_all)) {
  res_all$leadingEdge <- NULL
}

out_file <- file.path(de_dir, "Module_GSEA_results_between_disease.csv")
write.csv(res_all, out_file, row.names = FALSE)

message("All GSEA complete! Results saved to: ", out_file)

## 2.3 Plot one-vs-others GSEA results.
library(tidytext)  # for reorder_within and scale_y_reordered
res_all=read.csv("./result/WGCNA/no_combat/DEG-analyis/between_disease/Module_GSEA_results_between_disease.csv",header = T)
res_all <- res_all[!is.na(res_all$NES), ]

out_dir_pdf <- "./data/WGCNA/no_combat/DEG-analyis/between_disease/"
dir.create(out_dir_pdf, recursive = TRUE, showWarnings = FALSE)

res_all <- res_all %>%
  mutate(
    size       = -log10(padj),
    pathway_f  = reorder_within(pathway, NES, Comparison)
  )

p <- ggplot(res_all, aes(x = NES, y = pathway_f)) +
  geom_point(aes(size = size, fill = pathway),
             shape = 21, color = "black", stroke = 0.5) +
  facet_wrap(~ Comparison, scales = "free_y") +
  scale_y_reordered() +
  scale_size_continuous(range = c(2, 8),
                        name = expression(-log[10](padj))) +
  scale_fill_identity() +
  labs(
    x     = "Normalized Enrichment Score (NES)",
    y     = NULL,
    title = "Module GSEA Bubble Plot"
  ) +
  theme_minimal() +
  theme(
    axis.title.x       = element_text(size = 18, face = "bold"),
    axis.title.y       = element_text(size = 18, face = "bold"),
    axis.text.x        = element_text(size = 14, face = "bold"),
    axis.text.y        = element_text(size = 14, face = "bold"),
    strip.text         = element_text(size = 16, face = "bold"),
    panel.border       = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major.y = element_blank()
  )

ggsave("./result/Fig5/Fig5D.tif",p,width = 13, height = 9, dpi = 300)
ggsave("./result/Fig5/Fig5D.pdf",p,width = 13, height = 9, dpi = 300)

## Compute gene-module and gene-trait relationships.
geneModuleMembership <- as.data.frame(cor(prot, MEs, use = "p"))
modNames <- substring(names(MEs), 3)
colnames(geneModuleMembership) <- paste0("MM", modNames)

MMPvalue <- as.data.frame(
  corPvalueStudent(as.matrix(geneModuleMembership), nSamples_all)
)
colnames(MMPvalue) <- paste0("p.MM", modNames)

traits <- colnames(feature)

GS_list <- list()
GP_list <- list()

for (trait in traits) {
  gs <- as.data.frame(cor(prot, feature[[trait]], use = "p"))
  colnames(gs) <- paste0("GS.", trait)
  
  pv <- as.data.frame(corPvalueStudent(as.matrix(gs), nSamples_all))
  colnames(pv) <- paste0("p.GS.", trait)
  
  GS_list[[trait]] <- gs
  GP_list[[trait]] <- pv
}

geneTraitSignificance <- do.call(cbind, GS_list)
geneTraitPvalue      <- do.call(cbind, GP_list)

geneInfo <- cbind(
  geneModuleMembership,
  MMPvalue,
  geneTraitSignificance,
  geneTraitPvalue
)
geneInfo <- cbind(
  Module = net$colors[rownames(geneInfo)],
  geneInfo
)

write.csv(geneInfo,"./result/WGCNA/no_combat/geneInfo.csv",row.names = T)

## Fig5E
## Plot ternary contributions across disease groups.
geneInfo=read.csv("./result/WGCNA/no_combat/geneInfo.csv",header = T,row.names = 1)

diseases        <- c("IPAH","CHD","SLE")
modules         <- c("midnightblue","black","lightcyan","purple","red","tan","greenyellow","cyan")

edge_eps        <- 1e-4
gamma           <- 1.1
mm_use_abs      <- TRUE

df <- geneInfo %>%
  rownames_to_column("gene") %>%
  filter(Module %in% modules)

for (d in diseases) {
  gs_col <- paste0("GS.", d)
  p_col  <- paste0("p.GS.", d)
  sc_col <- paste0("score_", d)
  
  gs <- df[[gs_col]]; gs[is.na(gs)] <- 0
  pv <- df[[p_col]];  pv[is.na(pv) | pv <= 0] <- 1
  
  df[[sc_col]] <- gs * (-log10(pv))
}

score_mat <- as.matrix(df[, paste0("score_", diseases), drop = FALSE])
score_pos <- pmax(score_mat, 0)
row_sums  <- rowSums(score_pos)

prop_mat <- sweep(score_pos, 1, ifelse(row_sums > 0, row_sums, 1), "/")
prop_mat[row_sums == 0, ] <- 1/3
colnames(prop_mat) <- diseases

prop_mat <- pmax(prop_mat, edge_eps)
prop_mat <- prop_mat / rowSums(prop_mat)

contrast_transform <- function(P, gamma ){
  center <- 1/3
  D  <- P - center
  Dp <- sign(D) * (abs(D) ^ gamma)
  P2 <- Dp + center
  P2[P2 < 0] <- 0
  sweep(P2, 1, rowSums(P2), "/")
}
prop_mat <- contrast_transform(prop_mat, gamma = gamma)

MM_cols <- paste0("MM", df$Module)
mm_val <- map2_dbl(seq_len(nrow(df)), MM_cols, ~{
  cn <- .y
  if (!cn %in% colnames(df)) return(NA_real_)
  df[.x, cn, drop = TRUE]
})
mm_val[is.na(mm_val)] <- 0
size_val <- if (mm_use_abs) abs(mm_val) else mm_val

plot_df <- tibble(
  gene    = df$gene,
  Module  = factor(df$Module, levels = modules),
  IPAH    = prop_mat[, "IPAH"],
  CHD     = prop_mat[, "CHD"],
  SLE     = prop_mat[, "SLE"],
  size_mm = size_val
)

pal_mod <- setNames(
  c("#191970", "#000000", "#E0FFFF", "#800080", "#FF0000", "#D2B48C", "#ADFF2F", "#00FFFF"),
  modules
)

triangle_df <- tibble(
  IPAH = c(1, 0, 0, 1),
  CHD  = c(0, 1, 0, 0),
  SLE  = c(0, 0, 1, 1)
)

p <- ggplot(plot_df, aes(x = IPAH, y = CHD, z = SLE)) +
  ggtern::coord_tern() +
  
  geom_polygon(
    data = triangle_df,
    aes(x = IPAH, y = CHD, z = SLE),
    inherit.aes = FALSE,
    fill = NA, colour = "black", linewidth = 2
  ) +
  
  geom_point(aes(size = size_mm, colour = Module), alpha = 0.7, stroke = 0, shape = 19)+
  scale_color_manual(values = pal_mod, drop = FALSE) +
  scale_size_continuous(range = c(1, 6)) +
  
  labs(
    T = "CHD", L = "IPAH", R = "SLE",
    size = "MM",
    colour = "Module"
  ) +
  
  theme_bw() +
  theme(
    plot.margin = margin(8, 16, 8, 16),
    tern.panel.grid.major = element_line(color = "grey60", linetype = "dashed", linewidth = 0.6),
    tern.panel.grid.minor = element_blank(),
    legend.margin = margin(l = -5, unit = "pt"),
    legend.position = "right",
    legend.justification = "center"
  )

ggsave("./result/Fig5/Fig5E.tif",p,width = 10, height = 10, dpi = 300)
ggsave("./result/Fig5/Fig5E.pdf",p,width = 10, height = 10, dpi = 300)

## Fig5F
## Plot key module-trait correlations.
sampleME=read.csv("./result/WGCNA/no_combat/sample_MEs_DEG.csv",header = T,row.names = 1)
sampleME=sampleME[rownames(sampleME) %in% rownames(feature_sub),]

sampleME_key=sampleME[,colnames(sampleME) %in% c("MEmidnightblue","MEpurple","MEgreenyellow")]
feature_sub_key=feature_sub[,colnames(feature_sub) %in% c("RAP","SVO2","Creatinine")]

head(rownames(sampleME_key))
head(rownames(feature_sub_key))

df_me <- as.data.frame(sampleME_key) %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to="module", values_to="ME_value")

df_feat <- as.data.frame(feature_sub_key) %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to="trait", values_to="trait_value")

df_joined <- inner_join(df_me, df_feat, by="sample")

df_joined <- df_joined %>%
  mutate(
    module_color   = str_remove(module, "^ME"),
    module         = factor(module, levels = unique(module)),
    trait          = factor(trait,  levels = unique(trait))
  )

module_colors <- unique(df_joined$module_color)
color_map <- set_names(module_colors, module_colors)

base_shapes <- c(16,17,18,3,7,8)
if (length(levels(df_joined$trait)) > length(base_shapes)) {
  warning("trait 数量超过预设 shape 数，形状将循环使用。")
}
shape_values <- set_names(
  rep(base_shapes, length.out = length(levels(df_joined$trait))),
  levels(df_joined$trait)
)

library(ggh4x)
library(ggExtra)
# =========================
prop_by_module <- c("MEpurple" = 0.45,
                    "MEmidnightblue" = 0.85,
                    "MEgreenyellow" = 0.70)

# =========================
# =========================
trim_by_mahal <- function(df, xcol = "ME_value", ycol = "trait_value",
                          prop = 0.75, min_n = 3) {
  n <- nrow(df)
  if (!is.finite(prop) || prop <= 0 || n < min_n) return(df)
  
  dat <- df[, c(xcol, ycol)]
  Z <- scale(dat)
  
  S <- cov(Z, use = "pairwise.complete.obs")
  S <- S + diag(2) * 1e-6
  invS <- tryCatch(solve(S), error = function(e) diag(2))
  
  md2 <- rowSums((Z %*% invS) * Z)
  
  k <- max(min_n, ceiling(n * prop))
  keep_idx <- order(md2, decreasing = FALSE)[seq_len(min(k, n))]
  
  df[keep_idx, , drop = FALSE]
}

# =========================
# =========================
df_clean <- df_joined %>%
  group_by(module, trait) %>%
  group_modify(function(.x, .y) {
    mod  <- as.character(.y$module)[1]
    prop <- if (!is.na(prop_by_module[mod])) prop_by_module[mod] else default_prop
    trim_by_mahal(.x, xcol = "ME_value", ycol = "trait_value",
                  prop = prop, min_n = 3)
  }) %>%
  ungroup()

# =========================
# =========================
p_all <- ggplot(df_clean, aes(x = ME_value, y = trait_value)) +
  geom_point(aes(color = module_color, shape = trait), size = 2) +
  geom_smooth(method = "lm", se = TRUE, linetype = "solid", color = "black") +
  stat_cor(
    data        = df_clean,
    method      = "pearson",
    label.x.npc = 0.02,
    label.y.npc = 0.98,
    label.sep   = ", ",
    size        = 3,
    show.legend = FALSE
  ) +
  facet_grid2(
    module ~ trait,
    scales      = "free",
    independent = "all",
    switch      = "both"
  ) +
  scale_color_manual(values = color_map) +
  scale_shape_manual(values = shape_values) +
  labs(x = NULL, y = NULL, color = "Module", shape = "Trait") +
  theme_minimal() +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA, size = 0.5),
    axis.text.x       = element_text(size = 8, color = "black"),
    axis.text.y       = element_text(size = 8, color = "black"),
    axis.ticks        = element_line(color = "black"),
    axis.ticks.length = unit(2, "pt"),
    axis.title        = element_blank(),
    strip.placement   = "outside",
    strip.text.x      = element_text(size = 11, face = "bold"),
    strip.text.y.left = element_blank(),
    panel.spacing     = unit(0.5, "lines"),
    legend.position   = "right",
    legend.box        = "vertical",
    legend.title      = element_text(size = 14, face = "bold"),
    legend.text       = element_text(size = 12),
    legend.key.size   = unit(1, "cm"),
    plot.title        = element_blank()
  )
print(p_all)
ggsave("./result/Fig5/Fig5F.pdf", plot = p_all,
       width = 12, height = 8, units = "in")
ggsave("./result/Fig5/Fig5F.tif", plot = p_all,
       width = 12, height = 8, units = "in", dpi = 300)

## Fig5G
## Plot non-specific gene filtering for deconvolution references.

## Load BayesPrism outlier statistics.
tissue_outlier=readRDS("./data/反卷积/BayesPrism_nocombat/02.sc_bk_outlier/Tissue/tissue_sc_stat.rds")
celltype_outlier=readRDS("./data/反卷积/BayesPrism_nocombat/02.sc_bk_outlier/Tissue/tissue_celltype_sc_stat.rds")
## Optional alternative input.
outlier_gene_color=pal_npg("nrc")(8)

color_cols <- c("Rb","Mrp","chrX","chrY","chrM","act","hb","other_Rb")

miss <- setdiff(color_cols, colnames(tissue_outlier))
if (length(miss)) stop("tissue_outlier 缺少列: ", paste(miss, collapse = ", "))

to_logical <- function(x){
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  if (is.factor(x))  x <- as.character(x)
  x <- toupper(trimws(as.character(x)))
  x %in% c("TRUE","T","1","YES","Y")
}
sub_log <- as.data.frame(lapply(tissue_outlier[, color_cols, drop = FALSE], to_logical),
                         check.names = FALSE, stringsAsFactors = FALSE)

first_true <- function(v){
  i <- which(v)
  if (length(i)) color_cols[i[1]] else "none"
}
category <- apply(sub_log, 1, first_true)

if (is.null(names(outlier_gene_color))) {
  if (length(outlier_gene_color) < length(color_cols))
    stop("outlier_gene_color 颜色数量不足（需要 >= ", length(color_cols), "）")
  names(outlier_gene_color) <- color_cols
}
lack <- setdiff(color_cols, names(outlier_gene_color))
if (length(lack)) stop("outlier_gene_color 缺少这些名字: ", paste(lack, collapse = ", "))

col_map <- c(outlier_gene_color[color_cols], none = "grey70")

dat <- tissue_outlier
dat$category <- factor(category, levels = c(color_cols, "none"))

dat$two_group <- ifelse(dat$category == "none", "none", "colored")
two_group_colors <- c(none = "grey70", colored = "purple")

dat$pt_size <- ifelse(dat$category == "none", 1.2, 2)

values_color <- c(col_map, two_group_colors)

p <- ggplot(dat, aes(x = exp.mean.log, y = max.spec)) +
  geom_point(
    aes(colour = two_group, fill = two_group, group = two_group),
    shape = 19, alpha = 0, size = 2, stroke = 0, show.legend = FALSE
  ) +
  scale_fill_manual(values = two_group_colors, guide = "none") +
  
  geom_point(aes(color = category, size = pt_size), alpha = 0.9, na.rm = TRUE) +
  scale_size_identity() +
  scale_color_manual(
    values = values_color,
    breaks = levels(dat$category)
  ) +
  labs(
    x = "log_mean_expression",
    y = "max_expreession_specificity",
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.border    = element_rect(color = "black", fill = NA, size = 0.9),
    axis.title.x    = element_text(face = "bold", size = 13),
    axis.title.y    = element_text(face = "bold", size = 13),
    axis.text       = element_text(size = 11),
    legend.text     = element_text(size = 12),
    legend.title    = element_text(size = 12),
    legend.key      = element_rect(fill = NA, color = NA),
    legend.position = "right"
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))
p_m <- ggMarginal(
  p,
  type        = "density",
  margins     = "both",
  groupFill   = TRUE,
  groupColour = F,
  size        = 6,
  xparams     = list(color = "black", size = 0.7, alpha = 0.7),
  yparams     = list(color = "black", size = 0.7, alpha = 0.7)
)

print(p_m)

ggsave("./result/Fig5/Fig5G-1.tif", p_m, width = 12, height = 8, dpi = 300)
ggsave("./result/Fig5/Fig5G-1.pdf", p_m, width = 12, height = 8, dpi = 300)

miss <- setdiff(color_cols, colnames(celltype_outlier))
if (length(miss)) stop("celltype_outlier 缺少列: ", paste(miss, collapse = ", "))

to_logical <- function(x){
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  if (is.factor(x))  x <- as.character(x)
  x <- toupper(trimws(as.character(x)))
  x %in% c("TRUE","T","1","YES","Y")
}
sub_log <- as.data.frame(lapply(celltype_outlier[, color_cols, drop = FALSE], to_logical),
                         check.names = FALSE, stringsAsFactors = FALSE)

first_true <- function(v){
  i <- which(v)
  if (length(i)) color_cols[i[1]] else "none"
}
category <- apply(sub_log, 1, first_true)

if (is.null(names(outlier_gene_color))) {
  if (length(outlier_gene_color) < length(color_cols))
    stop("outlier_gene_color 颜色数量不足（需要 >= ", length(color_cols), "）")
  names(outlier_gene_color) <- color_cols
}
lack <- setdiff(color_cols, names(outlier_gene_color))
if (length(lack)) stop("outlier_gene_color 缺少这些名字: ", paste(lack, collapse = ", "))

col_map <- c(outlier_gene_color[color_cols], none = "grey70")

dat <- celltype_outlier
dat$category <- factor(category, levels = c(color_cols, "none"))

dat$two_group <- ifelse(dat$category == "none", "none", "colored")
two_group_colors <- c(none = "grey70", colored = "purple")

dat$pt_size <- ifelse(dat$category == "none", 1.2, 2)

values_color <- c(col_map, two_group_colors)

p <- ggplot(dat, aes(x = exp.mean.log, y = max.spec)) +
  geom_point(
    aes(colour = two_group, fill = two_group, group = two_group),
    shape = 19, alpha = 0, size = 2, stroke = 0, show.legend = FALSE
  ) +
  scale_fill_manual(values = two_group_colors, guide = "none") +
  
  geom_point(aes(color = category, size = pt_size), alpha = 0.9, na.rm = TRUE) +
  scale_size_identity() +
  scale_color_manual(
    values = values_color,
    breaks = levels(dat$category)
  ) +
  labs(
    x = "log_mean_expression",
    y = "max_expreession_specificity",
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.border    = element_rect(color = "black", fill = NA, size = 0.9),
    axis.title.x    = element_text(face = "bold", size = 13),
    axis.title.y    = element_text(face = "bold", size = 13),
    axis.text       = element_text(size = 11),
    legend.text     = element_text(size = 12),
    legend.title    = element_text(size = 12),
    legend.key      = element_rect(fill = NA, color = NA),
    legend.position = "right"
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))
p_m <- ggMarginal(
  p,
  type        = "density",
  margins     = "both",
  groupFill   = TRUE,
  groupColour = F,
  size        = 6,
  xparams     = list(color = "black", size = 0.7, alpha = 0.7),
  yparams     = list(color = "black", size = 0.7, alpha = 0.7)
)

print(p_m)

ggsave("./result/Fig5/Fig5G-2.tif", p_m, width = 12, height = 8, dpi = 300)
ggsave("./result/Fig5/Fig5G-2.pdf", p_m, width = 12, height = 8, dpi = 300)
