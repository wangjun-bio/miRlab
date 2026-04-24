library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(BayesPrism)
library(arrow)
library(tidyr)
library(pheatmap)
library(reshape2)
library(scales)
library(RColorBrewer)
library(rstatix)
library(ggpubr)
library(patchwork)
library(BayesPrism)
library(ggalluvial)

# Environment setup
setwd("./data")
RESULT_DIR <- "../result"
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

# Upstream BayesPrism input loading
bk.dat=read.csv("./counts_no_combat_all_Deconvolution.csv",header = T)
tissue.labels<- readLines("./key_4tissue_celltype_label.txt")
sc.dat=read_feather("./key_celltype_4tissue_counts_Deconvolution.feather")

# Input preprocessing
sc.dat=as.data.frame(sc.dat)
rownames(sc.dat)=sc.dat$cell_id;sc.dat$cell_id=NULL
rownames(bk.dat)=bk.dat$X;bk.dat$X=NULL
colnames(bk.dat) <- gsub("\\.", "-", colnames(bk.dat))
colnames(sc.dat) <- gsub("\\.", "-", colnames(sc.dat))

n_common <- length(intersect(colnames(sc.dat), colnames(bk.dat)))
n_diff <- length(setdiff(colnames(sc.dat), colnames(bk.dat))) +
  length(setdiff(colnames(bk.dat), colnames(sc.dat)))

cat("Matched gene columns: ", n_common, "\n")
cat("Unmatched gene columns: ", n_diff, "\n")

sort(table(tissue.labels))

# Reference correlation QC
plot.cor.phi (input=sc.dat,
              input.labels=tissue.labels,
              title="cell state correlation",
              pdf.prefix=file.path(RESULT_DIR, "key_celltype_cor"),
              cexRow=1.2, cexCol=1.2,
              margins=c(2,2))

# Outlier-gene QC
sc.stat <- plot.scRNA.outlier(
  input=sc.dat,
  cell.type.labels=tissue.labels,
  species="hs",
  return.raw=TRUE,
  pdf.prefix=file.path(RESULT_DIR, "sc_data_key_CELLTYPE_outlier") )
saveRDS(sc.stat, file.path(RESULT_DIR, "key_celltype_sc_stat.rds"))

# Gene filtering
sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs",
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                  exp.cells=5)

dim(sc.dat)
dim(sc.dat.filtered)

plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = bk.dat,
                 pdf.prefix=file.path(RESULT_DIR, "sc_bulk_key_celltype_cor") )

sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered,
                                         gene.type = "protein_coding")

# Marker-gene subset for deconvolution
markers <- readLines("./4tissue_celltype_integration_marker_top100.txt")

common_genes <- intersect(markers, colnames(sc.dat.filtered))
length(common_genes)

sc.dat.pc.celltype <- sc.dat.filtered[, common_genes, drop = FALSE]

# BayesPrism model construction
myPrism <- new.prism(
  reference=sc.dat.pc.celltype,
  mixture=bk.dat,
  input.type="count.matrix",
  cell.type.labels = tissue.labels,
  cell.state.labels = NULL,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1,
)
saveRDS(myPrism, file.path(RESULT_DIR, "4tissue_celltype_Integration_marker_bayes.rds"))

# BayesPrism inference
bp.res <- run.prism(prism = myPrism, n.cores=100)
saveRDS(bp.res, file.path(RESULT_DIR, "result_bayes_4tissue_celltype_marker_100.rds"))

# Downstream analysis configuration
config <- list(
  sampleid_file    = "./WGCNA/after-combat/sampleid.rds",
  sampleMEs_file   = "./WGCNA/no_combat/sample_MEs_DEG.csv",
  bayes_res_file   = "./BayesPrism_nocombat/03.bayes_result/Celltype_Integration/result_bayes_4tissue_celltype_marker_100_new.rds",

  tissue_color     = pal_igv("default")(51)[c(26:45)],
  group_color      = c(NOR  = "#909EC6",IPAH = "#EC926B",CHD  = "#7DBFA6",SLE  = "#D98DBF"),
  bubble_colors    = c("#134B87","#74B0D2","#FFE3B7","#C6403D","#780522"),

  output_dirs = list(
    barplot          = file.path(RESULT_DIR, "Fig6"),
    boxplot          = file.path(RESULT_DIR, "Fig6"),
    heatmap          = file.path(RESULT_DIR, "Fig6"),
    pca              = file.path(RESULT_DIR, "Fig6"),
    bubble           = file.path(RESULT_DIR, "Fig6"),
    scatter_cor      = file.path(RESULT_DIR, "Fig6")
  )
)

for (dir in config$output_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

# Load downstream inputs
sampleid  <- readRDS(config$sampleid_file)
sampleMEs <- read.csv(config$sampleMEs_file, header=TRUE, row.names=1)
bp.res    <- readRDS(config$bayes_res_file)

theta_tissue   <- get.fraction(bp=bp.res, which.theta="final", state.or.type="type")
theta_celltype <- get.fraction(bp=bp.res, which.theta="final", state.or.type="state")
theta_cv       <- bp.res@posterior.theta_f@theta.cv

common_samples <- intersect(rownames(sampleMEs), rownames(theta_tissue))
theta_subset   <- theta_tissue[common_samples, , drop=FALSE] %>% as.data.frame()

# Fig6-SF4-1: sample-level stacked bar plot
sample_group <- data.frame(
  sample = rownames(theta_tissue),
  group  = sapply(rownames(theta_tissue), function(s) {
    grp <- names(sampleid)[vapply(sampleid, function(v) s %in% v, logical(1))]
    if(length(grp)==0) NA_character_ else grp
  }),
  stringsAsFactors = FALSE
)

library(tibble)
df_long <- theta_tissue %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(
    cols      = -sample,
    names_to  = "tissue",
    values_to = "fraction"
  ) %>%
  left_join(sample_group, by="sample")
df_long$group <- factor(df_long$group, levels = c("NOR","IPAH","CHD","SLE"))
df_long$tissue <- factor(df_long$tissue, levels=unique(df_long$tissue), ordered=TRUE)

dir.create(file.path(RESULT_DIR, "supplementary"), recursive = TRUE, showWarnings = FALSE)
write.csv(df_long, file.path(RESULT_DIR, "supplementary", "Table_S9.csv"), row.names = T)

p <- ggplot(df_long, aes(x=sample, y=fraction, fill=tissue)) +
  geom_col(width=1) +
  scale_fill_manual(values=config$tissue_color) +
  facet_grid(~ group, scales="free_x", space="free_x") +
  scale_y_continuous(expand=c(0,0), limits=c(0,1)) +
  theme_bw() +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.title.y     = element_text(size=14, face="bold"),
    panel.spacing    = unit(0.2, "lines"),
    strip.background = element_blank(),
    strip.text       = element_text(face="bold", size=14),
    legend.position  = "bottom",
    legend.title     = element_text(size=14),
    legend.text      = element_text(size=12)
  ) +
  labs(x=NULL, y="Estimated Fraction", fill="Tissue Type")

dir.create(config$output_dirs$barplot, recursive=TRUE, showWarnings=FALSE)
ggsave(
  filename = file.path(config$output_dirs$barplot, "Fig6-SF4-1.pdf"),
  plot     = p, width=12, height=8, units="in"
)
ggsave(
  filename = file.path(config$output_dirs$barplot, "Fig6-SF4-1.tif"),
  plot     = p, width=12, height=8, units="in", device="tiff", dpi=300
)

group_levels  <- c("NOR","IPAH","CHD","SLE")
tissue_levels <- levels(df_long$tissue) %||% unique(df_long$tissue)

# Fig6F-1: group-level mean fraction bar plot
df_group_mean <- df_long %>%
  filter(!is.na(group)) %>%
  group_by(group, tissue) %>%
  dplyr::summarize(mean_fraction = mean(fraction, na.rm = TRUE), .groups = "drop") %>%
  complete(group = unique(df_long$group), tissue = tissue_levels, fill = list(mean_fraction = 0)) %>%
  mutate(
    group  = factor(group,  levels = group_levels),
    tissue = factor(tissue, levels = tissue_levels, ordered = TRUE)
  ) %>%
  group_by(group) %>%
  mutate(mean_fraction = if (sum(mean_fraction, na.rm = TRUE) > 0)
    mean_fraction / sum(mean_fraction, na.rm = TRUE) else 0) %>%
  ungroup()

keep <- c("Lung","Heart","Blood","Vasculature","Skin")

df_group_mean_subset <- df_group_mean %>%
  mutate(
    tissue = as.character(tissue),
    tissue = ifelse(tissue %in% keep, tissue, "Other_tissue")
  ) %>%
  group_by(group, tissue) %>%
  dplyr::summarise(mean_fraction = sum(mean_fraction, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    group  = factor(group,  levels = levels(df_group_mean$group)),
    tissue = factor(tissue, levels = c(keep, "Other_tissue"), ordered = TRUE)
  ) %>%
  group_by(group) %>%
  mutate(mean_fraction = mean_fraction / sum(mean_fraction, na.rm = TRUE)) %>%
  ungroup()
tissue_color=config$tissue_color[c(1:24)]
names(tissue_color)=tissue_levels
if (!"Other_tissue" %in% names(tissue_color)) {
  tissue_color["Other_tissue"] <- "#BDBDBD"
}

nor_order <- df_group_mean %>%
  dplyr::filter(group == "NOR") %>%
  dplyr::arrange(dplyr::desc(mean_fraction)) %>%
  dplyr::pull(tissue)

df_group_mean <- df_group_mean %>%
  dplyr::mutate(
    tissue = factor(tissue, levels = rev(nor_order)),
    group  = forcats::fct_rev(group)
  )

p_group <- ggplot(df_group_mean, aes(x = group, y = mean_fraction, fill = tissue)) +
  geom_col(width = 0.95, color = NA, show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(
    values = tissue_color,
    drop   = FALSE,
    name   = "Tissue Type",
    breaks = rev(levels(df_group_mean$tissue))
  ) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  labs(x = NULL, y = NULL) +
  scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
  theme_bw() +
  theme(
    panel.grid         = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.title.y       = element_text(size = 14, face = "bold"),
    axis.text.x        = element_text(size = 12, face = "bold"),
    axis.text.y        = element_text(size = 12, face = "bold"),
    legend.position    = "bottom",
    legend.title       = element_text(size = 14),
    legend.text        = element_text(size = 12)
  )

legend_df <- data.frame(
  group         = factor(rep(levels(df_group_mean$group)[1],
                             length(levels(df_group_mean$tissue))),
                         levels = levels(df_group_mean$group)),
  mean_fraction = 0,
  tissue        = factor(levels(df_group_mean$tissue),
                         levels = rev(levels(df_group_mean$tissue)))
)

p_group <- p_group +
  geom_point(
    data = legend_df,
    aes(x = group, y = mean_fraction, fill = tissue),
    inherit.aes = FALSE,
    shape = 21, size = 4, alpha = 0, show.legend = TRUE
  ) +
  guides(
    fill = guide_legend(
      nrow = 3, byrow = TRUE,
      override.aes = list(shape = 21, size = 4, alpha = 1, colour = "grey20")
    )
  )

dir.create(config$output_dirs$barplot, recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = file.path(config$output_dirs$barplot, "Fig6F-1.pdf"),
  plot     = p_group, width = 14, height = 6, units = "in"
)
ggsave(
  filename = file.path(config$output_dirs$barplot, "Fig6F-1.tif"),
  plot     = p_group, width = 14, height = 6, units = "in",
  device   = "tiff", dpi = 300, compression = "lzw"
)

# Fig6-SF4-2: fraction boxplot across groups
new_order <- df_long %>%
  group_by(tissue) %>%
  dplyr::summarize(mean_frac = mean(fraction, na.rm = TRUE), .groups="drop") %>%
  arrange(desc(mean_frac)) %>%
  pull(tissue)
df_long$tissue <- factor(df_long$tissue, levels=new_order)

p <- ggplot(df_long, aes(x = tissue, y = fraction, fill = group)) +
  geom_boxplot(
    position      = position_dodge(width = 0.8),
    width         = 0.8,
    outlier.shape = NA,
    size          = 1
  ) +
  geom_jitter(
    aes(color = group),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    size  = 1, alpha = 0.8
  )  +
  scale_fill_manual(name = "Disease Group", values = config$group_color) +
  scale_color_manual(name = "Disease Group", values = config$group_color) +
  theme_classic() +
  theme(
    axis.title.x    = element_text(size = 14, face = "bold"),
    axis.title.y    = element_text(size = 14, face = "bold"),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title    = element_text(size = 14, face = "bold"),
    legend.text     = element_text(size = 12)
  ) +
  labs(
    x = "Tissue Type",
    y = "Estimated Fraction"
  )

dir.create(config$output_dirs$boxplot, recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = file.path(config$output_dirs$boxplot, "Fig6-SF4-2.pdf"),
  plot     = p, width = 19, height = 8, units = "in"
)
ggsave(
  filename = file.path(config$output_dirs$boxplot, "Fig6-SF4-2.tif"),
  plot     = p, width = 19, height = 8, units = "in",
  device   = "tiff", dpi = 300
)

# Log2 fold-change summary versus NOR
nor_means <- df_long %>%
  filter(group == "NOR") %>%
  group_by(tissue) %>%
  dplyr::summarize(mean_NOR = mean(fraction, na.rm = TRUE), .groups = "drop")

fc_df <- df_long %>%
  group_by(tissue, group) %>%
  dplyr::summarize(mean_frac = mean(fraction, na.rm = TRUE), .groups = "drop") %>%
  filter(group != "NOR") %>%
  left_join(nor_means, by = "tissue") %>%
  mutate(log2FC = log2(mean_frac / mean_NOR))

df_stats_fc <- df_long %>%
  group_by(tissue) %>%
  wilcox_test(fraction ~ group, ref.group = "NOR") %>%
  add_significance("p") %>%
  dplyr::rename(group = group2) %>%
  select(tissue, group, p.signif)

fc_df <- fc_df %>%
  left_join(df_stats_fc, by = c("tissue", "group"))

tissue_order <- fc_df %>%
  group_by(tissue) %>%
  dplyr::summarize(mean_log2FC = mean(log2FC, na.rm = TRUE)) %>%
  arrange(desc(mean_log2FC)) %>%
  pull(tissue)
fc_df$tissue <- factor(fc_df$tissue, levels = tissue_order)

p <- ggplot(fc_df, aes(x = tissue, y = log2FC, fill = group)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width    = 0.8,
    color    = "NA"
  ) +
  geom_text(
    aes(y = log2FC, label = p.signif),
    position = position_dodge(width = 0.8),
    vjust    = -0.3,
    size     = 5
  ) +
  scale_fill_manual(name = "Disease Group", values = config$group_color) +
  theme_classic() +
  theme(
    axis.title.x     = element_text(size = 14, face = "bold"),
    axis.title.y     = element_text(size = 14, face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y      = element_text(size = 12),
    legend.position  = "bottom",
    legend.title     = element_text(size = 14, face = "bold"),
    legend.text      = element_text(size = 12),
    panel.border     = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(
    x = "Tissue Type",
    y = "log2 Fold Change vs NOR"
  )

ggsave(
  filename = file.path(config$output_dirs$boxplot, "log2FC_celltype_fraction.pdf"),
  plot     = p, width = 19, height = 8, units = "in"
)
ggsave(
  filename = file.path(config$output_dirs$boxplot, "log2FC_celltype_fraction.tif"),
  plot     = p, width = 19, height = 8, units = "in",
  device   = "tiff", dpi = 300
)

new_order <- df_long %>%
  group_by(tissue) %>%
  dplyr::summarize(mean_frac = mean(fraction, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_frac)) %>%
  pull(tissue)
df_long <- df_long %>%
  mutate(tissue = factor(tissue, levels = new_order))

# Fig6F-2: focused cell-type fraction boxplot
df_other <- df_long %>%
  mutate(
    tissue = if_else(
      tissue %in% c("Heart_Fibroblast", "Vasculature_Muscle_cell", "Blood_Erythrocyte", "Skin_Stromal_cell", "Skin_Myeloid_cell"),
      tissue,
      "Other_Tissue"
    )
  )
ordered_levels <- c("Heart_Fibroblast", "Vasculature_Muscle_cell", "Blood_Erythrocyte", "Skin_Stromal_cell", "Skin_Myeloid_cell")

df_other <- df_other %>%
  mutate(tissue = factor(tissue, levels = ordered_levels))

df_stats <- df_other %>%
  group_by(tissue) %>%
  wilcox_test(fraction ~ group, ref.group = "NOR") %>%
  add_significance("p") %>%
  add_xy_position(x = "group", dodge = 0.8) %>%
  dplyr::rename(test_group = group2)

p <- ggplot(df_other, aes(x = group, y = fraction, fill = group)) +

  geom_boxplot(
    position      = position_dodge(width = 0.8),
    width         = 0.8,
    outlier.shape = NA,
    size          = 1
  ) +
  geom_jitter(
    aes(color = group),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    size  = 1, alpha = 0.8
  ) +

  geom_text(
    data = df_stats,
    aes(x = test_group, y = Inf, label = p.signif),
    inherit.aes = FALSE,
    position    = position_dodge(width = 0.8),
    size        = 6,
    vjust       = 2
  ) +

  scale_fill_manual(name = "Disease Group", values = config$group_color) +
  scale_color_manual(guide = "none", values = config$group_color) +

  facet_wrap(~ tissue, nrow = 1, ncol = 6, scales = "free_y") +

  theme_classic() +
  theme(
    panel.border    = element_rect(color = "black", fill = NA, size = 1),
    strip.text      = element_text(size = 12, face = "bold"),
    strip.background= element_blank(),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y     = element_text(size = 12),
    axis.title.x    = element_blank(),
    axis.title.y    = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.title    = element_text(size = 14, face = "bold"),
    legend.text     = element_text(size = 12),
    panel.spacing   = unit(0.5, "lines")
  ) +

  labs(y = "Estimated Fraction")

dir.create(config$output_dirs$boxplot, recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = file.path(config$output_dirs$boxplot, "Fig6F-2.pdf"),
  plot     = p,
  width    = 12,
  height   = 3,
  units    = "in"
)
ggsave(
  filename = file.path(config$output_dirs$boxplot, "Fig6F-2.tif"),
  plot     = p,
  width    = 12,
  height   = 3,
  units    = "in",
  device   = "tiff",
  dpi      = 300
)

# Fig6-SF5-2: clustering heatmap
group_levels <- c("NOR","IPAH","CHD","SLE")
sample_group$group <- factor(sample_group$group, levels = group_levels, ordered = TRUE)
sample_group <- sample_group[order(sample_group$group), ]

annotation_col <- data.frame(
  Group = sample_group$group,
  row.names = sample_group$sample
)
annotation_col$Group <- factor(annotation_col$Group, levels = group_levels)

tissue_levels <- colnames(theta_tissue)
annotation_row <- data.frame(
  Tissue = tissue_levels,
  row.names = tissue_levels
)
annotation_row$Tissue <- factor(annotation_row$Tissue, levels = tissue_levels)

group_colors_used =config$group_color; names(group_colors_used)=group_levels
tissue_colors_used=config$tissue_color;names(tissue_colors_used)=tissue_levels

ann_colors <- list(Group = group_colors_used, Tissue = tissue_colors_used)

ann_colors$Group  <- ann_colors$Group[!is.na(names(ann_colors$Group))]
ann_colors$Tissue <- ann_colors$Tissue[!is.na(names(ann_colors$Tissue))]

palette_cont <- colorRampPalette(
  c("#134B87", "#0f7ab0", "#fdf4af", "#f9b64b", "#a51a49")
)(100)
breaks_cont <- seq(-1, 1, length.out = length(palette_cont) + 1)

global_rng <- range(as.matrix(theta_tissue), na.rm = TRUE)
min_all    <- global_rng[1]
max_all    <- global_rng[2]
theta_norm_global <- theta_tissue %>%
  as.data.frame() %>%
  mutate(across(everything(), ~ {
    if (max_all - min_all == 0) rep(0, length(.x))
    else 2 * (.x - min_all) / (max_all - min_all) - 1
  }))
mat_global <- t(theta_norm_global)

p_global <- pheatmap(
  mat               = mat_global,
  annotation_col    = annotation_col,
  annotation_row    = annotation_row,
  annotation_colors = ann_colors,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  border_color      = NA,
  fontsize_col      = 10,
  color             = palette_cont,
  breaks            = breaks_cont
)

dir.create(config$output_dirs$heatmap, recursive = TRUE, showWarnings = FALSE)
pdf(file.path(config$output_dirs$heatmap, "cluster_propotion_sample_heatmap_global_norm.pdf"), width = 12, height = 8)
print(p_global)
dev.off()
tiff(file.path(config$output_dirs$heatmap, "cluster_propotion_sample_heatmap_global_norm.tif"), width = 12, height = 8, units = "in", res = 300)
print(p_global)
dev.off()

theta_norm_tissue <- theta_tissue %>%
  as.data.frame() %>%
  mutate(across(everything(), ~ {
    rng <- range(.x, na.rm = TRUE)
    if (diff(rng) == 0) rep(0, length(.x))
    else 2 * (.x - rng[1]) / (rng[2] - rng[1]) - 1
  }))

mat_tissue <- t(theta_norm_tissue)
p_tissue <- pheatmap(
  mat               = mat_tissue,
  annotation_col    = annotation_col,
  annotation_row    = annotation_row,
  annotation_colors = ann_colors,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  border_color      = NA,
  fontsize_col      = 10,
  color             = palette_cont,
  breaks            = breaks_cont
)

pdf(file.path(config$output_dirs$heatmap, "Fig6-SF5-2.pdf"), width = 12, height = 8)
print(p_tissue)
dev.off()
tiff(file.path(config$output_dirs$heatmap, "Fig6-SF5-2.tif"), width = 12, height = 8, units = "in", res = 300)
print(p_tissue)
dev.off()

group_levels <- c("NOR","IPAH","CHD","SLE")
sample_group$group <- factor(sample_group$group, levels=group_levels, ordered=TRUE)

combined_data <- data.frame(
  Group = sample_group$group,
  theta_tissue,
  row.names = rownames(theta_tissue)
)

library(dplyr)
df_mat <- combined_data %>% select(-Group) %>% mutate(across(everything(), as.numeric)) %>% as.matrix()
df_mat <- df_mat[, apply(df_mat,2,var, na.rm=TRUE)>0]

res_pca <- prcomp(df_mat, center=TRUE, scale.=TRUE)
pve     <- res_pca$sdev^2/sum(res_pca$sdev^2)
pca_df  <- as.data.frame(res_pca$x[,1:2,drop=FALSE])
pca_df$Group <- factor(combined_data$Group, levels=group_levels)
colnames(pca_df)[1:2] <- c("PC1","PC2")

library(scales)
x_lab <- paste0("PC1 (", percent(pve[1],accuracy=0.1),")")
y_lab <- paste0("PC2 (", percent(pve[2],accuracy=0.1),")")

p_pca <- ggplot(pca_df, aes(x=PC1,y=PC2,color=Group,fill=Group)) +
  stat_ellipse(type="norm",geom="polygon",alpha=0.2,color=NA, shape=19, stroke=0) +
  geom_point(size=4,alpha=0.7) +
  scale_color_manual(name="Group",values=config$group_color,breaks=group_levels) +
  scale_fill_manual(name="Group",values=config$group_color,breaks=group_levels) +
  labs(x=x_lab,y=y_lab) +
  theme_bw(base_size=14) +
  theme(
    legend.position  = "bottom",
    legend.title     = element_text(size=14,face="bold"),
    legend.text      = element_text(size=12),
    axis.title       = element_text(size=14,face="bold"),
    axis.text        = element_text(size=12),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color="black",fill=NA,size=1)
  )

dir.create(config$output_dirs$pca, recursive=TRUE, showWarnings=FALSE)
ggsave(file.path(config$output_dirs$pca,"Fig6G.pdf"),p_pca,width=8,height=8,units="in")
ggsave(file.path(config$output_dirs$pca,"Fig6G.tif"),p_pca,width=8,height=8,units="in",device="tiff",dpi=300)

sampleMEs <- read.csv(config$sampleMEs_file,header=TRUE,row.names=1)
common_samples <- intersect(rownames(sampleMEs),rownames(theta_tissue))
theta_subset <- theta_tissue[rownames(theta_tissue)%in%common_samples,] %>% as.data.frame()

# Fig6-SF7: module-fraction correlation bubble plot
cor_mat <- cor(theta_subset, sampleMEs, use="pairwise.complete.obs")
cor_df  <- melt(cor_mat, varnames=c("Tissue","Module"), value.name="Correlation")
colors  <- config$bubble_colors

p <- ggplot(cor_df, aes(x=Module,y=Tissue)) +
  geom_point(aes(size=abs(Correlation), color=Correlation)) +
  scale_color_gradientn(
    colors = colors,
    limits = c(-0.5, 0.5),
    oob    = scales::squish,
    name   = "Correlation"
  ) +
  scale_size_continuous(range=c(1,12),name="|Correlation|",limits = c(0, 1) ) +
  theme_bw(base_size=14) +
  theme(
    panel.border     = element_rect(color="black",fill=NA,size=1),
    panel.grid.major = element_line(color="grey80",size=0.5),
    panel.grid.minor = element_line(color="grey90",size=0.25),
    axis.text.x      = element_text(angle=45,hjust=1,size=12),
    axis.text.y      = element_text(size=12),
    axis.title       = element_blank()
  )

dir.create(config$output_dirs$bubble, recursive=TRUE, showWarnings=FALSE)
ggsave(file.path(config$output_dirs$bubble,"Fig6-SF7.pdf"),p,width=10,height=12,units="in")
ggsave(file.path(config$output_dirs$bubble,"Fig6-SF7.tif"),p,width=10,height=12,units="in",device="tiff",dpi=300)

# Fig6H: selected module-fraction scatter plots
sampleME_key=sampleMEs[,colnames(sampleMEs) %in% c("MEmidnightblue","MEpurple","MEgreenyellow")]
theta_sub_key=theta_subset[,colnames(theta_subset) %in% c("Heart_Fibroblast", "Blood_Erythrocyte",  "Skin_Myeloid_cell")]

head(rownames(sampleME_key))
head(rownames(theta_sub_key))

df_me <- as.data.frame(sampleME_key) %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to="module", values_to="ME_value")

df_feat <- as.data.frame(theta_sub_key) %>%
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

base_shapes <- c(3,7,8)
if (length(levels(df_joined$trait)) > length(base_shapes)) {
  warning("The number of traits exceeds the predefined shape set; shapes will be recycled.")
}
shape_values <- set_names(
  rep(base_shapes, length.out = length(levels(df_joined$trait))),
  levels(df_joined$trait)
)

library(ggh4x)
library(ggExtra)
p_all <- ggplot(df_joined, aes(x = ME_value, y = trait_value)) +
  geom_point(aes(color = module_color, shape = trait), size = 2) +
  geom_smooth(method = "lm", se = TRUE, linetype = "solid", color = "black") +
  stat_cor(
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
  labs(
    x     = NULL,
    y     = NULL,
    color = "Module",
    shape = "Tissue"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major    = element_blank(),
    panel.grid.minor    = element_blank(),
    panel.border        = element_rect(color = "black", fill = NA, size = 0.5),
    axis.text.x         = element_text(size = 8, color = "black"),
    axis.text.y         = element_text(size = 8, color = "black"),
    axis.ticks          = element_line(color = "black"),
    axis.ticks.length   = unit(2, "pt"),
    axis.title          = element_blank(),
    strip.placement     = "outside",
    strip.text.x        = element_text(size = 11, face = "bold"),
    strip.text.y.left   = element_blank(),
    panel.spacing       = unit(0.5, "lines"),
    legend.position     = "right",
    legend.box          = "vertical",
    legend.title        = element_text(size = 14, face = "bold"),
    legend.text         = element_text(size = 12),
    legend.key.size     = unit(1, "cm"),
    plot.title          = element_blank()
  )

out_dir <- config$output_dirs$scatter_cor

ggsave(
  filename = file.path(out_dir, "Fig6H.pdf"),
  plot     = p_all,
  width    = 10, height = 10, units = "in"
)
ggsave(
  filename = file.path(out_dir, "Fig6H.tif"),
  plot     = p_all,
  width    = 10, height = 10, units = "in",
  device   = "tiff", dpi = 300
)

# Individual scatter plots with marginal densities
modules <- unique(df_joined$module)
traits  <- unique(df_joined$trait)

plot_list <- list()
for (mod in modules) {
  for (tr in traits) {
    subdf <- df_joined %>% filter(module == mod, trait == tr)

    if (nrow(subdf) == 0) next

    p <- ggplot(subdf, aes(x = ME_value, y = trait_value)) +
      geom_point(aes(color = module_color, shape = trait), size = 2, alpha = 0.8) +
      geom_smooth(method = "lm", se = TRUE, color = "black") +
      stat_cor(
        method      = "pearson",
        label.x.npc = 0.02,
        label.y.npc = 0.98,
        label.sep   = ", ",
        size        = 3,
        show.legend = FALSE
      ) +
      scale_color_manual(values = color_map) +
      scale_shape_manual(values = shape_values) +
      theme_minimal(base_size = 12) +
      theme(
        axis.title      = element_blank(),
        panel.grid      = element_blank(),
        panel.border    = element_rect(color = "black", fill = NA),
        legend.position = "none"
      )

    p_marginal <- ggMarginal(
      p,
      type    = "density",
      margins = "both",
      fill    = alpha(color_map[subdf$module_color[1]], 0.8),
      size    = 5
    )

    key <- paste0(mod, "_vs_", tr)
    plot_list[[key]] <- p_marginal
  }
}

for (nm in names(plot_list)) {
  ggsave(
    filename = file.path(out_dir, paste0(nm, ".pdf")),
    plot     = plot_list[[nm]],
    width    = 6,
    height   = 6,
    units    = "in"
  )
}

