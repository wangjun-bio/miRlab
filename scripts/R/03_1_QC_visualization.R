# =============================================================================
# QC Visualization for PAH cfRNA Project
# =============================================================================
# This script generates:
#   - Violin plots of clean reads distribution across groups
#   - Stacked barplot of RNA composition
#   - Violin plots of RPKM distribution for each RNA type
#   - Correlation scatter plots between groups
#   - PCA plots for each RNA type
# =============================================================================

# ----------------------------- Libraries -------------------------------------
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpmisc)
library(ggExtra)
library(ggsci)
library(viridis)

# ----------------------------- User configuration ----------------------------
# Modify these paths according to your local setup
data_dir   <- "data"       # directory containing input RDS files
output_dir <- "figures"    # directory to save plots

# Input file names (relative to data_dir)
df_merge_file   <- file.path(data_dir, "df_RNA_ratio.rds")
sampleID_file   <- file.path(data_dir, "sampleID_filted300W.rds")
df_rpkm_file    <- file.path(data_dir, "df_rpkm_rpm.rds")
RNAname_file    <- file.path(data_dir, "RNAname.rds")

# Create output directory if not exists
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ----------------------------- Load data -------------------------------------
df_merge <- readRDS(df_merge_file)
# Keep only samples with >3 million clean reads
df_merge <- df_merge[which(df_merge$clean_reads > 3000000), ]

sampleID <- readRDS(sampleID_file)

# ----------------------------- Helper function -------------------------------
extract_and_label <- function(df, ids, group_name) {
  df_subset <- df[df$ID %in% ids, ]
  df_subset$Group <- group_name
  return(df_subset[, c("Group", "clean_reads")])
}

# ----------------------------- Extract groups (no T4PNK) --------------------
df_IPAH_not4pnk <- extract_and_label(df_merge, sampleID$IPAH_not4pnk, "IPAH_not4pnk")
df_CHD_not4pnk  <- extract_and_label(df_merge, sampleID$CHD_not4pnk,  "CHD_not4pnk")
df_SLE_not4pnk  <- extract_and_label(df_merge, sampleID$SLE_not4pnk,  "SLE_not4pnk")
df_NOR_not4pnk  <- extract_and_label(df_merge, sampleID$NOR_not4pnk,  "NOR_not4pnk")

# T4PNK groups
df_IPAH_t4pnk   <- extract_and_label(df_merge, sampleID$IPAH_t4pnk,   "IPAH_t4pnk")
df_NOR_t4pnk    <- extract_and_label(df_merge, sampleID$NOR_t4pnk,    "NOR_t4pnk")

# ----------------------------- Violin plot function --------------------------
plot_violin_reads <- function(data, title, output_file) {
  data$Reads_in_Millions <- data$clean_reads / 1e6
  p <- ggplot(data, aes(x = Group, y = Reads_in_Millions, fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.2, color = "black", alpha = 0.5, outlier.shape = NA) +
    scale_y_continuous(labels = scales::label_number(scale = 1, suffix = "M")) +
    labs(title = title, x = "Group", y = "Clean Reads (Millions)") +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", family = "sans"),
      axis.title = element_text(size = 14, family = "sans"),
      axis.text = element_text(size = 12, family = "sans"),
      text = element_text(family = "sans")
    )
  ggsave(p, file = output_file, width = 20, height = 10, units = "cm")
  return(p)
}

# Generate violin plots for non‑T4PNK groups
plot_violin_reads(
  data = rbind(df_IPAH_not4pnk, df_CHD_not4pnk, df_SLE_not4pnk, df_NOR_not4pnk),
  title = "Clean Reads Distribution Across Groups (No T4PNK)",
  output_file = file.path(output_dir, "clean_reads_violin_not4pnk.pdf")
)

plot_violin_reads(
  data = rbind(df_IPAH_t4pnk, df_NOR_t4pnk),
  title = "Clean Reads Distribution Across Groups (T4PNK)",
  output_file = file.path(output_dir, "clean_reads_violin_t4pnk.pdf")
)

# ----------------------------- RNA proportion stack bar ----------------------
df_merge[, 3] <- as.numeric(df_merge[, 3])
df_merge[, 4:16] <- lapply(df_merge[, 4:16], as.numeric)

# Calculate proportions
df_merge[, c(4:8, 12:16, 11)] <- t(apply(df_merge, 1, function(row) {
  as.numeric(row[c(4:8, 12:16, 11)]) / as.numeric(row[3])
}))

df_long <- df_merge[, c(1, 4:8, 12:16, 11)]
df_long <- melt(df_long, id.vars = "ID")
df_long$Group <- "NA"

df_long$Group[which(df_long$ID %in% sampleID$IPAH_not4pnk)] <- "IPAH"
df_long$Group[which(df_long$ID %in% sampleID$NOR_not4pnk)]  <- "NOR"
df_long$Group[which(df_long$ID %in% sampleID$CHD_not4pnk)]  <- "CHD"
df_long$Group[which(df_long$ID %in% sampleID$SLE_not4pnk)]  <- "SLE"
df_long$Group[which(df_long$ID %in% sampleID$NOR_t4pnk)]    <- "NORT4PNK"
df_long$Group[which(df_long$ID %in% sampleID$IPAH_t4pnk)]   <- "IPAHT4PNK"

df_long$Group <- factor(df_long$Group,
  levels = c("NOR", "IPAH", "CHD", "SLE", "NORT4PNK", "IPAHT4PNK"))

# Order by group and ID
df_long <- df_long[order(df_long$Group, df_long$ID), ]
df_long$ID <- factor(df_long$ID, levels = unique(df_long$ID))

# Helper to get sample ranges for annotation
get_group_range <- function(group_name) {
  group_indices <- which(df_long$Group == group_name)
  group_ids <- unique(df_long$ID[group_indices])
  group_start <- which(levels(df_long$ID) == group_ids[1])
  group_end   <- which(levels(df_long$ID) == group_ids[length(group_ids)])
  return(list(start = group_start, end = group_end))
}

nor_range     <- get_group_range("NOR")
ipah_range    <- get_group_range("IPAH")
chd_range     <- get_group_range("CHD")
sle_range     <- get_group_range("SLE")
nort4pnk_range <- get_group_range("NORT4PNK")
ipaht4pnk_range <- get_group_range("IPAHT4PNK")

p <- ggplot(df_long, aes(x = ID, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Stacked Bar Chart for Samples", x = "Group", y = "Proportion") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12, family = "sans"),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  # Group annotations
  geom_segment(aes(x = nor_range$start + 0.5, xend = nor_range$end - 0.5,
                   y = -0.05, yend = -0.05), color = "black", size = 0.6) +
  annotate("text", x = (nor_range$start + nor_range$end)/2, y = -0.2,
           label = "NOR", size = 4, hjust = 0.5, angle = 90) +
  geom_segment(aes(x = ipah_range$start + 0.5, xend = ipah_range$end - 0.5,
                   y = -0.05, yend = -0.05), color = "black", size = 0.6) +
  annotate("text", x = (ipah_range$start + ipah_range$end)/2, y = -0.2,
           label = "IPAH", size = 4, hjust = 0.5, angle = 90) +
  geom_segment(aes(x = chd_range$start + 0.5, xend = chd_range$end - 0.5,
                   y = -0.05, yend = -0.05), color = "black", size = 0.6) +
  annotate("text", x = (chd_range$start + chd_range$end)/2, y = -0.2,
           label = "CHD", size = 4, hjust = 0.5, angle = 90) +
  geom_segment(aes(x = sle_range$start + 0.5, xend = sle_range$end - 0.5,
                   y = -0.05, yend = -0.05), color = "black", size = 0.6) +
  annotate("text", x = (sle_range$start + sle_range$end)/2, y = -0.2,
           label = "SLE", size = 4, hjust = 0.5, angle = 90) +
  geom_segment(aes(x = nort4pnk_range$start + 0.5, xend = nort4pnk_range$end - 0.5,
                   y = -0.05, yend = -0.05), color = "black", size = 0.6) +
  annotate("text", x = (nort4pnk_range$start + nort4pnk_range$end)/2, y = -0.35,
           label = "NORT4PNK", size = 4, hjust = 0.5, angle = 90) +
  geom_segment(aes(x = ipaht4pnk_range$start + 0.5, xend = ipaht4pnk_range$end - 0.5,
                   y = -0.05, yend = -0.05), color = "black", size = 0.6) +
  annotate("text", x = (ipaht4pnk_range$start + ipaht4pnk_range$end)/2, y = -0.35,
           label = "IPAHT4PNK", size = 4, hjust = 0.5, angle = 90) +
  coord_cartesian(clip = "off")

ggsave(p, file = file.path(output_dir, "RNAProportion.pdf"),
       width = 30, height = 10, units = "cm")

# ----------------------------- RPKM violin plots ----------------------------
df_rpkm <- readRDS(df_rpkm_file)
RNAname <- readRDS(RNAname_file)

prepare_data <- function(df_rpkm, RNA_list, samples, group_name) {
  df_subset <- df_rpkm[which(rownames(df_rpkm) %in% RNA_list),
                       which(colnames(df_rpkm) %in% samples)]
  df_subset <- log2(df_subset + 1)
  df_subset <- melt(df_subset)
  colnames(df_subset)[1] <- "Group"
  df_subset$Group <- group_name
  df_subset <- df_subset[df_subset$value > 0, ]
  return(df_subset)
}

plot_violin_rpkm <- function(data, title, y_label) {
  ggplot(data, aes(x = Group, y = value, fill = Group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.2, color = "black", alpha = 0.5, outlier.shape = NA) +
    scale_y_continuous(labels = scales::label_number(scale = 1)) +
    labs(title = title, x = "Group", y = y_label) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", family = "sans"),
      axis.title = element_text(size = 14, family = "sans"),
      axis.text = element_text(size = 12, family = "sans"),
      text = element_text(family = "sans")
    )
}

process_and_plot <- function(rna_type) {
  groups <- list(
    IPAH    = sampleID$IPAH_not4pnk,
    NOR     = sampleID$NOR_not4pnk,
    CHD     = sampleID$CHD_not4pnk,
    SLE     = sampleID$SLE_not4pnk,
    IPAHT4PNK = sampleID$IPAH_t4pnk,
    NORT4PNK  = sampleID$NOR_t4pnk
  )
  data <- do.call(rbind, lapply(names(groups), function(grp) {
    prepare_data(df_rpkm, RNAname[[rna_type]], groups[[grp]], grp)
  }))
  data$Group <- factor(data$Group,
    levels = c("NOR", "IPAH", "CHD", "SLE", "NORT4PNK", "IPAHT4PNK"))
  p <- plot_violin_rpkm(data,
    title = paste(rna_type, "RPKM Distribution Across Groups"),
    y_label = paste0("cf-", rna_type, " Log2(rpkm+1)"))
  ggsave(p, file = file.path(output_dir, paste0("rpkm_cf-", rna_type, ".pdf")),
         width = 20, height = 10, units = "cm")
}

# Run for each RNA type
rna_types <- c("mRNA", "lncRNA", "snRNA", "snoRNA", "miRNA", "rsRNA", "ysRNA", "tRFs", "mttRNA")
for (rt in rna_types) process_and_plot(rt)

# ----------------------------- Correlation scatter plots --------------------
cor_scatter <- function(x, y, z) {
  # x : RNA type (e.g., "mRNA")
  # y : reference group key (e.g., "NOR_not4pnk")
  # z : comparison group key (e.g., "IPAH_not4pnk")
  df_rpkm_local <- readRDS(df_rpkm_file)
  df_rpkm_local <- df_rpkm_local[which(rownames(df_rpkm_local) %in% RNAname[[x]]), ]
  df_rpkm_local[is.na(df_rpkm_local)] <- 0
  df_rpkm_local <- log2(df_rpkm_local + 1)

  CON <- as.data.frame(rowMeans(df_rpkm_local[, sampleID[[y]]]))
  EXP <- as.data.frame(rowMeans(df_rpkm_local[, sampleID[[z]]]))
  cor_df <- cbind(CON, EXP)
  cor_df <- cor_df[which(apply(cor_df, 1, sum) > 0), ]
  colnames(cor_df) <- c('CON', 'EXP')

  fit <- lm(CON ~ EXP, data = cor_df)
  p <- ggplot(cor_df, aes(x = CON, y = EXP)) +
    geom_point(size = 0.9) +
    geom_smooth(method = 'lm', se = TRUE, color = 'blue') +
    stat_poly_eq() +
    theme(axis.text.x = element_text(size = 30, color = 'black'),
          axis.text.y = element_text(size = 30, color = 'black'),
          axis.title.x = element_text(face = 'bold', size = 30, color = 'black'),
          axis.title.y = element_text(face = 'bold', size = 30, color = 'black')) +
    geom_rug() +
    labs(title = paste0("Scatter plot of normalized expression for ", x),
         x = paste0(y, " (N=", length(sampleID[[y]]), ")"),
         y = paste0(z, " (N=", length(sampleID[[z]]), ")"))
  ggsave(file.path(output_dir, paste0("Scatter_plot_", x, "_", y, "_vs_", z, ".pdf")),
         p, width = 18, height = 18, units = 'cm', dpi = 1200)
}

# Example: IPAH vs NOR
cor_scatter("mRNA",   "NOR_not4pnk", "IPAH_not4pnk")
cor_scatter("lncRNA", "NOR_not4pnk", "IPAH_not4pnk")
cor_scatter("snRNA",  "NOR_not4pnk", "IPAH_not4pnk")
cor_scatter("snoRNA", "NOR_not4pnk", "IPAH_not4pnk")
cor_scatter("miRNA",  "NOR_not4pnk", "IPAH_not4pnk")
cor_scatter("tRFs",   "NOR_not4pnk", "IPAH_not4pnk")
cor_scatter("mttRNA", "NOR_not4pnk", "IPAH_not4pnk")
cor_scatter("rsRNA",  "NOR_not4pnk", "IPAH_not4pnk")
cor_scatter("ysRNA",  "NOR_not4pnk", "IPAH_not4pnk")

# Similar comparisons can be added for CHD, SLE, PH etc.
# (the original script contains many such calls; you can uncomment and adapt as needed)

# ----------------------------- PCA plots -------------------------------------
PCA_plot <- function(x, y, z) {
  # x : expression matrix (data frame)
  # y : RNA type (character)
  # z : subset ("all", "not4pnk", "t4pnk")
  df_tmp <- x[which(rownames(x) %in% RNAname[[y]]), ]
  df_tmp <- data.frame(t(df_tmp))
  df_tmp <- log2(df_tmp + 1)
  df_tmp <- df_tmp[apply(df_tmp, 1, sum) > 0, ]
  df_tmp <- df_tmp[, apply(df_tmp, 2, sum) > 0]

  df_tmp$group <- "NA"
  df_tmp$group[which(rownames(df_tmp) %in% sampleID[["CHD_not4pnk"]])]   <- 'CHD_not4pnk'
  df_tmp$group[which(rownames(df_tmp) %in% sampleID[["IPAH_not4pnk"]])]  <- 'IPAH_not4pnk'
  df_tmp$group[which(rownames(df_tmp) %in% sampleID[["SLE_not4pnk"]])]   <- 'SLE_not4pnk'
  df_tmp$group[which(rownames(df_tmp) %in% sampleID[["NOR_not4pnk"]])]   <- 'NOR_not4pnk'
  df_tmp$group[which(rownames(df_tmp) %in% sampleID[["IPAH_t4pnk"]])]    <- 'IPAH_t4pnk'
  df_tmp$group[which(rownames(df_tmp) %in% sampleID[["NOR_t4pnk"]])]     <- 'NOR_t4pnk'

  # Select top 20% most variable genes
  medians <- apply(df_tmp[, -ncol(df_tmp)], 2, median)
  top_n <- max(1, floor(0.2 * (ncol(df_tmp) - 1)))
  top_genes <- names(sort(medians, decreasing = TRUE))[1:top_n]
  df_top <- df_tmp[, c(top_genes, "group")]

  df_pca <- prcomp(df_top[, -ncol(df_top)], center = TRUE, scale. = TRUE)
  df_pcs <- data.frame(df_pca$x, Species = df_top$group)
  pca_var <- df_pca$sdev^2 / sum(df_pca$sdev^2)
  palette <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4")

  p <- ggplot(df_pcs, aes(x = PC1, y = PC2, color = Species)) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    stat_ellipse(aes(x = PC1, y = PC2), linetype = 2, size = 0.5, level = 0.95,
                 inherit.aes = FALSE) +
    theme_bw() +
    scale_color_manual(values = palette) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = paste0("PC1: ", signif(pca_var[1] * 100, 3), "%"),
         y = paste0("PC2: ", signif(pca_var[2] * 100, 3), "%"),
         title = paste0("PCA of ", y)) +
    theme(plot.title = element_text(hjust = 0.5))

  suffix <- switch(z,
                   "all" = "_all_PCA.pdf",
                   "not4pnk" = "_not4pnk_PCA.pdf",
                   "t4pnk" = "_t4pnk_PCA.pdf")
  ggsave(p, file = file.path(output_dir, paste0("rpkm_", y, suffix)),
         width = 18, height = 18, units = 'cm', dpi = 1200)
}

# All samples
for (rt in rna_types) PCA_plot(df_rpkm, rt, "all")

# Only non‑T4PNK samples
df_rpkm_not4pnk <- df_rpkm[, which(colnames(df_rpkm) %in%
                                   c(sampleID$CHD_not4pnk, sampleID$IPAH_not4pnk,
                                     sampleID$SLE_not4pnk, sampleID$NOR_not4pnk))]
for (rt in rna_types) PCA_plot(df_rpkm_not4pnk, rt, "not4pnk")

# Only T4PNK samples
df_rpkm_t4pnk <- df_rpkm[, which(colnames(df_rpkm) %in%
                                 c(sampleID$NOR_t4pnk, sampleID$IPAH_t4pnk))]
for (rt in rna_types) PCA_plot(df_rpkm_t4pnk, rt, "t4pnk")

# ----------------------------- Combined PCA (all RNA types) -----------------
create_combined_pca <- function(df_rpkm_sub, output_file) {
  rna_types_comb <- c('miRNA', 'rsRNA', 'snoRNA', 'snRNA', 'mRNA', 'lncRNA', 'tRFs', 'ysRNA')
  all_high_genes <- c()
  for (rt in rna_types_comb) {
    gene_list <- RNAname[[rt]]
    df_tmp <- df_rpkm_sub[which(rownames(df_rpkm_sub) %in% gene_list), ]
    df_tmp <- data.frame(t(df_tmp))
    df_tmp <- log2(df_tmp + 1)
    df_tmp <- df_tmp[apply(df_tmp, 1, sum) > 0, ]
    df_tmp <- df_tmp[, apply(df_tmp, 2, sum) > 0]
    medians <- apply(df_tmp, 2, median)
    top_genes <- names(sort(medians, decreasing = TRUE))[1:100]
    all_high_genes <- union(all_high_genes, top_genes)
  }
  combined_df <- df_rpkm_sub[which(rownames(df_rpkm_sub) %in% all_high_genes), ]
  combined_df <- data.frame(t(combined_df))
  combined_df <- log2(combined_df + 1)
  combined_df <- combined_df[apply(combined_df, 1, sum) > 0, ]
  combined_df <- combined_df[, apply(combined_df, 2, sum) > 0]

  combined_df$group <- "NA"
  combined_df$group[which(rownames(combined_df) %in% sampleID[["CHD_not4pnk"]])]   <- 'CHD_not4pnk'
  combined_df$group[which(rownames(combined_df) %in% sampleID[["IPAH_not4pnk"]])]  <- 'IPAH_not4pnk'
  combined_df$group[which(rownames(combined_df) %in% sampleID[["SLE_not4pnk"]])]   <- 'SLE_not4pnk'
  combined_df$group[which(rownames(combined_df) %in% sampleID[["NOR_not4pnk"]])]   <- 'NOR_not4pnk'
  combined_df$group[which(rownames(combined_df) %in% sampleID[["IPAH_t4pnk"]])]    <- 'IPAH_t4pnk'
  combined_df$group[which(rownames(combined_df) %in% sampleID[["NOR_t4pnk"]])]     <- 'NOR_t4pnk'

  df_pca <- prcomp(combined_df[, -ncol(combined_df)], center = TRUE, scale. = TRUE)
  df_pcs <- data.frame(df_pca$x, Species = combined_df$group)
  pca_var <- df_pca$sdev^2 / sum(df_pca$sdev^2)
  palette <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4")

  p <- ggplot(df_pcs, aes(x = PC1, y = PC2, color = Species)) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    stat_ellipse(aes(x = PC1, y = PC2), linetype = 2, size = 0.5, level = 0.95,
                 inherit.aes = FALSE) +
    theme_bw() +
    scale_color_manual(values = palette) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = paste0("PC1: ", signif(pca_var[1] * 100, 3), "%"),
         y = paste0("PC2: ", signif(pca_var[2] * 100, 3), "%"),
         title = "PCA of Combined High Expression Genes") +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(p, file = output_file, width = 18, height = 18, units = 'cm', dpi = 1200)
}

# Only non‑T4PNK samples for combined PCA
create_combined_pca(df_rpkm_not4pnk,
                    file.path(output_dir, "rpkm_combined_RNA_PCA.pdf"))

cat("All QC visualizations completed.\n")