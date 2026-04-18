# =============================================================================
# Differential Expression Analysis for PAH cfRNA Project
# =============================================================================
# This script performs:
#   - DESeq2 analysis for multiple RNA types and group comparisons (1vs1 and 1vs3)
#   - Visualization of DE genes using heatmaps, PCA and t-SNE
# =============================================================================

# ----------------------------- Libraries -------------------------------------
library(readxl)
library(stringr)
library(ggplot2)
library(tidyverse)
library(Rtsne)
library(scatterplot3d)
library(RColorBrewer)
library(reshape2)
library(DESeq2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)
library(dplyr)
library(tidyr)

# ----------------------------- User configuration ----------------------------
data_dir        <- "data"                    # directory containing input RDS files
output_de_dir   <- file.path("results", "DE")   # CSV results directory
output_fig_dir  <- file.path("figures", "DE")   # figures directory

# Create output directories if not exist
for (d in c(output_de_dir, output_fig_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Input file names (relative to data_dir)
RNAname_file      <- file.path(data_dir, "RNAname.rds")
sampleID_file     <- file.path(data_dir, "sampleID_filted300W.rds")
df_counts_file    <- file.path(data_dir, "df_counts.rds")
df_rpkm_file      <- file.path(data_dir, "df_rpkm_rpm.rds")
df_merge_file     <- file.path(data_dir, "df_RNA_ratio.rds")

# ----------------------------- Load data -------------------------------------
RNAname   <- readRDS(RNAname_file)
sampleID  <- readRDS(sampleID_file)
df_counts <- readRDS(df_counts_file)
df_rpkm   <- readRDS(df_rpkm_file)
df_merge  <- readRDS(df_merge_file)

# Keep only samples with >3 million clean reads
df_merge <- df_merge[df_merge$clean_reads > 3000000, ]
df_counts <- df_counts[, colnames(df_counts) %in% df_merge$ID]

# ----------------------------- DE analysis function -------------------------
DE_analyze <- function(x, y, z) {
  # x : RNA type (e.g., "miRNA")
  # y : control group key (e.g., "NOR_not4pnk")
  # z : treatment group key (e.g., "IPAH_not4pnk")
  RNA_type <- x
  output_y <- y
  output_z <- z
  
  y_samples <- sampleID[[y]][sampleID[[y]] %in% colnames(df_counts)]
  z_samples <- sampleID[[z]][sampleID[[z]] %in% colnames(df_counts)]
  
  gene_list <- RNAname[[x]]
  countData <- df_counts[which(rownames(df_counts) %in% gene_list), c(y_samples, z_samples)]
  
  # Ensure numeric
  countData[] <- lapply(countData, as.numeric)
  # Filter genes with median count > 2
  countData <- countData[apply(countData, 1, median) > 2, ]
  
  sample_type <- c(rep("Con", length(y_samples)), rep("Exp", length(z_samples)))
  annotDF <- data.frame(condition = sample_type)
  rownames(annotDF) <- c(y_samples, z_samples)
  
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = annotDF,
                                design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", "Exp", "Con"))
  resdata <- merge(as.data.frame(res),
                   as.data.frame(counts(dds, normalized = TRUE)),
                   by = "row.names", sort = FALSE)
  
  outputname <- file.path(output_de_dir,
                          paste0("DE_", x, "_", output_y, "_", output_z, ".csv"))
  write.csv(resdata, file = outputname, row.names = FALSE)
}

# ----------------------------- Run DE analyses -------------------------------
# 1. 1vs1 comparisons (each disease vs NOR, with/without T4PNK)
rna_types <- c("miRNA", "mRNA", "lncRNA", "rsRNA", "ysRNA", "tRFs", "mttRNA", "snRNA", "snoRNA")

group_pairs_1vs1 <- list(
  c("NOR_t4pnk", "IPAH_t4pnk"),
  c("NOR_not4pnk", "IPAH_not4pnk"),
  c("NOR_not4pnk", "CHD_not4pnk"),
  c("NOR_not4pnk", "SLE_not4pnk")
)

for (rna in rna_types) {
  for (pair in group_pairs_1vs1) {
    DE_analyze(rna, pair[1], pair[2])
  }
}

# 2. 1vs3 comparisons (each disease vs the other two combined)
focus_groups <- c("IPAH_not4pnk", "CHD_not4pnk", "SLE_not4pnk")

for (rna in rna_types) {
  for (focus in focus_groups) {
    other_groups <- setdiff(c("NOR_not4pnk", "IPAH_not4pnk", "CHD_not4pnk", "SLE_not4pnk"), focus)
    merged_group_name <- paste0("Others_not_", focus)
    merged_samples <- unlist(sampleID[other_groups], use.names = FALSE)
    # Add temporary group to sampleID
    sampleID[[merged_group_name]] <- merged_samples
    DE_analyze(rna, focus, merged_group_name)
    # Remove temporary group
    sampleID[[merged_group_name]] <- NULL
  }
}

# ----------------------------- Visualization of DE results ------------------
# Helper functions for processing DE files and plotting

process_de_files <- function(base_path, rna_types, sample_types, analysis_type, top_n) {
  rna_sample_row_names <- list()
  for (rna in rna_types) {
    for (sample in sample_types) {
      if (analysis_type == "1Vs1") {
        filename <- paste0("DE_", rna, "_NOR_not4pnk_", sample, "_not4pnk.csv")
      } else {
        filename <- paste0("DE_", rna, "_", sample, "_not4pnk_Others_not_", sample, "_not4pnk.csv")
      }
      file_path <- file.path(base_path, filename)
      if (file.exists(file_path)) {
        processed_rows <- read_csv(file_path, show_col_types = FALSE) %>%
          filter(padj < 0.05) %>%
          arrange(desc(abs(log2FoldChange))) %>%
          slice_head(n = top_n) %>%
          pull(Row.names)
        list_name <- paste0(rna, "_", sample)
        rna_sample_row_names[[list_name]] <- processed_rows
        cat("Processed:", filename, "->", length(processed_rows), "genes\n")
      } else {
        cat("Warning: file not found -", file_path, "\n")
      }
    }
  }
  return(rna_sample_row_names)
}

merge_rna_types <- function(rna_sample_row_names) {
  element_names <- names(rna_sample_row_names)
  rna_categories <- sapply(strsplit(element_names, "_"), `[`, 1)
  unique_rna <- unique(rna_categories)
  merged_rna <- list()
  for (rna in unique_rna) {
    indices <- which(rna_categories == rna)
    merged_values <- unique(unlist(rna_sample_row_names[indices]))
    merged_rna[[rna]] <- merged_values
    cat("Merged", rna, ":", length(merged_values), "unique genes\n")
  }
  return(merged_rna)
}

plot_heatmap <- function(filtered_rpkm, merged_rna, sampleID, top_n, output_file) {
  # Prepare row annotation (RNA type)
  filtered_genes <- rownames(filtered_rpkm)
  gene_rna_type <- character(length(filtered_genes))
  names(gene_rna_type) <- filtered_genes
  for (rna_type in names(merged_rna)) {
    genes_in_rna <- merged_rna[[rna_type]]
    gene_rna_type[genes_in_rna[genes_in_rna %in% filtered_genes]] <- rna_type
  }
  annotation_row <- data.frame(RNA_type = gene_rna_type, row.names = filtered_genes)
  
  # Prepare column annotation (sample type)
  filtered_samples <- colnames(filtered_rpkm)
  sample_type <- character(length(filtered_samples))
  names(sample_type) <- filtered_samples
  sampleID_sub <- sampleID[names(sampleID) != "PH_not4pnk"]
  for (type in names(sampleID_sub)) {
    samples_in_type <- sampleID_sub[[type]]
    sample_type[samples_in_type[samples_in_type %in% filtered_samples]] <- type
  }
  annotation_col <- data.frame(Sample_type = sample_type, row.names = filtered_samples)
  
  # Scale data
  filtered_rpkm <- log2(filtered_rpkm + 1)
  rpkm_scaled <- t(scale(t(filtered_rpkm)))
  rpkm_scaled <- rpkm_scaled[apply(rpkm_scaled, 1, function(x) all(abs(x) <= 3)), ]
  
  # Order samples
  desired_order <- c("NOR_not4pnk", "IPAH_not4pnk", "CHD_not4pnk", "SLE_not4pnk")
  ordered_sample_names <- unlist(lapply(desired_order, function(x) {
    colnames(rpkm_scaled)[sample_type[colnames(rpkm_scaled)] == x]
  }))
  rpkm_scaled <- rpkm_scaled[, ordered_sample_names]
  annotation_col <- annotation_col[ordered_sample_names, , drop = FALSE]
  
  # Colors
  rna_colors <- brewer.pal(n = length(unique(annotation_row$RNA_type)), name = "Set3")
  names(rna_colors) <- unique(annotation_row$RNA_type)
  sample_colors <- brewer.pal(n = length(unique(annotation_col$Sample_type)), name = "Set2")
  names(sample_colors) <- unique(annotation_col$Sample_type)
  annotation_colors <- list(RNA_type = rna_colors, Sample_type = sample_colors)
  
  p <- pheatmap(rpkm_scaled,
                annotation_row = annotation_row,
                annotation_col = annotation_col,
                annotation_colors = annotation_colors,
                show_rownames = FALSE,
                show_colnames = FALSE,
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                treeheight_row = 0,
                treeheight_col = 0,
                color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(50),
                main = paste0("Filtered RPKM Heatmap (Top ", top_n, " DE Genes)"),
                filename = output_file,
                width = 7, height = 7)
  return(p)
}

plot_pca <- function(filtered_rpkm, sampleID, top_n, output_file) {
  df_tmp <- t(filtered_rpkm)
  df_tmp <- as.data.frame(df_tmp)
  df_tmp <- df_tmp[apply(df_tmp, 1, sum) > 0, ]
  df_tmp <- df_tmp[, apply(df_tmp, 2, sum) > 0]
  
  df_tmp$group <- "NA"
  df_tmp$group[rownames(df_tmp) %in% sampleID[["CHD_not4pnk"]]]  <- "CHD_not4pnk"
  df_tmp$group[rownames(df_tmp) %in% sampleID[["IPAH_not4pnk"]]] <- "IPAH_not4pnk"
  df_tmp$group[rownames(df_tmp) %in% sampleID[["SLE_not4pnk"]]]  <- "SLE_not4pnk"
  df_tmp$group[rownames(df_tmp) %in% sampleID[["NOR_not4pnk"]]]  <- "NOR_not4pnk"
  
  df_pca_data <- df_tmp[, -ncol(df_tmp)]
  df_pca_res <- prcomp(df_pca_data, center = TRUE, scale. = TRUE)
  df_pcs <- data.frame(df_pca_res$x, Species = df_tmp$group)
  pca_var <- df_pca_res$sdev^2 / sum(df_pca_res$sdev^2)
  palette <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4")
  
  p <- ggplot(df_pcs, aes(x = PC1, y = PC2, color = Species)) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    stat_ellipse(aes(x = PC1, y = PC2), linetype = 2, size = 0.5, level = 0.95, inherit.aes = FALSE) +
    theme_bw() +
    scale_color_manual(values = palette) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = paste0("PC1: ", signif(pca_var[1]*100, 3), "%"),
         y = paste0("PC2: ", signif(pca_var[2]*100, 3), "%"),
         title = paste0("PCA of Filtered RPKM (Top ", top_n, " DE Genes)")) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(output_file, p, width = 18, height = 18, units = "cm", dpi = 1200)
  return(p)
}

plot_tsne <- function(filtered_rpkm, sampleID, top_n, output_file) {
  df_tmp <- t(filtered_rpkm)
  df_tmp <- as.data.frame(df_tmp)
  df_tmp <- df_tmp[apply(df_tmp, 1, sum) > 0, ]
  df_tmp <- df_tmp[, apply(df_tmp, 2, sum) > 0]
  
  df_tmp$group <- "NA"
  df_tmp$group[rownames(df_tmp) %in% sampleID[["CHD_not4pnk"]]]  <- "CHD_not4pnk"
  df_tmp$group[rownames(df_tmp) %in% sampleID[["IPAH_not4pnk"]]] <- "IPAH_not4pnk"
  df_tmp$group[rownames(df_tmp) %in% sampleID[["SLE_not4pnk"]]]  <- "SLE_not4pnk"
  df_tmp$group[rownames(df_tmp) %in% sampleID[["NOR_not4pnk"]]]  <- "NOR_not4pnk"
  
  df_tsne_data <- df_tmp[, -ncol(df_tmp)]
  set.seed(1234)
  tsne_res <- Rtsne(as.matrix(df_tsne_data), dims = 2, perplexity = 30, verbose = FALSE, max_iter = 1000)
  df_tsne <- as.data.frame(tsne_res$Y)
  colnames(df_tsne) <- c("tSNE1", "tSNE2")
  df_tsne$Species <- df_tmp$group
  palette <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4")
  
  p <- ggplot(df_tsne, aes(x = tSNE1, y = tSNE2, color = Species)) +
    geom_point(size = 3) +
    theme_bw() +
    scale_color_manual(values = palette) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = "t-SNE 1", y = "t-SNE 2",
         title = paste0("t-SNE of Filtered RPKM (Top ", top_n, " DE Genes)")) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(output_file, p, width = 18, height = 18, units = "cm", dpi = 1200)
  return(p)
}

# ----------------------------- Main visualization loop -----------------------
analysis_types <- c("1Vs1", "1Vs3")
top_n_genes_seq <- seq(100, 100, by = 50)   # can be extended
rna_types_viz <- c("lncRNA", "rsRNA", "ysRNA", "mRNA", "miRNA", "tRFs", "snRNA", "snoRNA")
sample_types_viz <- c("CHD", "SLE", "IPAH")

# Ensure df_rpkm only includes samples with >3M reads
df_rpkm <- df_rpkm[, colnames(df_rpkm) %in% df_merge$ID]

for (analysis_type in analysis_types) {
  for (top_n in top_n_genes_seq) {
    cat("\nProcessing:", analysis_type, "top_n =", top_n, "\n")
    # Process DE files
    rna_sample_row_names <- process_de_files(output_de_dir, rna_types_viz, sample_types_viz, analysis_type, top_n)
    merged_rna <- merge_rna_types(rna_sample_row_names)
    
    # Subset expression matrix
    df_rpkm_filtered <- df_rpkm[, colnames(df_rpkm) %in%
                                  c(sampleID$CHD_not4pnk, sampleID$IPAH_not4pnk,
                                    sampleID$SLE_not4pnk, sampleID$NOR_not4pnk)]
    all_merged_genes <- unique(unlist(merged_rna))
    filtered_rpkm <- df_rpkm_filtered[rownames(df_rpkm_filtered) %in% all_merged_genes, ]
    if (nrow(filtered_rpkm) == 0) {
      cat("No genes selected for top_n =", top_n, "- skipping\n")
      next
    }
    
    # Plot heatmap, PCA, t-SNE
    heatmap_file <- file.path(output_fig_dir, paste0("DE_", analysis_type, "_heatmap_top", top_n, ".pdf"))
    plot_heatmap(filtered_rpkm, merged_rna, sampleID, top_n, heatmap_file)
    
    pca_file <- file.path(output_fig_dir, paste0("DE_", analysis_type, "_pca_top", top_n, ".pdf"))
    plot_pca(filtered_rpkm, sampleID, top_n, pca_file)
    
    tsne_file <- file.path(output_fig_dir, paste0("DE_", analysis_type, "_tsne_top", top_n, ".pdf"))
    plot_tsne(filtered_rpkm, sampleID, top_n, tsne_file)
  }
}

cat("All DE analyses and visualizations completed.\n")