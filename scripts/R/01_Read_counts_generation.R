# =============================================================================
# Read Counts Generation and Normalization
# =============================================================================
# This script:
#   1. Reads sample metadata and filters relevant samples
#   2. Reads raw count files for various RNA types (rsRNA, ysRNA, tRFs, miRNA,
#      mttRNA, mRNA, lncRNA, snRNA, snoRNA) from featureCounts outputs
#   3. Merges count matrices and integrates re-sequenced (buce) samples
#   4. Performs RPKM normalization for mRNA/lncRNA/snRNA/snoRNA (using transcript lengths)
#   5. Performs RPM normalization for other RNA types (miRNA, tRFs, etc.)
#   6. Saves final count matrix, normalized matrix, and sample/RNA annotations
# =============================================================================

# ----------------------------- Libraries -------------------------------------
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(stringr)
library(openxlsx)
library(ggrepel)
library(DESeq2)
library(readxl)
library(Rtsne)
library(scatterplot3d)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)
library(rtracklayer)
library(GenomicRanges)
library(dplyr)

# ----------------------------- User configuration ----------------------------
# Project root directory (set to your project root; if using RStudio project, 
# this is automatically set. Otherwise, set manually.)
# For command line, you may need to run: setwd("/path/to/project")

# Input directories (relative to project root)
meta_file      <- file.path("data", "Sampl_meta_202409_withoutfunction_buce_20241107.xlsx")
buce_file      <- file.path("data", "buce_ID_20241107.xlsx")
result_final   <- file.path("data", "featureCounts_raw")  # raw count files from featureCounts
gtf_file       <- file.path("data", "ref", "gencode.v38.annotation.gtf")

# Output directories
output_dir     <- file.path("results", "expression_matrices")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ----------------------------- 1. Read sample metadata ----------------------
meta <- read_xlsx(meta_file, sheet = "临床信息统计")

# Samples without T4PNK treatment and sequenced
meta_not4pnk <- meta[which(meta$实验方法 == "No T4PNK treatment" & 
                             meta$`已完成建库+测序（2024-9）` == "TRUE"), ]
# Samples with T4PNK treatment and sequenced
meta_t4pnk <- meta[which(meta$实验方法 == "T4PNK treatment" & 
                           meta$`已完成建库+测序（2024-9）` == "TRUE"), ]
sample_ID <- rbind(meta_not4pnk, meta_t4pnk)

# Keep only relevant disease categories
sample_ID <- sample_ID[sample_ID$`分类...86` %in% c("IPAH/HPAH", "CHD-ASD-PH", 
                                                    "健康对照", "SLE-PH"), ]

# ----------------------------- 2. Helper function to read count files -------
read_rna_counts <- function(rna_type, pattern_suffix, has_header = TRUE, 
                            swap_columns = FALSE, keep_cols = NULL) {
  # rna_type: name of RNA (e.g., "rsRNA")
  # pattern_suffix: file name suffix after sample ID (e.g., "Homo-rsRNA_cal.txt")
  # has_header: does file have a header row?
  # swap_columns: whether to swap first two columns (for miRNA/mttRNA)
  # keep_cols: for featureCounts output, columns to keep (e.g., c(1,7) for Geneid and count)
  
  # Get sample IDs that have this file
  seq_name <- paste0(sample_ID$Library_ID, "_")
  file_paths <- file.path(result_final, paste0(seq_name, pattern_suffix))
  exists <- file.exists(file_paths)
  sample_filtered <- sample_ID[exists, ]
  seq_name_filtered <- seq_name[exists]
  
  mydata <- list()
  for (i in seq_along(seq_name_filtered)) {
    cat("Reading", rna_type, ":", seq_name_filtered[i], "\n")
    if (has_header) {
      dat <- read.delim(file_paths[i], comment.char = "#", header = TRUE, sep = "")
    } else {
      dat <- read.delim(file_paths[i], comment.char = "#", header = FALSE, sep = "")
    }
    if (swap_columns) {
      dat[, c(1:2)] <- dat[, c(2:1)]
    }
    if (!is.null(keep_cols)) {
      dat <- dat[, keep_cols]
    }
    if (has_header && "Geneid" %in% colnames(dat)) {
      colnames(dat)[2] <- seq_name_filtered[i]
    } else if (!has_header) {
      colnames(dat)[2] <- seq_name_filtered[i]
    } else {
      colnames(dat)[ncol(dat)] <- seq_name_filtered[i]
    }
    mydata[[i]] <- dat
  }
  
  # Merge all data frames by the first column (feature name)
  df <- Reduce(function(x, y) full_join(x, y, by = colnames(mydata[[1]])[1]), mydata)
  rownames(df) <- df[[1]]
  df <- df[, -1, drop = FALSE]
  return(df)
}

# ----------------------------- 3. Read each RNA type ------------------------
# rsRNA
df_rsRNA <- read_rna_counts("rsRNA", "Homo-rsRNA_cal.txt", has_header = FALSE, swap_columns = FALSE)

# ysRNA
df_ysRNA <- read_rna_counts("ysRNA", "Homo-ysRNA_cal.txt", has_header = FALSE, swap_columns = FALSE)

# tRFs
df_tRFs <- read_rna_counts("tRFs", "tRFs_cal.txt", has_header = FALSE, swap_columns = FALSE)

# miRNA (requires swapping columns)
df_miRNA <- read_rna_counts("miRNA", "miRNA_cal.txt", has_header = FALSE, swap_columns = TRUE)

# mttRNA (requires swapping columns)
df_mttRNA <- read_rna_counts("mttRNA", "mttRNAs_cal.txt", has_header = FALSE, swap_columns = TRUE)

# mRNA (featureCounts output: keep Geneid and 7th column (counts))
df_mRNA <- read_rna_counts("mRNA", "mRNA_S1_transcript_cal.txt", has_header = TRUE, 
                           swap_columns = FALSE, keep_cols = c(1, 7))

# lncRNA
df_lncRNA <- read_rna_counts("lncRNA", "lncRNA_S1_transcript_cal.txt", has_header = TRUE,
                             swap_columns = FALSE, keep_cols = c(1, 7))

# snRNA
df_snRNA <- read_rna_counts("snRNA", "snRNA_S1_transcript_cal.txt", has_header = TRUE,
                            swap_columns = FALSE, keep_cols = c(1, 7))

# snoRNA
df_snoRNA <- read_rna_counts("snoRNA", "snoRNA_S1_transcript_cal.txt", has_header = TRUE,
                             swap_columns = FALSE, keep_cols = c(1, 7))

# ----------------------------- 4. Integrate re-sequenced (buce) samples -----
buce <- read_xlsx(buce_file)

combine_buce_df <- function(df) {
  df[] <- lapply(df, as.numeric)
  cols_to_delete <- c()
  for (i in 1:nrow(buce)) {
    original_id <- paste0(as.character(buce[i, 1]), "_")
    index_id <- paste0(as.character(buce[i, 2]), "_")
    original_col <- which(colnames(df) == original_id)
    index_col <- which(colnames(df) == index_id)
    if (length(original_col) == 0) {
      cat(original_id, "not found in df\n")
    }
    if (length(original_col) > 0 && length(index_col) == 0) {
      cat(index_id, "not found in df\n")
    }
    if (length(original_col) > 0 && length(index_col) > 0) {
      df[, original_col] <- df[, original_col] + df[, index_col]
      cols_to_delete <- c(cols_to_delete, index_col)
    }
  }
  if (length(cols_to_delete) > 0) {
    df <- df[, -cols_to_delete]
  }
  return(df)
}

df_rsRNA   <- combine_buce_df(df_rsRNA)
df_ysRNA   <- combine_buce_df(df_ysRNA)
df_tRFs    <- combine_buce_df(df_tRFs)
df_miRNA   <- combine_buce_df(df_miRNA)
df_mttRNA  <- combine_buce_df(df_mttRNA)
df_mRNA    <- combine_buce_df(df_mRNA)
df_lncRNA  <- combine_buce_df(df_lncRNA)
df_snRNA   <- combine_buce_df(df_snRNA)
df_snoRNA  <- combine_buce_df(df_snoRNA)

# ----------------------------- 5. Create RNA name list and count matrix -----
RNAname <- list(
  rsRNA   = rownames(df_rsRNA),
  ysRNA   = rownames(df_ysRNA),
  tRFs    = rownames(df_tRFs),
  mttRNA  = rownames(df_mttRNA),
  miRNA   = rownames(df_miRNA),
  mRNA    = rownames(df_mRNA),
  lncRNA  = rownames(df_lncRNA),
  snRNA   = rownames(df_snRNA),
  snoRNA  = rownames(df_snoRNA)
)
saveRDS(RNAname, file = file.path(output_dir, "RNAname.rds"))

# Combine all count matrices (rows = features, columns = samples)
df_counts <- do.call(rbind, list(df_rsRNA, df_ysRNA, df_tRFs, df_mttRNA, df_miRNA,
                                 df_mRNA, df_lncRNA, df_snRNA, df_snoRNA))
df_counts[is.na(df_counts)] <- 0

# ----------------------------- 6. Create sample ID lists --------------------
# Re-filter sample_ID to those present in count matrix columns
sample_ID <- sample_ID[paste0(sample_ID$Library_ID, "_") %in% colnames(df_counts), ]
meta <- sample_ID

# Non-T4PNK groups
meta_not4pnk <- meta[meta$实验方法 == "No T4PNK treatment" & 
                       meta$`已完成建库+测序（2024-9）` == "TRUE", ]
CHD_not4pnk  <- meta_not4pnk[meta_not4pnk$`分类...86` == "CHD-ASD-PH", ]$Library_ID
IPAH_not4pnk <- meta_not4pnk[meta_not4pnk$`分类...86` == "IPAH/HPAH", ]$Library_ID
SLE_not4pnk  <- meta_not4pnk[meta_not4pnk$`分类...86` == "SLE-PH", ]$Library_ID
NOR_not4pnk  <- meta_not4pnk[meta_not4pnk$`分类...86` == "健康对照", ]$Library_ID

# T4PNK groups
meta_t4pnk <- meta[meta$实验方法 == "T4PNK treatment" & 
                     meta$`已完成建库+测序（2024-9）` == "TRUE", ]
CHD_t4pnk   <- meta_t4pnk[meta_t4pnk$`分类...86` == "CHD-ASD-PH", ]$Library_ID
IPAH_t4pnk  <- meta_t4pnk[meta_t4pnk$`分类...86` == "IPAH/HPAH", ]$Library_ID
SLE_t4pnk   <- meta_t4pnk[meta_t4pnk$`分类...86` == "SLE-PH", ]$Library_ID
NOR_t4pnk   <- meta_t4pnk[meta_t4pnk$`分类...86` == "健康对照", ]$Library_ID

sampleID <- list(
  CHD_not4pnk = CHD_not4pnk, IPAH_not4pnk = IPAH_not4pnk, SLE_not4pnk = SLE_not4pnk,
  NOR_not4pnk = NOR_not4pnk, CHD_t4pnk = CHD_t4pnk, IPAH_t4pnk = IPAH_t4pnk,
  SLE_t4pnk = SLE_t4pnk, NOR_t4pnk = NOR_t4pnk
)
saveRDS(sampleID, file = file.path(output_dir, "sampleID.rds"))

# ----------------------------- 7. Build final count matrix (samples as columns) --
# Remove trailing underscore from column names
colnames(df_counts) <- sub("_[^_]*$", "", colnames(df_counts))
df_counts <- df_counts[, colnames(df_counts) %in% unlist(sampleID)]

# Save count matrix
saveRDS(df_counts, file = file.path(output_dir, "df_counts.rds"))
write.csv(df_counts, file = file.path(output_dir, "df_counts.csv"))

# ----------------------------- 8. RPKM normalization for mRNA/lncRNA/snRNA/snoRNA --
# Load GTF and compute longest transcript length per gene
gtf_data <- import(gtf_file)
exons <- gtf_data[gtf_data$type == "exon"]
exon_lengths <- width(exons)
transcript_ids <- mcols(exons)$transcript_id
gene_ids <- mcols(exons)$gene_id
gene_names <- mcols(exons)$gene_name

exon_info <- data.frame(transcript_id = transcript_ids, gene_id = gene_ids,
                        gene_name = gene_names, exon_length = exon_lengths)
transcript_lengths <- exon_info %>%
  group_by(transcript_id, gene_name) %>%
  summarise(total_transcript_length = sum(exon_length), .groups = "drop")

# Keep longest transcript per gene
longest_transcripts <- transcript_lengths %>%
  group_by(gene_name) %>%
  arrange(desc(total_transcript_length)) %>%
  slice_head(n = 1) %>%
  ungroup()

rpkm_convert <- function(rna_type) {
  idx <- RNAname[[rna_type]]
  counts_matrix <- df_counts[rownames(df_counts) %in% idx, , drop = FALSE]
  counts_matrix[] <- lapply(counts_matrix, as.numeric)
  gene_lengths <- longest_transcripts$total_transcript_length[
    match(rownames(counts_matrix), longest_transcripts$gene_name)]
  if (any(is.na(gene_lengths))) {
    warning("Some genes missing length info; set to 1")
    gene_lengths[is.na(gene_lengths)] <- 1
  }
  # RPKM = (counts * 1e6) / (gene_lengths/1000 * total_reads)
  # Actually standard: RPK = counts / (gene_lengths/1000); then RPKM = RPK / (total_reads/1e6)
  rpk <- counts_matrix / (gene_lengths / 1000)
  total_reads <- colSums(counts_matrix)
  rpkm_matrix <- sweep(rpk, 2, total_reads / 1e6, FUN = "/")
  return(rpkm_matrix)
}

df_mRNA_rpkm   <- rpkm_convert("mRNA")
df_lncRNA_rpkm <- rpkm_convert("lncRNA")
df_snRNA_rpkm  <- rpkm_convert("snRNA")
df_snoRNA_rpkm <- rpkm_convert("snoRNA")

# ----------------------------- 9. RPM normalization for other RNA types ----
rpm_convert <- function(rna_type) {
  idx <- RNAname[[rna_type]]
  counts_matrix <- df_counts[rownames(df_counts) %in% idx, , drop = FALSE]
  counts_matrix[] <- lapply(counts_matrix, as.numeric)
  total_reads <- colSums(counts_matrix)
  rpm_matrix <- sweep(counts_matrix, 2, total_reads / 1e6, FUN = "/")
  return(rpm_matrix)
}

df_miRNA_rpm  <- rpm_convert("miRNA")
df_tRFs_rpm   <- rpm_convert("tRFs")
df_mttRNA_rpm <- rpm_convert("mttRNA")
df_ysRNA_rpm  <- rpm_convert("ysRNA")
df_rsRNA_rpm  <- rpm_convert("rsRNA")

# ----------------------------- 10. Combine and save normalized matrix ------
df_rpkm_rpm <- rbind(df_mRNA_rpkm, df_lncRNA_rpkm, df_snRNA_rpkm, df_snoRNA_rpkm,
                     df_miRNA_rpm, df_tRFs_rpm, df_mttRNA_rpm, df_ysRNA_rpm, df_rsRNA_rpm)

saveRDS(df_rpkm_rpm, file = file.path(output_dir, "df_rpkm_rpm.rds"))
write.csv(df_rpkm_rpm, file = file.path(output_dir, "df_rpkm_rpm.csv"))

cat("All expression matrices saved to:", output_dir, "\n")