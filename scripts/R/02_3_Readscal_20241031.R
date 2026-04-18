# =============================================================================
# Readscal Processing and RNA Ratio Calculation
# =============================================================================
# This script:
#   1. Reads individual readscal output files (*_cal.txt) from different libraries
#   2. Calculates read counts for various RNA types (rsRNA, ysRNA, tRFs, miRNA, etc.)
#   3. Reads mRNA, lncRNA, snRNA, snoRNA summary files from featureCounts
#   4. Merges all data and calculates unmapped/uncharacterized reads
#   5. Integrates re-sequenced (buce) samples by adding their read counts
#   6. Filters to keep only samples that were actually sequenced
#   7. Saves final data frame as CSV and RDS
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
library(gridExtra)

# ----------------------------- User configuration ----------------------------
# Project root directory (set this to your project root)
# If running from RStudio with project, this is automatic; otherwise set manually.
# For command line, you may need to set working directory first.

# Input directories (relative to project root)
readscal_dir   <- file.path("data", "readscal")           # contains *_cal.txt files
summary_dir    <- file.path("results", "featureCounts", "summaries")  # *_assigned_counts.tsv files
buce_file      <- file.path("data", "buce_ID_20241107.xlsx")
meta_file      <- file.path("data", "Sampl_meta_202409_withoutfunction.xlsx")

# Output directory
output_dir     <- file.path("results", "readscal")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ----------------------------- 1. Read readscal files -----------------------
# Find all *_cal.txt files in readscal directory
file_list <- list.files(path = readscal_dir, pattern = "_cal\\.txt$", full.names = TRUE)
if (length(file_list) == 0) stop("No _cal.txt files found in ", readscal_dir)

# Read first file with header
df <- read.csv(file_list[1], sep = "", header = TRUE)
# Add prefix to Barcode column
lib_date <- paste0(strsplit(basename(file_list[1]), "_cal.txt")[[1]], "_")
df$Barcode <- paste0(lib_date, df$Barcode)
col_names <- colnames(df)

# Read remaining files (no header)
df_list <- lapply(file_list[-1], function(file) {
  df_tmp <- read.csv(file, sep = "", header = FALSE)
  df_tmp <- df_tmp[-1, ]  # remove header row
  colnames(df_tmp) <- col_names
  lib_date <- paste0(strsplit(basename(file), "_cal.txt")[[1]], "_")
  df_tmp$Barcode <- paste0(lib_date, df_tmp$Barcode)
  return(df_tmp)
})

# Combine all
df_combined <- do.call(rbind, c(list(df), df_list))
df <- df_combined

# Convert fraction strings (e.g., "123/456") to numeric counts (first part)
df[, 2:10] <- sapply(df[, 2:10], function(x) sapply(strsplit(as.character(x), "/"), `[`, 1))
df[, 2:10] <- lapply(df[, 2:10], as.numeric)

# Calculate derived columns (differences between successive columns)
new_df <- df
for (i in 4:10) {
  new_df[[i]] <- as.numeric(df[[i-1]]) - as.numeric(df[[i]])
}
new_df[[11]] <- as.numeric(df[, 10])  # last column as is

colnames(new_df) <- c("ID", "raw_reads", "clean_reads", "rsRNA", "ysRNA", "tRFs",
                      "miRNA", "piRNA", "rm25nt", "hg38_reads", "unmapped_reads")

# ----------------------------- 2. Read mRNA/lncRNA/snRNA/snoRNA summaries ---
# These are TSV files with sample_id and assigned_reads (from 02_2 script)

# mRNA
mRNA_file <- file.path(summary_dir, "mRNA_assigned_counts.tsv")
if (file.exists(mRNA_file)) {
  df_mRNA <- read.csv(mRNA_file, sep = "\t", header = TRUE)
  colnames(df_mRNA) <- c("ID", "mRNA_reads")
  df_merge <- merge(new_df, df_mRNA, by = "ID", all.x = TRUE)
} else {
  warning("mRNA summary file not found: ", mRNA_file)
  df_merge <- new_df
}

# lncRNA
lncRNA_file <- file.path(summary_dir, "lncRNA_assigned_counts.tsv")
if (file.exists(lncRNA_file)) {
  df_lncRNA <- read.csv(lncRNA_file, sep = "\t", header = TRUE)
  colnames(df_lncRNA) <- c("ID", "lncRNA_reads")
  df_merge <- merge(df_merge, df_lncRNA, by = "ID", all.x = TRUE)
}

# snRNA
snRNA_file <- file.path(summary_dir, "snRNA_assigned_counts.tsv")
if (file.exists(snRNA_file)) {
  df_snRNA <- read.csv(snRNA_file, sep = "\t", header = TRUE)
  colnames(df_snRNA) <- c("ID", "snRNA_reads")
  df_merge <- merge(df_merge, df_snRNA, by = "ID", all.x = TRUE)
}

# snoRNA
snoRNA_file <- file.path(summary_dir, "snoRNA_assigned_counts.tsv")
if (file.exists(snoRNA_file)) {
  df_snoRNA <- read.csv(snoRNA_file, sep = "\t", header = TRUE)
  colnames(df_snoRNA) <- c("ID", "snoRNA_reads")
  df_merge <- merge(df_merge, df_snoRNA, by = "ID", all.x = TRUE)
}

# Calculate uncharacterized / multiple-mapped reads
# hg38_uncharacterized = total hg38_reads - sum of assigned RNA types
rna_cols <- c("mRNA_reads", "lncRNA_reads", "snRNA_reads", "snoRNA_reads")
# Convert to numeric in case they are character
df_merge[, rna_cols] <- lapply(df_merge[, rna_cols], as.numeric)
df_merge$hg38_uncharacterized_multiplemapped <- df_merge$hg38_reads -
  rowSums(df_merge[, rna_cols], na.rm = TRUE)

# ----------------------------- 3. Integrate re-sequenced (buce) samples -----
if (file.exists(buce_file)) {
  buce <- read_xlsx(buce_file)
  # Ensure numeric columns for accumulation
  cols_to_sum <- 2:16  # adjust if column indices differ
  df_merge[, cols_to_sum] <- lapply(df_merge[, cols_to_sum], as.numeric)
  
  rows_to_delete <- c()
  for (i in 1:nrow(buce)) {
    original_id <- as.character(buce[i, 1])
    index_id   <- as.character(buce[i, 2])
    
    original_row <- which(df_merge$ID == original_id)
    index_row    <- which(df_merge$ID == index_id)
    
    if (length(original_row) > 0 && length(index_row) > 0) {
      # Add index row counts to original row
      df_merge[original_row, cols_to_sum] <- df_merge[original_row, cols_to_sum] +
        df_merge[index_row, cols_to_sum]
      rows_to_delete <- c(rows_to_delete, index_row)
    } else {
      if (length(original_row) == 0) cat(original_id, "not found in df\n")
      if (length(index_row) == 0) cat(index_id, "not found in df\n")
    }
  }
  if (length(rows_to_delete) > 0) {
    df_merge <- df_merge[-rows_to_delete, ]
  }
} else {
  warning("Buce file not found: ", buce_file)
}

# ----------------------------- 4. Filter to actual sequenced samples -------
if (file.exists(meta_file)) {
  meta <- read_xlsx(meta_file, sheet = "临床信息统计")
  # Samples without T4PNK treatment and marked as sequenced
  meta_not4pnk <- meta[which(meta$实验方法 == "No T4PNK treatment" &
                               meta$`已完成建库+测序（2024-9）` == "TRUE"), ]
  # Samples with T4PNK treatment and sequenced
  meta_t4pnk <- meta[which(meta$实验方法 == "T4PNK treatment" &
                             meta$`已完成建库+测序（2024-9）` == "TRUE"), ]
  sample_ID <- rbind(meta_not4pnk, meta_t4pnk)
  # Keep only the four disease categories
  sample_ID <- sample_ID[sample_ID$`分类...86` %in% c("IPAH/HPAH", "CHD-ASD-PH", "健康对照", "SLE-PH"), ]
  
  df_merge <- df_merge[which(df_merge$ID %in% sample_ID$Library_ID), ]
} else {
  warning("Meta file not found: ", meta_file)
}

# ----------------------------- 5. Save output -------------------------------
output_csv <- file.path(output_dir, "df_RNA_ratio.csv")
output_rds <- file.path(output_dir, "df_RNA_ratio.rds")
write.csv(df_merge, file = output_csv, row.names = FALSE)
saveRDS(df_merge, file = output_rds)

cat("Readscal processing completed.\n")
cat("Output saved to:\n", output_csv, "\n", output_rds, "\n")