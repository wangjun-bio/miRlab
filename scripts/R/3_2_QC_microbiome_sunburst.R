# =============================================================================
# QC and Microbiome Analysis for PAH cfRNA Project
# =============================================================================
# This script generates:
#   - Sunburst plot of RNA composition
#   - Statistics of detected RNA counts and read fractions per group
#   - Barplots of detected genes and read fractions for each group
# =============================================================================

# ----------------------------- Libraries -------------------------------------
library(openxlsx)
library(sunburstR)
library(htmlwidgets)
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(cowplot)
library(Cairo)

# ----------------------------- User configuration ----------------------------
data_dir    <- "data"          # directory containing input files
output_dir  <- "figures"       # directory to save plots
result_dir  <- "results"       # directory to save intermediate tables

# Create directories if not exist
for (d in c(output_dir, result_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Input file names (relative to data_dir)
RNA_ratio_file <- file.path(data_dir, "RNA_ratio.xlsx")
sampleID_files <- list(
  CHD  = file.path(data_dir, "CHD_not4pnk_ID.txt"),
  IPAH = file.path(data_dir, "IPAH_not4pnk_ID.txt"),
  NOR  = file.path(data_dir, "NOR_not4pnk_ID.txt"),
  SLE  = file.path(data_dir, "SLE_not4pnk_ID.txt")
)

# RPKM and RNA name files (for downstream statistics)
df_rpkm_file <- file.path(data_dir, "df_rpkm_rpm.rds")
RNAname_file <- file.path(data_dir, "RNAname.rds")
sampleID_full_file <- file.path(data_dir, "sampleID.rds")

# ----------------------------- Sunburst plot ---------------------------------
# Read RNA ratio data
RNA_ratio <- read_excel(RNA_ratio_file, sheet = "reads_cal-原始")

# Read sample IDs for each group
CHD_ID  <- read.table(sampleID_files$CHD,  quote = "\"", comment.char = "")$V1
IPAH_ID <- read.table(sampleID_files$IPAH, quote = "\"", comment.char = "")$V1
NOR_ID  <- read.table(sampleID_files$NOR,  quote = "\"", comment.char = "")$V1
SLE_ID  <- read.table(sampleID_files$SLE,  quote = "\"", comment.char = "")$V1

# Filter to keep only samples used in the paper (452 samples)
df <- RNA_ratio[RNA_ratio$ID %in% c(CHD_ID, IPAH_ID, NOR_ID, SLE_ID), ]

# Save filtered data to results directory (optional)
wb <- loadWorkbook(RNA_ratio_file)
addWorksheet(wb, "保留文章中的452个样本")
writeData(wb, "保留文章中的452个样本", df)
saveWorkbook(wb, file.path(result_dir, "RNA_ratio_filtered.xlsx"), overwrite = TRUE)

# Build data frame for sunburst plot
df_plot <- data.frame(
  paths = c(
    "Human genome - cf_tsRNA",
    "Human genome - cf_rsRNA",
    "Human genome - cf_miRNA",
    "Human genome - cf_ysRNA",
    "Human genome - cf_mRNA",
    "Human genome - cf_lncRNA",
    "Human genome - cf_piRNA",
    "Human genome - cf_snRNA",
    "Human genome - cf_snoRNA",
    "Human genome - uncharacterized",
    "Unmapped - Discarded (<36nt)",
    "Unmapped - Microbe_classified (<36nt)",
    "Unmapped - Microbe_unclassified (<36nt)"
  ),
  values = c(
    mean(df$tRNA_ratio),
    mean(df$rsRNA_ratio),
    mean(df$miRNA_ratio),
    mean(df$ysRNA_ratio),
    mean(df$mRNA_ratio),
    mean(df$lncRNA_ratio),
    mean(df$piRNA_ratio),
    mean(df$snRNA_ratio),
    mean(df$snoRNA_ratio),
    mean(df$uncharacterized_ratio),
    mean(df$trim36nt后丢弃的_ratio),
    mean(df$Microbe_classified_ratio),
    mean(df$Microbe_unclassified_ratio)
  )
)

# Generate sunburst plot (HTML widget)
sunburst_plot <- sunburst(df_plot)

# Save as HTML file (can be converted to PDF manually)
html_file <- file.path(output_dir, "sunburst_RNA_composition.html")
saveWidget(sunburst_plot, file = html_file, selfcontained = TRUE)
cat("Sunburst plot saved as", html_file, "\n")

# ----------------------------- Statistics per group -------------------------
# Load RPKM and RNA name data
df_rpkm <- readRDS(df_rpkm_file)
RNAname <- readRDS(RNAname_file)
sampleID_full <- readRDS(sampleID_full_file)

# Define sample groups (only non-T4PNK groups)
group_list <- list(
  CHD  = CHD_ID,
  IPAH = IPAH_ID,
  NOR  = NOR_ID,
  SLE  = SLE_ID
)

# Initialize Excel workbooks for output
wb_counts <- createWorkbook()
wb_ratio  <- createWorkbook()

for (group_name in names(group_list)) {
  sample_ids <- group_list[[group_name]]
  
  # Extract RNA ratio data for this group
  tmp_ratio <- RNA_ratio[RNA_ratio$ID %in% sample_ids, c(1, 21:31)]
  addWorksheet(wb_ratio, sheetName = group_name)
  writeData(wb_ratio, sheet = group_name, x = tmp_ratio)
  
  # Count detected RNAs (RPKM > 1) for each RNA type
  rna_count_list <- list()
  for (rna_type in names(RNAname)) {
    rna_list <- RNAname[[rna_type]]
    tmp <- df_rpkm[rownames(df_rpkm) %in% rna_list,
                   colnames(df_rpkm) %in% sample_ids]
    non_zero_counts <- colSums(tmp > 1)
    rna_count_list[[rna_type]] <- non_zero_counts
  }
  
  rna_counts_df <- as.data.frame(rna_count_list)
  rna_counts_df$Sample <- rownames(rna_counts_df)
  rna_counts_df <- rna_counts_df[, c("Sample", setdiff(names(rna_counts_df), "Sample"))]
  
  addWorksheet(wb_counts, sheetName = group_name)
  writeData(wb_counts, sheet = group_name, x = rna_counts_df)
}

# Save Excel files to results directory
saveWorkbook(wb_counts,
             file = file.path(result_dir, "RNA_counts_by_group.xlsx"),
             overwrite = TRUE)
saveWorkbook(wb_ratio,
             file = file.path(result_dir, "RNA_ratio_by_group.xlsx"),
             overwrite = TRUE)

# ----------------------------- Barplots for each group ----------------------
process_group_data <- function(group_name, count_file, ratio_file) {
  count_df <- read_excel(count_file, sheet = group_name)
  ratio_df <- read_excel(ratio_file, sheet = group_name)
  
  df <- data.frame(
    RNA_type = c("cf-rsRNA", "cf-ysRNA", "cf-tsRNA", "cf-miRNA",
                 "cf-mRNA", "cf-lncRNA", "cf-snRNA", "cf-snoRNA"),
    Detected_Genes = c(mean(count_df$rsRNA), mean(count_df$ysRNA),
                       mean(count_df$tRFs), mean(count_df$miRNA),
                       mean(count_df$mRNA), mean(count_df$lncRNA),
                       mean(count_df$snRNA), mean(count_df$snoRNA)),
    Detected_Genes_sd = c(sd(count_df$rsRNA)/sqrt(length(count_df$rsRNA)),
                          sd(count_df$ysRNA)/sqrt(length(count_df$ysRNA)),
                          sd(count_df$tRFs)/sqrt(length(count_df$tRFs)),
                          sd(count_df$miRNA)/sqrt(length(count_df$miRNA)),
                          sd(count_df$mRNA)/sqrt(length(count_df$mRNA)),
                          sd(count_df$lncRNA)/sqrt(length(count_df$lncRNA)),
                          sd(count_df$snRNA)/sqrt(length(count_df$snRNA)),
                          sd(count_df$snoRNA)/sqrt(length(count_df$snoRNA))),
    Reads_Fraction = c(mean(ratio_df$rsRNA_ratio), mean(ratio_df$ysRNA_ratio),
                       mean(ratio_df$tRNA_ratio), mean(ratio_df$miRNA_ratio),
                       mean(ratio_df$mRNA_ratio), mean(ratio_df$lncRNA_ratio),
                       mean(ratio_df$snRNA_ratio), mean(ratio_df$snoRNA_ratio)),
    Reads_Fraction_sd = c(sd(ratio_df$rsRNA_ratio)/sqrt(length(ratio_df$rsRNA_ratio)),
                          sd(ratio_df$ysRNA_ratio)/sqrt(length(ratio_df$ysRNA_ratio)),
                          sd(ratio_df$tRNA_ratio)/sqrt(length(ratio_df$tRNA_ratio)),
                          sd(ratio_df$miRNA_ratio)/sqrt(length(ratio_df$miRNA_ratio)),
                          sd(ratio_df$mRNA_ratio)/sqrt(length(ratio_df$mRNA_ratio)),
                          sd(ratio_df$lncRNA_ratio)/sqrt(length(ratio_df$lncRNA_ratio)),
                          sd(ratio_df$snRNA_ratio)/sqrt(length(ratio_df$snRNA_ratio)),
                          sd(ratio_df$snoRNA_ratio)/sqrt(length(ratio_df$snoRNA_ratio)))
  )
  
  desired_order <- c("cf-tsRNA", "cf-rsRNA", "cf-miRNA", "cf-ysRNA",
                     "cf-mRNA", "cf-lncRNA", "cf-snRNA", "cf-snoRNA")
  df$RNA_type <- factor(df$RNA_type, levels = desired_order)
  df <- df[order(df$RNA_type), ]
  df$RNA_type <- factor(df$RNA_type, levels = rev(levels(df$RNA_type)))
  
  df <- df %>%
    mutate(log_Genes = log10(Detected_Genes),
           log_Genes_err = log10(Detected_Genes + Detected_Genes_sd))
  return(df)
}

# File paths for the generated Excel files
count_file <- file.path(result_dir, "RNA_counts_by_group.xlsx")
ratio_file <- file.path(result_dir, "RNA_ratio_by_group.xlsx")

groups <- c("CHD", "IPAH", "NOR", "SLE")
df_list <- lapply(groups, function(g) process_group_data(g, count_file, ratio_file))

# Define colors for RNA types (Nature-style)
rna_colors <- c(
  "cf-lncRNA" = "#9B3A4D",
  "cf-miRNA"  = "#E2AE79",
  "cf-mRNA"   = "#D0DCAA",
  "cf-rsRNA"  = "#F0EEBB",
  "cf-snoRNA" = "#8CBDA7",
  "cf-snRNA"  = "#566CA5",
  "cf-tsRNA"  = "#70A0AC",
  "cf-ysRNA"  = "#7d3f98"
)

plot_group <- function(df, group_label) {
  p1 <- ggplot(df, aes(x = RNA_type, y = -log_Genes, fill = RNA_type)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_segment(aes(x = RNA_type, xend = RNA_type,
                     y = -log_Genes, yend = -log_Genes_err), linewidth = 0.6) +
    geom_segment(aes(x = as.numeric(RNA_type) - 0.1,
                     xend = as.numeric(RNA_type) + 0.1,
                     y = -log_Genes_err, yend = -log_Genes_err), linewidth = 0.6) +
    coord_flip() +
    labs(x = NULL, y = "log10(Detected Genes)", title = group_label) +
    scale_fill_manual(values = rna_colors) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = "black"),
      plot.background = element_rect(fill = "white", color = NA),
      axis.text.y = element_text(size = 12, family = "Arial", color = "black", hjust = 1),
      axis.ticks.y = element_line(color = "black", linewidth = 0.8),
      axis.ticks.x = element_line(color = "black", linewidth = 0.8),
      axis.text.x = element_text(size = 10, family = "Arial", color = "black"),
      axis.title.x = element_text(size = 12, family = "Arial", color = "black"),
      plot.title = element_text(size = 14, family = "Arial", face = "bold"),
      legend.position = "none"
    )
  
  p2 <- ggplot(df, aes(x = RNA_type, y = Reads_Fraction, fill = RNA_type)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_segment(aes(x = RNA_type, xend = RNA_type,
                     y = Reads_Fraction, yend = Reads_Fraction + Reads_Fraction_sd), linewidth = 0.6) +
    geom_segment(aes(x = as.numeric(RNA_type) - 0.1,
                     xend = as.numeric(RNA_type) + 0.1,
                     y = Reads_Fraction + Reads_Fraction_sd,
                     yend = Reads_Fraction + Reads_Fraction_sd), linewidth = 0.6) +
    coord_flip() +
    labs(x = NULL, y = "Fraction of reads (%)") +
    scale_fill_manual(values = rna_colors) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.ticks.x = element_line(color = "black", linewidth = 0.8),
      panel.background = element_rect(fill = "white", color = "black"),
      plot.background = element_rect(fill = "white", color = NA),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 10, family = "Arial", color = "black"),
      axis.title.x = element_text(size = 12, family = "Arial", color = "black"),
      legend.position = "none"
    )
  
  return(p1 + p2)
}

# Generate plots for each group and save as single PDF
plot_list <- lapply(seq_along(df_list), function(i) {
  plot_group(df_list[[i]], groups[i])
})
final_plot <- wrap_plots(plotlist = plot_list, ncol = 1)

ggsave(file.path(output_dir, "RNA_detection_by_group.pdf"),
       plot = final_plot, width = 32, height = 20, units = "cm", dpi = 300)

# ----------------------------- Combined NOR vs PH plots ---------------------
combine_group_data <- function(group_names, count_file, ratio_file) {
  count_dfs <- lapply(group_names, function(g) read_excel(count_file, sheet = g))
  ratio_dfs <- lapply(group_names, function(g) read_excel(ratio_file, sheet = g))
  count_all <- bind_rows(count_dfs)
  ratio_all <- bind_rows(ratio_dfs)
  
  df <- data.frame(
    RNA_type = c("cf-rsRNA", "cf-ysRNA", "cf-tsRNA", "cf-miRNA",
                 "cf-mRNA", "cf-lncRNA", "cf-snRNA", "cf-snoRNA"),
    Detected_Genes = c(mean(count_all$rsRNA), mean(count_all$ysRNA),
                       mean(count_all$tRFs), mean(count_all$miRNA),
                       mean(count_all$mRNA), mean(count_all$lncRNA),
                       mean(count_all$snRNA), mean(count_all$snoRNA)),
    Detected_Genes_sd = c(sd(count_all$rsRNA)/sqrt(length(count_all$rsRNA)),
                          sd(count_all$ysRNA)/sqrt(length(count_all$ysRNA)),
                          sd(count_all$tRFs)/sqrt(length(count_all$tRFs)),
                          sd(count_all$miRNA)/sqrt(length(count_all$miRNA)),
                          sd(count_all$mRNA)/sqrt(length(count_all$mRNA)),
                          sd(count_all$lncRNA)/sqrt(length(count_all$lncRNA)),
                          sd(count_all$snRNA)/sqrt(length(count_all$snRNA)),
                          sd(count_all$snoRNA)/sqrt(length(count_all$snoRNA))),
    Reads_Fraction = c(mean(ratio_all$rsRNA_ratio), mean(ratio_all$ysRNA_ratio),
                       mean(ratio_all$tRNA_ratio), mean(ratio_all$miRNA_ratio),
                       mean(ratio_all$mRNA_ratio), mean(ratio_all$lncRNA_ratio),
                       mean(ratio_all$snRNA_ratio), mean(ratio_all$snoRNA_ratio)),
    Reads_Fraction_sd = c(sd(ratio_all$rsRNA_ratio)/sqrt(length(ratio_all$rsRNA_ratio)),
                          sd(ratio_all$ysRNA_ratio)/sqrt(length(ratio_all$ysRNA_ratio)),
                          sd(ratio_all$tRNA_ratio)/sqrt(length(ratio_all$tRNA_ratio)),
                          sd(ratio_all$miRNA_ratio)/sqrt(length(ratio_all$miRNA_ratio)),
                          sd(ratio_all$mRNA_ratio)/sqrt(length(ratio_all$mRNA_ratio)),
                          sd(ratio_all$lncRNA_ratio)/sqrt(length(ratio_all$lncRNA_ratio)),
                          sd(ratio_all$snRNA_ratio)/sqrt(length(ratio_all$snRNA_ratio)),
                          sd(ratio_all$snoRNA_ratio)/sqrt(length(ratio_all$snoRNA_ratio)))
  )
  desired_order <- c("cf-tsRNA", "cf-rsRNA", "cf-miRNA", "cf-ysRNA",
                     "cf-mRNA", "cf-lncRNA", "cf-snRNA", "cf-snoRNA")
  df$RNA_type <- factor(df$RNA_type, levels = desired_order)
  df <- df[order(df$RNA_type), ]
  df$RNA_type <- factor(df$RNA_type, levels = rev(levels(df$RNA_type)))
  df <- df %>%
    mutate(log_Genes = log10(Detected_Genes),
           log_Genes_err = log10(Detected_Genes + Detected_Genes_sd))
  return(df)
}

df_NOR <- process_group_data("NOR", count_file, ratio_file)
df_PH  <- combine_group_data(c("CHD", "IPAH", "SLE"), count_file, ratio_file)

plot_NOR <- plot_group(df_NOR, "NOR")
plot_PH  <- plot_group(df_PH, "PH")
final_combined <- wrap_plots(plot_NOR, plot_PH, ncol = 1)

ggsave(file.path(output_dir, "NOR_PH_comparison.pdf"),
       plot = final_combined, width = 32, height = 20, units = "cm", dpi = 300,
       device = cairo_pdf)

cat("All QC plots and tables generated successfully.\n")