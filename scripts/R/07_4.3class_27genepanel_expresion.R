# =============================================================================
# 27-gene Panel Expression Analysis for PAH Subtypes
# =============================================================================
# This script:
#   1. Reads RPKM expression matrix and sample group IDs
#   2. Extracts expression of a predefined 27-gene panel
#   3. Calculates mean and SD of log2(RPKM+1) for each gene per group
#   4. Creates a barplot with error bars (mean + SD) faceted by RNA type
# =============================================================================

# ----------------------------- Libraries -------------------------------------
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# ----------------------------- User configuration ----------------------------
# Input file paths (relative to project root)
rpkm_file     <- file.path("data", "df_rpkm_300W.csv")
sampleID_dir  <- "data"   # directory containing *_ID.txt files
rna_dir       <- "data"   # directory containing RNA list files

# Output directory for figure
fig_dir <- file.path("figures", "ML", "27gene_panel")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# ----------------------------- Helper function -------------------------------
read_id_file <- function(filename) {
  path <- file.path(sampleID_dir, filename)
  if (!file.exists(path)) {
    warning("File not found: ", path)
    return(NULL)
  }
  return(trimws(readLines(path, warn = FALSE)))
}

# ----------------------------- 1. Read sample IDs ----------------------------
sampleID <- list(
  CHD_not4pnk  = read_id_file("CHD_not4pnk_ID.txt"),
  SLE_not4pnk  = read_id_file("SLE_not4pnk_ID.txt"),
  NOR_not4pnk  = read_id_file("NOR_not4pnk_ID.txt"),
  IPAH_not4pnk = read_id_file("IPAH_not4pnk_ID.txt"),
  IPAH_t4pnk   = read_id_file("IPAH_t4pnk_ID.txt"),
  NOR_t4pnk    = read_id_file("NOR_t4pnk_ID.txt")
)

# ----------------------------- 2. Read RNA name lists ------------------------
read_rna_file <- function(filename) {
  path <- file.path(rna_dir, filename)
  if (!file.exists(path)) {
    warning("File not found: ", path)
    return(NULL)
  }
  return(trimws(readLines(path, warn = FALSE)))
}

RNAname <- list(
  rsRNA  = read_rna_file("rsRNA.txt"),
  ysRNA  = read_rna_file("ysRNA.txt"),
  tRFs   = read_rna_file("tRFs.txt"),
  miRNA  = read_rna_file("miRNA.txt"),
  mttRNA = read_rna_file("mttRNA.txt"),
  snRNA  = read_rna_file("snRNA.txt"),
  snoRNA = read_rna_file("snoRNA.txt"),
  mRNA   = read_rna_file("mRNA.txt"),
  lncRNA = read_rna_file("lncRNA.txt")
)

# ----------------------------- 3. Read RPKM matrix ---------------------------
rpkm <- read.csv(rpkm_file, row.names = 1, check.names = FALSE)
# Replace dots in column names with hyphens (common in sample IDs)
colnames(rpkm) <- gsub("\\.", "-", colnames(rpkm))

# ----------------------------- 4. Define 27-gene panel -----------------------
# The gene order (X-axis order)
ID_list <- c(
  "CTDSPL", "LPP", "CTDSP2", "VMP1", "CTD-2034I21.1", 
  "RP11-274H24.2", "DLEU2", "RPPH1", "SNORD95", "RP11-20B24.2", 
  "SNORD22", "SNORD26", "tRF-40-8L8NRS9NS334L2H1", "SNORD14E", 
  "hsa-miR-142-3p", "hsa-miR-142-5p", "hsa-miR-4728-3p", 
  "hsa-miR-149-3p", "tRF-25-1KYK37ZXK1", "tRF-35-HXDD1KYK37ZXK1", 
  "Homo-5.8S-27", "tRF-17-Y78ZZYJ", "RNY3-1011", "RNY4-1435", 
  "Homo-28S-1286", "Homo-5.8S-278", "WDR82"
)

# Subset to genes present in RPKM matrix
target_genes <- intersect(ID_list, rownames(rpkm))
cat("Found", length(target_genes), "out of", length(ID_list), "genes in expression matrix\n")

# ----------------------------- 5. Extract and transform expression ----------
# Keep only PAH subtype samples (IPAH, CHD, SLE) for this analysis
groups_of_interest <- c("IPAH_not4pnk", "CHD_not4pnk", "SLE_not4pnk")
samples_to_keep <- unlist(sampleID[groups_of_interest], use.names = FALSE)
rpkm_sub <- rpkm[target_genes, samples_to_keep, drop = FALSE]

# Convert to long format and log2 transform
long_data <- data.frame()
for (group in groups_of_interest) {
  group_samples <- intersect(sampleID[[group]], colnames(rpkm_sub))
  if (length(group_samples) > 0) {
    temp <- rpkm_sub[, group_samples, drop = FALSE] %>%
      as.data.frame() %>%
      mutate(GeneID = rownames(.)) %>%
      pivot_longer(cols = -GeneID, names_to = "Sample", values_to = "RPKM") %>%
      mutate(PAH_Type = group)
    long_data <- rbind(long_data, temp)
  }
}

# Calculate mean and SD of log2(RPKM+1)
summary_df <- long_data %>%
  mutate(log_val = log2(RPKM + 1)) %>%
  group_by(GeneID, PAH_Type) %>%
  summarise(
    mean_val = mean(log_val, na.rm = TRUE),
    sd_val = sd(log_val, na.rm = TRUE),
    .groups = "drop"
  )

# ----------------------------- 6. Map RNA types to genes --------------------
# Build mapping from gene to RNA type using RNAname lists
gene_to_rna <- data.frame()
for (rna_type in names(RNAname)) {
  genes <- RNAname[[rna_type]]
  if (!is.null(genes)) {
    temp <- data.frame(GeneID = genes, RNA_Type = rna_type, stringsAsFactors = FALSE)
    gene_to_rna <- rbind(gene_to_rna, temp)
  }
}

# Merge with summary data
plot_df <- summary_df %>%
  left_join(gene_to_rna, by = "GeneID") %>%
  filter(!is.na(RNA_Type)) %>%
  mutate(
    GeneID = factor(GeneID, levels = ID_list),   # preserve order
    PAH_Type = factor(PAH_Type, levels = groups_of_interest)
  )

# ----------------------------- 7. Define colors for groups ------------------
group_colors <- c(
  "IPAH_not4pnk" = "#FC8D62",
  "CHD_not4pnk"  = "#8DA0CB",
  "SLE_not4pnk"  = "#E78AC3"
)

# ----------------------------- 8. Plot --------------------------------------
p <- ggplot(plot_df, aes(x = GeneID, y = mean_val, fill = PAH_Type)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), color = "black", size = 0.2) +
  geom_errorbar(aes(ymin = mean_val, ymax = mean_val + sd_val),
                position = position_dodge(0.8), width = 0.3, size = 0.3) +
  facet_grid(. ~ RNA_Type, scales = "free_x", space = "free_x") +
  theme_bw() +
  scale_fill_manual(values = group_colors) +
  labs(
    title = "Differential Expression of cfRNA Markers Across PAH Subtypes",
    y = expression(log[2](RPKM + 1)),
    x = NULL,
    fill = "Subtype"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold", size = 9),
    panel.spacing = unit(0.3, "lines"),
    legend.position = "bottom"
  )

print(p)
ggsave(file.path(fig_dir, "PAH_cfRNA_27gene_expression.pdf"),
       width = 18, height = 6, units = "in", dpi = 300)

cat("27-gene panel expression plot saved to:", file.path(fig_dir, "PAH_cfRNA_27gene_expression.pdf"), "\n")