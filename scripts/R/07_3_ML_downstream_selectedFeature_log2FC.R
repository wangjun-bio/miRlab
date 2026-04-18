# =============================================================================
# Downstream Analysis: Selected Features vs Log2FC Barplot
# =============================================================================
# This script:
#   1. Reads feature frequency files (from LASSO stability selection)
#   2. Filters features with frequency >= 25
#   3. Merges with differential expression results (log2FoldChange)
#   4. Creates a side-by-side barplot of Frequency and log2FC
# =============================================================================

# ----------------------------- Libraries -------------------------------------
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)

# ----------------------------- User configuration ----------------------------
# Base directories (relative to project root)
base_dir   <- file.path("results", "ML", "Tables", "Feature_AUC_ACC",
                        "NOR_vs_PH_2025-08-25_Top100DE_baseMean_ge100_dynamic_threshold_class_weight_balanced_30iter")
de_dir     <- file.path("results", "DE")
fig_dir    <- file.path("figures", "ML", "selected_features")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# RNA types (same order as in original script)
rna_types <- c("snoRNA", "snRNA", "lncRNA", "mRNA", "ysRNA", "miRNA", "rsRNA", "tRFs")

# ----------------------------- 1. Read feature frequency files --------------
feature_data <- list()
for (rt in rna_types) {
  file_pattern <- paste0("features_frequency_", rt, "_\\d{4}-\\d{2}-\\d{2}\\.csv")
  file_path <- list.files(path = base_dir, pattern = file_pattern, full.names = TRUE)[1]
  if (!is.na(file_path)) {
    df <- read.csv(file_path)
    # Keep features with frequency >= 25
    df <- df[df$Frequency >= 25, , drop = FALSE]
    feature_data[[rt]] <- df
    cat("Loaded", rt, "with", nrow(df), "features\n")
  } else {
    cat("Warning: No feature frequency file found for", rt, "\n")
  }
}

# ----------------------------- 2. Read DE results ---------------------------
de_data <- list()
for (rt in rna_types) {
  file_path <- file.path(de_dir, paste0("DE_", rt, "_NOR_not4pnk_PH_not4pnk.csv"))
  if (file.exists(file_path)) {
    de_data[[rt]] <- read.csv(file_path)
  } else {
    cat("Warning: DE file not found for", rt, "\n")
  }
}

# ----------------------------- 3. Merge feature frequency and DE log2FC -----
combined_df <- data.frame(
  GeneID = character(),
  RNA_type = character(),
  Frequency = numeric(),
  Log2FC = numeric(),
  stringsAsFactors = FALSE
)

common_types <- intersect(names(feature_data), names(de_data))
for (rt in common_types) {
  feature_df <- feature_data[[rt]]
  de_df <- de_data[[rt]]
  
  # Match genes by Feature (from feature file) and Row.names (from DE file)
  temp_feature <- data.frame(
    GeneID = feature_df$Feature,
    RNA_type = rt,
    Frequency = feature_df$Frequency,
    stringsAsFactors = FALSE
  )
  temp_de <- data.frame(
    GeneID = de_df$Row.names,
    Log2FC = de_df$log2FoldChange,
    stringsAsFactors = FALSE
  )
  temp_combined <- merge(temp_feature, temp_de, by = "GeneID", all.x = TRUE)
  combined_df <- rbind(combined_df, temp_combined)
}

# ----------------------------- 4. Prepare data for plotting ----------------
# Sort by RNA_type and GeneID
rna_order <- c("snoRNA", "snRNA", "lncRNA", "mRNA", "ysRNA", "miRNA", "rsRNA", "tRFs")
combined_df$RNA_type <- factor(combined_df$RNA_type, levels = rna_order)
combined_df <- combined_df[order(combined_df$RNA_type, combined_df$GeneID), ]
combined_df$GeneID <- factor(combined_df$GeneID, levels = unique(combined_df$GeneID))

# Split into Frequency and Log2FC for side-by-side plots
freq_df <- combined_df[, c("GeneID", "RNA_type", "Frequency")]
log2fc_df <- combined_df[, c("GeneID", "RNA_type", "Log2FC")]

# RNA type colors (Nature-style)
rna_colors <- c(
  "tRFs"   = "#70A0AC",
  "rsRNA"  = "#F0EEBB",
  "miRNA"  = "#E2AE79",
  "ysRNA"  = "#7d3f98",
  "mRNA"   = "#D0DCAA",
  "lncRNA" = "#9B3A4D",
  "snRNA"  = "#566CA5",
  "snoRNA" = "#8CBDA7"
)

# Left plot: Frequency
p1 <- ggplot(freq_df, aes(x = GeneID, y = Frequency, fill = RNA_type)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_y_continuous(position = "left", trans = "reverse") +
  scale_x_discrete(position = "bottom") +
  ylab("Frequency") +
  xlab("GeneID") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_line(color = "black"),
    axis.text.y = element_text(hjust = 1, angle = 0),
    axis.text.x.top = element_text(angle = 0, vjust = 0, margin = margin(b = 5)),
    legend.position = "none"
  ) +
  scale_fill_manual(values = rna_colors)

# Right plot: Log2FC
p2 <- ggplot(log2fc_df, aes(x = GeneID, y = Log2FC, fill = RNA_type)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_y_continuous(position = "left") +
  ylab("Log2FC") +
  xlab("") +
  theme_minimal() +
  theme(
    axis.ticks.x = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_text()
  ) +
  scale_fill_manual(values = rna_colors, name = "RNA Type")

# Combine plots
combined_plot <- p1 + p2 + plot_layout(widths = c(1, 1), guides = "collect")
ggsave(file.path(fig_dir, "selected_feature_log2FC.pdf"),
       combined_plot, width = 24, height = 24, units = "cm", dpi = 1200)

cat("Plot saved to:", file.path(fig_dir, "selected_feature_log2FC.pdf"), "\n")