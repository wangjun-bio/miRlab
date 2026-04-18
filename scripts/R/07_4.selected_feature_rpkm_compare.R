# =============================================================================
# Downstream Analysis: Expression Boxplots for LASSO-Selected Features
# =============================================================================
# This script:
#   1. Reads feature frequency file for combined cfRNA (all_lasso_selected_feature)
#   2. Filters features with frequency >= 25
#   3. Finds the best iteration and model based on test AUC
#   4. Reads risk scores to split samples into train/test and NOR/PH
#   5. Extracts expression of selected features from RPKM matrix
#   6. Creates boxplots (train/test stratified) for each feature
# =============================================================================

# ----------------------------- Libraries -------------------------------------
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)
library(stringr)
library(readr)

# ----------------------------- User configuration ----------------------------
# Base directories (consistent with 07_1 and 07_2)
base_dir   <- file.path("results", "ML", "Tables", "Feature_AUC_ACC",
                        "NOR_vs_PH_2025-08-26_combined_lasso_selected_features_baseMean_ge100")
risk_dir   <- file.path("results", "ML", "Tables", "Risk_Score",
                        "NOR_vs_PH_2025-08-26_combined_lasso_selected_features_baseMean_ge100")
rpkm_file  <- file.path("data", "df_rpkm_300W.csv")
sampleID_file <- file.path("data", "sampleID_filted300W.rds")

# Output directory
fig_dir <- file.path("figures", "ML", "selected_feature_rpkm")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# ----------------------------- 1. Read feature frequency file --------------
feature_pattern <- "features_frequency_all_lasso_selected_feature_\\d{4}-\\d{2}-\\d{2}\\.csv"
feature_file <- list.files(path = base_dir, pattern = feature_pattern, full.names = TRUE)[1]
if (is.na(feature_file)) stop("Feature frequency file not found in ", base_dir)

feature_df <- read.csv(feature_file)
# Keep features with frequency >= 25 (out of 30 iterations)
selected_features <- feature_df[feature_df$Frequency >= 25, "Feature"]
cat("Number of selected features (freq >= 25):", length(selected_features), "\n")

# ----------------------------- 2. Find best iteration and model ------------
auc_pattern <- "AUC_test_results_all_lasso_selected_feature_\\d{4}-\\d{2}-\\d{2}\\.csv"
auc_file <- list.files(path = base_dir, pattern = auc_pattern, full.names = TRUE)[1]
if (is.na(auc_file)) stop("AUC test file not found in ", base_dir)

auc_df <- read.csv(auc_file)
if ("Iteration" %in% colnames(auc_df)) {
  rownames(auc_df) <- auc_df$Iteration
  auc_df$Iteration <- NULL
}
max_val <- max(auc_df, na.rm = TRUE)
max_pos <- which(auc_df == max_val, arr.ind = TRUE)
best_iter <- as.numeric(rownames(auc_df)[max_pos[1, "row"]])
best_model <- colnames(auc_df)[max_pos[1, "col"]]
cat("Best model:", best_model, "at iteration", best_iter, "with AUC =", max_val, "\n")

# ----------------------------- 3. Read risk scores to get train/test split -
risk_pattern <- "RiskScore_All_all_lasso_selected_feature_\\d{4}-\\d{2}-\\d{2}\\.csv"
risk_file <- list.files(path = risk_dir, pattern = risk_pattern, full.names = TRUE)[1]
if (is.na(risk_file)) stop("Risk score file not found in ", risk_dir)

risk_df <- read.csv(risk_file)
risk_df <- risk_df[risk_df$Iteration == best_iter & risk_df$Model == best_model, ]
# Extract sample IDs for train and test sets
# The risk file contains columns: Sample, Model, RNA_Type, Risk_Score, True_Label, Predicted_Label, Iteration
# We need to split samples into train/test based on the original split.
# Since the risk file only contains test set predictions? Actually from the original script,
# the RiskScore_All file contains both train and test? Let's check: In classification script,
# all_risk_df was saved with a column "Set" (train/test). But the provided risk file may not have that.
# Alternative: read sampleID from RDS and intersect with samples in risk file.
# For this script, we'll assume risk_df contains only test samples? Actually from 07_2, risk_df had
# all samples (train+test) with a "Set" column? The original 07_4 script reads risk file and then
# uses sampleID to split into train/test based on presence in sampleID$PH_not4pnk etc.
# Let's replicate that logic.

sampleID <- readRDS(sampleID_file)
# Get all PH and NOR samples
all_PH <- sampleID$PH_not4pnk
all_NOR <- sampleID$NOR_not4pnk

# The risk file contains all samples (both train and test) because it was saved from cross-validation?
# Actually from the original 07_4 script, they read risk file and then used the iteration's test samples
# from the risk file itself? The original code: 
#   test_PH = df$Sample[df$Sample %in% sampleID$PH_not4pnk]
#   test_NOR = df$Sample[df$Sample %in% sampleID$NOR_not4pnk]
#   train_PH = sampleID$PH_not4pnk[!sampleID$PH_not4pnk %in% test_PH]
#   train_NOR = sampleID$NOR_not4pnk[!sampleID$NOR_not4pnk %in% test_NOR]
# This assumes that the risk file contains only the test samples for that iteration.
# Let's follow that: risk_df contains test samples only.

test_PH <- risk_df$Sample[risk_df$Sample %in% all_PH]
test_NOR <- risk_df$Sample[risk_df$Sample %in% all_NOR]
train_PH <- setdiff(all_PH, test_PH)
train_NOR <- setdiff(all_NOR, test_NOR)

cat("Train - PH:", length(train_PH), "NOR:", length(train_NOR), "\n")
cat("Test  - PH:", length(test_PH), "NOR:", length(test_NOR), "\n")

# ----------------------------- 4. Read RPKM matrix and extract expression ---
rpkm <- read.csv(rpkm_file, row.names = 1, check.names = FALSE)
colnames(rpkm) <- gsub("\\.", "-", colnames(rpkm))  # standardize sample IDs
# Keep only samples present in train/test sets
all_samples <- c(train_PH, train_NOR, test_PH, test_NOR)
rpkm <- rpkm[, colnames(rpkm) %in% all_samples, drop = FALSE]
# Log2 transform
rpkm <- log2(rpkm + 1)

# Subset to selected features
rpkm_sub <- rpkm[rownames(rpkm) %in% selected_features, , drop = FALSE]
cat("Number of selected features found in RPKM:", nrow(rpkm_sub), "\n")

# ----------------------------- 5. Prepare long format with group info -------
group_info <- data.frame(
  sample = all_samples,
  group = factor(
    c(rep("train_NOR", length(train_NOR)), rep("train_PH", length(train_PH)),
      rep("test_NOR", length(test_NOR)), rep("test_PH", length(test_PH))),
    levels = c("train_NOR", "train_PH", "test_NOR", "test_PH")
  ),
  group_type = factor(
    ifelse(grepl("train", group), "train", "test"),
    levels = c("train", "test")
  ),
  ph_type = factor(
    ifelse(grepl("PH", group), "PH", "NOR"),
    levels = c("NOR", "PH")
  )
)

long_df <- rpkm_sub %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "expression") %>%
  left_join(group_info, by = "sample") %>%
  group_by(gene, group) %>%
  # Remove extreme outliers (3*IQR)
  mutate(
    q1 = quantile(expression, 0.25, na.rm = TRUE),
    q3 = quantile(expression, 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    lower_bound = q1 - 3 * iqr,
    upper_bound = q3 + 3 * iqr
  ) %>%
  filter(expression >= lower_bound & expression <= upper_bound) %>%
  ungroup() %>%
  select(-q1, -q3, -iqr, -lower_bound, -upper_bound)

# ----------------------------- 6. Order genes by LASSO coefficient ----------
# Read LASSO coefficients to order genes by absolute value
coef_pattern <- "LASSO_Coefs_All_all_lasso_selected_feature_\\d{4}-\\d{2}-\\d{2}\\.csv"
coef_file <- list.files(path = risk_dir, pattern = coef_pattern, full.names = TRUE)[1]
if (is.na(coef_file)) {
  warning("LASSO coefficient file not found; ordering by gene name")
  gene_order <- sort(unique(long_df$gene))
} else {
  coef_df <- read.csv(coef_file, row.names = 1)
  iter_col <- paste0("iter_", best_iter)
  if (iter_col %in% colnames(coef_df)) {
    coef_vec <- coef_df[[iter_col]]
    names(coef_vec) <- rownames(coef_df)
    # Keep only selected features
    coef_vec <- coef_vec[names(coef_vec) %in% selected_features]
    gene_order <- names(sort(abs(coef_vec), decreasing = TRUE))
  } else {
    gene_order <- sort(unique(long_df$gene))
  }
}
long_df$gene <- factor(long_df$gene, levels = gene_order)

# ----------------------------- 7. Boxplot with median lines ----------------
# Define colors
ph_colors <- c("NOR" = "#1f77b4", "PH" = "#ff7f0e")

p <- ggplot(long_df, aes(x = group_type, y = expression, fill = ph_type)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black",
    size = 0.7,
    linewidth = 0.7,
    outlier.shape = 16
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 95,
    aes(group = interaction(group_type, ph_type)),
    color = "black",
    size = 8,
    position = position_dodge(width = 0.8),
    show.legend = FALSE
  ) +
  scale_x_discrete(expand = c(0.5, 0.5)) +
  facet_wrap(~gene, ncol = 5, scales = "free_y") +
  scale_fill_manual(values = ph_colors) +
  labs(x = NULL, y = "Expression") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10, color = "black", face = "bold"),
    strip.text.x = element_text(size = 10),
    legend.position = "bottom"
  )

ggsave(file.path(fig_dir, "selected_features_expression_boxplots.pdf"),
       plot = p, width = 36, height = 24, units = "cm", dpi = 300)

cat("Boxplots saved to:", file.path(fig_dir, "selected_features_expression_boxplots.pdf"), "\n")