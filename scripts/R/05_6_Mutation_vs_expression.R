# =============================================================================
# Mutation vs Expression Association Analysis
# =============================================================================
# This script:
#   1. Reads mutation presence matrix and RPKM expression matrix
#   2. For each gene, compares expression between mutated vs unmutated samples,
#      and mutation frequency between Disease and Normal groups
#   3. Generates scatter plot of Odds Ratio vs log2FoldChange
#   4. Creates boxplots for significantly associated genes
# =============================================================================

# ----------------------------- Libraries -------------------------------------
library(dplyr)
library(ggplot2)

# ----------------------------- User configuration ----------------------------
# Input directories (relative to project root)
mut_matrix_file <- file.path("results", "mut_presence_matrix", "ALL_DP0_GT11_presence.rds")
exp_matrix_file <- file.path("data", "df_rpkm_rpm.rds")
sampleID_files <- list(
  IPAH = file.path("data", "IPAH_not4pnk_ID.txt"),
  NOR  = file.path("data", "NOR_not4pnk_ID.txt"),
  SLE  = file.path("data", "SLE_not4pnk_ID.txt"),
  CHD  = file.path("data", "CHD_not4pnk_ID.txt")
)

# Output directory for figures
fig_dir <- file.path("figures", "mutation_vs_expression")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# ----------------------------- Load data -------------------------------------
# Mutation presence matrix (rows = variants, columns = samples)
mut_matrix <- readRDS(mut_matrix_file)
# Keep only variant ID and sample columns (remove annotation columns if present)
mut_matrix <- mut_matrix[, c(2, 14:ncol(mut_matrix))]

# Sum mutations per gene (if multiple variants per gene)
mut_matrix <- mut_matrix %>%
  group_by(Hugo_Symbol) %>%
  summarise(across(everything(), sum), .groups = "drop")

# Expression matrix (RPKM)
exp_matrix <- readRDS(exp_matrix_file)

# Convert mut_matrix to data frame with gene names as rownames
mut_matrix <- as.data.frame(mut_matrix)
rownames(mut_matrix) <- mut_matrix$Hugo_Symbol
mut_matrix$Hugo_Symbol <- NULL

# Keep only samples present in both matrices
common_samples <- intersect(colnames(exp_matrix), colnames(mut_matrix))
exp_matrix <- exp_matrix[, common_samples]
mut_matrix <- mut_matrix[, common_samples]

# Find common genes
common_genes <- intersect(rownames(exp_matrix), rownames(mut_matrix))
cat("Common genes:", length(common_genes), "\n")
exp_matrix <- exp_matrix[common_genes, ]
mut_matrix <- mut_matrix[common_genes, ]

# Binarize mutation matrix (presence/absence)
mut_matrix[mut_matrix > 1] <- 1

# Ensure same sample order
exp_matrix <- exp_matrix[, match(colnames(mut_matrix), colnames(exp_matrix))]

# ----------------------------- Sample grouping -------------------------------
# Read sample IDs for each group
read_sample_ids <- function(file_path) {
  read.table(file_path, stringsAsFactors = FALSE)$V1
}
IPAH_idx <- read_sample_ids(sampleID_files$IPAH)
NOR_idx  <- read_sample_ids(sampleID_files$NOR)
SLE_idx  <- read_sample_ids(sampleID_files$SLE)
CHD_idx  <- read_sample_ids(sampleID_files$CHD)

group_info <- data.frame(
  Sample = colnames(exp_matrix),
  Group = case_when(
    colnames(exp_matrix) %in% c(SLE_idx, CHD_idx, IPAH_idx) ~ "Disease",
    colnames(exp_matrix) %in% NOR_idx ~ "Normal",
    TRUE ~ "Unknown"
  ),
  stringsAsFactors = FALSE
)

# Convert to matrices for faster computation
exp_matrix <- as.matrix(exp_matrix)
mut_matrix <- as.matrix(mut_matrix)

# ----------------------------- Per-gene analysis ----------------------------
results <- data.frame()
for (gene in rownames(exp_matrix)) {
  expr <- exp_matrix[gene, ]
  mut <- mut_matrix[gene, ]
  df <- data.frame(
    Sample = names(expr),
    Expression = expr,
    Mutation = mut,
    Group = group_info$Group[match(names(expr), group_info$Sample)]
  )
  
  # log2FoldChange (Disease vs Normal)
  mean_expr_disease <- mean(df$Expression[df$Group == "Disease"], na.rm = TRUE)
  mean_expr_normal <- mean(df$Expression[df$Group == "Normal"], na.rm = TRUE)
  log2FC <- log2((mean_expr_disease + 1) / (mean_expr_normal + 1))
  
  # Mutation frequency comparison (Fisher's exact test)
  tbl <- table(df$Mutation, df$Group)
  if (all(dim(tbl) == c(2, 2))) {
    fisher_test <- fisher.test(tbl)
    fisher_p <- fisher_test$p.value
    odds_ratio <- as.numeric(1 / fisher_test$estimate)
  } else {
    fisher_p <- NA
    odds_ratio <- NA
  }
  
  # Expression difference (Wilcoxon test)
  wilcox_p <- tryCatch(wilcox.test(Expression ~ Group, data = df)$p.value,
                       error = function(e) NA)
  
  results <- rbind(results, data.frame(
    Gene = gene,
    Mean_Expr_Disease = mean_expr_disease,
    Mean_Expr_Normal = mean_expr_normal,
    log2FoldChange = log2FC,
    Expression_p = wilcox_p,
    Odds_Ratio = odds_ratio,
    Mutation_p = fisher_p
  ))
}

# Add classification labels
results <- results %>%
  mutate(
    Mutation_Sig = ifelse(Mutation_p < 0.05, "Yes", "No"),
    Expression_Sig = ifelse(Expression_p < 0.05, "Yes", "No"),
    Class = case_when(
      Mutation_Sig == "Yes" & Expression_Sig == "No" ~ "Real Alteration",
      Mutation_Sig == "Yes" & Expression_Sig == "Yes" ~
        ifelse(
          (log2FoldChange > log2(1.12) & Odds_Ratio > 1.5) |
          (log2FoldChange < log2(0.88) & Odds_Ratio < 0.5),
          "Possibly Expression-Driven", "Real Alteration"
        ),
      TRUE ~ "Other"
    )
  )

# Keep only genes with significant mutation association
results_sig <- results[results$Mutation_p < 0.05, ]

# ----------------------------- Scatter plot ----------------------------------
p <- ggplot(results_sig, aes(x = Odds_Ratio, y = log2FoldChange, color = Class)) +
  annotate("rect", xmin = -Inf, xmax = 0.5,
           ymin = log2(0.88), ymax = log2(1.12), fill = "grey80", alpha = 0.5) +
  annotate("rect", xmin = 1.5, xmax = Inf,
           ymin = log2(0.88), ymax = log2(1.12), fill = "grey80", alpha = 0.5) +
  geom_point(size = 6, alpha = 0.7, shape = 16) +
  geom_vline(xintercept = c(0.5, 1.5), linetype = "dashed", col = "grey40") +
  geom_hline(yintercept = c(log2(0.88), log2(1.12)), linetype = "dashed", col = "grey40") +
  scale_color_manual(values = c(
    "Real Alteration" = "#EF767A",
    "Possibly Expression-Driven" = "#456990",
    "Other" = "#48C0AA"
  )) +
  labs(x = "Odds Ratio", y = "log2FoldChange",
       title = "Expression vs Mutation", color = "Class") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = unit(0.15, "cm"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

ggsave(p, filename = file.path(fig_dir, "expression_vs_mutation_scatter.pdf"),
       width = 10, height = 10)

# Identify significant genes for boxplot
sig_genes <- results_sig %>%
  filter(Class %in% c("Real Alteration", "Possibly Expression-Driven")) %>%
  pull(Gene)

# If you have mut_results from previous MAF comparison, you can add those genes
# mut_results may be loaded from a file; here we skip as it's optional.

# ----------------------------- Boxplots for significant genes ----------------
group_colors <- c("Disease" = "#D6AFB9", "Normal" = "#7E9BB7", "Unknown" = "#8491B4")

for (gene in sig_genes) {
  df <- data.frame(
    Expression = exp_matrix[gene, ],
    Group = factor(group_info$Group[match(colnames(exp_matrix), group_info$Sample)],
                   levels = c("Normal", "Disease", "Unknown")),
    stringsAsFactors = FALSE
  )
  df_test <- df[df$Group %in% c("Disease", "Normal"), ]
  
  # Wilcoxon test
  if (nrow(df_test) >= 3 & length(unique(df_test$Group)) == 2) {
    p_val <- wilcox.test(Expression ~ Group, data = df_test)$p.value
    significance <- case_when(
      p_val < 0.001 ~ "***",
      p_val < 0.01  ~ "**",
      p_val < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  } else {
    significance <- "ns"
  }
  
  max_expr <- max(df$Expression, na.rm = TRUE)
  y_pos <- max_expr + 0.1 * max_expr
  
  p <- ggplot(df, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(width = 0.6, outlier.size = 1.5) +
    scale_fill_manual(values = group_colors) +
    annotate("text", x = 1.5, y = y_pos, label = significance,
             size = 5, fontface = "bold") +
    geom_segment(aes(x = 1, y = y_pos - 0.02 * max_expr,
                     xend = 2, yend = y_pos - 0.02 * max_expr),
                 color = "black", linewidth = 0.5) +
    labs(title = paste("Expression of", gene), x = "Group", y = "Expression Level (RPKM)") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(0.15, "cm"),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 12, face = "bold"),
      plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
      legend.position = "none"
    )
  
  ggsave(filename = file.path(fig_dir, paste0(gene, "_expression_boxplot.pdf")),
         plot = p, width = 3, height = 5, dpi = 300)
}

cat("All mutation vs expression plots saved to:", fig_dir, "\n")