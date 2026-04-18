# =============================================================================
# Downstream Analysis: Risk Score Visualization
# =============================================================================
# This script generates:
#   1. Waterfall plots (bar plots) of Risk Scores for each RNA type
#   2. Violin plots comparing Risk Score distributions between groups
#   3. Scatter plots with group-colored bands
#   4. Multi-class (3-class) Risk Score probability distributions
# =============================================================================

# ----------------------------- Libraries -------------------------------------
library(stringr)
library(ggsci)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)

# ----------------------------- User configuration ----------------------------
# Base directories (relative to project root) - consistent with 07_1
base_dir <- file.path("results", "ML", "Tables")
fig_dir  <- file.path("figures", "ML", "riskscore")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Subdirectories for different analysis types (English names)
dir_2class   <- file.path(base_dir, "Feature_AUC_ACC", "NOR_vs_PH_2025-08-25_Top100DE_baseMean_ge100_dynamic_threshold_class_weight_balanced_30iter")
dir_combined <- file.path(base_dir, "Feature_AUC_ACC", "NOR_vs_PH_2025-08-26_combined_lasso_selected_features_baseMean_ge100")
dir_3class   <- file.path(base_dir, "Feature_AUC_ACC", "PAH_SLE_CHD_vs_each_2025-08-26_3class_noNOR_baseMean_ge100")

# Risk Score subdirectories (parallel to Feature_AUC_ACC)
risk_2class   <- file.path(base_dir, "Risk_Score", "NOR_vs_PH_2025-08-25_Top100DE_baseMean_ge100_dynamic_threshold_class_weight_balanced_30iter")
risk_combined <- file.path(base_dir, "Risk_Score", "NOR_vs_PH_2025-08-26_combined_lasso_selected_features_baseMean_ge100")
risk_3class   <- file.path(base_dir, "Risk_Score", "PAH_SLE_CHD_vs_each_2025-08-26_3class_noNOR_baseMean_ge100")

# Sample ID file
sampleID_file <- file.path("data", "sampleID_filted300W.rds")

# ----------------------------- Helper function ------------------------------
# Read sample IDs
sampleID <- readRDS(sampleID_file)

# ----------------------------- 1. PAH, CHD, SLE vs NOR (multi-group) -------
RNA_type <- c('miRNA','lncRNA','mRNA','snRNA','snoRNA','tRFs','rsRNA','ysRNA')

# 1a. Find best iteration and model for each RNA type (from AUC test results)
result_best <- data.frame()
for (RNA in RNA_type) {
  file_pattern <- paste0("AUC_test_results_", RNA, "_\\d{4}-\\d{1,2}-\\d{1,2}\\.csv")
  matching_files <- list.files(dir_2class, pattern = file_pattern, full.names = FALSE)
  if (length(matching_files) == 0) next
  for (file in matching_files) {
    df <- read.csv(file.path(dir_2class, file))
    if ("Iteration" %in% colnames(df)) df$Iteration <- NULL
    max_val <- max(df, na.rm = TRUE)
    max_pos <- which(df == max_val, arr.ind = TRUE)
    row_names <- rownames(df)[max_pos[, "row"]]
    col_names <- colnames(df)[max_pos[, "col"]]
    date_str <- str_extract(file, "\\d{4}-\\d{1,2}-\\d{1,2}")
    current_result <- data.frame(
      File_Date = date_str,
      Iteration = row_names,
      ML_method = col_names,
      Max_AUC = max_val,
      RNA_type = RNA,
      stringsAsFactors = FALSE
    )
    result_best <- rbind(result_best, current_result)
  }
}
rownames(result_best) <- NULL

# 1b. Read Risk Score files for the best iteration and model
RS_list <- list()
for (RNA in RNA_type) {
  file_pattern <- paste0("RiskScore_All_", RNA, "_\\d{4}-\\d{1,2}-\\d{1,2}\\.csv")
  matching_files <- list.files(risk_2class, pattern = file_pattern, full.names = FALSE)
  if (length(matching_files) == 0) next
  for (file in matching_files) {
    df <- read.csv(file.path(risk_2class, file))
    # Add group information
    df <- df %>%
      rowwise() %>%
      mutate(group = names(sampleID)[which(sapply(sampleID, function(x) Sample %in% x))]) %>%
      ungroup()
    Iter <- result_best[result_best$RNA_type == RNA, "Iteration"]
    Method <- result_best[result_best$RNA_type == RNA, "ML_method"]
    df <- df[df$Iteration == Iter & df$RNA_Type == RNA & df$Model == Method, ]
    RS_list[[RNA]] <- df
  }
}

# 1c. Waterfall plot function
plot_waterfall <- function(df, show_legend = TRUE) {
  title_text <- paste0(unique(df$RNA_Type), "_", unique(df$Model))
  df <- df %>% arrange(desc(Risk_Score)) %>% mutate(Sample = factor(Sample, levels = Sample))
  p1 <- ggplot(df, aes(x = Sample, y = Risk_Score, fill = group)) +
    geom_col(width = 0.8) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = title_text, x = "", y = "Risk Score", fill = "Group") +
    theme_classic() +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = ifelse(show_legend, "right", "none"),
          axis.title.y = element_text(size = 12))
  p2 <- ggplot(df, aes(x = Sample, y = 1, fill = group)) +
    geom_tile() +
    scale_fill_brewer(palette = "Set2") +
    theme_void() +
    theme(legend.position = "none", plot.margin = margin(t = -5))
  combined <- p1 / p2 + plot_layout(heights = c(4, 0.3))
  return(combined)
}

plot_list <- lapply(seq_along(RS_list), function(i) {
  plot_waterfall(RS_list[[i]], show_legend = (i == 1))
})
final_plot <- wrap_plots(plotlist = plot_list, ncol = 4, nrow = 2)
ggsave(final_plot, file = file.path(fig_dir, "RiskScore_PAH_CHD_SLE_vs_NOR.pdf"),
       width = 33.6, height = 12, units = "cm")

# ----------------------------- 2. PH vs NOR (binary) ------------------------
# 2a. Best iteration and model for combined cfRNA (all_lasso_selected_feature)
RNA_comb <- "all_lasso_selected_feature"
result_best_comb <- data.frame()
file_pattern <- paste0("AUC_test_results_", RNA_comb, "_\\d{4}-\\d{1,2}-\\d{1,2}\\.csv")
matching_files <- list.files(dir_combined, pattern = file_pattern, full.names = FALSE)
if (length(matching_files) > 0) {
  for (file in matching_files) {
    df <- read.csv(file.path(dir_combined, file))
    if ("Iteration" %in% colnames(df)) df$Iteration <- NULL
    max_val <- max(df, na.rm = TRUE)
    max_pos <- which(df == max_val, arr.ind = TRUE)
    row_names <- rownames(df)[max_pos[, "row"]]
    col_names <- colnames(df)[max_pos[, "col"]]
    date_str <- str_extract(file, "\\d{4}-\\d{1,2}-\\d{1,2}")
    current_result <- data.frame(
      File_Date = date_str, Iteration = row_names, ML_method = col_names,
      Max_AUC = max_val, RNA_type = RNA_comb, stringsAsFactors = FALSE
    )
    result_best_comb <- rbind(result_best_comb, current_result)
  }
}
# 2b. Read Risk Score for combined cfRNA
RS_comb_list <- list()
file_pattern <- paste0("RiskScore_All_", RNA_comb, "_\\d{4}-\\d{1,2}-\\d{1,2}\\.csv")
matching_files <- list.files(risk_combined, pattern = file_pattern, full.names = FALSE)
if (length(matching_files) > 0) {
  for (file in matching_files) {
    df <- read.csv(file.path(risk_combined, file))
    # Subset to PH vs NOR only (remove CHD/SLE if present)
    sampleID_PH <- sampleID$PH_not4pnk
    sampleID_NOR <- sampleID$NOR_not4pnk
    df <- df[df$Sample %in% c(sampleID_PH, sampleID_NOR), ]
    df <- df %>%
      rowwise() %>%
      mutate(group = ifelse(Sample %in% sampleID_PH, "PH", "NOR")) %>%
      ungroup()
    Iter <- result_best_comb$Iteration
    Method <- result_best_comb$ML_method
    df <- df[df$Iteration == Iter & df$RNA_Type == RNA_comb & df$Model == Method, ]
    RS_comb_list[[RNA_comb]] <- df
  }
}
# 2c. Waterfall plot for PH vs NOR
if (length(RS_comb_list) > 0) {
  plot_list_comb <- lapply(seq_along(RS_comb_list), function(i) {
    plot_waterfall(RS_comb_list[[i]], show_legend = (i == 1))
  })
  final_plot_comb <- wrap_plots(plotlist = plot_list_comb, ncol = 4, nrow = 2)
  ggsave(final_plot_comb, file = file.path(fig_dir, "RiskScore_PH_Vs_NOR.pdf"),
         width = 33.6, height = 12, units = "cm")
}

# ----------------------------- 3. Violin plots for Risk Scores --------------
plot_violin <- function(df, show_legend = TRUE) {
  title_text <- paste0(unique(df$RNA_Type), "_", unique(df$Model))
  # For binary PH vs NOR, groups are "PH" and "NOR"; for multi-group, adjust levels
  df <- df %>% mutate(group = factor(group, levels = c("NOR", "PH", "IPAH_not4pnk", "CHD_not4pnk", "SLE_not4pnk")))
  p <- ggplot(df, aes(x = group, y = Risk_Score, fill = group)) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 1, color = "black", size = 0.5) +
    geom_jitter(aes(color = group), width = 0.15, size = 2, alpha = 1, shape = 16) +
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Set2") +
    labs(title = title_text, x = "", y = "Risk Score", fill = "Group", color = "Group") +
    theme_classic() +
    theme(text = element_text(color = "black"),
          plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
          legend.position = ifelse(show_legend, "right", "none"),
          axis.title.y = element_text(color = "black", size = 12),
          axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(color = "black", size = 12))
  return(p)
}
# Apply to PH vs NOR combined data
if (exists("RS_comb_list") && length(RS_comb_list) > 0) {
  plot_list_violin <- lapply(seq_along(RS_comb_list), function(i) {
    plot_violin(RS_comb_list[[i]], show_legend = (i == 1))
  })
  final_violin <- wrap_plots(plotlist = plot_list_violin, ncol = 4, nrow = 2)
  ggsave(final_violin, file = file.path(fig_dir, "RiskScore_PH_Vs_NOR_violin.pdf"),
         width = 38, height = 24, units = "cm")
}

# ----------------------------- 4. Scatter plot with group bands --------------
plot_scatter_band <- function(df, show_legend = TRUE) {
  title_text <- paste0(unique(df$RNA_Type), "_", unique(df$Model))
  df <- df %>% arrange(desc(Risk_Score)) %>%
    mutate(Sample = factor(Sample, levels = Sample), Index = row_number())
  p <- ggplot(df, aes(x = Index, y = Risk_Score)) +
    geom_rect(data = df %>% group_by(group) %>%
                summarise(xmin = min(Index), xmax = max(Index)),
              aes(xmin = xmin - 0.5, xmax = xmax + 0.5, ymin = -Inf, ymax = Inf, fill = group),
              inherit.aes = FALSE, alpha = 0.15) +
    geom_point(aes(color = group), size = 1.8, alpha = 0.8) +
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Set2") +
    labs(title = title_text, x = "Samples (sorted by Risk Score)", y = "Risk Score",
         color = "Group", fill = "Group") +
    theme_classic() +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.position = ifelse(show_legend, "right", "none"),
          axis.title = element_text(size = 12))
  return(p)
}
if (exists("RS_comb_list") && length(RS_comb_list) > 0) {
  plot_list_band <- lapply(seq_along(RS_comb_list), function(i) {
    plot_scatter_band(RS_comb_list[[i]], show_legend = (i == 1))
  })
  final_band <- wrap_plots(plotlist = plot_list_band, ncol = 4, nrow = 2)
  ggsave(final_band, file = file.path(fig_dir, "RiskScore_PH_Vs_NOR_scatter_band.pdf"),
         width = 33.6, height = 12, units = "cm")
}

# ----------------------------- 5. Multi-class (3-class) Risk Scores ---------
# 5a. Find best iteration and model for 3-class combined cfRNA
# (assuming the combined cfRNA was used for 3-class classification)
# This part uses files from dir_3class and risk_3class
RNA_3c <- "all_lasso_selected_feature"  # or specific name used in 3-class
result_best_3c <- data.frame()
file_pattern <- "AUC_test_allcfRNA_combine_for_3class_\\d{4}-\\d{2}-\\d{2}\\.csv"
matching_files <- list.files(dir_3class, pattern = file_pattern, full.names = TRUE)
if (length(matching_files) > 0) {
  df_auc <- read.csv(matching_files[1])
  # Compute mean AUC per model across 3 classes
  model_names <- c('GLMNETRIDGE', 'GLMNETLASSO', 'SVMLIN', 'SVMRAD', 'RF', 'EXTRATREES',
                   'NNET', 'LDA', 'C5', 'KNN', 'NB', 'RPART', 'GLM', 'GBM')
  mean_auc <- sapply(model_names, function(m) {
    cols <- paste0(m, "_Class_", 1:3)
    mean(rowMeans(df_auc[, cols], na.rm = TRUE), na.rm = TRUE)
  })
  best_model <- names(which.max(mean_auc))
  best_iter <- which.max(rowMeans(df_auc[, paste0(best_model, "_Class_", 1:3)], na.rm = TRUE))
  result_best_3c <- data.frame(Iteration = best_iter, ML_method = best_model, stringsAsFactors = FALSE)
}

# 5b. Read Risk Score for 3-class
if (nrow(result_best_3c) > 0) {
  risk_file <- list.files(risk_3class, pattern = "3class_RiskScores_allModels_allcfRNA_combine_for_3class_.*\\.csv", full.names = TRUE)[1]
  if (!is.na(risk_file)) {
    RiskScores <- read.csv(risk_file)
    RiskScores <- RiskScores[RiskScores$Iteration == result_best_3c$Iteration & 
                               RiskScores$Model == result_best_3c$ML_method, ]
    # Add group information (IPAH, CHD, SLE)
    sampleID_3c <- sampleID[c("IPAH_not4pnk", "CHD_not4pnk", "SLE_not4pnk")]
    RiskScores$group <- NA
    for (gn in names(sampleID_3c)) {
      RiskScores$group[RiskScores$SampleID %in% sampleID_3c[[gn]]] <- gn
    }
    # Convert to long format for probability distribution
    plot_data <- RiskScores %>%
      pivot_longer(cols = contains("Class_"), names_to = "Class", values_to = "Probability")
    
    p_3class <- ggplot(plot_data, aes(x = group, y = Probability, fill = Class)) +
      geom_violin(alpha = 0.7, position = position_dodge(width = 0.6), trim = FALSE,
                  scale = "width", width = 0.6) +
      geom_boxplot(width = 0.1, position = position_dodge(width = 0.6), outlier.shape = NA) +
      scale_x_discrete(limits = c("IPAH_not4pnk", "CHD_not4pnk", "SLE_not4pnk"),
                       expand = expansion(add = 0.4)) +
      scale_fill_manual(values = c("Class_1_prob" = "#FC8D62", "Class_2_prob" = "#8DA0CB",
                                   "Class_3_prob" = "#E78AC3"),
                        labels = c("Class 1", "Class 2", "Class 3")) +
      theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
            plot.background = element_rect(fill = "white", colour = "black", linewidth = 1),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            axis.text = element_text(colour = "black", family = "sans"),
            axis.title = element_text(colour = "black", family = "sans"),
            text = element_text(family = "sans", colour = "black"),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            panel.grid = element_blank(),
            legend.position = "top") +
      labs(title = "Distribution of Class Probabilities by Group",
           x = "Group", y = "Probability", fill = "Class")
    
    ggsave(file.path(fig_dir, "3class_models_riskscore.pdf"),
           plot = p_3class, width = 15, height = 15, units = "cm", device = cairo_pdf)
  }
}

cat("All risk score plots generated in:", fig_dir, "\n")