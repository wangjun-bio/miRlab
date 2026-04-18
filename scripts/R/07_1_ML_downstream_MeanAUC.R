# =============================================================================
# Downstream Analysis: Mean AUC Visualization for Machine Learning Models
# =============================================================================
# This script generates:
#   1. Mean AUC plots (train vs test) for a single RNA type
#   2. Faceted mean AUC plots for all RNA types
#   3. Mean AUC plots for mutation data
#   4. Mean AUC plots for combined cfRNA features (Lasso-selected)
#   5. Multi-class (3-class) AUC processing and plotting (mean and per-class)
#   6. Comparison of AUC across cfRNA, mutation, and microbe modalities
# =============================================================================

# ----------------------------- Libraries -------------------------------------
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(openxlsx)

# ----------------------------- User configuration ----------------------------
# Base directories (relative to project root)
base_dir <- file.path("results", "ML", "Tables")
fig_dir  <- file.path("figures", "ML")

# Subdirectories for different analysis types (all in English)
dir_2class       <- file.path(base_dir, "Feature_AUC_ACC", "NOR_vs_PH_2025-08-25_Top100DE_baseMean_ge100_dynamic_threshold_class_weight_balanced_30iter")
dir_3class       <- file.path(base_dir, "Feature_AUC_ACC", "PAH_SLE_CHD_vs_each_2025-08-26_3class_noNOR_baseMean_ge100")
dir_combined     <- file.path(base_dir, "Feature_AUC_ACC", "NOR_vs_PH_2025-08-26_combined_lasso_selected_features_baseMean_ge100")
dir_mutation     <- file.path(base_dir, "..", "Gene_mutation", "Tables", "Feature_AUC_ACC", "mutation_lasso_stability_0.5_fixed_train_threshold_2025-08-06")
dir_microbe      <- file.path(base_dir, "..", "microbe", "Tables", "Feature_AUC_ACC")
dir_multimodal   <- file.path(base_dir, "Feature_AUC_ACC", "mutation_microbe_cfRNA_integrated_NOR_vs_PH_cfRNA_baseMean_ge100")

# Output figure directory
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# ----------------------------- Helper functions ------------------------------
plot_mean_auc <- function(df_train, df_test, output_file, ylim = c(0.6, 1.0), xlab_angle = 0) {
  # df_train, df_test: data frames with columns as models, rows as iterations
  df_train <- df_train %>% mutate(Run = 1:n())
  df_test  <- df_test %>% mutate(Run = 1:n())
  if ("Iteration" %in% colnames(df_train)) df_train$Iteration <- NULL
  if ("Iteration" %in% colnames(df_test)) df_test$Iteration <- NULL
  
  df_train_long <- pivot_longer(df_train, -Run, names_to = "Model", values_to = "AUC") %>% mutate(Type = "Train")
  df_test_long  <- pivot_longer(df_test,  -Run, names_to = "Model", values_to = "AUC") %>% mutate(Type = "Test")
  df_all <- bind_rows(df_train_long, df_test_long)
  
  df_avg <- df_all %>% group_by(Model, Type) %>% summarise(mean_AUC = mean(AUC), .groups = "drop")
  model_order <- df_avg %>% filter(Type == "Test") %>% arrange(mean_AUC) %>% pull(Model)
  df_avg$Model <- factor(df_avg$Model, levels = model_order)
  
  shape_map <- c("Train" = 1, "Test" = 2)
  color_map <- c("Train" = "#0000FF", "Test" = "#FF0000")
  
  p <- ggplot(df_avg, aes(x = Model, y = mean_AUC, shape = Type, color = Type)) +
    geom_point(size = 3, stroke = 1.2) +
    scale_shape_manual(values = shape_map) +
    scale_color_manual(values = color_map) +
    coord_flip() +
    scale_y_continuous(limits = ylim, breaks = seq(ylim[1], ylim[2], by = 0.2)) +
    labs(x = NULL, y = "Mean AUC", shape = "Dataset", color = "Dataset") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 12, color = "black"),
          axis.text.x = element_text(size = 12, color = "black", angle = xlab_angle),
          axis.title.x = element_text(size = 14),
          legend.position = "top")
  ggsave(output_file, p, width = 12, height = 24, units = "cm")
  return(p)
}

# ----------------------------- 1. Single RNA type (example: tRFs) -----------
# (Optional, can be commented out if not needed)
# df_train <- read_csv(file.path(dir_2class, "AUC_train_results_tRFs_2025-08-06.csv"))
# df_test  <- read_csv(file.path(dir_2class, "AUC_test_results_tRFs_2025-08-06.csv"))
# plot_mean_auc(df_train, df_test, file.path(fig_dir, "AUC_tRFs_test_train.pdf"))

# ----------------------------- 2. All RNA types (2-class, NOR vs PH) ---------
group_names <- c("rsRNA", "ysRNA", "miRNA", "snRNA", "mRNA", "tRFs", "lncRNA", "snoRNA")

all_data <- list()
for (group in group_names) {
  train_pattern <- paste0("AUC_train_results_", group, "_.+\\.csv$")
  test_pattern  <- paste0("AUC_test_results_", group, "_.+\\.csv$")
  train_file <- list.files(path = dir_2class, pattern = train_pattern, full.names = TRUE)[1]
  test_file  <- list.files(path = dir_2class, pattern = test_pattern, full.names = TRUE)[1]
  if (is.na(train_file) || is.na(test_file)) next
  
  df_train <- read_csv(train_file, show_col_types = FALSE) %>%
    select(-contains("Iteration")) %>%
    mutate(Run = 1:n(), Group = group, Type = "Train") %>%
    pivot_longer(-c(Run, Group, Type), names_to = "Model", values_to = "AUC")
  df_test <- read_csv(test_file, show_col_types = FALSE) %>%
    select(-contains("Iteration")) %>%
    mutate(Run = 1:n(), Group = group, Type = "Test") %>%
    pivot_longer(-c(Run, Group, Type), names_to = "Model", values_to = "AUC")
  all_data[[group]] <- bind_rows(df_train, df_test)
}
df_all <- bind_rows(all_data)

# Save raw data for supplementary table
write.xlsx(df_all, file = file.path(fig_dir, "Table_S5_cfRNA_AUC_raw.xlsx"))

df_avg <- df_all %>% group_by(Group, Model, Type) %>% summarise(mean_AUC = mean(AUC), .groups = "drop")
model_order <- df_avg %>% filter(Type == "Test") %>% group_by(Model) %>%
  summarise(mean_overall = mean(mean_AUC)) %>% arrange(desc(mean_overall)) %>% pull(Model)
df_avg$Model <- factor(df_avg$Model, levels = model_order)

shape_map <- c("Train" = 1, "Test" = 2)
color_map <- c("Train" = "#0000FF", "Test" = "#FF0000")

p <- ggplot(df_avg, aes(x = Model, y = mean_AUC, shape = Type, color = Type)) +
  geom_point(size = 3, stroke = 1.2) +
  facet_wrap(~Group, nrow = 1, ncol = 8) +
  coord_flip() +
  scale_y_continuous(limits = c(0.6, 1.0), breaks = c(0.6, 0.8, 1.0)) +
  scale_shape_manual(values = shape_map) +
  scale_color_manual(values = color_map) +
  labs(x = NULL, y = "AUC", shape = "Dataset", color = "Dataset") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12),
        legend.position = "top",
        panel.spacing = unit(0.5, "cm"))
ggsave(file.path(fig_dir, "AUC_allRNA_test_train.pdf"), p, width = 33.6, height = 12, units = "cm")

# ----------------------------- 3. Mutation data (2-class) -------------------
train_file <- list.files(path = dir_mutation, pattern = "AUC_train_results_.*\\.csv", full.names = TRUE)[1]
test_file  <- list.files(path = dir_mutation, pattern = "AUC_test_results_.*\\.csv", full.names = TRUE)[1]
if (!is.na(train_file) && !is.na(test_file)) {
  df_train <- read_csv(train_file, show_col_types = FALSE)
  df_test  <- read_csv(test_file, show_col_types = FALSE)
  # Save raw data for supplementary table
  df_all_mut <- bind_rows(
    df_train %>% mutate(Run = 1:n(), Type = "Train") %>% pivot_longer(-c(Run, Type), names_to = "Model", values_to = "AUC"),
    df_test  %>% mutate(Run = 1:n(), Type = "Test")  %>% pivot_longer(-c(Run, Type), names_to = "Model", values_to = "AUC")
  )
  write.xlsx(df_all_mut, file = file.path(fig_dir, "Table_S5_mutation_AUC_raw.xlsx"))
  plot_mean_auc(df_train, df_test, file.path(fig_dir, "AUC_mutation_test_train.pdf"), ylim = c(0.6, 1.0))
}

# ----------------------------- 4. Combined cfRNA (Lasso-selected features) --
train_file <- list.files(path = dir_combined, pattern = "AUC_train_results_all_lasso_selected_feature_.*\\.csv", full.names = TRUE)[1]
test_file  <- list.files(path = dir_combined, pattern = "AUC_test_results_all_lasso_selected_feature_.*\\.csv", full.names = TRUE)[1]
if (!is.na(train_file) && !is.na(test_file)) {
  df_train <- read_csv(train_file, show_col_types = FALSE)
  df_test  <- read_csv(test_file, show_col_types = FALSE)
  df_all_comb <- bind_rows(
    df_train %>% mutate(Run = 1:n(), Type = "Train") %>% pivot_longer(-c(Run, Type), names_to = "Model", values_to = "AUC"),
    df_test  %>% mutate(Run = 1:n(), Type = "Test")  %>% pivot_longer(-c(Run, Type), names_to = "Model", values_to = "AUC")
  )
  write.xlsx(df_all_comb, file = file.path(fig_dir, "Table_S6_combined_cfRNA_AUC_raw.xlsx"))
  plot_mean_auc(df_train, df_test, file.path(fig_dir, "AUC_RNA_combine_test_train.pdf"), ylim = c(0.6, 1.0))
}

# ----------------------------- 5. Multi-class (3-class) AUC processing ------
# 5a. Calculate average AUC across three classes for each model
process_3class_auc <- function(input_dir, output_dir_fig) {
  file_pattern <- "AUC_(train|test)_DE13_(miRNA|tRFs|mRNA|lncRNA|snRNA|snoRNA|rsRNA|ysRNA)_\\d{4}-\\d{2}-\\d{2}\\.csv"
  input_files <- list.files(path = input_dir, pattern = file_pattern, full.names = TRUE)
  model_names <- c('GLMNETRIDGE', 'GLMNETLASSO', 'SVMLIN', 'SVMRAD', 'RF', 'EXTRATREES',
                   'NNET', 'LDA', 'C5', 'KNN', 'NB', 'RPART', 'GLM', 'GBM')
  columns_to_drop <- c()
  for (m in model_names) {
    columns_to_drop <- c(columns_to_drop, paste0(m, "_Class_1"), paste0(m, "_Class_2"), paste0(m, "_Class_3"))
  }
  
  for (input_file in input_files) {
    df <- read.csv(input_file)
    # Calculate average AUC across three classes (equal weight)
    for (m in model_names) {
      cols <- paste0(m, "_Class_", 1:3)
      df[[m]] <- rowMeans(df[, cols], na.rm = TRUE)
    }
    df <- df[, !(names(df) %in% columns_to_drop)]
    # Save averaged file
    out_name <- gsub("\\.csv$", "_average.csv", basename(input_file))
    write.csv(df, file.path(input_dir, out_name), row.names = FALSE)
  }
}
process_3class_auc(dir_3class, fig_dir)

# 5b. Plot faceted mean AUC for 3-class (average across classes)
group_names_3c <- c("rsRNA", "ysRNA", "miRNA", "snRNA", "mRNA", "tRFs", "lncRNA", "snoRNA")
all_data_3c <- list()
for (group in group_names_3c) {
  train_pattern <- paste0("AUC_train_DE13_", group, "_average_.+\\.csv$")
  test_pattern  <- paste0("AUC_test_DE13_", group, "_average_.+\\.csv$")
  train_file <- list.files(path = dir_3class, pattern = train_pattern, full.names = TRUE)[1]
  test_file  <- list.files(path = dir_3class, pattern = test_pattern, full.names = TRUE)[1]
  if (is.na(train_file) || is.na(test_file)) next
  
  df_train <- read_csv(train_file, show_col_types = FALSE) %>%
    select(-contains("Iteration")) %>%
    mutate(Run = 1:n(), Group = group, Type = "Train") %>%
    pivot_longer(-c(Run, Group, Type), names_to = "Model", values_to = "AUC")
  df_test <- read_csv(test_file, show_col_types = FALSE) %>%
    select(-contains("Iteration")) %>%
    mutate(Run = 1:n(), Group = group, Type = "Test") %>%
    pivot_longer(-c(Run, Group, Type), names_to = "Model", values_to = "AUC")
  all_data_3c[[group]] <- bind_rows(df_train, df_test)
}
df_all_3c <- bind_rows(all_data_3c)
df_avg_3c <- df_all_3c %>% group_by(Group, Model, Type) %>% summarise(mean_AUC = mean(AUC), .groups = "drop")
model_order_3c <- df_avg_3c %>% filter(Type == "Test") %>% group_by(Model) %>%
  summarise(mean_overall = mean(mean_AUC)) %>% arrange(desc(mean_overall)) %>% pull(Model)
df_avg_3c$Model <- factor(df_avg_3c$Model, levels = model_order_3c)

p_3c <- ggplot(df_avg_3c, aes(x = Model, y = mean_AUC, shape = Type, color = Type)) +
  geom_point(size = 3, stroke = 1.2) +
  facet_wrap(~Group, nrow = 1, ncol = 8) +
  coord_flip() +
  scale_y_continuous(limits = c(0.5, 1.0), breaks = c(0.6, 0.8, 1.0)) +
  scale_shape_manual(values = shape_map) +
  scale_color_manual(values = color_map) +
  labs(x = NULL, y = "AUC", shape = "Dataset", color = "Dataset") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12),
        legend.position = "top",
        panel.spacing = unit(0.5, "cm"))
ggsave(file.path(fig_dir, "AUC_multiclass_allRNA_test_train.pdf"), p_3c, width = 33.6, height = 12, units = "cm")

# 5c. Per-class AUC plots (separate for Class 1,2,3)
for (cls in 1:3) {
  df_current <- df_all_3c %>%
    filter(str_detect(Model, "_Class_")) %>%
    separate(Model, into = c("Model", "Class"), sep = "_Class_") %>%
    filter(Class == cls) %>%
    select(-Class)
  df_avg_cls <- df_current %>% group_by(Group, Model, Type) %>% summarise(mean_AUC = mean(AUC), .groups = "drop")
  model_order_cls <- df_avg_cls %>% filter(Type == "Test") %>% group_by(Model) %>%
    summarise(mean_overall = mean(mean_AUC)) %>% arrange(desc(mean_overall)) %>% pull(Model)
  df_avg_cls$Model <- factor(df_avg_cls$Model, levels = model_order_cls)
  
  p_cls <- ggplot(df_avg_cls, aes(x = Model, y = mean_AUC, shape = Type, color = Type)) +
    geom_point(size = 3, stroke = 1.2) +
    facet_wrap(~Group, nrow = 1, ncol = 8) +
    coord_flip() +
    scale_y_continuous(limits = c(0.5, 1.1), breaks = seq(0.6, 1.1, 0.2)) +
    scale_color_manual(values = c("Train" = "#0000FF", "Test" = "#FF0000")) +
    scale_shape_manual(values = c("Train" = 1, "Test" = 2)) +
    labs(x = "Machine Learning Models", y = "Average AUC", shape = "Dataset", color = "Dataset",
         title = paste0("Class ", cls, " AUC Comparison")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          strip.text.x = element_text(size = 10, face = "bold"),
          axis.text.y = element_text(size = 8, color = "black"),
          axis.text.x = element_text(size = 10, color = "black"),
          axis.title.x = element_text(size = 12),
          legend.position = "top")
  ggsave(file.path(fig_dir, paste0("AUC_Class", cls, "_AllRNA.pdf")), p_cls, width = 33.6, height = 12, units = "cm")
}

# ----------------------------- 6. Multi-modal comparison (cfRNA, mutation, microbe) ----
multi_modalities <- c("mutation", "microbe", "all_lasso_selected_feature")
all_data_multi <- list()
for (mod in multi_modalities) {
  train_pattern <- paste0("AUC_train_results_", mod, "_.+\\.csv$")
  test_pattern  <- paste0("AUC_test_results_", mod, "_.+\\.csv$")
  train_file <- list.files(path = dir_multimodal, pattern = train_pattern, full.names = TRUE)[1]
  test_file  <- list.files(path = dir_multimodal, pattern = test_pattern, full.names = TRUE)[1]
  if (is.na(train_file) || is.na(test_file)) next
  
  df_train <- read_csv(train_file, show_col_types = FALSE) %>%
    select(-contains("Iteration")) %>%
    mutate(Run = 1:n(), Group = mod, Type = "Train") %>%
    pivot_longer(-c(Run, Group, Type), names_to = "Model", values_to = "AUC")
  df_test <- read_csv(test_file, show_col_types = FALSE) %>%
    select(-contains("Iteration")) %>%
    mutate(Run = 1:n(), Group = mod, Type = "Test") %>%
    pivot_longer(-c(Run, Group, Type), names_to = "Model", values_to = "AUC")
  all_data_multi[[mod]] <- bind_rows(df_train, df_test)
}
df_all_multi <- bind_rows(all_data_multi)
# Save microbe raw data for supplementary table
df_microbe <- df_all_multi %>% filter(Group == "microbe")
write.xlsx(df_microbe, file = file.path(fig_dir, "Table_S5_microbe_AUC_raw.xlsx"))

df_avg_multi <- df_all_multi %>% group_by(Group, Model, Type) %>% summarise(mean_AUC = mean(AUC), .groups = "drop")
model_order_multi <- df_avg_multi %>% filter(Type == "Test") %>% group_by(Model) %>%
  summarise(mean_overall = mean(mean_AUC)) %>% arrange(desc(mean_overall)) %>% pull(Model)
df_avg_multi$Model <- factor(df_avg_multi$Model, levels = model_order_multi)

p_multi <- ggplot(df_avg_multi, aes(x = Model, y = mean_AUC, shape = Type, color = Type)) +
  geom_point(size = 3, stroke = 1.2) +
  facet_wrap(~Group, nrow = 1, ncol = 3) +
  coord_flip() +
  scale_y_continuous(limits = c(0.6, 1.0), breaks = c(0.6, 0.8, 1.0)) +
  scale_shape_manual(values = shape_map) +
  scale_color_manual(values = color_map) +
  labs(x = NULL, y = "AUC", shape = "Dataset", color = "Dataset") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12),
        legend.position = "top",
        panel.spacing = unit(0.5, "cm"))
ggsave(file.path(fig_dir, "AUC_allRNA_microbe_mutation_test_train.pdf"), p_multi, width = 18, height = 12, units = "cm")

cat("All mean AUC plots generated successfully.\n")