# =============================================================================
# Survival Analysis: Lasso-Cox and Time-Dependent ROC
# =============================================================================
# This script:
#   1. Reads expression data (cfRNA, mutation, microbe) and clinical metadata
#   2. For each data type (cfRNA, variant, microbe), performs LOOCV to compute
#      risk scores using Lasso-Cox models
#   3. Compares time-dependent AUC (24, 36, 48, 60 months) across modalities
#      and clinical biomarkers (NTproBNP, 6MWT, CI, REVEAL)
#   4. Generates ROC curves and AUC comparison plots
# =============================================================================

# ----------------------------- Libraries -------------------------------------
library(survival)
library(glmnet)
library(timeROC)
library(dplyr)
library(readxl)
library(ggplot2)
library(readr)
library(tidyr)
set.seed(123)

# ----------------------------- User configuration ----------------------------
# Input directories (relative to project root)
data_dir       <- "data"
clinical_file  <- file.path(data_dir, "Meta_clinical_20241208.xlsx")
reveal_file    <- file.path(data_dir, "Meta_REVEAL_241223.csv")
rpkm_file      <- file.path(data_dir, "df_rpkm_300W_allRNA_microbe_variant.csv")

# Output directories
fig_dir        <- file.path("figures", "survival")
res_dir        <- file.path("results", "survival")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)

# ----------------------------- Helper functions ------------------------------
read_id_file <- function(filename) {
  path <- file.path(data_dir, filename)
  if (!file.exists(path)) {
    warning("File not found: ", path)
    return(NULL)
  }
  return(trimws(readLines(path, warn = FALSE)))
}

read_rna_file <- function(filename) {
  path <- file.path(data_dir, filename)
  if (!file.exists(path)) {
    warning("File not found: ", path)
    return(NULL)
  }
  return(trimws(readLines(path, warn = FALSE)))
}

# ----------------------------- 1. Load clinical data ------------------------
# Load survival data (sheet "20241223")
meta <- read_excel(clinical_file, sheet = "20241223")
meta <- meta[, which(colnames(meta) %in% c('Seq_ID', '生存时间2021', '随访结果2021（0=存活，1=死亡，2=失访）'))]
colnames(meta) <- c("SampleID", "Status", "Time")
meta <- meta[!is.na(meta$Time), ]
meta <- meta[!is.na(meta$SampleID), ]
meta$Status <- as.numeric(meta$Status)

# Load clinical biomarkers (NTproBNP, 6MWT, CI)
Meta_clinical <- read_excel(clinical_file, sheet = "Sheet1")
Meta_clinical <- Meta_clinical[, c("Seq_ID", "六mwt", "ntprobnp", "心脏指数")]
Meta_clinical <- na.omit(Meta_clinical)
colnames(Meta_clinical) <- c("SampleID", "SixMWT", "NTproBNP", "CI")
meta <- merge(meta, Meta_clinical, by = "SampleID")

# Load REVEAL scores
Meta_REVEAL <- read_csv(reveal_file, show_col_types = FALSE)
colnames(Meta_REVEAL)[1] <- "SampleID"
colnames(Meta_REVEAL)[13] <- "REVEAL"
Meta_REVEAL <- Meta_REVEAL[, c("SampleID", "REVEAL")]
meta <- merge(meta, Meta_REVEAL, by = "SampleID")

# ----------------------------- 2. Load sample IDs ---------------------------
sampleID <- list(
  CHD_not4pnk  = read_id_file("CHD_not4pnk_ID.txt"),
  SLE_not4pnk  = read_id_file("SLE_not4pnk_ID.txt"),
  NOR_not4pnk  = read_id_file("NOR_not4pnk_ID.txt"),
  IPAH_not4pnk = read_id_file("IPAH_not4pnk_ID.txt"),
  IPAH_t4pnk   = read_id_file("IPAH_t4pnk_ID.txt"),
  NOR_t4pnk    = read_id_file("NOR_t4pnk_ID.txt")
)

# ----------------------------- 3. Load RNA name lists -----------------------
RNAname <- list(
  rsRNA   = read_rna_file("rsRNA.txt"),
  ysRNA   = read_rna_file("ysRNA.txt"),
  tRFs    = read_rna_file("tRFs.txt"),
  miRNA   = read_rna_file("miRNA.txt"),
  mttRNA  = read_rna_file("mttRNA.txt"),
  snRNA   = read_rna_file("snRNA.txt"),
  snoRNA  = read_rna_file("snoRNA.txt"),
  mRNA    = read_rna_file("mRNA.txt"),
  lncRNA  = read_rna_file("lncRNA.txt"),
  variant = read_rna_file("variant.txt"),
  microbe = read_rna_file("microbe.txt")
)

# ----------------------------- 4. Load expression matrix --------------------
df_rpkm <- read.csv(rpkm_file, row.names = 1, check.names = FALSE)
colnames(df_rpkm) <- gsub("\\.", "-", colnames(df_rpkm))

# Keep only samples present in meta
common_samples <- intersect(colnames(df_rpkm), meta$SampleID)
df_rpkm <- df_rpkm[, common_samples]

# ----------------------------- 5. Process expression data per type ----------
process_expr <- function(df, rna_names, type) {
  if (type %in% c("variant", "microbe")) {
    expr <- df[rownames(df) %in% rna_names, ]
    return(expr)
  } else {
    expr <- df[rownames(df) %in% rna_names, ]
    # Filter rows with median >= 10 (optional, adjust as needed)
    row_medians <- apply(expr, 1, median, na.rm = TRUE)
    expr <- expr[row_medians >= 10, ]
    expr <- log2(expr + 1)
    return(expr)
  }
}

expr_mRNA   <- process_expr(df_rpkm, RNAname$mRNA,   "mRNA")
expr_lncRNA <- process_expr(df_rpkm, RNAname$lncRNA, "lncRNA")
expr_miRNA  <- process_expr(df_rpkm, RNAname$miRNA,  "miRNA")
expr_rsRNA  <- process_expr(df_rpkm, RNAname$rsRNA,  "rsRNA")
expr_ysRNA  <- process_expr(df_rpkm, RNAname$ysRNA,  "ysRNA")
expr_tRNA   <- process_expr(df_rpkm, RNAname$tRFs,   "tRFs")
expr_snRNA  <- process_expr(df_rpkm, RNAname$snRNA,  "snRNA")
expr_snoRNA <- process_expr(df_rpkm, RNAname$snoRNA, "snoRNA")
expr_variant<- process_expr(df_rpkm, RNAname$variant,"variant")
expr_microbe<- process_expr(df_rpkm, RNAname$microbe,"microbe")

# Combine cfRNA types
expr_pools <- list(
  cfRNA = rbind(expr_mRNA, expr_lncRNA, expr_miRNA, expr_rsRNA,
                expr_ysRNA, expr_tRNA, expr_snRNA, expr_snoRNA),
  variant = expr_variant,
  microbe = expr_microbe
)

# ----------------------------- 6. Feature selection function ----------------
select_features <- function(expr, meta_sub, top_n = 100, p_val_cutoff = 0.05) {
  common_samples <- intersect(colnames(expr), meta_sub$SampleID)
  expr <- expr[, common_samples]
  meta_sub <- meta_sub[match(common_samples, meta_sub$SampleID), ]
  
  p_values <- apply(expr, 1, function(x) {
    fit <- coxph(Surv(meta_sub$Time, meta_sub$Status) ~ x)
    summary(fit)$coefficients[, "Pr(>|z|)"]
  })
  significant_genes <- names(p_values)[p_values < p_val_cutoff]
  if (length(significant_genes) > 0) {
    sorted_p <- sort(p_values[significant_genes])
    selected <- names(sorted_p)[1:min(top_n, length(sorted_p))]
    return(expr[selected, ])
  } else {
    return(expr[FALSE, ])
  }
}

# ----------------------------- 7. LOOCV Lasso-Cox and timeROC --------------
bootstrap_ci <- function(data, marker, time_point, n_boot = 200) {
  aucs <- numeric(n_boot)
  for (i in 1:n_boot) {
    idx <- sample(seq_len(nrow(data)), replace = TRUE)
    boot_dat <- data[idx, ]
    roc_boot <- timeROC(T = boot_dat$Time, delta = boot_dat$Status,
                        marker = boot_dat[[marker]], cause = 1, times = time_point)
    aucs[i] <- roc_boot$AUC[which(roc_boot$times == time_point)]
  }
  ci <- quantile(aucs, c(0.025, 0.975), na.rm = TRUE)
  return(ci)
}

# Main analysis for a given top_n (here fixed at 30 for brevity, but can loop)
top_n <- 30
output_dir <- file.path(fig_dir, paste0("top_n_", top_n))
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

results <- list()
for (type in c("cfRNA", "variant", "microbe")) {
  cat("Processing type:", type, "\n")
  common_samples <- intersect(colnames(expr_pools[[type]]), meta$SampleID)
  meta_sub <- meta[match(common_samples, meta$SampleID), ]
  n_samples <- nrow(meta_sub)
  
  # LOOCV to get risk scores
  loocv_scores <- numeric(n_samples)
  names(loocv_scores) <- common_samples
  
  for (i in 1:n_samples) {
    train_idx <- (1:n_samples)[-i]
    test_idx <- i
    meta_train <- meta_sub[train_idx, ]
    
    if (type == "cfRNA") {
      s_mRNA   <- select_features(expr_mRNA[, train_idx], meta_train, top_n = top_n)
      s_lncRNA <- select_features(expr_lncRNA[, train_idx], meta_train, top_n = top_n)
      s_miRNA  <- select_features(expr_miRNA[, train_idx], meta_train, top_n = top_n)
      s_rsRNA  <- select_features(expr_rsRNA[, train_idx], meta_train, top_n = top_n)
      s_ysRNA  <- select_features(expr_ysRNA[, train_idx], meta_train, top_n = top_n)
      s_tRFs   <- select_features(expr_tRNA[, train_idx], meta_train, top_n = top_n)
      s_snRNA  <- select_features(expr_snRNA[, train_idx], meta_train, top_n = top_n)
      s_snoRNA <- select_features(expr_snoRNA[, train_idx], meta_train, top_n = top_n)
      expr_train <- rbind(s_mRNA, s_lncRNA, s_miRNA, s_rsRNA, s_ysRNA, s_tRFs, s_snRNA, s_snoRNA)
    } else {
      expr_train <- select_features(expr_pools[[type]][, train_idx], meta_train, top_n = top_n)
    }
    expr_train <- expr_train[complete.cases(expr_train), ]
    genes_loop <- rownames(expr_train)
    if (length(genes_loop) == 0) next
    
    x_train <- t(as.matrix(expr_train))
    y_train <- Surv(meta_train$Time, meta_train$Status)
    fit_loop <- cv.glmnet(x_train, y_train, family = "cox", alpha = 1, nfolds = 10)
    
    x_test <- t(as.matrix(expr_pools[[type]][genes_loop, test_idx, drop = FALSE]))
    loocv_scores[i] <- predict(fit_loop, newx = x_test, s = "lambda.min", type = "link")
    if (i %% 10 == 0) cat(i, " ")
  }
  cat("\n")
  
  # Build risk table
  risk_table <- data.frame(SampleID = names(loocv_scores), risk_score = loocv_scores)
  colnames(risk_table)[2] <- paste0("risk_", type)
  results[[type]] <- risk_table
}

# Merge all risk scores with meta
risk_table_final <- Reduce(function(x, y) merge(x, y, by = "SampleID", all = TRUE), results)
risk_table_final <- merge(risk_table_final, meta, by = "SampleID")
write.csv(risk_table_final, file = file.path(res_dir, "risk_score_table.csv"), row.names = FALSE)

# ----------------------------- 8. Time-dependent ROC ------------------------
time_points <- c(24, 36, 48, 60)
models <- c("risk_cfRNA", "risk_variant", "risk_microbe", "NTproBNP", "SixMWT", "CI", "REVEAL")
model_names <- c("cfRNA", "Variant", "Microbe", "NTproBNP", "6MWT", "CI", "REVEAL")
colors <- c("red", "blue", "green", "purple", "orange", "brown", "pink")

auc_results <- data.frame(Time = time_points)
for (tp in time_points) {
  pdf_file <- file.path(output_dir, paste0("TimeROC_", tp, "m.pdf"))
  pdf(pdf_file, width = 8, height = 8)
  plot(NULL, xlim = c(0,1), ylim = c(0,1),
       xlab = "1 - Specificity", ylab = "Sensitivity",
       main = paste0("Time-dependent ROC at ", tp, " months (top_n=", top_n, ")"))
  abline(0,1, lty=2, col="gray")
  
  legend_text <- c()
  for (i in seq_along(models)) {
    marker <- models[i]
    name <- model_names[i]
    roc_obj <- timeROC(T = risk_table_final$Time, delta = risk_table_final$Status,
                       marker = risk_table_final[[marker]], cause = 1, times = tp)
    lines(roc_obj$FP[,2], roc_obj$TP[,2], col = colors[i], lwd = 2)
    auc_val <- roc_obj$AUC[which(roc_obj$times == tp)]
    ci_val <- bootstrap_ci(risk_table_final, marker, time_point = tp, n_boot = 200)
    legend_text <- c(legend_text, paste0(name, ": AUC=", round(auc_val,3),
                                         " (95%CI ", round(ci_val[1],3), "-",
                                         round(ci_val[2],3), ")"))
    # Store AUC results
    if (tp == 24) auc_results[auc_results$Time == tp, paste0(name, "_AUC")] <- auc_val
    # (simplified: you can expand to store all time points)
  }
  legend("bottomright", legend = legend_text, col = colors, lty = 1, cex = 0.8, bty = "n")
  dev.off()
}

# Also create a combined AUC over time plot
auc_over_time <- data.frame()
for (tp in time_points) {
  for (i in seq_along(models)) {
    roc_obj <- timeROC(T = risk_table_final$Time, delta = risk_table_final$Status,
                       marker = risk_table_final[[models[i]]], cause = 1, times = tp)
    auc_val <- roc_obj$AUC[which(roc_obj$times == tp)]
    auc_over_time <- rbind(auc_over_time, data.frame(Time = tp, Model = model_names[i], AUC = auc_val))
  }
}
p_auc <- ggplot(auc_over_time, aes(x = Time, y = AUC, color = Model)) +
  geom_line(size = 1.2) + geom_point(size = 3) +
  scale_y_continuous(limits = c(0,1)) +
  labs(title = "AUC Comparison Over Time", x = "Time (months)", y = "AUC") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(file.path(output_dir, "AUC_comparison_over_time.pdf"), p_auc, width = 8, height = 6)

cat("Survival analysis completed. Results saved to:", output_dir, "\n")