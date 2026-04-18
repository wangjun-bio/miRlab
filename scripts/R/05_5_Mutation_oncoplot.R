# =============================================================================
# Mutation Oncoprint and Visualization
# =============================================================================
# This script:
#   1. Reads MAF files for each group (CHD, IPAH, NOR, SLE)
#   2. Merges MAF objects and creates oncoPrint (waterfall plot)
#   3. Compares mutation frequencies between PH (CHD+IPAH+SLE) and NOR
#   4. Generates forest plot, co-barplot, and MAF summary plots
# =============================================================================

# ----------------------------- Libraries -------------------------------------
library(maftools)
library(dplyr)
library(ComplexHeatmap)
library(grid)
library(readxl)
library(tidyr)
library(tibble)
library(circlize)

# ----------------------------- User configuration ----------------------------
# Input directories (relative to project root)
maf_dir       <- file.path("results", "maftools")     # contains CHD_DP0_GT11.maf etc.
clinical_file <- file.path("data", "Meta_clinical.xlsx")
sampleID_file <- file.path("data", "sampleID.rds")

# Output directory for figures
fig_dir <- file.path("figures", "mutation")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Mutation types to keep for non-synonymous variants (used in read.maf)
vc_nonSyn <- c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
               "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
               "Silent", "Splice_Site", "Translation_Start_Site")

# ----------------------------- Helper function -------------------------------
check_variant_types <- function(maf_obj, name) {
  cat("Checking variant types in", name, ":\n")
  print(table(maf_obj@data$Variant_Classification))
}

# ----------------------------- Load MAF files --------------------------------
# Read MAF for each group (using DP0_GT11 as an example; adjust if needed)
maf_chd  <- read.maf(maf = file.path(maf_dir, "CHD_DP0_GT11.maf"), vc_nonSyn = vc_nonSyn)
maf_ipah <- read.maf(maf = file.path(maf_dir, "IPAH_DP0_GT11.maf"), vc_nonSyn = vc_nonSyn)
maf_nor  <- read.maf(maf = file.path(maf_dir, "NOR_DP0_GT11.maf"), vc_nonSyn = vc_nonSyn)
maf_sle  <- read.maf(maf = file.path(maf_dir, "SLE_DP0_GT11.maf"), vc_nonSyn = vc_nonSyn)

# Add group labels to clinical data
maf_chd@clinical.data$Group  <- "CHD"
maf_ipah@clinical.data$Group <- "IPAH"
maf_nor@clinical.data$Group  <- "NOR"
maf_sle@clinical.data$Group  <- "SLE"

# Check variant type distributions
check_variant_types(maf_chd, "CHD")
check_variant_types(maf_ipah, "IPAH")
check_variant_types(maf_nor, "NOR")
check_variant_types(maf_sle, "SLE")

# Merge all MAF objects
combined_maf <- merge_mafs(mafs = list(maf_chd, maf_ipah, maf_nor, maf_sle), vc_nonSyn = vc_nonSyn)

# Merge PH (CHD+IPAH+SLE) for comparison with NOR
maf_ph <- merge_mafs(mafs = list(maf_chd, maf_ipah, maf_sle), vc_nonSyn = vc_nonSyn)

# Save combined MAF data as CSV (optional)
write.csv(as.data.frame(combined_maf@data),
          file = file.path(fig_dir, "TableS2_combined_mutations.csv"),
          row.names = FALSE)

# ----------------------------- Compare PH vs NOR -----------------------------
compare_results <- mafCompare(
  m1 = maf_ph,
  m2 = maf_nor,
  m1Name = "PH",
  m2Name = "NOR",
  minMut = 3
)

# View significantly different genes
mut_results <- compare_results$results
write.csv(mut_results, file = file.path(fig_dir, "mut_diff_PH_vs_NOR.csv"), row.names = FALSE)

# Forest plot of mutation odds ratios
pdf(file.path(fig_dir, "DE_forest_plot.pdf"), width = 10, height = 8)
forestPlot(
  mafCompareRes = compare_results,
  pVal = 0.05,
  color = c('red', 'blue'),
  geneFontSize = 0.8
)
dev.off()

# Co-barplot of top altered genes
coBarplot(
  m1 = maf_nor,
  m2 = maf_ph,
  m1Name = "NOR",
  m2Name = "PH",
  genes = compare_results$results[compare_results$results$pval < 0.05, ]$Hugo_Symbol
)

# ----------------------------- OncoPrint (waterfall plot) --------------------
# Helper function to create mutation matrix for oncoPrint
maf_to_oncomatrix <- function(maf, top_n = 20) {
  maf_df <- maf@data
  top_genes <- maf_df %>%
    count(Hugo_Symbol, name = "n") %>%
    arrange(desc(n)) %>%
    slice_head(n = top_n) %>%
    pull(Hugo_Symbol)
  
  mut_matrix <- maf_df %>%
    filter(Hugo_Symbol %in% top_genes) %>%
    dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification) %>%
    group_by(Hugo_Symbol, Tumor_Sample_Barcode) %>%
    summarise(Mutation = paste(unique(as.character(Variant_Classification)), collapse = ";"),
              .groups = "drop") %>%
    pivot_wider(names_from = Tumor_Sample_Barcode, values_from = Mutation, values_fill = "") %>%
    tibble::column_to_rownames("Hugo_Symbol")
  return(mut_matrix)
}

# Generate mutation matrix for top 20 genes
mut_matrix <- maf_to_oncomatrix(combined_maf, top_n = 20)

# Prepare sample annotation (clinical data)
sampleID <- readRDS(sampleID_file)
meta_clinical <- read_excel(clinical_file)

# Function to get sample type (CHD/IPAH/NOR/SLE)
get_sample_type <- function(col_name) {
  for (type in names(sampleID)) {
    if (col_name %in% sampleID[[type]]) return(type)
  }
  return(NA)
}

# Function to get clinical info (age, gender, NYHA, etc.)
get_clinical_info <- function(col_name, info_type) {
  match_idx <- which(meta_clinical$Seq_ID == col_name)
  if (length(match_idx) > 0) return(meta_clinical[[info_type]][match_idx[1]])
  return(NA)
}

# Build annotation data frame
sample_anno <- data.frame(
  Group = sapply(colnames(mut_matrix), get_sample_type),
  Gender = sapply(colnames(mut_matrix), get_clinical_info, info_type = "性别"),
  NYHA_Class = sapply(colnames(mut_matrix), get_clinical_info, info_type = "心功能分级"),
  Age = as.numeric(sapply(colnames(mut_matrix), get_clinical_info, info_type = "年龄")),
  SixMWT = as.numeric(sapply(colnames(mut_matrix), get_clinical_info, info_type = "六mwt")),
  mPAP = as.numeric(sapply(colnames(mut_matrix), get_clinical_info, info_type = "肺动脉平均压")),
  CI = as.numeric(sapply(colnames(mut_matrix), get_clinical_info, info_type = "心脏指数")),
  Ntprobnp = as.numeric(sapply(colnames(mut_matrix), get_clinical_info, info_type = "ntprobnp")),
  row.names = colnames(mut_matrix),
  check.names = FALSE
)

# Define annotation colors
discrete_colors <- list(
  Group = c("CHD_not4pnk" = "#E41A1C", "IPAH_not4pnk" = "#377EB8",
            "NOR_not4pnk" = "#4DAF4A", "SLE_not4pnk" = "#984EA3"),
  Gender = c("1" = "#66C2A5", "2" = "#FC8D62"),
  NYHA_Class = c("1" = "#8DA0CB", "2" = "#FD8D3C", "3" = "#FC4E2A", "4" = "#E41A1C")
)

continuous_colors <- list(
  Age = colorRamp2(c(min(sample_anno$Age, na.rm = TRUE),
                     median(sample_anno$Age, na.rm = TRUE),
                     max(sample_anno$Age, na.rm = TRUE)),
                   c("#F2F0F7", "#9E9AC8", "#4A1486")),
  SixMWT = colorRamp2(c(min(sample_anno$SixMWT, na.rm = TRUE),
                        median(sample_anno$SixMWT, na.rm = TRUE),
                        max(sample_anno$SixMWT, na.rm = TRUE)),
                      c("#FEE0D2", "#FC9272", "#DE2D26")),
  mPAP = colorRamp2(c(min(sample_anno$mPAP, na.rm = TRUE),
                      median(sample_anno$mPAP, na.rm = TRUE),
                      max(sample_anno$mPAP, na.rm = TRUE)),
                    c("#E5F5E0", "#A1D99B", "#006D2C")),
  CI = colorRamp2(c(min(sample_anno$CI, na.rm = TRUE),
                    median(sample_anno$CI, na.rm = TRUE),
                    max(sample_anno$CI, na.rm = TRUE)),
                  c("#DEEBF7", "#9ECAE1", "#3182BD")),
  Ntprobnp = colorRamp2(c(min(sample_anno$Ntprobnp, na.rm = TRUE),
                          median(sample_anno$Ntprobnp, na.rm = TRUE),
                          max(sample_anno$Ntprobnp, na.rm = TRUE)),
                        c("#FFF5F0", "#FC8D59", "#D73027"))
)
anno_colors <- c(discrete_colors, continuous_colors)

# Define alter_fun for oncoPrint
actual_types <- unique(unlist(strsplit(as.matrix(mut_matrix), ";")))
actual_types <- actual_types[actual_types != ""]
col <- setNames(c("#336499", "#CA9832", "#C90030", "#984EA3", "#EE8BB9",
                  "#FFFF33", "#A65628", "#F781BF", "#00969B", "#66C2A5"),
                c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Ins",
                  "Frame_Shift_Del", "Splice_Site", "In_Frame_Ins", "In_Frame_Del",
                  "Nonstop_Mutation", "Translation_Start_Site", "Other_Mutation"))
for (atype in actual_types) {
  if (!atype %in% names(col)) col[atype] <- "#CCCCCC"
}

alter_fun <- list(background = function(x, y, w, h)
  grid.rect(x, y, w, h, gp = gpar(fill = "#F0F0F0", col = NA)))
for (atype in actual_types) {
  alter_fun[[atype]] <- local({
    this_color <- col[[atype]]
    function(x, y, w, h) {
      grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = this_color, col = NA))
    }
  })
}

# Top annotation
ha_top <- HeatmapAnnotation(
  df = sample_anno,
  col = anno_colors,
  which = "column",
  show_legend = TRUE,
  annotation_name_side = "left"
)

# Draw oncoPrint
pdf(file.path(fig_dir, "oncoprint.pdf"), width = 36, height = 10)
p <- oncoPrint(
  as.matrix(mut_matrix),
  alter_fun = alter_fun,
  col = col,
  top_annotation = ha_top,
  remove_empty_columns = TRUE,
  remove_empty_rows = TRUE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  column_title = "Mutation Landscape"
)
print(p)
dev.off()

# ----------------------------- Additional MAF summary plots -----------------
# MAF summary (variant counts per sample, classification)
pdf(file.path(fig_dir, "maf_summary.pdf"), width = 14, height = 10)
plotmafSummary(maf = combined_maf, rmOutlier = TRUE, addStat = 'median',
               dashboard = TRUE, titvRaw = FALSE)
dev.off()

# Ti/Tv ratio
pdf(file.path(fig_dir, "titv_plot.pdf"), width = 12, height = 10)
laml.titv <- titv(maf = combined_maf, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)
dev.off()

# Somatic interactions (co-occurrence/exclusivity)
pdf(file.path(fig_dir, "somatic_interactions.pdf"), width = 12, height = 10)
somaticInteractions(maf = combined_maf, top = 25, pvalue = c(0.05, 0.1))
dev.off()

# Extract dbSNP IDs for significantly different genes
top_genes <- mut_results[mut_results$pval < 0.05, ]$Hugo_Symbol
result_dbSNP <- lapply(top_genes, function(gene) {
  unique(combined_maf@data[combined_maf@data$Hugo_Symbol == gene, "dbSNP"])
})
names(result_dbSNP) <- top_genes
print(head(result_dbSNP, 31))

cat("All mutation oncoplot and summary plots saved to:", fig_dir, "\n")