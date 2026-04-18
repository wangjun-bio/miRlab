# =============================================================================
# Microbiome QC and Diversity Analysis
# =============================================================================
# This script:
#   1. Reads KrakenUniq genus-level reports
#   2. Removes contaminants (pre-identified)
#   3. Performs beta diversity (PCoA with Bray-Curtis)
#   4. Performs alpha diversity (Shannon, Simpson, evenness)
#   5. Creates stacked barplot of top genera
#   6. Differential abundance analysis (Wilcoxon test, Maaslin2)
# =============================================================================

# ----------------------------- Libraries -------------------------------------
library(dplyr)
library(data.table)
library(pheatmap)
library(vegan)
library(ggplot2)
library(tidyverse)
library(phyloseq)
library(permute)
library(lattice)
library(ggpubr)
library(readxl)
library(sva)
library(ggrepel)
library(magrittr)
library(reshape2)
library(ggsci)
library(ggprism)
library(rstatix)
library(ggpmisc)
library(limma)
library(tidyr)
library(edgeR)
library(metagenomeSeq)
library(Maaslin2)
library(scales)

# ----------------------------- User configuration ----------------------------
# Input directories (relative to project root)
kraken_dir    <- file.path("data", "krakenuniq", "20250707")
contam_file   <- file.path("results", "microbiome", "Microbe_contamination_df_0816.txt")
meta_file     <- file.path("data", "Meta_clinical_20241208.xlsx")
sampleID_file <- file.path("data", "sampleID_filted300W.rds")

# Output directories
fig_dir   <- file.path("figures", "microbiome")
res_dir   <- file.path("results", "microbiome")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE)

# ----------------------------- 1. Read KrakenUniq reports --------------------
files <- list.files(kraken_dir, pattern = "_reportfile.tsv", full.names = TRUE)

# Load contamination list (precomputed from decontamination step)
contamination <- read.delim(contam_file, header = TRUE)$rowname

result <- lapply(files, function(x) {
  tmpname <- gsub("_reportfile.tsv", "", basename(x))
  tmp <- read.table(x, header = TRUE, sep = "\t")
  tmp <- tmp[tmp$rank == "genus", ]
  tmp$taxName <- gsub(" ", "", tmp$taxName)
  tmp <- tmp[!(tmp$taxName %in% contamination), ]
  tmp <- tmp[, c(9, 2)]
  colnames(tmp)[2] <- tmpname
  return(tmp)
})

df_raw_microbeRNA <- Reduce(function(x, y) merge(x, y, by = "taxName", all = TRUE), result)
rownames(df_raw_microbeRNA) <- df_raw_microbeRNA[, 1]
df_raw_microbeRNA <- df_raw_microbeRNA[, -1]
df_raw_microbeRNA <- apply(df_raw_microbeRNA, 2, function(row) { row[is.na(row)] <- 0; row })
df_raw_microbeRNA <- as.data.frame(df_raw_microbeRNA)

# ----------------------------- 2. Read sample metadata ----------------------
files_meta <- list.files(file.path("data", "sampleID_txt"), pattern = "not4pnk_ID.txt", full.names = TRUE)
meta_list <- lapply(files_meta, function(x) {
  tmp <- fread(x, header = FALSE)
  group <- gsub("_not4pnk_ID.txt", "", basename(x))
  tmp$group <- group
  return(tmp)
})
meta_df <- do.call(bind_rows, meta_list)

# Keep only samples present in microbiome matrix
df_raw_microbeRNA <- df_raw_microbeRNA[, colnames(df_raw_microbeRNA) %in% meta_df$V1]

# Load clinical metadata
Meta <- read_excel(meta_file)
sampleID <- readRDS(sampleID_file)
Meta_filtered <- Meta[Meta$Seq_ID %in% sampleID$PH_not4pnk | Meta$Seq_ID %in% sampleID$NOR_not4pnk, ]

# Batch information for ComBat
Meta_filtered <- Meta_filtered %>%
  mutate(batch = case_when(
    `出库批次` == "2022.5.12" ~ 1,
    `出库批次` == "2022.11.11" ~ 2,
    `出库批次` == "2023.9.25" ~ 3,
    TRUE ~ NA_real_
  ))

# ----------------------------- 3. Beta diversity (PCoA) ----------------------
rel_abundance <- apply(df_raw_microbeRNA, 2, function(x) x / sum(x))
sample_ids <- colnames(rel_abundance)
batch_info <- Meta_filtered$batch[match(sample_ids, Meta_filtered$Seq_ID)]
batch_factor <- factor(batch_info)

# Remove batch effect using ComBat
rel_abundance_rmbatch <- ComBat(dat = rel_abundance, batch = batch_info)

# Calculate Bray-Curtis distance
RS_RF <- t(rel_abundance_rmbatch) %>% as.data.frame()
dataT <- na.omit(RS_RF)
dataT <- dataT[rowSums(dataT, na.rm = TRUE) > 0, , drop = FALSE]
dist <- vegdist(dataT, method = "bray")
dist <- as.matrix(dist)

# PCoA
pcoa <- cmdscale(dist, k = 3, eig = TRUE)
pc12 <- pcoa$points[, 1:2]
pc <- round(pcoa$eig / sum(pcoa$eig) * 100, digits = 2)
pc12 <- as.data.frame(pc12)
colnames(pc12) <- c("pc_x", "pc_y")
pc12$sample <- rownames(pc12)

# Map groups
sd <- data.frame(sample = rownames(RS_RF))
index <- match(sd$sample, meta_df$V1)
sd$group <- meta_df$group[index]
pc12 <- merge(pc12, sd, by = "sample")
pc12$group <- factor(pc12$group, levels = unique(meta_df$group))

# Plot
p_beta <- ggplot(pc12, aes(x = pc_x, y = pc_y, color = group)) +
  geom_point(size = 0.75) +
  geom_segment(aes(x = mean(pc_x), y = mean(pc_y), xend = pc_x, yend = pc_y, color = group),
               data = pc12 %>% group_by(group) %>% mutate(x_mean = mean(pc_x), y_mean = mean(pc_y))) +
  stat_ellipse(geom = "polygon", level = 0.9, linetype = 2, size = 0.5,
               aes(fill = group), alpha = 0.1, show.legend = TRUE) +
  coord_fixed(ratio = 1.25) +
  ylab(paste0("PCoA2 (", round(pc[2], 2), "%)")) +
  xlab(paste0("PCoA1 (", round(pc[1], 2), "%)")) +
  scale_fill_manual(values = c("#9B3A4D", "#E2AE79", "#D0DCAA", "#70A0AC", "#8CBDA7", "#566CA5", "#F0EEBB", "#7d3f98")) +
  scale_color_manual(values = c("#9B3A4D", "#E2AE79", "#D0DCAA", "#70A0AC", "#8CBDA7", "#566CA5", "#F0EEBB", "#7d3f98")) +
  theme_bw() +
  guides(color = "none", fill = guide_legend(override.aes = list(alpha = 1))) +
  labs(title = "PCoA of Relative Abundance") +
  theme(legend.title = element_blank(),
        legend.position = "right",
        panel.background = element_blank(),
        plot.title = element_text(size = 15, color = "black", hjust = 0.5, face = "bold"))

ggsave(plot = p_beta, filename = file.path(fig_dir, "beta_diversity_PCoA.pdf"),
       height = 5.5, width = 5.5)

# ----------------------------- 4. Top genera stacked barplot ----------------
# Calculate mean relative abundance per group
groups <- c("CHD", "IPAH", "SLE", "NOR")
rel_list <- list()
for (g in groups) {
  samples <- sampleID[[paste0(g, "_not4pnk")]]
  rel_sub <- rel_abundance[, colnames(rel_abundance) %in% samples, drop = FALSE]
  rel_list[[g]] <- rowMeans(rel_sub)
}
df_ratio <- as.data.frame(rel_list)

# Get top 10 genera per group
top_rows <- unique(unlist(lapply(df_ratio, function(x) rownames(df_ratio)[order(x, decreasing = TRUE)][1:10])))
df_plot <- df_ratio[top_rows, ]
other_row <- 1 - colSums(df_plot)
df_plot <- rbind(df_plot, other_row)
rownames(df_plot)[nrow(df_plot)] <- "others"

# Long format for plotting
df_long <- df_plot %>%
  rownames_to_column("genus") %>%
  pivot_longer(cols = -genus, names_to = "disease", values_to = "proportion")
df_long$genus <- factor(df_long$genus, levels = c(setdiff(unique(df_long$genus), "others"), "others"))
df_long$disease <- factor(df_long$disease, levels = c("NOR", "IPAH", "CHD", "SLE"))

# Colors
custom_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
  "#525252", "#F4CA26", "#9b43a8", "#63ACBE", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#009E73", "#56B4E9",
  "#E41A1C", "#6B8E6B", "#777777"
)

p_stacked <- ggplot(df_long, aes(x = disease, y = proportion, fill = genus)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Disease", y = "Proportion", fill = "Genus") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 10),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(p_stacked, file = file.path(fig_dir, "genus_proportion.pdf"),
       width = 17, height = 20, units = "cm")

# ----------------------------- 5. Alpha diversity ----------------------------
# Transpose for diversity calculation
df_t <- as.data.frame(t(df_raw_microbeRNA))
shannon <- diversity(df_t, index = "shannon")
simpson <- diversity(df_t, index = "simpson")
S <- specnumber(df_t)
evenness <- shannon / log(S)
richness <- estimateR(df_t)[1, ]

diversity_df <- data.frame(shannon, simpson, evenness, richness)
diversity_df$ID <- rownames(diversity_df)

# Add group info
diversity_df$Group <- meta_df$group[match(diversity_df$ID, meta_df$V1)]
diversity_df$Group <- factor(diversity_df$Group, levels = c("NOR", "IPAH", "CHD", "SLE"))

# Plot Shannon diversity
p_alpha <- ggplot(diversity_df, aes(x = Group, y = shannon, fill = Group)) +
  geom_boxplot(width = 0.35) +
  scale_fill_manual(values = c("#9B3A4D", "#E2AE79", "#D0DCAA", "#F0EEBB")) +
  stat_compare_means(method = "anova", label.y = max(diversity_df$shannon, na.rm = TRUE) + 1) +
  stat_compare_means(label = "p.signif", ref.group = "NOR") +
  labs(x = "", y = "Shannon Index", title = "α Diversity", fill = "") +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")

ggsave(plot = p_alpha, filename = file.path(fig_dir, "alpha_diversity_PH.pdf"),
       height = 3.5, width = 3.5)

# ----------------------------- 6. Differential abundance (Wilcoxon) ---------
# Prepare data for PH vs NOR comparison
group_nor <- sampleID$NOR_not4pnk
group_ph <- sampleID$PH_not4pnk
expr_filtered <- df_raw_microbeRNA[, intersect(c(group_nor, group_ph), colnames(df_raw_microbeRNA))]

wilcox_results <- apply(expr_filtered, 1, function(gene_expr) {
  test <- wilcox.test(x = gene_expr[group_nor], y = gene_expr[group_ph],
                      alternative = "two.sided", exact = FALSE)
  data.frame(p.value = test$p.value, statistic = test$statistic)
}) %>% bind_rows(.id = "gene")

wilcox_results$log2FoldChange <- apply(expr_filtered, 1, function(x) {
  log2((mean(x[group_ph], na.rm = TRUE) + 1e-6) / (mean(x[group_nor], na.rm = TRUE) + 1e-6))
})
wilcox_results$padj <- p.adjust(wilcox_results$p.value, method = "fdr")
write.csv(wilcox_results, file = file.path(res_dir, "Wilcoxon_test_results.csv"), row.names = FALSE)

# Volcano plot
volcano_data <- wilcox_results %>%
  mutate(log10_padj = -log10(padj),
         significance = case_when(
           padj < 0.05 & log2FoldChange > 0.5 ~ "Upregulated",
           padj < 0.05 & log2FoldChange < -0.5 ~ "Downregulated",
           TRUE ~ "Not significant"
         )) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  mutate(label = ifelse(row_number() <= 10 & padj < 0.05, gene, ""))

colors <- c("Upregulated" = "#FF6B6B", "Downregulated" = "#4E79A7", "Not significant" = "gray80")

p_volcano <- ggplot(volcano_data, aes(x = log2FoldChange, y = log10_padj, color = significance)) +
  geom_point(size = 2, alpha = 1, stroke = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_text_repel(aes(label = label), size = 3, box.padding = 0.3, point.padding = 0.1,
                  max.overlaps = 20, segment.color = "grey40", segment.size = 0.2) +
  scale_color_manual(values = colors) +
  scale_x_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 1)) +
  labs(x = expression(Log[2]~fold~change~"(PH vs NOR)"),
       y = expression(-Log[10]~adjusted~italic(P)),
       color = NULL) +
  theme_classic(base_size = 10) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        legend.position = "top")

ggsave(p_volcano, filename = file.path(fig_dir, "Volcano_Plot_microbe.pdf"),
       width = 4, height = 5, device = cairo_pdf, dpi = 600)

# ----------------------------- 7. Maaslin2 analysis -------------------------
# Prepare matrix and metadata
matrix <- df_raw_microbeRNA
metas <- meta_df
colnames(metas) <- c("SampleID", "group")
rownames(metas) <- metas$SampleID
metas <- metas[colnames(matrix), , drop = FALSE]

maaslin_dir <- file.path(res_dir, "maaslin2")
if (!dir.exists(maaslin_dir)) dir.create(maaslin_dir)

Maaslin2(
  input_data = matrix,
  input_metadata = metas,
  output = maaslin_dir,
  analysis_method = "LM",
  correction = "BH",
  normalization = "TMM",
  plot_heatmap = TRUE,
  plot_scatter = FALSE,
  heatmap_first_n = 50,
  min_prevalence = 0,
  max_significance = 1,
  fixed_effects = c("group"),
  reference = c("group", "NOR")
)

# Read results and create heatmap
res_row <- read.csv(file.path(maaslin_dir, "all_results.tsv"), sep = "\t")
res_filter <- res_row %>%
  group_by(value) %>%
  arrange(qval) %>%
  slice_head(n = 15)
res_filter$sig <- -log(res_filter$qval) * sign(res_filter$coef)
res_filter$label <- ifelse(res_filter$coef > 0, "+", "-")

# Clustering
mat <- res_filter %>%
  pivot_wider(id_cols = feature, names_from = value, values_from = sig) %>%
  column_to_rownames("feature") %>%
  as.matrix()
row_hc <- hclust(dist(mat, method = "euclidean"), method = "complete")
col_hc <- hclust(dist(t(mat), method = "euclidean"), method = "complete")
row_order <- row_hc$labels[row_hc$order]
col_order <- col_hc$labels[col_hc$order]

masslin2_clust <- res_filter %>%
  mutate(feature = factor(feature, levels = row_order),
         value = factor(value, levels = col_order))

p_maaslin <- ggplot(masslin2_clust, aes(x = value, y = feature, fill = sig)) +
  geom_tile(color = "grey80", size = 0.3) +
  scale_fill_gradient2(low = "#004b79", mid = "white", high = "#7f181b",
                       midpoint = 0, limits = c(-15, 15), oob = squish) +
  geom_text(aes(label = label), size = 3, color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(color = "black", size = 8),
        legend.title = element_text(angle = 90)) +
  labs(fill = "-log(FDR)*sign(coef)")

ggsave(p_maaslin, filename = file.path(fig_dir, "maaslin2_feature_heatmap.pdf"),
       height = 12, width = 3.5)

cat("All microbiome QC and diversity analyses completed.\n")