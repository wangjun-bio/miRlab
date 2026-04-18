# =============================================================================
# Microbiome Decontamination using KrakenUniq Outputs
# =============================================================================
# This script:
#   1. Reads KrakenUniq report files for water controls and patient samples
#   2. Identifies contaminants using:
#      - decontam (prevalence method)
#      - Wilcoxon test comparing water controls vs samples
#   3. Also incorporates known contaminants from:
#      - Skin-derived microbes
#      - Non-human viruses (eukaryotic hosts other than human)
#      - Published laboratory contamination lists
#   4. Outputs a filtered genus-level abundance matrix for downstream analysis
# =============================================================================

# ----------------------------- Libraries -------------------------------------
library(dplyr)
library(data.table)
library(pheatmap)
library(tibble)
library(here)
library(openxlsx)
library(decontam)
library(phyloseq)

# ----------------------------- User configuration ----------------------------
# Input directories (relative to project root)
kraken_dir      <- file.path("data", "krakenuniq")           # contains report files
h20_dir         <- file.path(kraken_dir, "h20_all", "krakenuniq_output")
patient_dir     <- file.path(kraken_dir, "20250707")
skin_meta_file  <- file.path(kraken_dir, "mmc2.xlsx")
skin_data_file  <- file.path(kraken_dir, "mmc3.xlsx")
virus_host_file <- file.path(kraken_dir, "virushostdb.daily.tsv")

# Output directory
output_dir <- file.path("results", "microbiome")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ----------------------------- 1. Water control samples ---------------------
files_h20 <- list.files(h20_dir, pattern = "_reportfile.tsv", full.names = TRUE)

result_h20 <- lapply(files_h20, function(x) {
  tmpname <- gsub("_reportfile.tsv", "", basename(x))
  tmp <- read.table(x, header = TRUE, sep = "\t")
  tmp <- tmp[tmp$rank == "genus", ]
  tmp$taxName <- gsub(" ", "", tmp$taxName)
  tmp <- tmp[, c(9, 2)]
  colnames(tmp)[2] <- tmpname
  return(tmp)
})

h20_all <- Reduce(function(x, y) merge(x, y, by = "taxName", all = TRUE), result_h20)
rownames(h20_all) <- h20_all[, 1]
h20_all <- h20_all[, -1]
h20_all <- apply(h20_all, 2, function(row) { row[is.na(row)] <- 0; row })
colnames(h20_all) <- paste0(colnames(h20_all), "_water")

# ----------------------------- 2. Skin-derived microbes ---------------------
skin_meta <- read.xlsx(skin_meta_file, check.names = FALSE)
colnames(skin_meta) <- skin_meta[1, ]
skin_meta <- skin_meta[-1, ]
skin_meta <- skin_meta[skin_meta$Timepoint == "First", ]

skin <- read.xlsx(skin_data_file, check.names = FALSE)
colnames(skin) <- skin[1, ]
skin <- skin[-1, ]
skin <- skin[6:nrow(skin), ]
skin <- skin %>%
  mutate(genus = sapply(strsplit(Taxa, ";"), function(x) x[length(x) - 1]))
skin[, 1] <- skin$genus
skin <- skin %>% select(-genus)
skin <- as.data.frame(skin)
skin[, -1] <- lapply(skin[, -1], as.numeric)
skin <- skin %>% group_by(Taxa) %>% summarise(across(everything(), sum), .groups = "drop")
skin <- skin[-1, ]
rownames(skin) <- skin[, 1]
skin <- skin[, -1]
skin <- skin[, colnames(skin) %in% skin_meta$SAMPLE]

# Keep genera present in at least 10% of skin samples
min_sample_count <- ceiling(0.1 * ncol(skin))
keep_rows <- rowSums(skin > 0) >= min_sample_count
skin <- skin[keep_rows, ]
filtered_skin <- rownames(skin)
filtered_skin_df <- data.frame(skin_derived = rep(1, length(filtered_skin)),
                               row.names = filtered_skin)

# ----------------------------- 3. Non-human viruses -------------------------
virus <- fread(virus_host_file)
virus <- virus %>%
  mutate(superkingdom = sapply(strsplit(`host lineage`, ";"), function(x) x[2]))
virus <- virus[virus$superkingdom == " Eukaryota", ]
virus <- virus[virus$`host name` != "Homo sapiens", ]
virus <- virus %>%
  mutate(genus = sapply(strsplit(`virus lineage`, ";"), function(x) x[length(x) - 1]))
virus <- virus[!grepl("unclassified", virus$`virus lineage`), ]
virus$genus <- gsub(" ", "", virus$genus)
virus <- virus$genus[!duplicated(virus$genus)]
virus_df <- data.frame(virus = rep(1, length(virus)), row.names = virus)

# ----------------------------- 4. Published lab contaminants ----------------
lab_contaminants <- c(
  'Afipia', 'Aquabacterium', 'Asticcacaulis', 'Aurantimonas', 'Beijerinckia', 'Bosea',
  'Bradyrhizobium', 'Brevundimonas', 'Caulobacter', 'Craurococcus', 'Devosia', 'Hoeflea',
  'Mesorhizobium', 'Methylobacterium', 'Novosphingobium', 'Ochrobactrum', 'Paracoccus',
  'Pedomicrobium', 'Phyllobacterium', 'Rhizobium', 'Roseomonas', 'Sphingobium', 'Sphingomonas','Sphingopyxis',
  'Acidovorax', 'Azoarcus', 'Azospira', 'Burkholderia', 'Comamonas',
  'Cupriavidus', 'Curvibacter', 'Delftia', 'Duganella', 'Herbaspirillum', 'Janthinobacterium', 'Kingella',
  'Leptothrix', 'Limnobacter', 'Massilia', 'Methylophilus', 'Methyloversatilis', 'Oxalobacter', 'Pelomonas',
  'Polaromonas', 'Ralstonia','Schlegelella', 'Sulfuritalea', 'Undibacterium', 'Variovorax',
  'Acinetobacter','Enhydrobacter', 'Enterobacter', 'Escherichia' ,'Nevskia', 'Pseudomonas', 'Pseudoxanthomonas', 'Psychrobacter',
  'Stenotrophomonas','Xanthomonas',
  'Aeromicrobium', 'Arthrobacter', 'Beutenbergia', 'Brevibacterium', 'Corynebacterium', 'Curtobacterium', 'Dietzia',
  'Geodermatophilus', 'Janibacter', 'Kocuria', 'Microbacterium', 'Micrococcus', 'Microlunatus', 'Patulibacter', 'Propionibacterium',
  'Rhodococcus', 'Tsukamurella',
  'Abiotrophia', 'Bacillus', 'Brevibacillus', 'Brochothrix', 'Facklamia', 'Paenibacillus', 'Streptococcus',
  'Chryseobacterium', 'Dyadobacter', 'Flavobacterium', 'Hydrotalea', 'Niastella', 'Olivibacter', 'Pedobacter', 'Wautersiella',
  'Deinococcus', 'Homo'
)

lab_contam_df <- data.frame(lab_contamination = rep(1, length(lab_contaminants)),
                            row.names = lab_contaminants)

# ----------------------------- 5. Read patient samples ----------------------
files_patient <- list.files(patient_dir, pattern = "_reportfile.tsv", full.names = TRUE)

result_patient <- lapply(files_patient, function(x) {
  tmpname <- gsub("_reportfile.tsv", "", basename(x))
  tmp <- read.table(x, header = TRUE, sep = "\t")
  tmp <- tmp[tmp$rank == "genus", ]
  tmp$taxName <- gsub(" ", "", tmp$taxName)
  tmp <- tmp[, c(9, 2)]
  colnames(tmp)[2] <- tmpname
  return(tmp)
})

PAN <- Reduce(function(x, y) merge(x, y, by = "taxName", all = TRUE), result_patient)
rownames(PAN) <- PAN[, 1]
PAN <- PAN[, -1]
PAN <- apply(PAN, 2, function(row) { row[is.na(row)] <- 0; row })
PAN <- as.data.frame(PAN)

# ----------------------------- 6. Combine all samples for decontam ----------
# Merge patient samples and water controls
ALL <- list(PAN, h20_all)
ALL <- lapply(ALL, function(x) { x <- rownames_to_column(x); return(x) })
ALL <- Reduce(function(x, y) merge(x, y, by = "rowname", all = TRUE), ALL)
ALL[is.na(ALL)] <- 0
rownames(ALL) <- ALL[, 1]
ALL <- ALL[, -1]

# Create sample metadata
sample_data <- data.frame(sample = colnames(ALL))
sample_data$is.neg <- sapply(sample_data$sample, function(x) {
  if (grepl("water", x)) "WaterControl" else "Sample"
})
sample_data$is.neg <- sample_data$is.neg == "WaterControl"
rownames(sample_data) <- sample_data$sample
sample_data <- sample_data[, "is.neg", drop = FALSE]

# Phyloseq object
OTU <- otu_table(as.matrix(t(ALL)), taxa_are_rows = FALSE)
SAM <- sample_data(sample_data)
ps <- phyloseq(OTU, SAM)

# decontam prevalence method
contam_df <- isContaminant(ps, method = "prevalence", neg = "is.neg")
decontam_contaminants <- rownames(contam_df)[contam_df$contaminant == TRUE]

# ----------------------------- 7. Wilcoxon test for contamination -----------
samples <- rownames(sample_data)[!sample_data$is.neg]
h20_samples <- rownames(sample_data)[sample_data$is.neg]

# Calculate relative abundance
zero_sum_cols <- names(which(colSums(ALL) == 0))
if (length(zero_sum_cols) > 0) {
  ALL_ratio <- apply(ALL[, !names(ALL) %in% zero_sum_cols, drop = FALSE], 2,
                     function(x) x / sum(x))
} else {
  ALL_ratio <- apply(ALL, 2, function(x) x / sum(x))
}

ALL_wilcox <- lapply(rownames(ALL_ratio), function(x) {
  tmp <- ALL_ratio[x, ]
  tmp_df <- data.frame(sample = names(tmp), ratio = tmp)
  tmp_df$group <- ifelse(grepl("water", tmp_df$sample), "WaterControl", "samples")
  tmp_df$group <- factor(tmp_df$group, levels = c("WaterControl", "samples"))
  if (length(unique(tmp_df$group)) < 2) {
    return(data.frame(taxa = x, p.value = NA, padj = NA, significant = "Not Significant",
                      mean_ratio = NA, mean_ratio_h20 = NA, log2FoldChange = NA, increase = NA))
  }
  wilcox_result <- wilcox.test(ratio ~ group, data = tmp_df)
  y <- data.frame(taxa = x, p.value = wilcox_result$p.value)
  y$padj <- p.adjust(y$p.value, method = "BH")
  y$significant <- ifelse(y$padj < 0.01, "Significant", "Not Significant")
  y$mean_ratio <- mean(tmp_df$ratio[tmp_df$group == "samples"])
  y$mean_ratio_h20 <- mean(tmp_df$ratio[tmp_df$group == "WaterControl"])
  y$log2FoldChange <- log2(y$mean_ratio / y$mean_ratio_h20)
  y$increase <- ifelse(y$log2FoldChange > 0, "Increase", "Decrease")
  return(y)
})
ALL_wilcox <- do.call(bind_rows, ALL_wilcox)
wilcox_contaminants <- ALL_wilcox %>%
  filter(significant == "Significant", increase == "Decrease") %>%
  pull(taxa)

# ----------------------------- 8. Combine all contaminant lists -------------
all_contamination <- unique(c(lab_contaminants, filtered_skin, virus, decontam_contaminants, wilcox_contaminants))
write.table(data.frame(contaminant = all_contamination),
            file = file.path(output_dir, "Microbe_contamination.txt"),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

# ----------------------------- 9. Filter patient samples --------------------
PAN_filtered <- PAN[!(rownames(PAN) %in% all_contamination), ]

# Save filtered matrix
write.csv(PAN_filtered, file = file.path(output_dir, "PAH_microbeRNA_filtered.csv"))

# ----------------------------- 10. Subset to final 452 samples --------------
# Load sample IDs (if available)
sampleID_file <- file.path("data", "sampleID_filted300W.rds")
if (file.exists(sampleID_file)) {
  sampleID <- readRDS(sampleID_file)
  sampleID <- unlist(sampleID[!names(sampleID) %in% c("CHD_t4pnk", "IPAH_t4pnk", "SLE_t4pnk", "NOR_t4pnk", "PH_not4pnk")])
  PAN_filtered <- PAN_filtered[, colnames(PAN_filtered) %in% sampleID]
  write.csv(PAN_filtered, file = file.path(output_dir, "PAH_microbeRNA_452samples.csv"))
}

cat("Decontamination completed. Output saved to:", output_dir, "\n")