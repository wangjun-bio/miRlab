# =============================================================================
# Mutation Detection and MAF Conversion from ANNOVAR Outputs
# =============================================================================
# This script:
#   1. Reads ANNOVAR CSV files for each sample under different filtering conditions
#   2. Filters variants (functional impact, population frequency, conservation)
#   3. Converts to MAF format
#   4. Builds mutation presence matrices (sparse for large datasets)
#   5. Calculates mutation frequencies per group and compares with NOR
# =============================================================================

# ----------------------------- Libraries -------------------------------------
library(tidyverse)
library(maftools)
library(Matrix)  # for sparse matrix

# ----------------------------- User configuration ----------------------------
# Define input/output base paths (all relative to project root)
base_paths <- list(
  # Input: ANNOVAR CSV files (should be placed in data/annovar/csv_file/)
  annovar_csv = file.path("data", "annovar", "csv_file"),
  
  # Output directories (will be created automatically)
  maf_output      = file.path("results", "maftools"),
  matrix_output   = file.path("results", "mut_presence_matrix"),
  result_output   = file.path("results", "mutation_freq"),
  
  # Input metadata files (RDS)
  sample_id_rds   = file.path("data", "sampleID.rds"),
  df_ratio_rds    = file.path("data", "df_RNA_ratio.rds")
)

# Create output directories if not exist
for (dir in c(base_paths$maf_output, base_paths$matrix_output, base_paths$result_output)) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# Filter conditions to process (substrings of ANNOVAR output file names)
filter_conditions <- c("DP0_GT01", "DP0_GT11", "DP30_GT11")

# Sample groups (must match names in sampleID.rds)
sample_groups <- c("CHD_not4pnk", "IPAH_not4pnk", "NOR_not4pnk", "SLE_not4pnk")

# ----------------------------- Helper functions ------------------------------

# Get file paths for samples in a given group and filter condition
get_sample_files <- function(condition, sample_group) {
  sampleID <- readRDS(base_paths$sample_id_rds)
  df_merge <- readRDS(base_paths$df_ratio_rds) %>%
    filter(clean_reads > 3000000) %>%
    pull(ID)
  
  group_samples <- sampleID[[sample_group]]
  valid_samples <- intersect(group_samples, df_merge)
  
  file.path(base_paths$annovar_csv,
            paste0(valid_samples, "_filtered_", condition, ".hg38_multianno.csv"))
}

# Convert ANNOVAR CSV to MAF format (with extensive filtering)
convert_to_maf <- function(csv_file) {
  sample_id <- strsplit(tools::file_path_sans_ext(basename(csv_file)), "_filtered")[[1]][1]
  
  # Handle empty file
  if (file.info(csv_file)$size == 0) {
    message(paste("Empty file:", sample_id))
    return(tibble(
      Hugo_Symbol = NA_character_, Chromosome = NA_character_,
      Start_Position = NA_integer_, End_Position = NA_integer_,
      Reference_Allele = NA_character_, Tumor_Seq_Allele1 = NA_character_,
      Tumor_Seq_Allele2 = NA_character_, AAChange = NA_character_,
      dbSNP = NA_character_, Depth = NA_integer_,
      Variant_Classification = NA_character_, Variant_Type = NA_character_,
      Genomic_Location = NA_character_, Tumor_Sample_Barcode = sample_id
    ))
  }
  
  annovar_data <- read.csv(csv_file, stringsAsFactors = FALSE) %>%
    filter(avsnp151 != ".") %>%
    # Keep exonic and splicing variants
    filter(Func.refGeneWithVer %in% c("exonic", "splicing")) %>%
    # Remove synonymous SNVs
    filter(ExonicFunc.refGeneWithVer != "synonymous SNV") %>%
    # Filter missense: require at least two tools predicting damaging
    filter(
      (ExonicFunc.refGeneWithVer == "nonsynonymous SNV" &
         ((REVEL_score > 0.5 & SIFT_pred == "D") |
          (REVEL_score > 0.5 & Polyphen2_HVAR_pred == "D") |
          (SIFT_pred == "D" & Polyphen2_HVAR_pred == "D"))) |
        ExonicFunc.refGeneWithVer %in% c("stopgain", "stoploss", "startloss", "startgain",
                                         "frameshift insertion", "frameshift deletion",
                                         "nonframeshift deletion", "frameshift substitution",
                                         "nonframeshift substitution", "nonframeshift insertion") |
        Func.refGeneWithVer == "splicing"
    ) %>%
    # Population frequency filtering (gnomAD, controls, East Asian)
    mutate(
      AF = as.numeric(AF), AF_popmax = as.numeric(AF_popmax),
      non_topmed_AF_popmax = as.numeric(non_topmed_AF_popmax),
      controls_AF_popmax = as.numeric(controls_AF_popmax),
      AF_eas = as.numeric(AF_eas)
    ) %>%
    filter(
      is.na(AF) | AF < 0.01,
      is.na(AF_popmax) | AF_popmax < 0.01,
      is.na(non_topmed_AF_popmax) | non_topmed_AF_popmax < 0.01,
      is.na(controls_AF_popmax) | controls_AF_popmax < 0.01,
      is.na(AF_eas) | AF_eas < 0.01
    ) %>%
    # Conservation filters
    mutate(GERP_RS = as.numeric(GERP.._RS), CADD_PHRED = as.numeric(CADD_phred)) %>%
    filter(
      (is.na(GERP_RS) | GERP_RS > 4) &
      (is.na(CADD_PHRED) | CADD_PHRED > 20)
    )
  
  if (nrow(annovar_data) == 0) {
    message(paste("No variants after filtering:", sample_id))
    return(tibble(
      Hugo_Symbol = NA_character_, Chromosome = NA_character_,
      Start_Position = NA_integer_, End_Position = NA_integer_,
      Reference_Allele = NA_character_, Tumor_Seq_Allele1 = NA_character_,
      Tumor_Seq_Allele2 = NA_character_, AAChange = NA_character_,
      dbSNP = NA_character_, Depth = NA_integer_,
      Variant_Classification = NA_character_, Variant_Type = NA_character_,
      Genomic_Location = NA_character_, Tumor_Sample_Barcode = sample_id
    ))
  }
  
  # Add depth information from dp_info.txt
  dp_file <- file.path(base_paths$annovar_csv, paste0(sample_id, "_dp_info.txt"))
  if (file.exists(dp_file)) {
    dp_info <- read.table(dp_file, header = FALSE, col.names = c("CHROM", "POS", "DP"))
    annovar_data <- merge(annovar_data, dp_info, by.x = c("Chr", "Start"), by.y = c("CHROM", "POS"), all.x = TRUE)
  } else {
    annovar_data$DP <- NA
  }
  
  # Convert to MAF format
  annovar_data %>%
    mutate(
      Hugo_Symbol = Gene.refGeneWithVer,
      Chromosome = Chr,
      Start_Position = Start,
      End_Position = End,
      Reference_Allele = Ref,
      Tumor_Seq_Allele1 = Ref,
      Tumor_Seq_Allele2 = Alt,
      AAChange = AAChange.refGeneWithVer,
      dbSNP = avsnp151,
      Depth = DP,
      Genomic_Location = Func.refGeneWithVer,
      Variant_Classification = case_when(
        ExonicFunc.refGeneWithVer == "nonsynonymous SNV" ~ "Missense_Mutation",
        ExonicFunc.refGeneWithVer == "stopgain" ~ "Nonsense_Mutation",
        ExonicFunc.refGeneWithVer == "stoploss" ~ "Nonstop_Mutation",
        ExonicFunc.refGeneWithVer == "startloss" ~ "Translation_Start_Site",
        ExonicFunc.refGeneWithVer == "startgain" ~ "Translation_Start_Site",
        ExonicFunc.refGeneWithVer == "frameshift insertion" ~ "Frame_Shift_Ins",
        ExonicFunc.refGeneWithVer == "frameshift deletion" ~ "Frame_Shift_Del",
        ExonicFunc.refGeneWithVer == "nonframeshift deletion" ~ "In_Frame_Del",
        ExonicFunc.refGeneWithVer == "frameshift substitution" ~ "Frame_Shift_Del",
        ExonicFunc.refGeneWithVer == "nonframeshift substitution" ~ "In_Frame_Del",
        ExonicFunc.refGeneWithVer == "nonframeshift insertion" ~ "In_Frame_Ins",
        Func.refGeneWithVer == "splicing" ~ "Splice_Site",
        TRUE ~ "Other_Mutation"
      ),
      Variant_Type = as.character(ifelse(nchar(coalesce(Ref, "")) == nchar(coalesce(Alt, "")), "SNP",
                                         ifelse(nchar(coalesce(Ref, "")) < nchar(coalesce(Alt, "")), "INS", "DEL"))),
      Tumor_Sample_Barcode = sample_id
    ) %>%
    dplyr::select(Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele,
                  Tumor_Seq_Allele1, Tumor_Seq_Allele2, AAChange, dbSNP, Depth,
                  Variant_Classification, Variant_Type, Tumor_Sample_Barcode, Genomic_Location)
}

# Build mutation presence matrix (0/1 for each variant per sample)
build_mut_matrix <- function(maf_data, use_sparse = FALSE) {
  df <- maf_data %>%
    mutate(variant = paste0(Chromosome, ":", Start_Position, "-", End_Position,
                            ":", Reference_Allele, ">", Tumor_Seq_Allele1, "+", Tumor_Seq_Allele2))
  
  annotation_df <- df %>%
    group_by(variant) %>%
    summarise(
      Hugo_Symbol = first(Hugo_Symbol),
      Chromosome = first(Chromosome),
      Start_Position = first(Start_Position),
      End_Position = first(End_Position),
      Reference_Allele = first(Reference_Allele),
      Tumor_Seq_Allele1 = first(Tumor_Seq_Allele1),
      Tumor_Seq_Allele2 = first(Tumor_Seq_Allele2),
      AAChange = first(AAChange),
      dbSNP = first(dbSNP),
      Variant_Classification = first(Variant_Classification),
      Variant_Type = first(Variant_Type),
      Depth_Mean = mean(as.numeric(Depth), na.rm = TRUE),
      .groups = "drop"
    )
  
  if (use_sparse) {
    row_idx <- as.integer(factor(df$variant))
    col_idx <- as.integer(factor(df$Tumor_Sample_Barcode))
    presence_mat <- sparseMatrix(
      i = row_idx, j = col_idx, x = 1,
      dims = c(nlevels(factor(df$variant)), nlevels(factor(df$Tumor_Sample_Barcode))),
      dimnames = list(levels(factor(df$variant)), levels(factor(df$Tumor_Sample_Barcode)))
    )
    presence_df <- as.data.frame(as.matrix(presence_mat)) %>%
      mutate(variant = rownames(presence_mat))
  } else {
    presence_df <- df %>%
      dplyr::select(variant, Tumor_Sample_Barcode) %>%
      mutate(present = 1) %>%
      pivot_wider(names_from = Tumor_Sample_Barcode, values_from = present, values_fill = 0)
  }
  
  left_join(annotation_df, presence_df, by = "variant")
}

# Calculate mutation frequencies across groups and merge
calculate_and_merge_freq <- function(condition, maf_groups) {
  freq_list <- map2(maf_groups, names(maf_groups), function(maf, group) {
    temp_maf <- file.path(base_paths$maf_output, paste0(str_remove(group, "_not4pnk"), "_", condition, ".maf"))
    # write MAF to temp file for read.maf (but maf is already a data frame, we can directly compute)
    # Instead, compute frequencies directly from maf data frame
    variant_data <- maf
    
    variant_data %>%
      group_by(dbSNP, Hugo_Symbol, Chromosome, Start_Position, End_Position) %>%
      summarise(
        Mutated_Samples = n_distinct(Tumor_Sample_Barcode),
        Total_Samples = length(unique(variant_data$Tumor_Sample_Barcode)),
        Avg_Depth = mean(Depth, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        Frequency = Mutated_Samples / Total_Samples,
        group = str_remove(group, "_not4pnk")
      ) %>%
      rename_with(~paste0(., "_", str_remove(group, "_not4pnk")),
                  c(Mutated_Samples, Total_Samples, Avg_Depth, Frequency)) %>%
      dplyr::select(-group)
  })
  
  merged_freq <- freq_list %>%
    map(~mutate(., key = paste(Chromosome, Start_Position, End_Position, sep = "_"))) %>%
    reduce(full_join, by = "key")
  
  merged_freq <- merged_freq %>%
    rowwise() %>%
    mutate(
      Hugo_Symbol = first(na.omit(c_across(contains("Hugo_Symbol")))),
      Chromosome = first(na.omit(c_across(contains("Chromosome")))),
      Start_Position = first(na.omit(c_across(contains("Start_Position")))),
      End_Position = first(na.omit(c_across(contains("End_Position")))),
      dbSNP = first(na.omit(c_across(contains("dbSNP"))))
    ) %>%
    ungroup() %>%
    dplyr::select(-contains(c(".x", ".y")))
  
  merged_freq %>%
    mutate(
      diffmut_NOR_IPAH = Frequency_IPAH - Frequency_NOR,
      diffmut_NOR_CHD = Frequency_CHD - Frequency_NOR,
      diffmut_NOR_SLE = Frequency_SLE - Frequency_NOR,
      diffdepth_NOR_IPAH = (Avg_Depth_IPAH + 1) / (Avg_Depth_NOR + 1),
      diffdepth_NOR_CHD = (Avg_Depth_CHD + 1) / (Avg_Depth_NOR + 1),
      diffdepth_NOR_SLE = (Avg_Depth_SLE + 1) / (Avg_Depth_NOR + 1)
    ) %>%
    dplyr::select(-key)
}

# ----------------------------- Main processing function ----------------------
process_condition <- function(condition) {
  message("===== Processing filter condition: ", condition, " =====")
  
  # Get sample files for each group
  group_files <- map(sample_groups, ~get_sample_files(condition, .x)) %>%
    set_names(sample_groups)
  
  # Convert ANNOVAR to MAF for each group
  maf_list <- map(group_files, ~bind_rows(lapply(.x, convert_to_maf)))
  
  # Save individual group MAF files
  walk2(maf_list, names(maf_list), function(maf, group) {
    maf_path <- file.path(base_paths$maf_output, paste0(str_remove(group, "_not4pnk"), "_", condition, ".maf"))
    write.table(maf, maf_path, sep = "\t", quote = FALSE, row.names = FALSE)
  })
  all_maf <- bind_rows(maf_list)
  write.table(all_maf, file.path(base_paths$maf_output, paste0("ALL_", condition, ".maf")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Build mutation presence matrices
  use_sparse_all <- (condition == "DP0_GT11")
  mat_list <- map(maf_list, ~build_mut_matrix(.x, use_sparse = FALSE))
  all_mat <- build_mut_matrix(all_maf, use_sparse = use_sparse_all)
  
  # Save presence matrices
  walk2(mat_list, names(mat_list), function(mat, group) {
    mat_path <- file.path(base_paths$matrix_output, paste0(str_remove(group, "_not4pnk"), "_", condition, "_presence.rds"))
    saveRDS(mat, mat_path)
  })
  saveRDS(all_mat, file.path(base_paths$matrix_output, paste0("ALL_", condition, "_presence.rds")))
  
  # Calculate and merge frequencies
  merged_freq <- calculate_and_merge_freq(condition, maf_list)
  write.csv(merged_freq, file.path(base_paths$result_output, paste0("Mutation_call_", condition, ".csv")),
            row.names = FALSE)
  
  message("===== Filter condition ", condition, " completed =====")
}

# ----------------------------- Run all conditions ---------------------------
walk(filter_conditions, process_condition)

message("All mutation detection steps completed.")