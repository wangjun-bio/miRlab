# =============================================================================
# REVEAL 2.0 and COMPERA 2.0 Risk Stratification
# =============================================================================
# This script:
#   1. Reads clinical data from CSV (Meta_REVEL_241211.csv)
#   2. Calculates REVEAL 2.0 risk scores and subgroups (low/moderate/high)
#   3. Calculates COMPERA 2.0 risk strata (low/intermediate/high)
#   4. Outputs sample IDs for low and moderate REVEAL groups
# =============================================================================

# ----------------------------- Libraries -------------------------------------
library(readxl)
library(dplyr)

# ----------------------------- User configuration ----------------------------
# Input file (relative to project root)
input_file <- file.path("data", "Meta_REVEL_241211.csv")

# Output files (saved to results/REVEAL/)
output_dir <- file.path("results", "REVEAL")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ----------------------------- Read clinical data ----------------------------
clinical_info <- read.csv(input_file, fileEncoding = "GBK")
table(clinical_info$分类)

# ----------------------------- Calculate REVEAL 2.0 -------------------------
revel_info <- clinical_info[, c("Seq_ID", "姓名...8", "分类", "性别", "年龄", "心功能分级",
                                "右心导管sbp", "右心导管hr", "六mwt", "ntprobnp", "右心房压", "肺血管阻力")]
colnames(revel_info) <- c("Seq_ID", "name", "WHO", "sex", "age", "NYHA",
                          "SBP", "HR", "sixmwt", "ntprobnp", "mRAP", "PVR")

# Replace NA with 0
revel_info <- revel_info %>%
  mutate(across(-Seq_ID, ~ ifelse(is.na(.), 0, .)))

# Filter out samples with missing key variables
revel_info <- revel_info %>%
  filter(!(NYHA == 0 | sixmwt == 0 | ntprobnp == 0))

# REVEAL 2.0 scoring (based on original code)
revel_info <- revel_info %>% 
  mutate(REVEL = case_when(
      WHO == "SLE-PH" ~ 1,
      WHO == "HPAH" ~ 2,
      TRUE ~ 0) +
    case_when(sex == 1 & age > 60 ~ 2, TRUE ~ 0) +
    case_when(NYHA == 1 ~ -1, NYHA == 3 ~ 1, NYHA == 4 ~ 2, TRUE ~ 0) +
    case_when(SBP < 110 ~ 1, TRUE ~ 0) +
    case_when(HR > 96 ~ 1, TRUE ~ 0) +
    case_when(sixmwt >= 440 ~ -2, sixmwt >= 320 & sixmwt < 440 ~ -1, sixmwt < 165 ~ 1, TRUE ~ 0) +
    case_when(ntprobnp < 300 ~ -2, ntprobnp >= 300 & ntprobnp < 1100 ~ 1, ntprobnp >= 1100 ~ 2, TRUE ~ 0) +
    case_when(mRAP > 20 ~ 1, TRUE ~ 0) +
    case_when(PVR < 5 ~ -1, TRUE ~ 0) + 6) %>%
  mutate(sub_groups = case_when(
    REVEL <= 6 ~ "low",
    REVEL > 6 & REVEL <= 10 ~ "moderate",
    REVEL > 10 ~ "high"
  ))

# Remove invalid Seq_ID
revel_info <- revel_info[which(revel_info$Seq_ID != "#N/A"), ]
colnames(revel_info)[14] <- "REVEAL_sub_groups"

# ----------------------------- Calculate COMPERA 2.0 ------------------------
calculate_compera_2.0 <- function(data) {
  data %>%
    mutate(
      high_risk_WHO = ifelse(NYHA == 4, TRUE, FALSE),
      high_risk_6MWD = ifelse(sixmwt < 165, TRUE, FALSE),
      high_risk_BNP = ifelse(ntprobnp > 1100, TRUE, FALSE),
      low_risk_WHO = ifelse(NYHA %in% c(1, 2), TRUE, FALSE),
      low_risk_6MWD = ifelse(sixmwt > 440, TRUE, FALSE),
      low_risk_BNP = ifelse(ntprobnp < 300, TRUE, FALSE),
      COMPERA_2.0 = case_when(
        low_risk_WHO & low_risk_6MWD & low_risk_BNP ~ "low",
        high_risk_WHO | high_risk_6MWD | high_risk_BNP ~ "high",
        TRUE ~ "intermediate"
      )
    )
}

tmp <- calculate_compera_2.0(revel_info)
tmp <- tmp[which(tmp$Seq_ID != "#N/A"), ]
compera_info <- tmp[, c("Seq_ID", "COMPERA_2.0")]

# Merge REVEAL and COMPERA
stratification <- merge(revel_info, compera_info, by = "Seq_ID")

# Save combined stratification table
write.csv(stratification, file = file.path(output_dir, "Meta_REVEAL_241223.csv"),
          row.names = FALSE, fileEncoding = "GBK")

# ----------------------------- Extract sample IDs for low and moderate REVEAL groups
low_seq_ids <- revel_info$Seq_ID[revel_info$sub_groups == "low"]
moderate_seq_ids <- revel_info$Seq_ID[revel_info$sub_groups == "moderate"]

write.table(low_seq_ids, file = file.path(output_dir, "REVEAL_low.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(moderate_seq_ids, file = file.path(output_dir, "REVEAL_moderate.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("REVEAL and COMPERA stratification completed.\n")
cat("Output files saved to:", output_dir, "\n")