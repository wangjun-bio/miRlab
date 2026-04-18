# =============================================================================
# Combine Differential Expression Results Across RNA Types
# =============================================================================
# This script merges DE results (|log2FC| > 0.5 & padj < 0.05) from all RNA
# types for each disease group (SLE, CHD, IPAH) into a single CSV file.
# =============================================================================

# ----------------------------- Libraries -------------------------------------
library(dplyr)

# ----------------------------- User configuration ----------------------------
# Directory containing DE CSV files (relative to project root)
de_dir <- file.path("results", "DE")

# RNA types to process (order can be adjusted)
rna_types <- c("snoRNA", "snRNA", "tRFs", "ysRNA", "rsRNA",
               "lncRNA", "mRNA", "miRNA")

# Disease groups to process
target_groups <- c("SLE", "CHD", "IPAH")

# ----------------------------- Helper function -------------------------------
process_group_de <- function(group_name, de_dir, rna_types) {
  cat("\n--- Processing group:", group_name, "---\n")
  
  # Read DE files for each RNA type (NOR_not4pnk vs group_name)
  de_list <- lapply(rna_types, function(type) {
    file_name <- file.path(de_dir,
                           paste0("DE_", type, "_NOR_not4pnk_", group_name, "_not4pnk.csv"))
    if (file.exists(file_name)) {
      df <- read.csv(file_name, stringsAsFactors = FALSE)
      return(df)
    } else {
      warning(paste("File not found:", file_name))
      return(NULL)
    }
  })
  names(de_list) <- rna_types
  de_list <- de_list[!sapply(de_list, is.null)]
  
  # Filter by |log2FC| > 0.5 and padj < 0.05, add RNA_Type column
  de_filtered_list <- lapply(names(de_list), function(type) {
    df <- de_list[[type]]
    filtered <- subset(df, abs(log2FoldChange) > 0.5 & padj < 0.05)
    if (nrow(filtered) > 0) {
      filtered$RNA_Type <- type
    }
    return(filtered)
  })
  names(de_filtered_list) <- names(de_list)
  
  # Combine all RNA types into one data frame
  final_df <- do.call(rbind, de_filtered_list)
  
  if (!is.null(final_df) && nrow(final_df) > 0) {
    # Move RNA_Type to first column
    col_order <- c("RNA_Type", colnames(final_df)[colnames(final_df) != "RNA_Type"])
    final_df <- final_df[, col_order]
    
    # Save combined results
    output_file <- file.path(de_dir,
                             paste0("PAH_cfRNA_", group_name, "_Differential_Expression_Total.csv"))
    write.csv(final_df, output_file, row.names = FALSE)
    
    cat(group_name, "completed! Total significant genes:", nrow(final_df), "\n")
    print(table(final_df$RNA_Type))
  } else {
    cat(group_name, "has no significant genes meeting criteria.\n")
  }
  
  return(final_df)
}

# ----------------------------- Run analysis ----------------------------------
all_results <- lapply(target_groups, function(g) {
  process_group_de(group_name = g, de_dir = de_dir, rna_types = rna_types)
})
names(all_results) <- target_groups

cat("\nAll group DE summaries have been saved to:\n", de_dir, "\n")