library(scplotter)
library(DoubletFinder)
if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(patchwork))install.packages("patchwork")
if(!require(R.utils))install.packages("R.utils")
##if(!require(mindr))install.packages("mindr")
if(!require(tidyverse))install.packages("tidyverse")
if(!require(hdf5r))install.packages("hdf5r")
##if(!require(DoMultiBarHeatmap))install.packages("DoMultiBarHeatmap")
if(!require(ggplot2))install.packages("ggplot2")
library(ggpubr)
library(viridis)
library(patchwork) 
library(harmony)
library(scales)
library(devtools)
library(presto)
library(ggsci)
library(monocle3)
library(ggrepel)
library(tidydr)
library(reshape2)
library(ComplexHeatmap)
library(rcartocolor)
library(ggh4x)
library(cowplot)
library(ggsignif)
library(DoubletFinder)
library("FactoMineR")
library(hexbin)
library(RColorBrewer)
library(GSVA)
library(AUCell)
library(UCell)
library("factoextra")
ROOT_DIR <- normalizePath(".", mustWork = FALSE)
DATA_DIR <- file.path(ROOT_DIR, "data")
RESULT_DIR <- file.path(ROOT_DIR, "result")
setwd(DATA_DIR)

dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULT_DIR, "Fig6"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULT_DIR, "Fig6-SF"), recursive = TRUE, showWarnings = FALSE)



# Fig6 scoring analysis
# Inputs and gene sets

# Load single-cell object
sc_data=readRDS(file.path(DATA_DIR, "GSE201333_RAW", "GSM6058681_TabulaSapiens.h5ad", "recluster_key_celltype_tissue.rds"))
# Check identities
unique(Idents(sc_data))


# Load module genes
geneInfo=read.csv(file.path(DATA_DIR, "WGCNA", "no_combat", "geneInfo.csv"),row.names = 1,header = T)

IPAH_genes <- rownames(geneInfo)[ geneInfo$Module %in% c("midnightblue","black","lightcyan") ]
CHD_genes <- rownames(geneInfo)[ geneInfo$Module  %in% c( "purple","red") ]
SLE_genes   <- rownames(geneInfo)[ geneInfo$Module  %in% c( "greenyellow","tan","cyan") ]

# Random control genes
# Set seed
set.seed(42)
random_select <- sample(rownames(geneInfo), size = 150, replace = FALSE)


# Check gene-set sizes
length(IPAH_genes)
length(CHD_genes)
length(SLE_genes)
length(random_select)


# Harmonize gene symbols
IPAH_genes=gsub("\\.", "-", IPAH_genes)
CHD_genes=gsub("\\.", "-", CHD_genes)
SLE_genes=gsub("\\.", "-", SLE_genes)
random_select=gsub("\\.", "-", random_select)

# Module scoring
# AddModuleScore

sc_data <- AddModuleScore(
  object   = sc_data,
  features = list(IPAH_genes),
  name     = "Add_IPAH_score",
  ctrl     = 300,
  seed     = 42
)

sc_data <- AddModuleScore(
  object   = sc_data,
  features = list(SLE_genes),
  name     = "Add_SLE_score",
  ctrl     = 300,
  seed     = 42
)

sc_data <- AddModuleScore(
  object   = sc_data,
  features = list(CHD_genes),
  name     = "Add_CHD_score",
  ctrl     = 300,
  seed     = 42
)

sc_data <- AddModuleScore(
  object   = sc_data,
  features = list(random_select),
  name     = "Add_random_score",
  ctrl     = 300,
  seed     = 42
)


# Module scoring

# UCell scoring
# Build expression matrix
expr_mat <- GetAssayData(sc_data, slot = "data")

# Score signatures with UCell
features_list <- list(
  UCELL_IPAH = IPAH_genes,
  UCELL_SLE  = SLE_genes,
  UCELL_CHD  = CHD_genes,
  UCELL_random  = random_select
)
ucell_scores <- ScoreSignatures_UCell(
  expr_mat,
  features = features_list
)

# Save UCell scores
saveRDS(ucell_scores, file = file.path(RESULT_DIR, "ucell_scores_celltype.rds"))


# Module scoring

# AUCell scoring
# Build expression matrix
expr_mat <- as.matrix(GetAssayData(sc_data, assay = "RNA", slot = "data"))

# Build rankings
cells_rankings <- AUCell_buildRankings(
  expr_mat,
  nCores    = 1,
  plotStats = FALSE
)

# Score signatures with AUCell
features_list <- list(
  AUCell_IPAH = IPAH_genes,
  AUCell_SLE  = SLE_genes,
  AUCell_CHD  = CHD_genes,
  AUCell_random=random_select
)
auc_res <- AUCell_calcAUC(
  geneSets   = features_list,
  rankings   = cells_rankings,
  aucMaxRank = nrow(expr_mat) * 0.05
)

# Extract score matrix
auc_mat <- getAUC(auc_res)

# Save AUCell scores
saveRDS(auc_mat, file = file.path(RESULT_DIR, "aucell_scores_celltype.rds"))

# Module scoring


# Merge scores into metadata
ucell_scores=readRDS(file.path(RESULT_DIR, "ucell_scores_celltype.rds"))
aucell_scores=readRDS(file.path(RESULT_DIR, "aucell_scores_celltype.rds"))
aucell_scores=t(aucell_scores)

sc_data <- AddMetaData(sc_data, metadata = as.data.frame(ucell_scores))
sc_data <- AddMetaData(sc_data, metadata = as.data.frame(aucell_scores))


# Build long-format tables
# IPAH
IPAH_long <- sc_data@meta.data %>%
  rownames_to_column(var = "barcode") %>%
  select(
    barcode,
    donor,
    tissue_celltype,
    organ_tissue,
    celltype_integration,
    contains("IPAH")
  ) %>%
  pivot_longer(
    cols      = contains("IPAH"),
    names_to  = "Algorithm",
    values_to = "Score"
  )


# SLE
SLE_long <- sc_data@meta.data %>%
  rownames_to_column(var = "barcode") %>%
  select(
    barcode,
    donor,
    tissue_celltype,
    organ_tissue,
    celltype_integration,
    contains("SLE")
  ) %>%
  pivot_longer(
    cols      = contains("SLE"),
    names_to  = "Algorithm",
    values_to = "Score"
  )


# CHD
CHD_long <- sc_data@meta.data %>%
  rownames_to_column(var = "barcode") %>%
  select(
    barcode,
    donor,
    tissue_celltype,
    organ_tissue,
    celltype_integration,
    contains("CHD")
  ) %>%
  pivot_longer(
    cols      = contains("CHD"),
    names_to  = "Algorithm",
    values_to = "Score"
  )

# Random
random_long <- sc_data@meta.data %>%
  rownames_to_column(var = "barcode") %>%
  select(
    barcode,
    donor,
    tissue_celltype,
    organ_tissue,
    celltype_integration,
    contains("random")
  ) %>%
  pivot_longer(
    cols      = contains("random"),
    names_to  = "Algorithm",
    values_to = "Score"
  )

# Load cached long-format tables
IPAH_long=readRDS(file.path(DATA_DIR, "IPAH_long_celltype.rds"))
SLE_long=readRDS(file.path(DATA_DIR, "SLE_long_celltype.rds"))
CHD_long=readRDS(file.path(DATA_DIR, "CHD_long_celltype.rds"))
random_long=readRDS(file.path(DATA_DIR, "random_long_celltype.rds"))


# Module scoring






# Boxplots across methods

# Plot configuration
# =========================================================================================
config <- list(
  data = list(
    IPAH_long = IPAH_long,
    CHD_long  = CHD_long,
    SLE_long  = SLE_long,
    random_long = random_long,
    sc_data   = sc_data
  ),
  alg = list(
    levels = c(
      "Add_IPAH_score1", "UCELL_IPAH_UCell", "AUCell_IPAH",
      "Add_SLE_score1",  "UCELL_SLE_UCell",  "AUCell_SLE",
      "Add_CHD_score1", "UCELL_CHD_UCell", "AUCell_CHD"
    ),
    labels = rep(c("AddModuleScore", "UCell", "AUCell"), 3),
    colors = c(
      AddModuleScore = "#ea5b57",
      UCell          = "#fac03d",
      AUCell         = "#0c3c5f"
    ),
    group_colors=c(IPAH =  "#008B45FF",CHD="#EE0000FF",SLE="#631879FF")
  ),
  tissue_color = pal_igv("default")(24),
  celltype_integration_color=pal_igv("default")(51)[c(26:45)],
  jitter_prop = 0.1,
  out = list(
    IPAH = list(
      box_ct = file.path(RESULT_DIR, "Fig6-SF"),
      box_ot = file.path(RESULT_DIR, "Fig6-SF"),
      umap   = file.path(RESULT_DIR, "Fig6")
    ),
    SLE = list(
      box_ct = file.path(RESULT_DIR, "Fig6-SF"),
      box_ot = file.path(RESULT_DIR, "Fig6-SF"),
      umap   = file.path(RESULT_DIR, "Fig6")
    ),
    CHD = list(
      box_ct = file.path(RESULT_DIR, "Fig6-SF"),
      box_ot = file.path(RESULT_DIR, "Fig6-SF"),
      umap   = file.path(RESULT_DIR, "Fig6")
    )
  ),
  dims = list(
    box_ct = c(width = 24, height = 8),
    box_ot = c(width = 24, height = 8),
    umap   = c(width = 16, height = 8)
  )
)

for (mod in names(config$out)) {
  for (sub in names(config$out[[mod]])) {
    dir.create(config$out[[mod]][[sub]], recursive = TRUE, showWarnings = FALSE)
  }
}

# Helper functions
# Build plotting data
# Helper functions
# Get algorithm prefix
get_prefix <- function(x) {
  sub("_.*", "", as.character(x))
}

make_plot_df <- function(long_df, random_long, config) {
  long_df2 <- long_df %>%
    mutate(prefix = get_prefix(Algorithm))
  random_long2 <- random_long %>%
    mutate(prefix = get_prefix(Algorithm))
  
  merged <- merge(
    long_df2,
    random_long2[, c("barcode", "prefix", "Score")],
    by = c("barcode", "prefix"),
    all.x = TRUE,
    suffixes = c("", "_rand")
  )
  if ("Score_rand" %in% names(merged)) {
    merged$random_score <- merged$Score_rand
    merged$Score_rand <- NULL
  } else if ("Score" %in% names(random_long2)) {
    # fallback: if merge didn't add suffix (unlikely), try to copy
    merged$random_score <- merged$Score
    warning("make_plot_df: fallback copying Score to random_score (check merge logic).")
  } else {
    stop("make_plot_df: could not find random score column after merge.")
  }
  
  missing_rand <- sum(is.na(merged$random_score))
  if (missing_rand > 0) {
    warning(sprintf(
      "make_plot_df: %d/%d rows have no matching random_score (barcode+prefix).",
      missing_rand, nrow(merged)
    ))
  }
  
  result <- merged %>%
    mutate(
      Algorithm = factor(
        Algorithm,
        levels = config$alg$levels,
        labels = config$alg$labels
      )
    ) %>%
    group_by(tissue_celltype,Algorithm) %>%
    mutate(.m = mean(Score, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(tissue_celltype = fct_reorder(tissue_celltype, .m, .desc = TRUE)) %>%
    select(-.m) %>%
    group_by(organ_tissue, Algorithm) %>%
    mutate(.o = mean(Score, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(organ_tissue = fct_reorder(organ_tissue, .o, .desc = TRUE)) %>%
    select(-.o, -prefix)
  
  return(result)
}

config$plots <- list(
  IPAH_plot = make_plot_df(config$data$IPAH_long, config$data$random_long, config),
  SLE_plot  = make_plot_df(config$data$SLE_long,  config$data$random_long, config),
  CHD_plot  = make_plot_df(config$data$CHD_long,  config$data$random_long, config)
)
plots=config$plots
saveRDS(plots,file.path(RESULT_DIR, "cell_score_plots_key_celltype.rds"))

# Helper functions
# Load cached plotting data
# Helper functions
plots=readRDS(file.path(DATA_DIR, "cell_score_plots_tissue.rds"))
plots=readRDS(file.path(RESULT_DIR, "cell_score_plots_key_celltype.rds"))

# =========================================================================================



# Grouped score plots
# Use cached tissue plots or key-celltype plots as needed

color_map <- c(config$alg$colors, Random = "#888888")
random_color <- "#888888"
signature_summary_color <- "purple"

# Keep UCell and AUCell only
plots$IPAH_plot= plots$IPAH_plot[ plots$IPAH_plot$Algorithm != "AddModuleScore",]
plots$SLE_plot= plots$SLE_plot[ plots$SLE_plot$Algorithm != "AddModuleScore",]
plots$CHD_plot= plots$CHD_plot[ plots$CHD_plot$Algorithm != "AddModuleScore",]


# Disease-specific plotting tables
plot_list <- list(
  IPAH = plots$IPAH_plot,
  SLE  = plots$SLE_plot,
  CHD  = plots$CHD_plot
)

# Plot signature vs random by group

plot_sig_random_by_group <- function(plot_list, config, group_var = "celltype_integration") {
  stopifnot(is.list(plot_list), is.character(group_var), length(group_var) == 1)
  group_var_sym <- rlang::sym(group_var)
  `%||%` <- function(x, y) if (!is.null(x)) x else y
  
  for (disease in names(plot_list)) {
    plot_df <- plot_list[[disease]]
    if (!group_var %in% names(plot_df)) {
      warning(sprintf("Skip '%s': grouping column '%s' not found.", disease, group_var))
      next
    }
    
    algorithms <- unique(plot_df$Algorithm)
    out_dir    <- config$out[[disease]]$box_ot
    dims       <- config$dims$box_ot
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    message(sprintf("[START SUMMARY AGG] %s signature vs random means per '%s' across algorithms",
                    disease, group_var))
    
    summary_means <- plot_df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c(group_var, "Algorithm")))) %>%
      dplyr::summarise(
        signature_score = mean(Score, na.rm = TRUE),
        no_sig_score    = mean(random_score, na.rm = TRUE),
        .groups = "drop"
      )
    utils::write.csv(
      summary_means,
      file.path(out_dir, paste0(disease, "_Algorithm_mean.csv")),
      row.names = TRUE
    )
    
    organ_order_df <- summary_means %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_var))) %>%
      dplyr::summarise(mean_sig_all = mean(signature_score, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(mean_sig_all))
    organ_order <- organ_order_df[[group_var]]
    message(sprintf("  %s order for summary %s: %s",
                    group_var, disease, paste(utils::head(organ_order, 5), collapse = ", ")))
    
    summary_long <- summary_means %>%
      tidyr::pivot_longer(
        cols = c(signature_score, no_sig_score),
        names_to = "Type",
        values_to = "Value"
      ) %>%
      dplyr::mutate(
        Type = dplyr::recode(Type, signature_score = "Signature", no_sig_score = "Random"),
        Type = factor(Type, levels = c("Signature", "Random")),
        !!group_var_sym := factor(.data[[group_var]], levels = organ_order)
      )
    
    p_summary_agg <- ggplot2::ggplot(
      summary_long,
      ggplot2::aes(x = .data[[group_var]], y = Value, fill = Type)
    ) +
      ggplot2::geom_boxplot(
        position = ggplot2::position_dodge(0.8),
        outlier.shape = NA,
        alpha = 0.8,
        size = 0.8
      ) +
      ggplot2::geom_jitter(
        ggplot2::aes(color = Type),
        position = ggplot2::position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
        size = 1.5,
        alpha = 0.9,
        show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(values = c(Signature = signature_summary_color, Random = random_color)) +
      ggplot2::scale_color_manual(values = c(Signature = signature_summary_color, Random = random_color)) +
      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold", size = 10),
        axis.text.y = ggplot2::element_text(size = 12, face = "bold"),
        legend.title = ggplot2::element_blank(),
        legend.position = "right"
      ) +
      ggplot2::labs(
        x = "Celltype Type",
        y = "Mean Score",
        title = paste0(disease, " — Aggregated Signature vs Random by Celltype"),
        subtitle = "Each box summarizes three algorithm-specific means; tissues ordered by signature"
      )
    
    fname_base_summary <- paste0("MODULE_Celltype_", disease, "_aggregated_sig_vs_random_by_tissue")
    ggplot2::ggsave(
      plot = p_summary_agg,
      filename = file.path(out_dir, paste0(fname_base_summary, ".pdf")),
      width = dims["width"], height = dims["height"], units = "in"
    )
    ggplot2::ggsave(
      plot = p_summary_agg,
      filename = file.path(out_dir, paste0(fname_base_summary, ".tif")),
      width = dims["width"], height = dims["height"], units = "in",
      device = "tiff", dpi = 300
    )
    message(sprintf("[DONE SUMMARY AGG] %s saved to %s", disease, out_dir))
    
    for (alg in algorithms) {
      message(sprintf("[START] %s - algorithm %s per-%s signature vs random",
                      disease, alg, group_var))
      alg_df <- dplyr::filter(plot_df, Algorithm == alg)
      if (nrow(alg_df) == 0) next
      
      organ_order_alg <- alg_df %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group_var))) %>%
        dplyr::summarise(mean_sig = mean(Score, na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(mean_sig)) %>%
        dplyr::pull(!!group_var_sym)
      
      long_plot_df <- alg_df %>%
        dplyr::mutate(!!group_var_sym := factor(.data[[group_var]], levels = organ_order_alg)) %>%
        tidyr::pivot_longer(
          cols = c(Score, random_score),
          names_to  = "Type",
          values_to = "Value"
        ) %>%
        dplyr::mutate(
          Type = dplyr::recode(Type, Score = "Signature", random_score = "Random"),
          Type = factor(Type, levels = c("Signature", "Random"))
        )
      
      sig_color <- config$alg$colors[[as.character(alg)]] %||% "#000000"
      fill_map  <- c(Signature = sig_color, Random = random_color)
      
      p_alg <- ggplot2::ggplot(
        long_plot_df,
        ggplot2::aes(x = .data[[group_var]], y = Value, fill = Type)
      ) +
        ggplot2::geom_boxplot(
          position      = ggplot2::position_dodge(0.8),
          outlier.shape = NA,
          alpha         = 0.7,
          size          = 0.8,
          coef          = 0.5
        ) +
        ggplot2::scale_fill_manual(values = fill_map) +
        ggplot2::scale_y_continuous(
          limits = c(0, 0.1),
          breaks = seq(0, 0.1, by = 0.02)
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          panel.grid    = ggplot2::element_blank(),
          axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1, face = "bold", size = 10),
          axis.text.y   = ggplot2::element_text(size = 12, face = "bold"),
          legend.title  = ggplot2::element_blank(),
          legend.position = "right"
        ) +
        ggplot2::labs(
          x     = "Celltype",
          y     = "Score",
          title = paste0(disease, " — Algorithm: ", alg),
          subtitle = "Signature vs Random per tissue"
        )
      
      fname_base_alg <- paste0("MODULE_Celltype_", disease, "_", alg, "_sig_vs_random_by_tissue")
      ggplot2::ggsave(
        plot     = p_alg,
        filename = file.path(out_dir, paste0(fname_base_alg, ".pdf")),
        width    = dims["width"], height = dims["height"], units = "in"
      )
      ggplot2::ggsave(
        plot     = p_alg,
        filename = file.path(out_dir, paste0(fname_base_alg, ".tif")),
        width    = dims["width"], height = dims["height"], units = "in",
        device   = "tiff", dpi = 300
      )
      message(sprintf("[DONE] %s - algorithm %s saved to %s", disease, alg, out_dir))
    }
  }
  invisible(NULL)
}


# Plot by integrated cell type
# Rebuild input tables for UMAP plotting
plot_sig_random_by_group(plot_list, config, group_var = "celltype_integration")





# UMAP preparation
# Prepare UMAP data using AUCell scores
# Extract AUCell results
plot_IPAH=plots$IPAH_plot[plots$IPAH_plot$Algorithm =="AUCell",]
plot_SLE=plots$SLE_plot[plots$SLE_plot$Algorithm =="AUCell",]
plot_CHD=plots$CHD[plots$CHD_plot$Algorithm =="AUCell",]


# Merge the three score tables
# Helper to find barcode and score columns
pick_barcode_and_score <- function(df, group_tag) {
  nm_raw  <- names(df)
  nm_trim <- str_trim(nm_raw)
  
  bc_idx <- which(str_to_lower(nm_trim) == "barcode")
  if (length(bc_idx) != 1) {
    stop(sprintf("[%s] Cannot find a unique barcode column; candidates: %s",
                 group_tag,
                 paste(nm_raw[str_detect(nm_trim, "barcode", ignore_case=TRUE)], collapse=", ")))
  }
  bc_col <- nm_raw[bc_idx]
  
  sc_idx <- which(str_to_lower(nm_trim) == "score")
  if (length(sc_idx) == 0) {
    cand <- which(
      str_detect(nm_trim, regex("\\bscore\\b", ignore_case = TRUE)) &
        !str_detect(nm_trim, regex("random|z|std|scale", ignore_case = TRUE))
    )
    if (length(cand) == 0)
      stop(sprintf("[%s] Cannot find a Score column; current columns: %s", group_tag, paste(nm_raw, collapse=", ")))
    sc_idx <- cand[1]
  }
  sc_col <- nm_raw[sc_idx]
  
  df %>%
    select(
      barcode = all_of(bc_col),
      !!paste0("Score_", group_tag) := .data[[sc_col]]
    )
}

Score_all <- pick_barcode_and_score(plot_IPAH, "IPAH") %>%
  full_join(pick_barcode_and_score(plot_CHD, "CHD"), by = "barcode") %>%
  full_join(pick_barcode_and_score(plot_SLE, "SLE"), by = "barcode")


# Join UMAP coordinates and metadata
# UMAP coordinates
emb <- config$data$sc_data@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column("barcode")

# Metadata
meta <- FetchData(config$data$sc_data, vars = c("donor", "organ_tissue","celltype_integration")) %>%
  as.data.frame() %>%
  rownames_to_column("barcode")

# Merge all plotting inputs
umap_df <- emb %>%
  left_join(meta, by = "barcode") %>%
  inner_join(Score_all %>% select(barcode, Score_IPAH, Score_CHD, Score_SLE),
             by = "barcode")





# Combined score UMAP
# Harmonize UMAP column names
if ("UMAP_1" %in% names(umap_df) && !"umap_1" %in% names(umap_df)) umap_df <- rename(umap_df, umap_1 = UMAP_1)
if ("UMAP_2" %in% names(umap_df) && !"umap_2" %in% names(umap_df)) umap_df <- rename(umap_df, umap_2 = UMAP_2)

# Highlight fractions
top_ipah <- 0.35
top_chd  <- 0.15
top_sle  <- 0.45

# Thresholds
pos_ipah <- pmax(umap_df$Score_IPAH, 0)
pos_chd  <- pmax(umap_df$Score_CHD,  0)
pos_sle  <- pmax(umap_df$Score_SLE,  0)

q_ipah <- as.numeric(quantile(pos_ipah, probs = 1 - top_ipah, na.rm = TRUE))
q_chd  <- as.numeric(quantile(pos_chd,  probs = 1 - top_chd,  na.rm = TRUE))
q_sle  <- as.numeric(quantile(pos_sle,  probs = 1 - top_sle,  na.rm = TRUE))

mx_ipah <- max(pos_ipah, na.rm = TRUE); den_ipah <- max(1e-8, mx_ipah - q_ipah)
mx_chd  <- max(pos_chd,  na.rm = TRUE); den_chd  <- max(1e-8, mx_chd  - q_chd )
mx_sle  <- max(pos_sle,  na.rm = TRUE); den_sle  <- max(1e-8, mx_sle  - q_sle )

# Layered score plot
p_mix <- ggplot(umap_df, aes(x = umap_1, y = umap_2)) +
  stat_bin_hex(bins = 100, na.rm = TRUE, fill = "grey85", color = NA) +
  
  stat_summary_hex(
    aes(
      z     = pmax(Score_IPAH, 0),
      alpha = after_stat(pmin(pmax((pmax(..value.., 0) - q_ipah) / den_ipah, 0), 1))
    ),
    fun = mean, bins = 100, na.rm = TRUE,
    fill = "#EC926B", color = NA
  ) +
  stat_summary_hex(
    aes(
      z     = pmax(Score_CHD, 0),
      alpha = after_stat(pmin(pmax((pmax(..value.., 0) - q_chd) / den_chd, 0), 1))
    ),
    fun = mean, bins = 100, na.rm = TRUE,
    fill =  "#7DBFA6", color = NA
  ) +
  stat_summary_hex(
    aes(
      z     = pmax(Score_SLE, 0),
      alpha = after_stat(pmin(pmax((pmax(..value.., 0) - q_sle) / den_sle, 0), 1))
    ),
    fun = mean, bins = 100, na.rm = TRUE,
    fill = "#D98DBF", color = NA
  ) +
  scale_alpha(range = c(0, 1), guide = "none") +
  guides(alpha = "none")+
  labs(x = "UMAP1", y = "UMAP2") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.box = "vertical",
    legend.spacing.y = unit(4, "pt")
  )

# Legend scaffolds
cx <- mean(umap_df$umap_1, na.rm = TRUE)
cy <- mean(umap_df$umap_2, na.rm = TRUE)
leg_vals <- data.frame(x = cx, y = cy, v = seq(0, 1, length.out = 50))

# IPAH legend
p_mix <- p_mix +
  ggnewscale::new_scale_fill() +
  geom_point(
    data = leg_vals,
    aes(x = x, y = y, fill = v),
    inherit.aes = FALSE, size = 0, alpha = 0, show.legend = TRUE
  ) +
  scale_fill_gradient(
    low = "grey85", high = "#EC926B",
    name = "Score_IPAH",
    limits = c(0, 1), breaks = c(0, 1), labels = c("Low", "High"),
    guide = "colourbar"
  )

# CHD legend
p_mix <- p_mix +
  ggnewscale::new_scale_colour() +
  geom_point(
    data = leg_vals,
    aes(x = x, y = y, colour = v),
    inherit.aes = FALSE, size = 0, alpha = 0, show.legend = TRUE
  ) +
  scale_colour_gradient(
    low = "grey85", high =  "#7DBFA6",
    name = "Score_CHD",
    limits = c(0, 1), breaks = c(0, 1), labels = c("Low", "High"),
    guide = "colourbar"
  )

# SLE legend
p_mix <- p_mix +
  ggnewscale::new_scale_colour() +
  geom_point(
    data = leg_vals,
    aes(x = x, y = y, colour = v),
    inherit.aes = FALSE, size = 0, alpha = 0, show.legend = TRUE
  ) +
  scale_colour_gradient(
    low = "grey85", high = "#D98DBF",
    name = "Score_SLE",
    limits = c(0, 1), breaks = c(0, 1), labels = c("Low", "High"),
    guide = "colourbar"
  ) +
  theme(legend.title.position = "top")+theme_dr(
    xlength = 0.3, ylength = 0.3,
    arrow   = arrow(length = unit(0.2,"inches"), type = "closed"))+theme(panel.grid   = element_blank())


# Tissue/cell-type UMAP and final layout

umap_df$organ_tissue=factor(umap_df$organ_tissue,levels = unique(sc_data$organ_tissue),ordered = T)
umap_df$celltype_integration=factor(umap_df$celltype_integration,levels = unique(sc_data$celltype_integration),ordered = T)

p_tissue <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = celltype_integration)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_manual(values = config$celltype_integration_color) +
  guides(color = guide_legend(
    nrow = 5, ncol = 5,
    override.aes = list(size = 6, alpha = 1)
  )) +
  theme_dr(
    xlength = 0.3, ylength = 0.3,
    arrow   = arrow(length = unit(0.2,"inches"), type = "closed")
  ) +
  labs(x = "UMAP1", y = "UMAP2", title = "Tissue") +
  theme(
    panel.grid   = element_blank(),
    legend.title = element_blank(),
    legend.text  = element_text(size = 15),
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5)
  )

# Combine panels
combined <- p_mix + p_tissue +
  plot_layout(ncol = 2, widths = c(1,1), guides = "collect") &
  theme(legend.position = "bottom")


# Save figures
ggsave(file.path(RESULT_DIR, "Fig6", "Fig6E.tiff"), combined, width = 16, height = 9, dpi = 300, compression = "lzw")
ggsave(file.path(RESULT_DIR, "Fig6", "Fig6E.pdf"), combined, width = 16, height = 9, dpi = 300)



