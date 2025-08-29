library(ggrepel)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(data.table)
library(ggpubr)
library(reshape2)
library(ggsci)
library ("phyloseq")
library("ggplot2")
library(permute)
library(lattice)
library(vegan) 
library(dplyr)
library(gapminder)
library(ggprism)
library(rstatix)
library(ggpmisc)
library(limma)
library(tidyr)
setwd('D:/R_project/Microbe cfRNA/MicrobeRNA_R')


# Main plot and supplement plots
load('../image/2.Multi train_CA.Rdata')
save.image('../image/3.Discribed picture')
load('../image/3.Discribed picture')

#### PAN VS NOR CA volcano ####
res_dds_pan = fread('/mnt/data3/yiyonghao/MicroRNA/process_file/0819/dds_result/PAN_vs_NOR_.csv')
res_dds_RS = fread('/mnt/data3/yiyonghao/MicroRNA/process_file/0819/dds_result/Other_vs_NOR_.csv')
dds_list = list(CA = res_dds_pan,RS = res_dds_RS)
for (i in names(dds_list)) {
  res_dds_pan = dds_list[[i]]
  
  res_dds_pan = res_dds_pan %>% data.frame %>% na.omit() %>% mutate(change = case_when(
    padj < 0.05 & log2FoldChange > .5 ~ "UP",
    padj < 0.05 & log2FoldChange < -.5 ~ "DOWN",
    TRUE ~ "Not_sig"
  ))
  print(i)
  print(table(res_dds_pan$change))
  }

for (i in names(dds_list)) {
  res_dds_pan = dds_list[[i]]
  
  res_dds_pan = res_dds_pan %>% data.frame %>% na.omit() %>% mutate(change = case_when(
    pvalue < 0.05 & log2FoldChange > .5 ~ "UP",
    pvalue < 0.05 & log2FoldChange < -.5 ~ "DOWN",
    TRUE ~ "Not_sig"
  ))
  res_dds_pan$change = factor(res_dds_pan$change,level = c('DOWN','Not_sig','UP'))
  res_dds_pan %$% change %>% table()
  res_dds_pan$change %>% table()
  # if(i == 'RS'){
  #   res_dds_pan = res_dds_pan[res_dds_pan$V1 != 'Picrophilus',]
  # }
  rownames(res_dds_pan) = res_dds_pan[,1];res_dds_pan = res_dds_pan[,-1]
  
  # top6 significant features
  sig = res_dds_pan %>% na.omit() %>% filter(change != "Not_sig")
  sig %<>% mutate(abs_FC = abs(log2FoldChange))
  sig = sig[order(sig$abs_FC, decreasing = T),]
  sig_head = sig %>% head(6)
  sig_head$change = factor(sig_head$change,level = c('DOWN','UP'))
  g = ggplot(res_dds_pan, aes(x = log2FoldChange, y = -log10(pvalue), color = change)) +
    geom_point(data = res_dds_pan %>% filter(change == "Not_sig"), color = "grey90") +
    geom_point(data = res_dds_pan %>% filter(change != "Not_sig"), aes(color = change, alpha = 0.85)) +
    geom_point(data = sig_head, aes(x = log2FoldChange, y = -log10(pvalue), fill = change), color = "black",shape = 21, size = 2.5) + 
    scale_color_manual(values = c("#003366", "#990033")) +
    geom_vline(xintercept = 0.5, color = "#990033") + 
    geom_vline(xintercept = -0.5, color = "#003366") +
    geom_hline(yintercept = -log10(0.05), color = "grey50") +
    scale_fill_manual(values = c("DOWN"="#003366","UP"="#990033"))+


    geom_text_repel(data = sig_head, aes(x = log2FoldChange,
                                         y = -log10(pvalue),label=rownames(sig_head),
                                         # size=0.5,
                                         # family="italic",
                                         fontface="italic",
    ),
    box.padding=unit(0.5, "lines"),
    point.padding=unit(0.1, "lines"),
    segment.size = 0.1,
    size = 3,
    arrow = arrow(length=unit(0.01, "npc")),
    force = 1.5, max.iter = 3e3,
    color = "black")+
    theme_bw()+
    theme(legend.position = "none", 
          axis.title= element_text( size = 12),
          plot.background = element_blank(), 
          # panel.grid.major = element_blank(), 
          # panel.grid.minor = element_blank(), 
          plot.title =element_text(size = 15, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 10,color = 'black'))+
    if (i == 'RS') {
      labs(title = "Resp. vs. NOR", x = "log2(FoldChange)", y = "-log10(FDR)") 
    }else{
      labs(title = "Pan vs. NOR", x = "log2(FoldChange)", y = "-log10(FDR)") 
    }
    if (i == 'RS') {
      ggsave(plot = g,filename = '/mnt/data3/yiyonghao/MicroRNA/process_file/0819/Plot/Resp vs Nor vol.pdf',width = 5.5,height = 5.5)
    }else{
      ggsave(plot = g,filename = '/mnt/data3/yiyonghao/MicroRNA/process_file/0819/Plot/Pan vs Nor vol.pdf',width = 5.5,height = 5.5)
    }
  
}


#### stack top 6 ###
CA_cohort_feature = readRDS('/mnt/data3/yiyonghao/MicroRNA/process_file/0819/CA_cohort_0820.rds')
RS_cohort_feature = read_rds('/mnt/data3/yiyonghao/MicroRNA/process_file/0819/RS_cohort_0820.rds') %>%
  select(classification_name,contains('PN'))
rownames(RS_cohort_feature) = RS_cohort_feature[,1];RS_cohort_feature = RS_cohort_feature[,-1]
df.all = CA_cohort_feature
rownames(df.all) = df.all[,1];df.all = df.all[,-1]
df.all = bind_cols(df.all,RS_cohort_feature)

df_all.ordered <- df.all[order(rowSums(df.all), decreasing = T),]
df_all.ordered_feature = rownames(df_all.ordered[1:7,])

df.all_top <- df_all.ordered %>%
  rownames_to_column("genus") %>%
  pivot_longer(cols = -genus, names_to = "variable", values_to = "value")
df.all_top$disease = lapply(df.all_top$variable,function(x) strsplit(x,'_')[[1]][[1]])
df.all_top$disease = gsub(' ','',df.all_top$disease)

df.all_top.sum <- df.all_top %>% 
  dplyr::select(1,4,3) %>% 
  group_by(genus,disease) %>% 
  summarise(sum=sum(value))

df.all_top.sum$genus = df.all_top.sum$genus  %>% as.character()
df.all_top.sum$label = ifelse(df.all_top.sum$genus %in% df_all.ordered_feature,df.all_top.sum$genus,'others')
df.all_top.sum$label = factor(df.all_top.sum$label,level = c(df_all.ordered_feature,'others'))
df.all_top.sum$disease = factor(df.all_top.sum$disease,levels = c('BRC','CRC','GC','HCC','MEN','LC','NOR','PN'))
p = df.all_top.sum %>% ggplot(data = ., aes(x = disease, y = sum, fill = label)) + 
  geom_col(position = "fill", width = 0.85) +
  scale_fill_manual(values = c("#9B3A4D", "#E2AE79", "#D0DCAA", "#F0EEBB", "#8CBDA7", "#566CA5", "#70A0AC",'lightgrey')) +
  theme_classic() +
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.text = element_text(color = 'black',size = 10))+
  labs(y = 'Ratio %')
ggsave(p,filename = "/mnt/data3/yiyonghao/MicroRNA/process_file/0819/Plot/stack_plot_top6.pdf", width = 7.5, height = 4.5)
getwd()



#### beta diversity ###
RS_RF = fread('/mnt/data3/yiyonghao/MicroRNA/process_file/0819/RS_selected_Masslin2.csv') %>% as.data.frame()
CA_RF = fread('/mnt/data3/yiyonghao/MicroRNA/process_file/0819/CA_selected_Masslin2.csv') %>% as.data.frame()
rownames(RS_RF) = RS_RF[,1];RS_RF = RS_RF[,-1]
rownames(CA_RF) = CA_RF[,1];CA_RF = CA_RF[,-1]
library(vegan)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(egg)
beta_list = list(CA = CA_RF, RS = RS_RF)

# ===== NEW STEP 1: 先预计算两个队列的全局 PCoA1/2 范围，用于统一坐标 =====
get_pcoa_xy <- function(mat) {
  mat <- mat[, sapply(mat, is.numeric), drop = FALSE]
  mat <- na.omit(mat)
  if (nrow(mat) < 3) return(list(x=numeric(0), y=numeric(0)))
  d <- vegdist(mat, method = "bray")
  pcoa <- cmdscale(d, k = 3, eig = TRUE)
  list(x = pcoa$points[,1], y = pcoa$points[,2])
}

xy1 <- get_pcoa_xy(beta_list$CA)
xy2 <- get_pcoa_xy(beta_list$RS)

x_all <- c(xy1$x, xy2$x)
y_all <- c(xy1$y, xy2$y)
stopifnot(length(x_all) > 0 && length(y_all) > 0)

xlim_glob <- range(x_all, na.rm = TRUE)
ylim_glob <- range(y_all, na.rm = TRUE)

# 给一点边距，避免点贴边
pad <- 0.05
xspan <- diff(xlim_glob); yspan <- diff(ylim_glob)
xlim_glob <- c(xlim_glob[1] - pad * xspan, xlim_glob[2] + pad * xspan)
ylim_glob <- c(ylim_glob[1] - pad * yspan, ylim_glob[2] + pad * yspan)

# ===== 循环开始 =====
for (i in names(beta_list)) {
  mat <- beta_list[[i]]
  mat <- mat[, sapply(mat, is.numeric), drop = FALSE]
  mat <- na.omit(mat)

  # 距离矩阵（Bray-Curtis）
  d <- vegdist(mat, method = "bray")

  # 样本分组
  sd <- data.frame(sample = rownames(mat))
  sd$group <- vapply(sd$sample, function(x) strsplit(x, "_")[[1]][1], character(1))

  # PCoA
  set.seed(42)
  pcoa <- cmdscale(d, k = 3, eig = TRUE)
  pc12 <- as.data.frame(pcoa$points[, 1:2])
  colnames(pc12) <- c("pc_x", "pc_y")
  pc_var <- round(pcoa$eig / sum(pcoa$eig) * 100, 2)
  pc12$sample <- rownames(pc12)
  pc12 <- merge(pc12, sd, by = "sample")

  if (i == "CA") {
    pc12$group <- factor(pc12$group, level = c('NOR','LC','CRC','GC','HCC','BRC','MEN','PN'))
  } else {
    pc12$group <- factor(pc12$group, level = c('NOR','LC','PN'))
  }

  # === 统计检验 ===
  set.seed(42)
  ad <- adonis2(d ~ group, data = sd, permutations = 999, by = "margin")
  R2 <- ad$R2[1]
  p_permanova <- ad[1, "Pr(>F)"]

  bd <- betadisper(d, sd$group)
  set.seed(42)
  bd_perm <- permutest(bd, permutations = 999)
  p_betadisper <- bd_perm$tab[1, "Pr(>F)"]

  # out_prefix <- if (i == "CA") "../plot/NC" else "../plot/RS"
  # write.csv(as.data.frame(ad), paste0(out_prefix, "_PERMANOVA.csv"), quote = FALSE)
  # capture.output(anova(bd), file = paste0(out_prefix, "_betadisper_anova.txt"))
  # capture.output(bd_perm, file = paste0(out_prefix, "_betadisper_permutest.txt"))

  pc12 <- pc12 %>%
    group_by(group) %>%
    mutate(x_mean = mean(pc_x), y_mean = mean(pc_y))

  p_pcoa <- ggplot(pc12, aes(x = pc_x, y = pc_y, color = group)) +
    geom_point(size = .75) +
    geom_segment(aes(x = x_mean, y = y_mean, xend = pc_x, yend = pc_y, color = group)) +
    stat_ellipse(geom = "polygon", level = 0.9, linetype = 2, size = 0.5,
                 aes(fill = group), alpha = 0.1, show.legend = TRUE) +
    xlim(xlim_glob) + ylim(ylim_glob) +
    coord_fixed(ratio = 1.25, clip = "on") +
    xlab(paste0("PCoA1 (", pc_var[1], "%)")) +
    ylab(paste0("PCoA2 (", pc_var[2], "%)")) +
    scale_fill_manual(values = c("#9B3A4D", "#E2AE79", "#D0DCAA", "#F0EEBB", "#8CBDA7", "#566CA5", "#70A0AC", "#7d3f98")) +
    scale_color_manual(values = c("#9B3A4D", "#E2AE79", "#D0DCAA", "#F0EEBB", "#8CBDA7", "#566CA5", "#70A0AC", "#7d3f98")) +
    theme_bw() +
    guides(color = "none", fill = guide_legend(override.aes = list(alpha = 1))) +
    labs(
      title   = "PCoA",
      caption = sprintf("R² = %.3f, p = %.3g", R2, p_permanova,p_betadisper)
    ) +
    theme(
      legend.title  = element_blank(),
      legend.position = "right",
      panel.background = element_blank(),
      plot.title    = element_text(size = 15, color = "black", hjust = 0.5, face = "bold"),
      plot.margin   = margin(5.5, 5.5, 5.5, 5.5)
    )

  p_pcoa_fixed <- egg::set_panel_size(
    p_pcoa,
    width  = grid::unit(7, "cm"),
    height = grid::unit(7, "cm")
  )

  ggsave(filename = ifelse(i == "CA", "/mnt/data3/yiyonghao/MicroRNA/process_file/0819/Plot/NC_PCoA.pdf", "/mnt/data3/yiyonghao/MicroRNA/process_file/0819/Plot/RS_PCoA.pdf"),
         plot = p_pcoa_fixed,
         width = 12, height = 12, units = "cm", dpi = 300)

  bd_df <- data.frame(sample = names(bd$distances),
                      group = sd$group,
                      dist_to_centroid = bd$distances)

  p_disp <- ggplot(bd_df, aes(x = group, y = dist_to_centroid, fill = group)) +
    geom_boxplot(width = .35) +
    scale_fill_manual(values = c("#9B3A4D", "#E2AE79", "#D0DCAA", "#F0EEBB", "#8CBDA7", "#566CA5", "#70A0AC", "#7d3f98")) +
    labs(x = "", y = "Distance to centroid", title = "β-diversity dispersion (Bray–Curtis)",
         caption = sprintf("betadisper (permutation test) p = %.3g", p_betadisper)) +
    theme_bw() +
    theme(axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = 0.5),
          plot.margin = margin(5.5, 5.5, 5.5, 5.5))

  p_disp_fixed <- egg::set_panel_size(p_disp,
                                      width = grid::unit(10, "cm"),
                                      height = grid::unit(5.5, "cm"))

  ggsave(filename = ifelse(i == "CA", "/mnt/data3/yiyonghao/MicroRNA/plot/NC_betadisper.pdf", "/mnt/data3/yiyonghao/MicroRNA/plot/RS_betadisper.pdf"),
         plot = p_disp_fixed,
         width = 15, height = 12, units = "cm", dpi = 300)
}




# drawing feature PASS heatmap
library(scales)
feature_CA = read.csv('/mnt/data3/yiyonghao/MicroRNA/process_file/0819/CA_Multi_featureImp.csv')
rownames(feature_CA) = feature_CA[,1];feature_CA = feature_CA[,-1]
row_names_list <- lapply(feature_CA, function(col) {
  rownames(feature_CA)[!is.na(col)]
})
features <- Reduce(intersect, row_names_list)
feature_RS = read.csv('/mnt/data3/yiyonghao/MicroRNA/process_file/0819/RS_Multi_featureImp.csv')
rownames(feature_RS) = feature_RS[,1];feature_RS = feature_RS[,-1]
row_names_list <- lapply(feature_RS, function(col) {
  rownames(feature_RS)[!is.na(col)]
})
tmp <- Reduce(intersect, row_names_list)
features = c(features,tmp)  %>%  unique()

masslin2_CA = fread('/mnt/data3/yiyonghao/MicroRNA/process_file/0819/CA_masslin2/Maaslin2/all_results.tsv')
masslin2_RS = fread('/mnt/data3/yiyonghao/MicroRNA/process_file/0819/RS_masslin2/Maaslin2/all_results.tsv')
masslin2_RS = masslin2_RS[masslin2_RS$value == 'PN',]

masslin2 = bind_rows(masslin2_RS,masslin2_CA)
masslin2 = masslin2[masslin2$feature %in% features,]
masslin2$sig = -log(masslin2$qval)*sign(masslin2$coef)
masslin2$label = ifelse(masslin2$coef>0,'+','-')
masslin2$feature %>% unique() %>% length()

top50_feature = masslin2 %>% 
  group_by(feature) %>%
  summarise(mean = mean(abs(sig), na.rm = TRUE)) %>%
  arrange(desc(mean)) %>%
  slice_head(n = 50)
masslin2 = masslin2 %>%
  filter(feature %in% top50_feature$feature)
# Cluster
mat <- masslin2 %>%
  pivot_wider(
    id_cols     = feature,
    names_from  = value,
    values_from = sig
  ) %>%
  column_to_rownames("feature") %>%
  as.matrix()

row_hc <- hclust(dist(mat, method = "euclidean"), method = "complete")
col_hc <- hclust(dist(t(mat), method = "euclidean"), method = "complete")

row_order <- row_hc$labels[row_hc$order]
col_order <- col_hc$labels[col_hc$order]
masslin2_clust <- masslin2 %>%
  mutate(
    feature = factor(feature, levels = row_order),
    value   = factor(value,   levels = col_order)
  )
# filter features not exitsed in all phenotypes
masslin2_clust = masslin2_clust %>%
  group_by(feature) %>%
  filter(n() >= 7) %>%
  ungroup()

p = ggplot(masslin2_clust, aes(x = value, y = feature, fill = sig)) +
  geom_tile(color = 'grey80',size = 0.3) +
  scale_fill_gradient2(
    low      = "#004b79",
    mid      = "white",
    high     = "#7f181b",
    midpoint = 0,
    limits   = c(-15, 15),
    oob      = squish
  ) +theme_minimal()+
  geom_text(aes(label = label),
            size  = 3,
            color = "black")+  
  theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  panel.grid = element_blank(),
  axis.title = element_blank(),
  axis.text = element_text(color = 'black',size = 8),
  legend.title = element_text(angle = 90)
)+
  labs(fill = "-log(FDR)*sign(coef)")

ggsave(plot = p,
       filename = '/mnt/data3/yiyonghao/MicroRNA/process_file/0819/Plot/feature_heatmap.pdf',
       height = 8,
       width = 3.5)


CA_feature = fread('/mnt/data3/yiyonghao/MicroRNA/process_file/0819/Table/selected_features_CA.csv',header = T)
CA_feature = do.call(c,CA_feature)  %>% unique()
RS_feature = fread('/mnt/data3/yiyonghao/MicroRNA/process_file/0819/Table/selected_features_RS.csv',header = T)
RS_feature = do.call(c,RS_feature)  %>% unique()
features = c(CA_feature,RS_feature) %>% unique()
