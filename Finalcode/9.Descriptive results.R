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
res_dds_pan = fread('D:/R_project/Microbe cfRNA/MicrobeRNA_R/dds_result/PAN_vs_NOR_.csv')
res_dds_RS = fread('D:/R_project/Microbe cfRNA/MicrobeRNA_R/dds_result/Other_vs_NOR_.csv')
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
  
  res_dds_pan %$% change %>% table()
  res_dds_pan$change %>% table()
  if(i == 'RS'){
    res_dds_pan = res_dds_pan[res_dds_pan$V1 != 'Picrophilus',]
  }
  rownames(res_dds_pan) = res_dds_pan[,1];res_dds_pan = res_dds_pan[,-1]
  
  # top10 significant features
  sig = res_dds_pan %>% na.omit() %>% filter(change != "Not_sig")
  sig %<>% mutate(abs_FC = abs(log2FoldChange))
  sig = sig[order(sig$abs_FC, decreasing = T),]
  sig_head = sig %>% head(10)
  
  g = ggplot(res_dds_pan, aes(x = log2FoldChange, y = -log10(pvalue), color = change)) +
    geom_point(data = res_dds_pan %>% filter(change == "Not_sig"), color = "grey90") +
    geom_point(data = res_dds_pan %>% filter(change != "Not_sig"), aes(color = change, alpha = 0.85)) +
    geom_point(data = sig_head, aes(x = log2FoldChange, y = -log10(pvalue), fill = change), color = "black",shape = 21, size = 2.5) + 
    scale_color_manual(values = c("#003366", "#990033")) +
    scale_fill_manual(values = c("#003366", "#990033")) +
    geom_vline(xintercept = 0.5, color = "#990033") + 
    geom_vline(xintercept = -0.5, color = "#003366") +
    geom_hline(yintercept = -log10(0.05), color = "grey50") +



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
      labs(title = "Resp. Disease vs. NOR", x = "log2(FoldChange)", y = "-log10(FDR)") 
    }else{
      labs(title = "Pan vs. NOR", x = "log2(FoldChange)", y = "-log10(FDR)") 
    }
    if (i == 'RS') {
      ggsave(plot = g,filename = '../plot/Resp vs Nor vol.pdf',width = 5.5,height = 5.5)
    }else{
      ggsave(plot = g,filename = '../plot/Pan vs Nor vol.pdf',width = 5.5,height = 5.5)
    }
  
}


#### stack top 6 ###
CA_cohort_feature = readRDS('../process_file/1.generate the exp matrix/CA_cohort_0801.rds')
RS_cohort_feature = read_rds('../process_file/1.generate the exp matrix/RS_cohort_0801.rds') %>%
  select(classification_name,contains('PN'))
rownames(RS_cohort_feature) = RS_cohort_feature[,1];RS_cohort_feature = RS_cohort_feature[,-1]
df.all = CA_cohort_feature
rownames(df.all) = df.all[,1];df.all = df.all[,-1]
df.all = bind_cols(df.all,RS_cohort_feature)

df_all.ordered <- df.all[order(rowSums(df.all), decreasing = T),]
df_all.ordered = df_all.ordered[1:7,]

df.all_top <- df_all.ordered %>%
  rownames_to_column("genus") %>%
  pivot_longer(cols = -genus, names_to = "variable", values_to = "value")
df.all_top$disease = lapply(df.all_top$variable,function(x) strsplit(x,'_')[[1]][[1]])
df.all_top$disease = gsub(' ','',df.all_top$disease)

df.all_top.sum <- df.all_top %>% 
  dplyr::select(1,4,3) %>% 
  group_by(genus,disease) %>% 
  summarise(sum=sum(value))

df.all_top.sum$disease %<>% factor(.)
df.all_top.sum$genus = factor(df.all_top.sum$genus,level = unique(rownames(df_all.ordered)))
df.all_top.sum$disease = factor(df.all_top.sum$disease,levels = c('BRC','CRC','GC','HCC','MEN','LC','NOR','PN'))
df.all_top.sum %>% ggplot(data = ., aes(x = disease, y = sum, fill = genus)) + 
  geom_col(position = "fill", width = 0.85) +
  scale_fill_manual(values = c("#9B3A4D", "#E2AE79", "#D0DCAA", "#F0EEBB", "#8CBDA7", "#566CA5", "#70A0AC")) +
  theme_classic() +
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.text = element_text(color = 'black',size = 10))+
  labs(y = 'Ratio %')
ggsave(filename = "../plot/stack_plot_top6.pdf", width = 7.5, height = 4.5)
getwd()



#### beta diversity ###
RS_RF = fread('../process_file/2.Multi training/RS_cohort_0715/RS_selected_RF.csv') %>% as.data.frame()
CA_RF = fread('../process_file/2.Multi training/CA_cohort_0715/CA_selected_RF.csv') %>% as.data.frame()
rownames(RS_RF) = RS_RF[,1];RS_RF = RS_RF[,-1]
rownames(CA_RF) = CA_RF[,1];CA_RF = CA_RF[,-1]
library(vegan)
library(ggpubr)
library(dplyr)
library(ggplot2)

beta_list = list(CA = CA_RF, RS = RS_RF)

for (i in names(beta_list)) {
  mat <- beta_list[[i]]

  # delete the NA and retain the numeric value
  mat <- mat[, sapply(mat, is.numeric), drop = FALSE]
  mat <- na.omit(mat)

  # Normalized to relative richness
  mat <- decostand(mat, method = "total")

  # --- Bray-Curtis ---
  d <- vegdist(mat, method = "bray")

  # --- sample group ---
  sd <- data.frame(sample = rownames(mat))
  sd$group <- vapply(sd$sample, function(x) strsplit(x, '_')[[1]][1], character(1))
  sd$group <- factor(sd$group, levels = sd$group[!duplicated(sd$group)])

  # --- PCoA（cmdscale） ---
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

  # ===================== Statistical test =====================
  # PERMANOVA
  set.seed(42)
  ad <- adonis2(d ~ group, data = sd, permutations = 999, by = "margin")
  R2 <- ad$R2[1]
  p_permanova <- ad$`Pr(>F)`[1]

  # betadisper
  bd <- betadisper(d, sd$group)
  set.seed(42)
  bd_perm <- permutest(bd, permutations = 999)
  p_betadisper <- bd_perm$tab[1, "Pr(>F)"]

  # summary output
  out_prefix <- if (i == "CA") "../plot/NC" else "../plot/RS"
  write.csv(as.data.frame(ad), paste0(out_prefix, "_PERMANOVA.csv"), quote = FALSE)
  capture.output(anova(bd), file = paste0(out_prefix, "_betadisper_anova.txt"))
  capture.output(bd_perm, file = paste0(out_prefix, "_betadisper_permutest.txt"))

  # ===================== PCoA =====================
  pc12 <- pc12 %>%
    group_by(group) %>%
    mutate(x_mean = mean(pc_x), y_mean = mean(pc_y))

  p_pcoa <- ggplot(pc12, aes(x = pc_x, y = pc_y, color = group)) +
    geom_point(size = .75) +
    geom_segment(aes(x = x_mean, y = y_mean, xend = pc_x, yend = pc_y, color = group)) +
    stat_ellipse(geom = "polygon", level = 0.9, linetype = 2, size = 0.5,
                 aes(fill = group), alpha = 0.1, show.legend = TRUE) +
    coord_fixed(ratio = 1.25) +
    xlab(paste0("PCoA1 (", pc_var[1], "%)")) +
    ylab(paste0("PCoA2 (", pc_var[2], "%)")) +
    scale_fill_manual(values = c("#9B3A4D", "#E2AE79", "#D0DCAA", "#F0EEBB", "#8CBDA7", "#566CA5", "#70A0AC", "#7d3f98")) +
    scale_color_manual(values = c("#9B3A4D", "#E2AE79", "#D0DCAA", "#F0EEBB", "#8CBDA7", "#566CA5", "#70A0AC", "#7d3f98")) +
    theme_bw() +
    guides(color = "none", fill = guide_legend(override.aes = list(alpha = 1))) +
    labs(
      title = "PCoA (Bray–Curtis)",
      caption = sprintf("PERMANOVA: R² = %.3f, p = %.3g;  betadisper p = %.3g", R2, p_permanova, p_betadisper)
    ) +
    theme(
      legend.title = element_blank(),
      legend.position = "right",
      panel.background = element_blank(),
      plot.title = element_text(size = 15, color = "black", hjust = 0.5, face = "bold")
    )


  # ===================== boxplot =====================
  bd_df <- data.frame(sample = names(bd$distances),
                      group = sd$group,
                      dist_to_centroid = bd$distances)

  p_disp <- ggplot(bd_df, aes(x = group, y = dist_to_centroid, fill = group)) +
    geom_boxplot(width = .35) +
    scale_fill_manual(values = c("#9B3A4D", "#E2AE79", "#D0DCAA", "#F0EEBB", "#8CBDA7", "#566CA5", "#70A0AC", "#7d3f98")) +
    labs(x = "", y = "Distance to centroid", title = "β-diversity dispersion (Bray–Curtis)") +
    theme_bw() +
    theme(axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(caption = sprintf("betadisper (permutation test) p = %.3g", p_betadisper))

  # 保存
  if (i == "CA") {
    ggsave(plot = p_pcoa, filename = "../plot/NC_PCoA.pdf", height = 5.5, width = 5.5)
    ggsave(plot = p_disp, filename = "../plot/NC_betadisper.pdf", height = 3.5, width = 6)
  } else {
    ggsave(plot = p_pcoa, filename = "../plot/RS_PCoA.pdf", height = 5.5, width = 5.5)
    ggsave(plot = p_disp, filename = "../plot/RS_betadisper.pdf", height = 3.5, width = 6)
  }
}

#### α diversity ####
library ("phyloseq")
library("ggplot2")
library(permute)
library(lattice)
library(vegan) 
library(tidyverse)  
library(ggpubr)
RS_RF = fread('../process_file/2.Multi training/RS_cohort_0715/RS_selected_RF.csv') %>% as.data.frame()
rownames(RS_RF) = RS_RF[,1];RS_RF = RS_RF[,-1]
CA_RF = fread('../process_file/2.Multi training/CA_cohort_0715/CA_selected_RF.csv') %>% as.data.frame()
rownames(CA_RF) = CA_RF[,1];CA_RF = CA_RF[,-1]

cohort_list = list(RS = RS_RF,CA = CA_RF)
for (i in names(cohort_list)) {
  divdata=cohort_list[[i]]
  RS_RF = cohort_list[[i]]
  mode(divdata)
  divdata<-as.data.frame(divdata)
  diversity(divdata, index="shannon")
  shannon_diversity=diversity(divdata,index="shannon") 
  simpson_diversity=diversity(divdata,index="simpson")
  S<-specnumber(divdata)
  evenness<-shannon_diversity/log(S)
  diversity<-cbind(shannon_diversity,simpson_diversity,evenness)
  print(diversity)
  
  
  richness <- estimateR(RS_RF)[1,]
  print(richness)
  diversity_richness<-cbind(shannon_diversity,simpson_diversity,evenness,richness)
  print(diversity_richness)
  if (i == 'RS') {
    write.table(diversity_richness,"RS_diversity_all.txt", sep="\t")
  }else{
    write.table(diversity_richness,"CA_diversity_all.txt", sep="\t")
  }


  RS_meta_shannon = data.frame(ID = rownames(RS_RF))
  RS_meta_shannon$group = lapply(RS_meta_shannon$ID,function(x) strsplit(x,'_')[[1]][[1]]) %>% as.character()
  colnames(RS_meta_shannon) = c('ID','Group')
  
  
  generate_vs_reference <- function(factor_variable, reference) {
    lvls <- levels(as.factor(factor_variable))
    others <- setdiff(lvls, reference)
    lapply(others, function(x) c(reference, x))
  }
  comparisons <- generate_vs_reference(RS_meta_shannon$Group, "NOR")
  comparisons
  
  plotdata = diversity_richness %>% data.frame() %>% mutate(ID = rownames(diversity_richness)) %>% 
    merge(., RS_meta_shannon, by = "ID") %>% 
    column_to_rownames("ID") %>% melt() %>% 
    filter(variable == "shannon_diversity")
  if (i == 'RS') {
    plotdata$Group = factor(plotdata$Group,level = c('NOR','LC','PN'))
  }else{
    plotdata$Group = factor(plotdata$Group,level = c('NOR','LC','CRC','GC','HCC','BRC','MEN'))
  }

  
  
  p = ggplot(data = plotdata, aes(x = Group, y = value, fill = Group)) +
    geom_boxplot(width = .35) +
    scale_fill_manual(values = c("#9B3A4D", "#E2AE79", "#D0DCAA", "#F0EEBB", "#8CBDA7", "#566CA5", "#70A0AC",'#7d3f98'))+
    stat_compare_means(method="anova",label.y = max(plotdata$value)+1) + 
    stat_compare_means(label="p.signif", ref.group = "NOR") +
    labs(x = "", y = "Shannon Index", title = "α Diversity",
         fill = '') +
    theme_bw() +
    theme(axis.text = element_text(color = "black",size = 12), plot.title = element_text(hjust = 0.5))+
    guides(color = 'none',
           fill = guide_legend(override.aes = list(alpha = 1)))
  
  if (i == 'RS') {
    ggsave(plot = p,
           filename = '../plot/alpha_diversity_RS.pdf',
           height = 3.5,
           width = 6)
  }else{
    ggsave(plot = p,
           filename = '../plot/alpha_diversity_NC.pdf',
           height = 3.5,
           width = 6)
  }

  
}


# drawing feature PASS heatmap
library(scales)
feature_CA = read.csv('../process_file/2.Multi training/CA_cohort/CA_RF_top50_features.csv')
feature_RS = read.csv('../process_file/2.Multi training/RS_cohort/RS_RF_top50_features.csv')
features = c(feature_CA$X,feature_RS$X)
features = features[!duplicated(features)]

masslin2_CA = fread('../process_file/2.Multi training/CA_cohort_0715/Maaslin2/all_results.tsv')
masslin2_RS = fread('../process_file/2.Multi training/RS_cohort_0715/Maaslin2/all_results.tsv')
masslin2_RS = masslin2_RS[masslin2_RS$value == 'PN',]

masslin2 = bind_rows(masslin2_RS,masslin2_CA)
masslin2 = masslin2[masslin2$feature %in% features,]
masslin2$sig = -log(masslin2$qval)*sign(masslin2$coef)
masslin2$label = ifelse(masslin2$coef>0,'+','-')

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
       filename = '../plot/feature_heatmap.pdf',
       height = 12,
       width = 3.5)
