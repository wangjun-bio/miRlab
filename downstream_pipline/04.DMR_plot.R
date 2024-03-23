args=commandArgs(T)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
library(GenomicAlignments)
library(ggplot2)
library(plyr)
library(GenomicRanges)
library(getopt)
library(Rsamtools)
library(devtools)
library(Homo.sapiens)
library(biovizBase)
library(BSgenome.Hsapiens.UCSC.hg19)
if(!require("annotatr")) BiocManager::install("annotatr",ask = F,update = F)
if(!require("rtracklayer")) BiocManager::install("rtracklayer",ask = F,update = F)
if(!require("SummarizedExperiment")) BiocManager::install("SummarizedExperiment",ask = F,update = F)
library(stringr)
library(BSgenome)
library(MEDIPS)
library(Biostrings)
library(tidyverse)
library(parallel)
library(data.table)
if(!require("Rmisc")) BiocManager::install("Rmisc",ask = F,update = F)
library(Rmisc)
library(parallel)
library('DESeq2')
library(pheatmap)
library(cowplot)

padj_val = args[2]
padj_val = as.numeric(padj_val)
cat ('padj_value is ', padj_val,'\n')
log2fold_val = args[3]
log2fold_val = as.numeric(log2fold_val)
cat ('log2foldChange is ', log2fold_val,'\n')

Rdata = args[1] #02.DMR_generagtion.R运行完成后产生的文件
Rdata = as.character(Rdata)
save_plot_name_combine = paste0(strsplit(Rdata,'_')[[1]][1],'_',strsplit(Rdata,'_')[[1]][2],'_',Sys.Date(),'_DMRs_ratio_combine','_padj_',padj_val,'_log2FC_',log2fold_val,'_median.pdf')
cat ('Rdata file is ',Rdata,'\n')
load(Rdata)
remove(Rdata)
gc()

fragment_to_bin_info_total = function(x){
  sample_name = x
  load(paste0('Rdata/frags_',sample_name,'_Input.Rdata'))
  length = width(frags)
  df_length = as.data.frame(length)
  df_length$sample=paste0(sample_name,'_Input_total')
  Input_100_150=length(which(df_length$length<151 & df_length$length>99))
  Input_151_220=length(which(df_length$length<221 & df_length$length>150))
  Input_ratio = Input_100_150/Input_151_220
  cat(paste0(sample_name, '_Input frag ratio is ',Input_ratio,'...','\n'))
  remove(length,df_length,Input_151_220,Input_100_150)
  gc()
  load(paste0('Rdata/frags_',sample_name,'_Ip.Rdata'))
  cat(paste0(sample_name, '_Ip is processing at ',date(),'...'), sep = "\n")
  bin_fragment=data.frame(frag_size=width(frags[queryHits(findOverlaps(frags,detected_bin_total))]),
                          bin_ID=detected_bin_total[subjectHits(findOverlaps(frags,detected_bin_total))]$bin_ID,
                          methy_info=detected_bin_total[subjectHits(findOverlaps(frags,detected_bin_total))]$Methy_info,
                          chr_info = detected_bin_total[subjectHits(findOverlaps(frags,detected_bin_total))]$Chr)
  total_frags=as.numeric(length(frags))
  count_name=paste(sample_name,'_counts',sep='')
  size_ratio_name=paste(sample_name,'_size_ratio',sep='')
  size_ratio_name_corrected = paste(sample_name,'_size_ratio_corrected',sep='')
  bin_info=bin_fragment %>% group_by(bin_ID) %>% summarise(count = n(),
                                                           size_ratio = length(frag_size[which(frag_size <= 150 & frag_size >= 100)])/length(frag_size[which(frag_size <= 220 & frag_size > 150)]),
                                                           size_ratio_corrected = (length(frag_size[which(frag_size <= 150 & frag_size >= 100)])/length(frag_size[which(frag_size <= 220 & frag_size > 150)]))/Input_ratio,
                                                           methy_stat = methy_info[1],
                                                           seqnames = chr_info[1]
                                                           )
  colnames(bin_info)=c('bin_ID',count_name,size_ratio_name,size_ratio_name_corrected,'methy_stat','seqnames')
  return(bin_info)
}

for (i in rev(seq(0.5,log2fold_val,by=0.1))){
  idx_hyper_sig0.05_lf1 = which(res$padj<padj_val & res$log2FoldChange>i & res$baseMean>10)
  idx_hypo_sig0.05_lf1 = which(res$padj<padj_val & res$log2FoldChange<(-i) & res$baseMean>10)
  if (bin_size == '300bp') bins_IP_hyper=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=300, cut.last.tile.in.chrom=T)[idx_hyper_sig0.05_lf1]
  if (bin_size == '500bp') bins_IP_hyper=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=500, cut.last.tile.in.chrom=T)[idx_hyper_sig0.05_lf1]
  if (bin_size == '1kb') bins_IP_hyper=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=1000, cut.last.tile.in.chrom=T)[idx_hyper_sig0.05_lf1]
  if (bin_size == '2kb') bins_IP_hyper=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=2000, cut.last.tile.in.chrom=T)[idx_hyper_sig0.05_lf1]
  if (bin_size == '5kb') bins_IP_hyper=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=5000, cut.last.tile.in.chrom=T)[idx_hyper_sig0.05_lf1]
  if (bin_size == '10kb') bins_IP_hyper=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=10000, cut.last.tile.in.chrom=T)[idx_hyper_sig0.05_lf1]
  
  if (bin_size == '300bp') bins_IP_hypo=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=300, cut.last.tile.in.chrom=T)[idx_hypo_sig0.05_lf1]
  if (bin_size == '500bp') bins_IP_hypo=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=500, cut.last.tile.in.chrom=T)[idx_hypo_sig0.05_lf1]
  if (bin_size == '1kb') bins_IP_hypo=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=1000, cut.last.tile.in.chrom=T)[idx_hypo_sig0.05_lf1]
  if (bin_size == '2kb') bins_IP_hypo=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=2000, cut.last.tile.in.chrom=T)[idx_hypo_sig0.05_lf1]
  if (bin_size == '5kb') bins_IP_hypo=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=5000, cut.last.tile.in.chrom=T)[idx_hypo_sig0.05_lf1]
  if (bin_size == '10kb') bins_IP_hypo=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=10000, cut.last.tile.in.chrom=T)[idx_hypo_sig0.05_lf1]
  bins_IP_hyper$Methy_info = 'Hyper'
  bins_IP_hypo$Methy_info = 'Hypo'
  bins_IP_hyper$bin_ID = idx_hyper_sig0.05_lf1
  bins_IP_hypo$bin_ID = idx_hypo_sig0.05_lf1
  bins_IP_total = c(bins_IP_hyper,bins_IP_hypo)
  detected_bin_total = bins_IP_total
  detected_bin_total$Chr=as.character(seqnames(detected_bin_total))
  cl=makeCluster(60)
  clusterExport(cl,'detected_bin_total')
  clusterExport(cl,'bin_size')
  clusterExport(cl,'fragment_to_bin_info_total')
  clusterEvalQ(cl,library(GenomicAlignments))
  clusterEvalQ(cl,library(plyr))
  clusterEvalQ(cl,library(GenomicRanges))
  clusterEvalQ(cl,library(stringr))
  clusterEvalQ(cl,library(Biostrings))
  clusterEvalQ(cl,library(tidyverse))
  clusterEvalQ(cl,library(parallel))

################################################################################################################################################################################################
##画图
  Results_expr=parLapply(cl,sample_expr,function(x) fragment_to_bin_info_total(x))
  Total_bin_info_expr=Reduce(function(x,y) merge(x = x, y = y, by = "bin_ID",all=TRUE), Results_expr)
  
  Results_control=parLapply(cl,sample_control,function(x) fragment_to_bin_info_total(x))
  Total_bin_info_control=Reduce(function(x,y) merge(x = x, y = y, by = "bin_ID",all=TRUE), Results_control)
  
  stopCluster(cl)
  
  idx_expr = which(rowSums(Total_bin_info_expr[,(5*seq(1,(dim(Total_bin_info_expr)[2]-1)/5))-3] >= 20) == length(sample_expr))
  idx_control = which(rowSums(Total_bin_info_control[,(5*seq(1,(dim(Total_bin_info_control)[2]-1)/5))-3] >= 20) == length(sample_control))
  idx_counts = intersect(idx_expr,idx_control)
  
  
  idx_below_five_expr = which(rowSums(Total_bin_info_expr[,(5*seq(1,(dim(Total_bin_info_expr)[2]-1)/5))-1] > 10) != 0)
  idx_below_five_control = which(rowSums(Total_bin_info_control[,(5*seq(1,(dim(Total_bin_info_control)[2]-1)/5))-1] > 10) != 0)
  
  
  idx_below_five = union(idx_below_five_expr,idx_below_five_control)
  
  idx = setdiff(idx_counts,idx_below_five)
  
  expr_ratio = Total_bin_info_expr[,(5*seq(1,(dim(Total_bin_info_expr)[2]-1)/5))-1]

  for (y in 1:nrow(expr_ratio)){
    expr_ratio$median_expr[y] = median(as.numeric(expr_ratio[y,]))
  }

  control_ratio = Total_bin_info_control[,(5*seq(1,(dim(Total_bin_info_control)[2]-1)/5))-1]

  for (z in 1:nrow(control_ratio)){
    control_ratio$median_control[z] = median(as.numeric(control_ratio[z,]))
  }

  rownames(expr_ratio) = Total_bin_info_expr$bin_ID
  rownames(control_ratio) = Total_bin_info_control$bin_ID
  
  control_ratio$methy_stat = Total_bin_info_control$methy_stat.x
  control_ratio$seqnames = Total_bin_info_control$seqnames.x 
  
  expr_ratio = expr_ratio[idx,]
  control_ratio = control_ratio[idx,]
  df1 = as.data.frame(cbind(expr_ratio$median_expr,control_ratio$median_control))
  rownames(df1) = rownames(control_ratio)
  colnames(df1) = c('Breast','Healthy')
  df1$methy_stat = control_ratio$methy_stat
  df1$seqnames = control_ratio$seqnames
  bin_tmp = rep(rownames(df1),2)
  df1 = melt(df1)
  colnames(df1)= c('methy_stat', 'Chr','group', 'size_ratio_median')
  df1$bin_id=as.character(bin_tmp)

  p_name = paste0('p_',i)
  p=ggplot(df1,aes(x=bin_id,y=size_ratio_median,group=group, colour=group))+
    scale_color_manual(values = c('purple','black')) +
    geom_line(size=0.1,linetype='solid')+
    facet_wrap(~methy_stat)+
    scale_y_continuous(breaks=seq(-10, 10, 2))+
    annotate('text',x = df1$bin_id[500],y=10,label = paste0('padj=',padj_val,',','\t','log2foldchange=',i),size = 4)+
    theme(
      legend.text = element_text(size = 10, face = "bold"),
      axis.line.x=element_line(size = 0.5,color="black"),
      axis.line.y=element_line(size = 0.5,color="black"),###显示x,y轴
      axis.ticks.x=element_blank(), ###显示x轴刻度线
      axis.ticks.y=element_line(colour = 'black',size = 0.2), ###显示y轴刻度线
      axis.ticks.length=unit(0.5,"lines"),##设置X轴上的刻度上的标尺
      ###刻度标签设置，以及坐标轴titile
      axis.text.x=element_blank(),
      axis.text.y=element_text(size = 20, face= 'plain', color = 'black',vjust = 0.5, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      ##取消边框背景设置
      panel.background = element_rect(fill = "transparent",colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.background = element_rect(fill = "transparent",colour = NA),
      text=element_text(size=20)
      )
  assign (p_name,p)
}
  
if (length((rev(seq(0.5,log2fold_val,by=0.1)))) == 1) p = plot_grid(p_0.5,ncol = 1,align = 'v')
if (length((rev(seq(0.5,log2fold_val,by=0.1)))) == 2) p = plot_grid(p_0.6,p_0.5,ncol = 1,align = 'v')
if (length((rev(seq(0.5,log2fold_val,by=0.1)))) == 3) p = plot_grid(p_0.7,p_0.6,p_0.5,ncol = 1,align = 'v')
if (length((rev(seq(0.5,log2fold_val,by=0.1)))) == 4) p = plot_grid(p_0.8,p_0.7,p_0.6,p_0.5,ncol = 1,align = 'v')
if (length((rev(seq(0.5,log2fold_val,by=0.1)))) == 5) p = plot_grid(p_0.9,p_0.8,p_0.7,p_0.6,p_0.5,ncol = 1,align = 'v')
if (length((rev(seq(0.5,log2fold_val,by=0.1)))) == 6) p = plot_grid(p_1,p_0.9,p_0.8,p_0.7,p_0.6,p_0.5,ncol = 1,align = 'v')
if (length((rev(seq(0.5,log2fold_val,by=0.1)))) == 7) p = plot_grid(p_1.1,p_1,p_0.9,p_0.8,p_0.7,p_0.6,p_0.5,ncol = 1,align = 'v')
if (length((rev(seq(0.5,log2fold_val,by=0.1)))) == 8) p = plot_grid(p_1.2,p_1.1,p_1,p_0.9,p_0.8,p_0.7,p_0.6,p_0.5,ncol = 1,align = 'v')
if (length((rev(seq(0.5,log2fold_val,by=0.1)))) == 9) p = plot_grid(p_1.3,p_1.2,p_1.1,p_1,p_0.9,p_0.8,p_0.7,p_0.6,p_0.5,ncol = 1,align = 'v')
if (length((rev(seq(0.5,log2fold_val,by=0.1)))) == 10) p = plot_grid(p_1.4,p_1.3,p_1.2,p_1.1,p_1,p_0.9,p_0.8,p_0.7,p_0.6,p_0.5,ncol = 1,align = 'v')
if (length((rev(seq(0.5,log2fold_val,by=0.1)))) == 11) p = plot_grid(p_1.5,p_1.4,p_1.3,p_1.2,p_1.1,p_1,p_0.9,p_0.8,p_0.7,p_0.6,p_0.5,ncol = 1,align = 'v')
if (length((rev(seq(0.5,log2fold_val,by=0.1)))) == 12) p = plot_grid(p_1.6,p_1.5,p_1.4,p_1.3,p_1.2,p_1.1,p_1,p_0.9,p_0.8,p_0.7,p_0.6,p_0.5,ncol = 1,align = 'v')
if (length((rev(seq(0.5,log2fold_val,by=0.1)))) == 13) p = plot_grid(p_1.7,p_1.6,p_1.5,p_1.4,p_1.3,p_1.2,p_1.1,p_1,p_0.9,p_0.8,p_0.7,p_0.6,p_0.5,ncol = 1,align = 'v')
if (length((rev(seq(0.5,log2fold_val,by=0.1)))) == 14) p = plot_grid(p_1.8,p_1.7,p_1.6,p_1.5,p_1.4,p_1.3,p_1.2,p_1.1,p_1,p_0.9,p_0.8,p_0.7,p_0.6,p_0.5,ncol = 1,align = 'v')
if (length((rev(seq(0.5,log2fold_val,by=0.1)))) == 15) p = plot_grid(p_1.9,p_1.8,p_1.7,p_1.6,p_1.5,p_1.4,p_1.3,p_1.2,p_1.1,p_1,p_0.9,p_0.8,p_0.7,p_0.6,p_0.5,ncol = 1,align = 'v')
if (length((rev(seq(0.5,log2fold_val,by=0.1)))) == 16) p = plot_grid(p_2,p_1.9,p_1.8,p_1.7,p_1.6,p_1.5,p_1.4,p_1.3,p_1.2,p_1.1,p_1,p_0.9,p_0.8,p_0.7,p_0.6,p_0.5,ncol = 1,align = 'v')
if (length((rev(seq(0.5,log2fold_val,by=0.1)))) == 17) p = plot_grid(P_2.1,p_2,p_1.9,p_1.8,p_1.7,p_1.6,p_1.5,p_1.4,p_1.3,p_1.2,p_1.1,p_1,p_0.9,p_0.8,p_0.7,p_0.6,p_0.5,,ncol = 1,align = 'v')
ggsave(save_plot_name_combine,p,width=20,height=4*length((rev(seq(0.5,log2fold_val,by=0.1)))),dpi=300)

