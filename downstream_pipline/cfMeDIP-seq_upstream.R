## -----------------------------------------------
# bam file path
bam_name = list.files(path='./mappeddata_bwa/',pattern='*_uniqmap_dedup.bam',full.names=TRUE)

# load R package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager") 
library(GenomicAlignments)
library(GenomicRanges)
library(getopt)
library(ggplot2)
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

# 1、设置读取bam文件的参数
param = ScanBamParam(flag = scanBamFlag(isDuplicate = NA, isUnmappedQuery = NA, isSecondaryAlignment = NA),
                     mapqFilter = 13)

# 2、获取参考基因组上1-22号染色体上，所有CG位点的位置，生成一个gragnse --- CG位点数据分析时使用
chrs <- names(Hsapiens)[1:22]
cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
cpgr <- do.call(c, lapply(1:22, function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 2))))

# 3、获取Homo参考基因组各组分的位置信息
 annots = c('hg19_genes_promoters','hg19_genes_cds','hg19_genes_5UTRs','hg19_genes_3UTRs',
           'hg19_genes_exons','hg19_genes_introns','hg19_genes_intergenic','hg19_cpg_islands','hg19_cpgs')
# 构建注释
annotations = build_annotations(genome = 'hg19', annotations = annots)

# 4、计算CpG enrichment score 变量准备
BSgenome="BSgenome.Hsapiens.UCSC.hg19"
dataset = get(ls(paste("package:", BSgenome, sep = ""))) 
# 构建输出文件框架
CpG_er_score=matrix(nrow=14,ncol=length(bam_name))  # 有多少个样品，就准备多少行nrow
CpG_er_score=as.data.frame(CpG_er_score)
CpG_er_score_rownames=c('region CpG number','region C number','region G number','region total length','region relH',
                  'region GoGe','genome CpG number','genome C number','genome G number','genome total length',
                  'genome relH','genome GoGe','enrichment score relH', 'enrichment score GoGe' )
rownames(CpG_er_score)=CpG_er_score_rownames

# 5、存储处理后frags中reads数量的data.frame
number_frags = matrix(nrow = 1, ncol=length(bam_name))
number_frags = as.data.frame(number_frags)
row.names(number_frags)[1] = 'total_number_of_frags'

# 6、存储不同基因组组分上fragments的数量
number_fragments_genome_component = list()

# 7、一个样品中每条染色体上fragments的数量
number_fragments_chromosome = matrix(nrow = 22, ncol=length(bam_name))
number_fragments_chromosome = as.data.frame(number_fragments_chromosome)
row.names(number_fragments_chromosome) = c(paste('chr',1:22,sep=''))

# 8、计算样品的所有fragments在全基因组不同大小bin中分布的数量
Autosomes=paste0('chr',c(1:22))
bins_300bp=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=300, cut.last.tile.in.chrom=T)
bins_500bp=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=500, cut.last.tile.in.chrom=T)
bins_1kb=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=1000, cut.last.tile.in.chrom=T)
bins_5kb=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=5000, cut.last.tile.in.chrom=T)
bins_10kb=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=10000, cut.last.tile.in.chrom=T)
bins_2kb=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=2000, cut.last.tile.in.chrom=T)
cat('Genome bins have been created...','\n')
uniq=0
extend=0
shift=0

# 9、样品中CG end偏好性的data.frame
CG_preference = matrix(nrow = 1, ncol=length(bam_name))
CG_preference = as.data.frame(CG_preference)
row.names(CG_preference)[1] = 'ratio_CG_end_preference'
dir.create('plots')
dir.create('Rdata')
# 设置小循环起始变量
n = 0
m = 0
## ----------------------------------
# 开始批量处理样品的bam文件
## ----------------------------------
for (a in bam_name){
  n = n + 1  
  # 设置变量名
  sample_name = strsplit(basename(a), "_")[[1]][1]
  # 读入bam文件
  galp = readGAlignmentPairs(a, use.names=TRUE, param = param)
  frags = granges(keepSeqlevels(galp, paste0("chr", 1:22), pruning.mode="coarse"),
                   on.discordant.seqnames="drop")
  w.all = width(frags)
  frags = frags[which(w.all>=50 & w.all<=650)]  # 对fragments长度进行过滤
  number_frags[1,n] = length(frags)  # 将处理后每个frags的数量保存！！！
  colnames(number_frags)[n] = sample_name
  cat('the operation of read in bam file is done','\n')
  cat('\n')
  # 获取样品中每条染色体上fragments的数量
  fragment_number_chr = as.data.frame(table(seqnames(frags)))
  fragment_number_chr$ratio = fragment_number_chr$Freq / sum(fragment_number_chr$Freq)  # 每条染色体上fragments数量占总数的比例
  number_fragments_chromosome[,n] = fragment_number_chr$Freq
  colnames(number_fragments_chromosome)[n] = sample_name
  
  p=ggplot(fragment_number_chr, aes(x=Var1, y=ratio)) + # 基础画图
    geom_bar(stat="identity",position="identity") +
    labs(x='Chromosome', y='Ratio of fragments number', title=sample_name) +
    geom_text(aes(label = round(ratio,digits=4), vjust = -0.8, hjust = 0.5), show.legend = FALSE) +
    theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), # 坐标轴标题大小
          axis.text=element_text(size=14))  # 坐标轴刻度大小
  ggsave(paste0('plots/barplot_fragments_number_on_chr_',sample_name,'.png'),p,width=40, height=20, dpi=300)
  cat('barplot of fragment number on chromosome is done','\n')
  cat('\n')
  ## barplot: the number of fragments in different components  样品的fragments在基因组不同组分中的分布情况
  ## ---------------------------------------------------
  # 获取不同基因组组分上覆盖的fragments数量
  component_count = countOverlaps(annotations, frags)
  annotations$frag_count = component_count
  annotations2 = as.data.frame(annotations)
  frag_on_components = annotations2[,c('type','frag_count')]
  # 不同基因组组分分组求fragments的数量之和
  frag_on_components2 = as.data.frame(tapply(frag_on_components$frag_count, frag_on_components$type, sum))
  names(frag_on_components2)[1] = 'number_of_fragments'
  frag_on_components2$type = row.names(frag_on_components2)
  # 将不同组分的fragments数量保存到list中
  number_fragments_genome_component[[sample_name]] = frag_on_components2
  # 画图---每种基因组分上fragments的数量
  p = ggplot(frag_on_components2, aes(x=type, y=number_of_fragments)) + # 基础画图
    geom_bar(stat="identity",position="identity") +
    labs(x='Genome components', y='Number of fragments', title=sample_name) +
    theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), # 坐标轴标题大小
          axis.text=element_text(angle = 30, size=14))  # 坐标轴刻度大小
  ggsave(paste0('plots/barplot_fragments_number_on_genome_components_',sample_name,'.png'),p, width=18, height=10, dpi=300)
  # 画图---每种基因组分上fragments数量占总数的比例
  frag_on_components2$ratio = frag_on_components2$number_of_fragments / sum(frag_on_components2$number_of_fragments)
  p = ggplot(frag_on_components2, aes(x=type, y=ratio)) + # 基础画图
    geom_bar(stat="identity",position="identity") +
    labs(x='Genome components', y='Ratio of fragments number', title=sample_name) +
    geom_text(aes(label = round(ratio,digits=4), vjust = -0.8, hjust = 0.5), show.legend = FALSE) +
    theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), # 坐标轴标题大小
          axis.text=element_text(angle = 30, size=14))  # 坐标轴刻度大小
  ggsave(paste0('plots/barplot_fragments_number_on_genome_components2_',sample_name,'.png'),p,width=18, height=10, dpi=300)
  cat('barplot of fragment number on different genome component is done','\n')
  cat('\n')
  ## --------------------------------
  ## 计算CpG enrichment score
  ## --------------------------------
  dataset = get(ls(paste("package:", BSgenome, sep = ""))) 
  y=getSeq(dataset, frags)
  regions.CG=sum(as.numeric(vcountPattern("CG",y)))
  regions.C=sum(as.numeric(vcountPattern("C",y)))
  regions.G=sum(as.numeric(vcountPattern("G",y)))
  regions.length=sum(as.numeric(width(y)))
  regions.relH=as.numeric(regions.CG)/as.numeric(regions.length)*100
  regions.GoGe=(as.numeric(regions.CG)*as.numeric(regions.length))/(as.numeric(regions.C)*as.numeric(regions.G)) 
  genome.CG=28217007
  genome.C=585012752
  genome.G=585358256
  genome.length=3095677412
  genome.relH=genome.CG/genome.length*100
  genome.GoGe=(genome.CG*genome.length)/(genome.C*genome.G)
  enrichment.score.relH=regions.relH/genome.relH  
  enrichment.score.GoGe=regions.GoGe/genome.GoGe
  colnames(CpG_er_score)[n]=sample_name
  CpG_er_score[1,n]=regions.CG
  CpG_er_score[2,n]=regions.C
  CpG_er_score[3,n]=regions.G
  CpG_er_score[4,n]=regions.length
  CpG_er_score[5,n]=regions.relH
  CpG_er_score[6,n]=regions.GoGe
  CpG_er_score[7,n]=genome.CG
  CpG_er_score[8,n]=genome.C
  CpG_er_score[9,n]=genome.G
  CpG_er_score[10,n]=genome.length
  CpG_er_score[11,n]=genome.relH
  CpG_er_score[12,n]=genome.GoGe
  CpG_er_score[13,n]=enrichment.score.relH
  CpG_er_score[14,n]=enrichment.score.GoGe
  cat(regions.CG,'CpG have been found in the',sample_name, 'bam file','\n')
  cat(regions.C,'C have been found in the',sample_name,'bam file','\n')
  cat(regions.G,'G have been found in the',sample_name,'bam file','\n')
  cat(regions.length,'base pair have been found in the',sample_name,'bam file','\n')
  cat(genome.CG,'CpG in total have been found in the Hg19 genome','\n')
  cat(genome.C,'C in total have been found in the Hg19 genome','\n')
  cat(genome.G,'G in total have been found in the Hg19 genome','\n')
  cat(genome.length,'base pair in total have been found in the Hg19 genome','\n')
  cat('CpG enrichment score for',sample_name,'is',enrichment.score.relH,'\n')
  cat('\n')
  remove(galp,y)
  gc()  
  ## ------------------------------------------------
  ## barplot: fragments上含有的CG数量
  ## ------------------------------------------------
  # 统计每条fragments上出现的 CG 位点的数量
  # frags = frags[1:2000]   # 临时
  CpG_count = queryHits(findOverlaps(frags,cpgr))  # 查看fragments与CG的overlap情况，从而获得每条fragments上CG的数量
  CpG_count2 = CpG_count
  # 统计每个fragments含有CG位点的数量
  temp = matrix(nrow=length(frags),ncol=1)
  temp = as.data.frame(temp)
  temp$V1 = row.names(temp)
  colnames(temp)[1] = 'CpG_count'
  CpG_count3 = as.data.frame(table(CpG_count))
  CpG_count3$CpG_count = as.factor(CpG_count3$CpG_count)
  # CpG_count3$Freq = as.factor(CpG_count3$Freq)
  temp2 = merge(temp, CpG_count3, by='CpG_count', all=T)
  temp2 = temp2[order(as.numeric(temp2$CpG_count)),]
  temp2$Freq = as.numeric(temp2$Freq)
  temp2[is.na(temp2)] = 0
  frags$CpG_count = temp2$Freq  # 在每条fragments后面加上CG位点的数量
  save(frags,file=paste('Rdata/frags_',sample_name,'.Rdata',sep=''))
  # 计算末端CG偏好性
  CpG_end_preference = (sum(frags$CpG_count) - regions.CG)/length(frags)
  CG_preference[1,n] = CpG_end_preference
  colnames(CG_preference)[n] = sample_name
  # 获取样品中每条fragments中出现的CG位点的数量，包括0，生成data.frame
  fragment_CpG_distribution = frags$CpG_count
  fragment_CpG_distribution = as.data.frame(table(fragment_CpG_distribution))
  fragment_CpG_distribution$ratio = fragment_CpG_distribution$Freq / sum(fragment_CpG_distribution$Freq)
  # 画图
  p = ggplot(fragment_CpG_distribution, aes(x=fragment_CpG_distribution, y=ratio)) + # 基础画图
    geom_bar(stat="identity",position="identity") +
    labs(x='Number of CG', y='Ratio of fragments number', title=sample_name) +
    geom_text(aes(label = round(ratio,digits=4), vjust = -0.8, hjust = 0.5), show.legend = FALSE) +
    theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), # 坐标轴标题大小
          axis.text=element_text(size=14))  # 坐标轴刻度大小
  ggsave(paste0('plots/barplot_fragments_number_about_CG_',sample_name,'.png'),p, width=18, height=10, dpi=400)
  cat('barplot of fragment number about CG is done','\n')
  cat('\n')
  gc()  # 释放内存
   
  ## --------------------------------------------
  ## 计算样品的所有fragments在全基因组不同大小bin中分布的数量
  ## --------------------------------------------
  CpG_count3 = unique(CpG_count2)
  frags = frags[CpG_count3]  # 保留含有CpG位点的frags
  output=paste('plots/SaturationPlot_',sample_name,'.pdf',sep='')
  sr_name=paste('sr_',sample_name,sep='')
  pdf(output)
  sr=MEDIPS.saturation(file=a, BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq, window_size = 300, nit=10, nrit=1, empty_bins=TRUE, rank=FALSE, chr.select=Autosomes)
  cat('Finished calclulating Saturation for',sample_name,'...','\n')
  assign(sr_name,sr)
  MEDIPS.plotSaturation(sr)
  dev.off()
  cat('Finished drawing Plot for',sample_name,'...','\n')
  # 计算样品与不同size的bin中有多少个fragments
  counts_300bp=countOverlaps(bins_300bp,frags)
  counts_500bp=countOverlaps(bins_500bp,frags)
  counts_1kb=countOverlaps(bins_1kb,frags)
  counts_5kb=countOverlaps(bins_5kb,frags)
  counts_10kb=countOverlaps(bins_10kb,frags)
  counts_2kb=countOverlaps(bins_2kb,frags)
  # 计算RPKM值
  rpkm_300bp=round((counts_300bp*10^9)/(300*length(frags)),digits=5)
  rpkm_500bp=round((counts_500bp*10^9)/(500*length(frags)),digits=5)
  rpkm_1kb=round((counts_1kb*10^9)/(1000*length(frags)),digits=5)
  rpkm_5kb=round((counts_5kb*10^9)/(5000*length(frags)),digits=5)
  rpkm_10kb=round((counts_10kb*10^9)/(10000*length(frags)),digits=5)
  rpkm_2kb=round((counts_2kb*10^9)/(2000*length(frags)),digits=5)
  # 计算bin内部RPKM和count
  count_name=paste(sample_name,'_count',sep='')
  rpkm_name=paste(sample_name,'_rpkm',sep='')
  bins_300bp$count=counts_300bp
  bins_500bp$count=counts_500bp
  bins_1kb$count=counts_1kb
  bins_5kb$count=counts_5kb
  bins_10kb$count=counts_10kb
  bins_2kb$count=counts_2kb
  m = m + 1
  colnames(mcols(bins_300bp))[m]=count_name
  colnames(mcols(bins_500bp))[m]=count_name
  colnames(mcols(bins_1kb))[m]=count_name
  colnames(mcols(bins_5kb))[m]=count_name
  colnames(mcols(bins_10kb))[m]=count_name
  colnames(mcols(bins_2kb))[m]=count_name
  bins_300bp$rpkm=rpkm_300bp
  bins_500bp$rpkm=rpkm_500bp
  bins_1kb$rpkm=rpkm_1kb
  bins_5kb$rpkm=rpkm_5kb
  bins_10kb$rpkm=rpkm_10kb
  bins_2kb$rpkm=rpkm_2kb
  m = m + 1
  colnames(mcols(bins_300bp))[m]=rpkm_name
  colnames(mcols(bins_500bp))[m]=rpkm_name
  colnames(mcols(bins_1kb))[m]=rpkm_name
  colnames(mcols(bins_5kb))[m]=rpkm_name
  colnames(mcols(bins_10kb))[m]=rpkm_name
  colnames(mcols(bins_2kb))[m]=rpkm_name
  remove(counts_300bp,counts_500bp,counts_1kb,counts_5kb,counts_10kb,counts_2kb,rpkm_300bp,rpkm_500bp,rpkm_1kb,rpkm_5kb,rpkm_10kb,rpkm_2kb,sr)
  cat('Finished Processing for',sample_name,'...','\n')
  
  
}
bins_300bp_df = data.frame(bins_300bp)
bins_500bp_df = data.frame(bins_500bp)
bins_1kb_df = data.frame(bins_1kb)
bins_2kb_df = data.frame(bins_2kb)
bins_5kb_df = data.frame(bins_5kb)
bins_10kb_df = data.frame(bins_10kb)
bins_300bp_df = bins_300bp_df[,-5]
bins_500bp_df = bins_500bp_df[,-5]
bins_1kb_df = bins_1kb_df[,-5]
bins_2kb_df = bins_2kb_df[,-5]
bins_5kb_df = bins_5kb_df[,-5]
bins_10kb_df = bins_10kb_df[,-5]
write.table(bins_300bp_df, file='bins_300bp_df.txt')
write.table(bins_500bp_df, file='bins_500bp_df.txt')
write.table(bins_1kb_df, file='bins_1kb_df.txt')
write.table(bins_2kb_df, file='bins_2kb_df.txt')
write.table(bins_5kb_df, file='bins_5kb_df.txt')
write.table(bins_10kb_df, file='bins_10kb_df.txt')
rm(frags,cpgr,annotations,annotations2,frag_on_components,cgs,CpG_count,
   bins_10kb_df,bins_1kb_df,bins_2kb_df,bins_300bp_df,bins_500bp_df,bins_5kb_df,
   dataset,frag_on_components2,fragment_CpG_distribution,fragment_number_chr,
   p,param,temp,temp2)
gc()
save.image('total_message_samples.Rdata')
# ----------------------------------
# 绘制参考基因组的，每条染色体碱基数量占全基因组的比例
# ----------------------------------
genome = BSgenome.Hsapiens.UCSC.hg19
genome2 = as.data.frame(seqlengths(genome))
genome2$seq_name = row.names(genome2)
genome3 = genome2[1:22,]
names(genome3)[1] = 'seq_length'
genome3$seq_name2 = factor(genome3$seq_name, levels=c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8',
                                                      'chr9','chr10','chr11','chr12','chr13','chr14','chr15',
                                                      'chr16','chr17','chr18','chr19','chr20','chr21','chr22')) ## 设置柱条的顺序
genome3$ratio = genome3$seq_length / sum(genome3$seq_length)  # 每条染色体上fragments数量占总数的比例
# 绘制柱状图
p = ggplot(genome3, aes(x=seq_name2, y=ratio)) + # 基础画图
  geom_bar(stat="identity",position="identity") +
  labs(x='Chromosome', y='Ratio', title='chromosome length in whole genome') +
  theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), # 坐标轴标题大小
        axis.text=element_text(size=14))  # 坐标轴刻度大小
ggsave(file='plots/barplot_chr_lenght_in_whole_genome.png',p, width=18, height=10, dpi=300)

# 两个样品一起画图，进行比较
if(F){
  aa = fragment_number_chr
  aa$sample = c('sf1')
  bb$sample = c('sf2')
  cc = rbind(aa,bb)
  ggplot(cc, aes(x=Var1, y=Freq, fill=sample)) + # 基础画图
    geom_bar(stat="identity",position="dodge") +
    labs(x='Chromosome', y='Number of fragments') +
    theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), # 坐标轴标题大小
          axis.text=element_text(size=14))  # 坐标轴刻度大小
}
