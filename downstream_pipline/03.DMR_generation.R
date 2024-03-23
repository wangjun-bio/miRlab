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

#1、先从csv文件中读取要寻找DMR的两组样品信息；
#2、读取01.frags_generation产生的frags文件；
#3、最后在分别在Input和IP的frag文件中找到DMR。

# 获取参考基因组上1-22号染色体上，所有CG位点的位置，生成一个gragnse --- CG位点数据分析时使用
chrs <- names(Hsapiens)[1:22]
cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))
cpgr <- do.call(c, lapply(1:22, function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 2))))

bin_size = args[2]
#生成需要的bin
if (bin_size == '300bp') bins_IP=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=300, cut.last.tile.in.chrom=T)
if (bin_size == '500bp') bins_IP=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=500, cut.last.tile.in.chrom=T)
if (bin_size == '1kb') bins_IP=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=1000, cut.last.tile.in.chrom=T)
if (bin_size == '2kb') bins_IP=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=2000, cut.last.tile.in.chrom=T)
if (bin_size == '5kb') bins_IP=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=5000, cut.last.tile.in.chrom=T)
if (bin_size == '10kb') bins_IP=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=10000, cut.last.tile.in.chrom=T)

if (bin_size == '300bp') bins_Input=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=300, cut.last.tile.in.chrom=T)
if (bin_size == '500bp') bins_Input=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=500, cut.last.tile.in.chrom=T)
if (bin_size == '1kb') bins_Input=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=1000, cut.last.tile.in.chrom=T)
if (bin_size == '2kb') bins_Input=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=2000, cut.last.tile.in.chrom=T)
if (bin_size == '5kb') bins_Input=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=5000, cut.last.tile.in.chrom=T)
if (bin_size == '10kb') bins_Input=tileGenome(seqlengths(Hsapiens)[1:22],tilewidth=10000, cut.last.tile.in.chrom=T)

cat('Genome bins have been created...','\n')

#分别针对IP样品与Input样品进行读取，保存frags，并返回所有样品在不同bin中的counts数量
#Rdata/ 为
counts_table_IP = function(x){
	samples = strsplit(x,'_')[[1]]
	n=0
	for (s in samples){
		n=n+1
		if (file.exists(paste0('Rdata/','frags_',s,'_Ip.Rdata'))) {
			cat ('reading ',s,'Ip from Rdata dir...','\n')
			load(paste0('Rdata/','frags_',s,'_Ip.Rdata'))
			bins_IP$count = countOverlaps(bins_IP,frags)
			colnames(mcols(bins_IP))[n]=as.character(s)
		}
		else {
			cat ('reading ',s,'Ip from bam file...','\n')
			s_path = paste0('/home/yangming/data/cfMEDIP/fragmentation_profile/bam/',s,'_Ip_uniqmap_dedup_sorted.bam')
			param = ScanBamParam(flag = scanBamFlag(isDuplicate = NA, isUnmappedQuery = NA, isSecondaryAlignment = NA),
                     mapqFilter = 13)
			galp = readGAlignmentPairs(s_path, use.names=TRUE, param = param)
			frags = granges(keepSeqlevels(galp, paste0("chr", 1:22), pruning.mode="coarse"),
                   on.discordant.seqnames="drop")
			w.all = width(frags)
			frags = frags[which(w.all>=50 & w.all<=650)]
			CpG_count = queryHits(findOverlaps(frags,cpgr))
			frags = frags[unique(CpG_count)]
			save(frags,file=paste0('Rdata/','frags_',s,'_Ip.Rdata'))
			bins_IP$count = countOverlaps(bins_IP,frags)
			colnames(mcols(bins_IP))[n]=as.character(s)
		}			
	}
	return(bins_IP)
}

counts_table_Input = function(x){
	samples = strsplit(x,'_')[[1]]
	m=0
	for (s in samples){
		m=m+1
		if (file.exists(paste0('Rdata/','frags_',s,'_Input.Rdata'))) {
			cat ('reading ',s,'Input from Rdata dir...','\n')
			load(paste0('Rdata/','frags_',s,'_Input.Rdata'))
			bins_Input$count = countOverlaps(bins_Input,frags)
			colnames(mcols(bins_Input))[m]=as.character(s)
		}
		else {
			cat ('reading ',s,'Input from bam file...','\n')
			s_path = paste0('/home/yangming/data/cfMEDIP/fragmentation_profile/bam/',s,'_Input_uniqmap_dedup_sorted.bam')
			param = ScanBamParam(flag = scanBamFlag(isDuplicate = NA, isUnmappedQuery = NA, isSecondaryAlignment = NA),
                     mapqFilter = 13)
			galp = readGAlignmentPairs(s_path, use.names=TRUE, param = param)
			frags = granges(keepSeqlevels(galp, paste0("chr", 1:22), pruning.mode="coarse"),
                   on.discordant.seqnames="drop")
			w.all = width(frags)
			frags = frags[which(w.all>=50 & w.all<=650)]
			save(frags,file=paste0('Rdata/','frags_',s,'_Input.Rdata'))
			bins_Input$count = countOverlaps(bins_Input,frags)
			colnames(mcols(bins_Input))[m]=as.character(s)
		}
	}
	return(bins_Input)
}
#将输入样品进行分类，分为control，expr组
sample_group = args[1]
sample_group = as.character(sample_group)
output_name = strsplit(sample_group,'.csv')[[1]]
sample_group = read.csv(sample_group)
sample_control = as.character(sample_group$control)
sample_control = sample_control[nchar(sample_control)>0]	#去除掉样品名为空的字符串
sample_expr = as.character(sample_group$sample) 
sample_expr = sample_expr[nchar(sample_expr)>0] #去除掉样品名为空的字符串
samples = c(sample_control,sample_expr) 


#将所有sample一次性写入一个名字中，随后多核读取，并构建完成所有样品在不同bin中的counts数量
n=0
for (i in samples){
  n=n+1
  if (n == 1) {
    T_samples = i
  }
  else {
    T_samples = paste0(T_samples,'_',i)
  }
}
remove(sample_group,cgs)
gc()
# 读取多个样品，计算不同bin中的counts 数量
Results_IP = counts_table_IP(T_samples)
Results_Input = counts_table_Input(T_samples)
remove(T_samples)
gc()

#生成用于构建差异分析的dataframe
bins_IP = data.frame(Results_IP)
bins_Input = data.frame(Results_Input)

bins_IP = bins_IP[,-5]
bins_Input = bins_Input[,-5]

IP_control = bins_IP[,c('seqnames','start','end',sample_control)]
IP_expr = bins_IP[,c('seqnames','start','end',sample_expr)]
Total_IP = cbind (IP_control,IP_expr)
idx = 4+length(sample_control)
Total_IP = Total_IP[,-c(idx,idx+1,idx+2)]
Total_IP = Total_IP[,4:(4+length(samples)-1)]

Input_control = bins_Input[,c('seqnames','start','end',sample_control)]
Input_expr = bins_Input[,c('seqnames','start','end',sample_expr)]
Total_Input = cbind (Input_control,Input_expr)
idx = 4+length(sample_control)
Total_Input = Total_Input[,-c(idx,idx+1,idx+2)]
Total_Input = Total_Input[,4:(4+length(samples)-1)]

#进行差异分析
cat('starting DMR analysis based on DESeq2...')
annotDF=data.frame(condition=c(rep('control',length(sample_control)),rep(('expr'),length(sample_expr))))
rownames(annotDF)=colnames(Total_IP)
dds = DESeqDataSetFromMatrix(countData=Total_IP, colData=annotDF, design= ~ condition)
dds = DESeq(dds)
res = results(dds,contrast=c("condition","expr","control"))
remove(dds,annotDF)
gc()
save.image(paste0(output_name,'_DMR_',Sys.Date(),'.Rdata'))

#找出差异显著的基因组bin区域，并计算这些区域内包含cfDNA 的长度
idx_hyper_sig0.05_lf1 = which(res$pvalue<0.05 & res$log2FoldChange>1 & res$baseMean>10)
idx_hypo_sig0.05_lf1 = which(res$pvalue<0.05 & res$log2FoldChange<(-1) & res$baseMean>10)

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

#生成每个不同样品在IP，Hyper，Hypo，Input中cfDNA长度的dataframe
cfDNA_length_Ip = function(x){
	load(paste0('Rdata/','frags_',x,'_Ip.Rdata'))
	length = width(frags)
	df_length = as.data.frame(length)
	df_length$sample=paste0(x,'_Ip_total')
	return(df_length)
}
cfDNA_length_Hyper = function(x){
	load(paste0('Rdata/','frags_',x,'_Ip.Rdata'))
	idx_hyper_frags = subjectHits(findOverlaps(bins_IP_hyper,frags))
	frags_hyper = frags[idx_hyper_frags]
	length = width(frags_hyper)
	df_length = as.data.frame(length)
	df_length$sample = paste0(x,'_IP_hyper')
	return(df_length)
}

cfDNA_length_Hypo = function(x){
	load(paste0('Rdata/','frags_',x,'_Ip.Rdata'))
	idx_hypo_frags = subjectHits(findOverlaps(bins_IP_hypo,frags))
	frags_hypo = frags[idx_hypo_frags]
	length = width(frags_hypo)
	df_length = as.data.frame(length)
	df_length$sample = paste0(x,'_IP_hypo')
	return(df_length)
}
cfDNA_length_Input = function(x){
	load(paste0('Rdata/','frags_',x,'_Input.Rdata'))
	length = width(frags)
	df_length = as.data.frame(length)
	df_length$sample=paste0(x,'_Input_total')
	return(df_length)
}

#构建一个输出的表格，用于存放不同IP，Hyper，Hypo，Input中的ratio
output_ratio = matrix(nrow = length(samples),ncol = 12)
output_ratio = as.data.frame(output_ratio)
colnames(output_ratio)=c('IP_short_counts','IP_long_counts','IP_ratio',
						  'Hyper_short_counts','Hyper_long_counts','Hyper_ratio',
						  'Hypo_short_counts','Hypo_long_counts','Hypo_ratio',
						  'Input_short_counts','Input_long_counts','Input_ratio')


#逐个读取样品，计算每个样品在IP，Input，Hyper，Hypo中的cfDNA长度，100-150/151-220的比值
cat ('generaing fragments ratio in IP、Hyper region、Hypo region、Input...')
o=0
dir.create('plot')
for (s in samples){
	o=o+1
	results_Ip = cfDNA_length_Ip(s)
	results_Hyper = cfDNA_length_Hyper(s)
	results_Hypo = cfDNA_length_Hypo(s)
	results_Input = cfDNA_length_Input(s)

	Ip_100_150=length(which(results_Ip$length<151 & results_Ip$length>99))
	Ip_151_220=length(which(results_Ip$length<221 & results_Ip$length>150))

	Hyper_100_150=length(which(results_Hyper$length<151 & results_Hyper$length>99))
	Hyper_151_220=length(which(results_Hyper$length<221 & results_Hyper$length>150))

	Hypo_100_150=length(which(results_Hypo$length<151 & results_Hypo$length>99))
	Hypo_151_220=length(which(results_Hypo$length<221 & results_Hypo$length>150))

	Input_100_150=length(which(results_Input$length<151 & results_Input$length>99))
	Input_151_220=length(which(results_Input$length<221 & results_Input$length>150))
	
	rownames(output_ratio)[o] = as.character(s)
	output_ratio[o,1] = Ip_100_150
	output_ratio[o,2] = Ip_151_220
	output_ratio[o,3] = Ip_100_150/Ip_151_220
	output_ratio[o,4] = Hyper_100_150
	output_ratio[o,5] = Hyper_151_220
	output_ratio[o,6] = Hyper_100_150/Hyper_151_220
	output_ratio[o,7] = Hypo_100_150
	output_ratio[o,8] = Hypo_151_220
	output_ratio[o,9] = Hypo_100_150/Hypo_151_220
	output_ratio[o,10] = Input_100_150
	output_ratio[o,11] = Input_151_220
	output_ratio[o,12] = Input_100_150/Input_151_220
	cat ('Short(100-150)/Long(151-220) ratio in ',s,' is:','\n',
		'IP:',Ip_100_150,'/',Ip_151_220,'=',Ip_100_150/Ip_151_220,'\n',
		'Hyper:',Hyper_100_150,'/',Hyper_151_220,'=',Hyper_100_150/Hyper_151_220,'\n',
		'Hypo:',Hypo_100_150,'/',Hypo_151_220,'=',Hypo_100_150/Hypo_151_220,'\n',
		'Input:',Input_100_150,'/',Input_151_220,'=',Input_100_150/Input_151_220,'\n')
	results = rbind(results_Ip,results_Hyper,results_Hypo,results_Input)
	dir.create(paste0('plot/',s,'_plot'))
	#针对每一个样品在IP，Hyper，Hypo，Input中的 ratio进行画图
  	p = ggplot(results,aes(x=length,group=factor(sample),color=factor(sample))) +
  	scale_color_manual(values = c('#3399FF','#990000','#009933','#FF9900')) +  #定义颜色，浏览器中的网页保存了不同颜色的16进制的颜色对照表，可以在火狐浏览器中查看
  	geom_density(size=1,linetype='solid') +  #定义图片类型
  	scale_x_continuous(limits = c(50,300),breaks = c(50,100,167,200,300)) +       #定义X轴的长度，注意X轴的长度会影响最终density图的样子，不同x轴取值影响y值
  	geom_vline(xintercept=c(120,160),color = 'darkgray',size =1, linetype = 'longdash') +     #在X轴上添加两条竖线
  	theme(
    ##设置图例文字大小
    	legend.text = element_text(size = 10, face = "bold"),                                   #定义图例
    #legend.position = 'none',                                                              #可以用此段代码去除图例
    #设置图例文字旁边的方框大小
    #legend.key.size=unit(0.05,'cm'),
    	axis.line=element_line(size = 1,color="black"),###显示x,y轴
    	axis.ticks.x=element_line(colour = 'black',size = 1), ###显示x轴刻度线
    	axis.ticks.y=element_line(colour = 'black',size = 1), ###显示y轴刻度线
    	axis.ticks.length=unit(0.5,"lines"),##设置X轴上的刻度上的标尺
    ###刻度标签设置，以及坐标轴titile
    	axis.text.x=element_text(size = 20,face = "bold",color = 'black' ,vjust = 0.5, hjust = 0.5),
    	axis.text.y=element_text(size = 20,face = "bold", color = 'black',vjust = 0.5, hjust = 0.5),
    	axis.title.x = element_text(size = 20,face = "bold", color = 'black',vjust = 0.5, hjust = 0.5),
    	axis.title.y = element_text(size = 20,face = "bold", color = 'black',vjust = 0.5, hjust = 0.5),
    	##取消边框背景设置
    	panel.background = element_rect(fill = "transparent",colour = NA),
    	panel.grid.minor = element_blank(),
    	panel.grid.major = element_blank(),
    	plot.background = element_rect(fill = "transparent",colour = NA))
  	ggsave(paste0('plot/',s,'_plot/',s,'_50-300.pdf'),p,width=12,height=6,dpi=300)
}
write.table(output_ratio, file = paste0('ratio_in_IP_Hyper_Hypo_Input in',output_name,'_',Sys.Date(),'.csv'))
remove(p,results,results_Ip,results_Hyper,results_Hypo,results_Input,
	   Ip_100_150,Ip_151_220,Hyper_100_150,Hyper_151_220,
	   Hypo_100_150,Hypo_151_220,Input_100_150,Input_151_220,s,output_ratio)
gc()

cat ('generating volcano plot...')
#画Volcano图
vol_df = as.data.frame(res)
vol_df = vol_df[which(vol_df$baseMean>10),]
#计算差异显著的bin
vol_df$Group = 'not-significant'
vol_df$Group[which((vol_df$padj < 0.05) & (vol_df$log2FoldChange > 1))] = 'up-regulated'
vol_df$Group[which((vol_df$padj < 0.05) & (vol_df$log2FoldChange < (-1)))] = 'down-regulated'

p=ggplot(vol_df, aes(x = log2FoldChange, y = -log10(padj),color = factor(Group)))+
  scale_color_manual(values = c("#2f5688","#BBBBBB","#CC0000"))+    #定义颜色
  geom_point(alpha=1,size =2,pch=16)+        #定义图片类型
  scale_x_continuous(expand = c(0, 0),limits=c(-3,3),breaks=c(0,-1,-2,1,2))+         #定义x轴
  scale_y_continuous(limits=c(0,5),breaks = seq(0,5,by = 1))+                        #定义y轴
  geom_vline(xintercept = c(-1,1), size = 1,lty = 2)+             #在x轴的虚线
  geom_hline(yintercept = -log10(0.05), size = 1,lty = 2)+        #在Y轴的虚线
  theme_bw()+ #清楚背景
  theme(
    legend.position = 'none',
    #legend.text = element_text(size = 10, face = "bold"), 
    axis.line = element_blank(),
    axis.ticks.x=element_line(colour = 'black',size = 1), ###显示x轴刻度线
    axis.ticks.y=element_line(colour = 'black',size = 1), ###显示y轴刻度线
    axis.ticks.length=unit(0.3,"lines"),##设置X轴上的刻度上的标尺
    axis.text.x = element_text(size = 12,color = 'black' ,vjust = 0.5, hjust = 0.5),
    axis.text.y=element_text(size = 12, color = 'black',vjust = 0.5, hjust = 0.5),
    axis.title.x = element_text(size = 14, color = 'black',vjust = 0, hjust = 0.5),
    axis.title.y = element_text(size = 14, color = 'black',vjust = 2, hjust = 0.5),
    legend.title = element_text(face ="bold"),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background = element_rect(colour = 'black',size = 1.5),
    panel.border = element_blank())+
  
  labs(x = "Log2FoldChange",
       y = "-log10(Adjust P-value)")

ggsave(paste0('plot/',output_name,'_volcano.pdf'),p,width=8,height=8,dpi=300)
remove(p,vol_df)
gc()

cat ('generatin heatmap...')
#画heatmap
for (i in seq(length(Total_IP))){
  if (bin_size == '300bp') tmp = 300
  if (bin_size == '500bp') tmp = 500
  if (bin_size == '1kb') tmp=1000
  if (bin_size == '2kb') tmp =2000
  if (bin_size == '5kb') tmp = 5000
  if (bin_size == '10kb') tmp =10000
  Total_IP[,i] = round(((Total_IP[,i])*10^9)/(tmp*sum(Total_IP[,i])),digits=5)
}
Hyper = Total_IP[idx_hyper_sig0.05_lf1,]
Hypo = Total_IP[idx_hypo_sig0.05_lf1,]
heat_df = rbind(Hyper,Hypo)
annotation_col = data.frame(CancerType = c(c(rep('H',length(sample_control))),c(rep('C',length(sample_expr)))))
rownames(annotation_col) = colnames(Total_IP)
ann_colors = list(CancerType = c( H = '#00FF00',C = '#000099'))
pheatmap(heat_df,scale = 'row',clustering_distance_rows = "correlation",
         #color = colorRampPalette(c("blue","white","red"))(100),
         clustering_method = "average",border=FALSE,
         show_rownames=F,show_colnames=F,cluster_row = FALSE,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         cellwidth = 15,cellheight = 0.05,
         width = 10,height = 20,
         filename = paste0('plot/',output_name,'_heatmap.pdf'))
remove(heat_df,annotation_col,ann_colors,Hyper,Hypo)
gc()








