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

#定义提取不同样品fragment长度的函数，生成一个保存不同样本fragment长度的文档
cfDNA_length = function(x){
  sample_name=as.character(x)
  frag_name=paste0('Rdata/','frags_',x,'.Rdata')
  load(frag_name)
  length=width(frags)
  df_length=as.data.frame(length)
  df_length$sample=sample_name
  return(df_length)
}

#输入要画图的文件，第一列为样品名称，第二列为样品的IP的frags名称，第三列为样品的Input的frags名称
info = args[1]
info = as.character(info)
group = read.csv(info)


#对每一个行的IP与Input样品画图
for (i in seq(1,dim(group)[1])){
  plot_group = c(as.character(group[i,2]),as.character(group[i,3]))
  results = lapply(plot_group,cfDNA_length)
  results = reduce(results,rbind)
  sample_name = as.character(group[i,1])
  cat (sample_name,'is processing...','\n')
  dir.create(paste0(sample_name,'_plot'))
  p = ggplot(results,aes(x=length,group=factor(sample),color=factor(sample))) + 
  scale_color_manual(values = c('#3399FF','#990000')) +                                     #定义颜色，浏览器中的网页保存了不同颜色的16进制的颜色对照表，可以在火狐浏览器中查看
  geom_density(size=1,linetype='solid') + 
  scale_x_continuous(limits = c(50,650),breaks = c(50,100,167,200,300,400,500,600)) +       #定义X轴的长度，注意X轴的长度会影响最终density图的样子，不同x轴取值影响y值
  geom_vline(xintercept=c(100,150),color = 'darkgray',size =1, linetype = 'longdash') +     #在X轴上添加两条竖线
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
  ggsave(paste0(sample_name,'_plot/',sample_name,'_50-650.pdf'),p,width=12,height=6,dpi=300)
  rm(p)
  gc()

  #画x轴50-300
  p = ggplot(results,aes(x=length,group=factor(sample),color=factor(sample))) + 
  scale_color_manual(values = c('#3399FF','#990000')) + 
  geom_density(size=1,linetype='solid') + 
  scale_x_continuous(limits = c(50,300),breaks = c(50,100,167,200,300)) + 
  geom_vline(xintercept=c(100,150),color = 'darkgray',size =1, linetype = 'longdash') +
  theme(
    ##设置图例文字大小
    legend.text = element_text(size = 10, face = "bold"),
    #legend.position = 'none',
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
  ggsave(paste0(sample_name,'_plot/',sample_name,'_50-300.pdf'),p,width=12,height=6,dpi=300)
  rm(p)
  gc()


    #画X轴100-220
  p = ggplot(results,aes(x=length,group=factor(sample),color=factor(sample))) + 
  scale_color_manual(values = c('#3399FF','#990000')) + 
  geom_density(size=1,linetype='solid') + 
  scale_x_continuous(limits = c(100,220),breaks = c(50,100,167,200)) + 
  geom_vline(xintercept=c(100,150),color = 'darkgray',size =1, linetype = 'longdash') +
  theme(
    ##设置图例文字大小
    legend.text = element_text(size = 10, face = "bold"),
    #legend.position = 'none',
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
  ggsave(paste0(sample_name,'_plot/',sample_name,'_100-220.pdf'),p,width=12,height=6,dpi=300)
}
