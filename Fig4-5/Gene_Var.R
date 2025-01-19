
data <- read.csv('../fitermRNA.csv',row.names = 1,check.names = F)
group <- read.csv('../group_list.csv',row.names = 1)
group$id <- rownames(group)
colnames(group) <- c('Group','Tumor_Sample_Barcode')
#group$Tumor_Sample_Barcode <- gsub('-01A',' ',group$Tumor_Sample_Barcode)
#group$Tumor_Sample_Barcode <- gsub('-11A',' ',group$Tumor_Sample_Barcode)

library(tidyverse)
library(maftools)
library(data.table)
mafFilePath <- dir(path = "./SNV/",
                   pattern="masked.maf.gz$",
                   full.names = T,
                   recursive=T)
head(mafFilePath)
maf.list <- lapply(mafFilePath, fread, 
                   skip = 7, 
                   sep = "\t", 
                   header = T)
maf.merge <- do.call(rbind,maf.list)
maf.merge$Tumor_Sample_Barcode <- gsub('01A.*','01A',maf.merge$Tumor_Sample_Barcode)
maf <- read.maf(maf = maf.merge,clinicalData = group)

#Shows sample summry.
getSampleSummary(maf)
#Shows gene summary.
getGeneSummary(maf)
#Shows all fields in MAF
getFields(maf)
#shows clinical data associated with samples
getClinicalData(maf)

#提取所需要的分组
maf1 = subsetMaf(maf = maf,tsb = group$Tumor_Sample_Barcode,mafObj = F)
maf1_list = read.maf(maf1,clinicalData = group)
getClinicalData(maf1_list)
output <- file.path('D:/乳腺癌/var/')

png(paste(output,'summary.png',sep = '/'),width = 10*300,height = 10*300,res = 300)
pdf(paste(output,'summary.pdf',sep = '/'),width = 10,height = 10)
plotmafSummary(maf = maf1_list, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()
#oncoplot for top ten mutated genes.
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)


png(paste(output,'oncoplot.png',sep = '/'),width = 10*300,height = 10*300,res = 300)
pdf(paste(output,'oncoplot.pdf',sep = '/'),width = 10,height = 10)
oncoplot(maf = maf1_list,
         clinicalFeatures = c("Group"),
         top = 10,
         colors = vc_cols,
         sortByAnnotation=T,borderCol=NULL
)
dev.off()


laml.titv = titv(maf = maf1_list, plot = FALSE, useSyn = TRUE)
#plot titv summary
png(paste(output,'titv.png',sep = '/'),width = 10*300,height = 10*300,res = 300)
pdf(paste(output,'titv.pdf',sep = '/'),width = 10,height = 10)
plotTiTv(res = laml.titv)
dev.off()

png(paste(output,'PIK3CA.png',sep = '/'),width = 8*300,height = 8*300,res = 300)
pdf(paste(output,'PIK3CA.pdf',sep = '/'),width = 8,height = 8)
lollipopPlot(
  maf = maf1_list,
  gene = 'PIK3CA',
  AACol = 'HGVSp_Short',
  showMutationRate = T
)
dev.off()

plotProtein(gene = "TP53", refSeqID = "NM_001126113")

my_data = data.frame(pos = sample.int(912, 15, replace = TRUE), count = sample.int(30, 15, replace = TRUE))
head(my_data)

png(paste(output,'TP53.png',sep = '/'),width = 8*300,height = 8*300,res = 300)
pdf(paste(output,'TP53.pdf',sep = '/'),width = 8,height = 8)
lollipopPlot(
  maf = maf1_list,
  gene = 'TP53',
  AACol = 'HGVSp_Short',
  showMutationRate = T,
  refSeqID = 'NM_000546'
)
dev.off()

png(paste(output,'VAF.png',sep = '/'),width = 10*300,height = 10*300,res = 300)
pdf(paste(output,'VAF.pdf',sep = '/'),width = 10,height = 10)
plotVaf(maf = maf1_list, vafCol = NULL)
dev.off()

maf.sig = oncodrive(maf = maf1_list, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = maf.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)



##########TMB###############

get_TMB <- function(file) {
  # 需要用到的列
  use_cols <- c(
    "Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode", 
    "HGVSc", "t_depth", "t_alt_count"
  )
  # 删除这些突变类型
  mut_type <- c(
    "5'UTR", "3'UTR", "3'Flank", "5'Flank", "Intron", "IGR",
    "Silent", "RNA", "Splice_Region"
  )
  # 读取文件
  df <- file
  data <- df %>% subset(!Variant_Classification %in% mut_type) %>%
    # 计算 VAF
    mutate(vaf = t_alt_count / t_depth) %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(mut_num = n(), TMB = mut_num / 30, MaxVAF = max(vaf))
  return(data)
}
#计算
tmb_table_wt_log = tmb(maf = maf1_list)
quantile(tmb_table_wt_log$total_perMB)



paad.filtered <- get_TMB(maf1)
max(paad.filtered$TMB)
#计算tmb值
paad.filtered <- as.data.frame(paad.filtered)
rownames(paad.filtered) <- paad.filtered$Tumor_Sample_Barcode
clin <- maf1_list@clinical.data
id = intersect(clin$Tumor_Sample_Barcode,paad.filtered$Tumor_Sample_Barcode)
clin <- as.data.frame(clin)
rownames(clin) <- clin$Tumor_Sample_Barcode
clin <- clin[id,]
paad.filtered <- paad.filtered[id,]
paad.filtered <- inner_join(paad.filtered,clin,by = 'Tumor_Sample_Barcode')

library(ggplot2)
library(ggpubr)

paad.filtered <- paad.filtered %>% 
  filter(TMB != max(TMB))
max(paad.filtered$TMB)

ggboxplot(paad.filtered, x = "Group", y = "TMB",
          fill = "Group") +
  stat_compare_means(comparisons = list(
    c("High", "Low")
  )) +
  stat_compare_means(label.y = 4.5) 








