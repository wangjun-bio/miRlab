# PS: this script is used to calculate the length distribution of sequences.
# modification: 2025-08-11
##### Calculate the length distribution ######
setwd('/mnt/data3/yiyonghao/NC_paper_rawdata/')
paths = c('upload','PN','Brain')
# length
lengths = list()
for (i in paths) {
   files = list.files(paste0('/mnt/data3/yiyonghao/NC_paper_rawdata/',i,'/decompressed/length'),full.names = T)
   result = lapply(files,function(x) {
        tmpname = basename(x)
        tmpname = gsub('.tsv','',tmpname)
        tmp = fread(x)
        return(tmp)
    })
    result = bind_rows(result)
    colnames(result) = c('reads','length')
    lengths[[i]] = result
}
lengths = do.call(bind_rows,lengths)

# classification 
classes = list()
for (i in paths) {
   files = list.files(paste0('/mnt/data3/yiyonghao/NC_paper_rawdata/',i,'/krakenuniq_output'),full.names = T,pattern = '_readclassification.tsv')
   result = lapply(files,function(x) {
        tmpname = basename(x)
        tmpname = gsub('_readclassification.tsv','',tmpname)
        tmp = fread(x)
        tmp = tmp[,c(2,1)]
        return(tmp)
    })
    result = bind_rows(result)
    colnames(result) = c('reads','classification')
    classes[[i]] = result
}
classes = do.call(bind_rows,classes)
length_all = merge(classes,length,by = 'reads',all = T) %>% na.omit()

# randomly select 1M reads to exhibit the length distribution
set.seed(123)
result = length_all
result_sub = result[sample(nrow(result),1000000),]
result_sub$classification = ifelse(is.na(result_sub$classification) & result_sub$length < 36,'Dropped',result_sub$classification)
result_sub = result_sub[!is.na(result_sub$length),]

result_sub$classification = gsub('U','Unclassified',result_sub$classification)
result_sub$classification = gsub('C','Classified',result_sub$classification)
result_sub$classification = factor(result_sub$classification,levels = c('Dropped','Unclassified','Classified'))
result_sub = result_sub[!is.na(result_sub$classification),]
p = ggplot(data = result_sub,aes(x = length))+
  geom_histogram(bins = 60,binwidth = 1,aes(fill = classification))+
  geom_vline(xintercept = 36, linetype = "dashed", color = "black", size = 0.5)+
  theme_classic()+
  scale_fill_manual(values = c('lightgrey','#005670','#84bd00'))+
  labs(fill = '',
       x = 'Length(nt)',
       y = 'Reads')+
  theme(axis.title = element_text(color = 'black',size = 12),
        axis.text = element_text(color = 'black',size = 10))
ggsave(plot = p,
       width = 6.5,
       height = 3.5,
       filename = './length_distribution.pdf')




## draw classification pie plot
NC_C = fread('./upload/summary_C_U.txt')
PN_C = fread('./PN/summary_C_U.txt')
Brain_C = fread('./Brain/summary_C_U.txt')
C = bind_rows(NC_C,PN_C,Brain_C)
df <- data.frame(
  category = c("C",'UC'),
  count    = c(sum(C$C_count),sum(C$U_count))
)


p = ggplot(df, aes(x = "", y = count, fill = category)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = scales::percent(count / sum(count))),
    position = position_stack(vjust = 0.5)
  ) +
  scale_fill_brewer(palette = "Pastel1") +
  scale_fill_manual(values = c('#84bd00','#005670'))+
  theme_void() +  # 去掉坐标轴和背景网格
  labs(fill = "") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
ggsave(plot = p,
       width = 6.5,
       height = 6.5,
       filename = './classification_distribution.pdf')