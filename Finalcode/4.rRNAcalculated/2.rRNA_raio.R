#### calculate the rRNA ratio
library(data.table)
library(ggplot2)
library(ggsci)
library(languageserver)
setwd('/mnt/data3/yiyonghao/NC_paper_rawdata')

# 读取mapping_summarySSU.txt文件
SSU = fread('/mnt/data3/yiyonghao/NC_paper_rawdata/rRNA_region/rRNA_calculate/mapping_summary_SSU.txt')
colnames(SSU) = c('sample','mapped','unmapped')
LSU = fread('/mnt/data3/yiyonghao/NC_paper_rawdata/rRNA_region/rRNA_calculate/mapping_summary_LSU.txt')
colnames(LSU) = c('sample','mapped','unmapped')

rRNA = data.frame(
  sample = SSU$sample,
  mapped = SSU$mapped + LSU$mapped,
  unmapped = SSU$unmapped
)

df <- data.frame(
  category = c("mapped",'unmapped'),
  count    = c(sum(rRNA$mapped),sum(rRNA$unmapped))
)



p = ggplot(df, aes(x = "", y = count, fill = category)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = scales::percent(count / sum(count))),
    position = position_stack(vjust = 0.5)
  ) +
  scale_fill_brewer(palette = "Pastel1") +
  scale_fill_manual(values = c('#4298b5','#ff9933'))+
  theme_void() +  # 去掉坐标轴和背景网格
  labs(fill = "") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
p
ggsave(plot = p,
       width = 6.5,
       height = 6.5,
       filename = './rRNA_distribution.pdf')
write.table(rRNA,'/mnt/data3/yiyonghao/NC_paper_rawdata/rRNA_region/rRNA_ratio.txt',col.names = T,row.names = F,quote = F,sep = '\t')
