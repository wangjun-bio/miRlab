##########################################################################################################################
#柱形图
setwd("C:/Users/树金/Desktop/代码/fig-1")
getwd()
library(ggplot2)

file_list <- list.files(pattern = "_filtered-统计count.csv")
data <- list()
for (file in file_list) {
  data[[file]] <- read.csv(file = file, header = T, sep = ",")
}

N_CIRI2_new_with_median_filtered <- data[[1]]
N_DCC_new_with_median_filtered <- data[[2]]
P_CIRI2_new_with_median_filtered <- data[[3]]
P_DCC_new_with_median_filtered <- data[[4]]
T_CIRI2_new_with_median_filtered <- data[[5]]
T_DCC_new_with_median_filtered <- data[[6]]

data <- N_CIRI2_new_with_median_filtered
colnames(data)[1] <- "grop"
rownames(data) <- data$grop
data <- data[,-1]
mat <- as.matrix(data)

barplot(mat,                                        
        col = c("#CCCCCC", "#E0E0E0", "#1b98e0", "#0e5a8a"), 
        border = "NA",ylab = "Number of CircRNAs",xlab = "Samples", 
        main = "CIRI2_RAW",,cex.axis = 1,cex.names = 1,
        axisnames = T,xlim = c(0,32),ylim = c(0,2000),xpd = T,width = 0.7)

legend("topright", 
       legend = c("1 reads", "2 to 5 readds", "6 to 10 readds", ">10 readds"),
       fill = c("#CCCCCC", "#E0E0E0", "#1b98e0", "#0e5a8a"),y.intersp = 0.9) 
##########################################################################################################################
#箱型图

file_list <- list.files(pattern = "NPT.csv")
data <- lapply(file_list, read.csv)
box_CIRI2_1 <- data[[1]]
box_DCC_1 <- data[[2]]
box_CIRI2_2 <- data[[3]]
box_DCC_2 <- data[[4]]

box <- box_CIRI2_1
p <- ggplot(box, aes(x = grop, y = counts, fill = grop)) +
  geom_boxplot(aes(alpha = 0.5)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75), size = 3, shape = 21, alpha = 0.7) +
  scale_fill_manual(values = c("#d9291d", "#68af31", "#eb9c23")) +  # 设置箱线图的颜色
  labs(title = "Expression Level of 3 Types of Tissue",
       x = "Tissue Type", y = "Counts") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_text(color = "red"),
        panel.grid = element_blank()) +  # 去除背景网格线
  theme_bw() +
  coord_cartesian(xlim = c(0, 4)) +  # 设置 x 轴范围从0到4
  coord_cartesian(ylim = c(0, 2000))

p

pdf_file <- "筛选后-box-CIRI2.pdf"
ggsave(pdf_file, plot = p, device = "pdf", width = 4, height = 6, units = "in")


