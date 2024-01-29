library(ggplot2)
library(ggalt)

setwd("C:/Users/树金/Desktop/代码/fig-6")
getwd()

file_list <- list.files(pattern = ".csv")
data <- list()
for (file in file_list) {
  data[[file]] <- read.csv(file = file, header = T, sep = ",")
}

verify.GSE235850 <- data[[2]]
verify.GSE221240 <- data[[1]]

m <- verify.GSE235850[,1:3]
colnames(m)[2] <- "myself"
colnames(m)[3] <- "GEO"

K <- verify.GSE221240[,1:3]
colnames(K)[2] <- "myself"
colnames(K)[3] <- "GEO"

p1 <- ggplot(aes(x = myself, xend = GEO, y = ID), data = m) +
             geom_dumbbell(colour_x = "#FFB6C1", colour_xend = "#4169E1", size_x = 10, size_xend = 10, size = 2, color = "gray") +
             theme_light() +
             theme(panel.grid.minor.x = element_blank(),
                   text = element_text(size = 12),  # 调整主要文本的字体大小
                   axis.title.x = element_text(size = 14),  # 调整x轴标题的字体大小
                   axis.title.y = element_text(size = 14),  # 调整y轴标题的字体大小
                   axis.text.x = element_text(size = 10),   # 调整x轴标签的字体大小
                   axis.text.y = element_text(size = 10)) +    # 调整y轴标签的字体大小
             labs(title = paste("verify—GSE235850"))+
             xlab("circRNAs") +
             ylab("logFC")
p1
ggsave("GEO-verify_GSE235850.pdf", plot = p1, width = 17, height = 10)

p2 <- ggplot(aes(x = myself, xend = GEO, y = ID), data = K) +
  geom_dumbbell(colour_x = "#FFB6C1", colour_xend = "#4169E1", size_x = 10, size_xend = 10, size = 2, color = "gray") +
  theme_light() +
  theme(panel.grid.minor.x = element_blank(),
        text = element_text(size = 12),  # 调整主要文本的字体大小
        axis.title.x = element_text(size = 14),  # 调整x轴标题的字体大小
        axis.title.y = element_text(size = 14),  # 调整y轴标题的字体大小
        axis.text.x = element_text(size = 10),   # 调整x轴标签的字体大小
        axis.text.y = element_text(size = 10)) +    # 调整y轴标签的字体大小
  labs(title = paste("verify—GSE221240"))+
  xlab("circRNAs") +
  ylab("logFC")
p2
ggsave("GEO-verify_GSE221240.pdf", plot = p2, width = 17, height = 10)


