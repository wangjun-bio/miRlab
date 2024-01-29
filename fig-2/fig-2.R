#韦恩图
library(VennDiagram)
setwd("C:/Users/树金/Desktop/代码/fig-2/A")
getwd()

file_list <- list.files(pattern = "_filtered.csv")
data <- list()
for (file in file_list) {
  data[[file]] <- read.csv(file = file, header = T, sep = ",")
}

N_CIRI2_new_with_median_filtered <- data[[1]]
N_DCC_new_with_median_filtered <- data[[2]]

venn_list <- list(group1 = N_CIRI2_new_with_median_filtered$ID, 
                  group2 = N_DCC_new_with_median_filtered$ID)
venn.diagram(venn_list, filename = 'venn1.png', imagetype = 'png', 
             fill = c('red', 'blue'), alpha = 0.50, cat.col = rep('black', 2), 
             col = 'black', cex = 1.5, fontfamily = 'serif', 
             cat.cex = 1.5, cat.fontfamily = 'serif')
#########################################################################################################################
file_list <- list.files(pattern = "_join.csv")
data <- list()
for (file in file_list) {
  data[[file]] <- read.csv(file = file, header = T, sep = ",")
}

N_inner_join <- data[[1]]
P_inner_join <- data[[2]]
T_inner_join <- data[[3]]

venn_list <- list(group1 = N_inner_join$ID, 
                  group2 = GRH38.circRNA$Genomic.position, 
                  group3 = P_inner_join$ID, 
                  group4 = T_inner_join$ID)
venn.diagram(venn_list, filename = 'venn2.png', imagetype = 'png', 
             fill = c('red', 'blue', 'green', 'orange'), alpha = 0.50, 
             cat.col = c('red', 'blue', 'green', 'orange'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('red', 'blue', 'green', 'orange'), cex = 1.5, fontfamily = 'serif')
#########################################################################################################################
#相关性分析
library(ggpubr)
library(ggplot2)
setwd("C:/Users/树金/Desktop/代码/fig-2/B")
getwd()

file_list <- list.files(pattern = "average.csv")
data <- list()
for (file in file_list) {
  data[[file]] <- read.csv(file = file, header = T, sep = ",")
}

average <- data[[1]]
average[is.na(average)] <- 0
colnames(average)[1] <- "CIRI2_N"
colnames(average)[2] <- "DCC_N"
colnames(average)[4] <- "CIRI2_P"
colnames(average)[5] <- "DCC_P"
colnames(average)[7] <- "CIRI2_T"
colnames(average)[8] <- "DCC_T"
m <- average
m_matrix_1 <- m[, c("CIRI2_N", "DCC_N")]
m_matrix_2 <- m[, c("CIRI2_P", "DCC_P")]
m_matrix_3 <- m[, c("CIRI2_T", "DCC_T")]

p <- ggscatter(data = m_matrix_1, x = "CIRI2_N", y = "DCC_N", 
               color = "red", fill = "lightgray",
               add = "reg.line", conf.int = TRUE, 
               add.params = list(color = "black", fill = "blue"),
               cor.coef = TRUE,
               cor.method = "pearson",
               cor.coef.coord = c(-3, 1)) + 
  ggtitle("A") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 50)) +
  annotate("text", x = 20, y = 50, 
           label = paste("Pearson Correlation:", round(cor(m$CIRI2_N, m$DCC_N), 2)))

p
ggsave("A.pdf", plot = p, width = 7, height = 7)



