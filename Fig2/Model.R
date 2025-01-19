library(mlr)
library(caret)
library(mlr3)
library(mlr3verse)
library(tidyverse)
library(randomForest)
library(pROC)
data1 <- read.csv('deg_tRNA.csv',row.names = 1)
data1 <- as.data.frame(t(data1))
colnames(data1) <- make.names(colnames(data1))
library(readxl)
group <- read.csv('F:\\R_script\\20230919 乳腺癌\\expr\\group2.csv')
target <- group$group
target <- as.factor(target)
data1$target <- target

train_rf_100_times <- function(data, target_column) {
  # 准备存储AUC值的矩阵
  auc_matrix <- matrix(NA, nrow=100, ncol=2)
  colnames(auc_matrix) <- c("Train_AUC", "Test_AUC")
  
  best_model <- NULL
  best_test_auc <- 0
  best_roc_curve <- NULL
  best_tpr <- NULL
  best_fpr <- NULL
  
  for (i in 1:100) {
    set.seed(i)
    
    # 划分训练集和测试集
    trainIndex <- createDataPartition(data[[target_column]], p = .5, 
                                      list = FALSE, 
                                      times = 1)
    trainData <- data[ trainIndex,]
    testData  <- data[-trainIndex,]
    
    # 训练随机森林模型
    rf_model <- randomForest(as.formula(paste(target_column, "~ .")), data=trainData)
    
    # 计算训练集AUC
    train_pred <- predict(rf_model, trainData, type="prob")[,2]
    train_roc <- roc(trainData[[target_column]], train_pred)
    train_auc <- auc(train_roc)
    
    # 计算测试集AUC
    test_pred <- predict(rf_model, testData, type="prob")[,2]
    test_roc <- roc(testData[[target_column]], test_pred)
    test_auc <- auc(test_roc)
    
    # 存储AUC值
    auc_matrix[i, ] <- c(train_auc, test_auc)
    
    # 检查是否是最优模型
    if (test_auc > best_test_auc) {
      best_test_auc <- test_auc
      best_model <- rf_model
      best_roc_curve <- test_roc
      best_tpr <- test_roc$sensitivities
      best_fpr <- 1 - test_roc$specificities
    }
  }
  
  # 绘制最优模型的ROC曲线
  plot(best_roc_curve, main="Best Model ROC Curve", col="#1c61b6", lwd=2)
  abline(a=0, b=1, lty=2, col="gray")
  text(0.6, 0.4, paste("AUC =", round(best_test_auc, 3)), col="#1c61b6")
  
  # 返回AUC矩阵、最优模型和最优模型的TPR/FPR矩阵
  return(list(AUC_Matrix = auc_matrix, Best_Model = best_model, TPR = best_tpr, FPR = best_fpr,best_roc = best_roc_curve))
}

result <- train_rf_100_times(data1, 'target')

library(reshape2)

create_roc_value_dataframe <- function(auc_matrix) {
  roc_value <- melt(data.frame(auc_matrix))
  roc_value$group <- rep(c('Train', 'Test'), each = 100)
  roc_value$group <- factor(roc_value$group, levels = c('Train', 'Test'))
  
  return(roc_value)
}

ys_roc_value <- create_roc_value_dataframe(result$AUC_Matrix)
rs_roc_value <- create_roc_value_dataframe(result$AUC_Matrix)
t_roc_value <- create_roc_value_dataframe(result$AUC_Matrix)
m_roc_value <- create_roc_value_dataframe(result$AUC_Matrix)
mi_roc_value <- create_roc_value_dataframe(result$AUC_Matrix)
pi_roc_value <- create_roc_value_dataframe(result$AUC_Matrix)
lnc_roc_value <- create_roc_value_dataframe(result$AUC_Matrix)

library(ggplot2)
top.mar=0.2
right.mar=0.2
bottm.mar=0.2
left.mar=0.2
mytheme <- theme_classic()+
  theme(text = element_text(family = 'sans',colour = 'black',size = 11),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.6,colour = 'black'),
        axis.ticks.length = unit(1.5,units = 'mm'),
        plot.margin = unit(x=c(top.mar,right.mar,left.mar,bottm.mar),units = 'inches'))
P <- ggplot(t_roc_value, aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  geom_point(aes(color = group), position = position_jitter(width = 0.2, height = 0), alpha = 0.5,size = 3) +
  scale_fill_manual(values = c("Train" = "lightblue", "Test" = "#DC143C")) +
  scale_color_manual(values = c("Train" = "blue", "Test" = "red")) +
  labs(title = "tRNA",
       x = NULL,
       y = "AUC Value") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.25),plot.title = element_text(hjust = 0.5)
  ) +theme_bw()
P
ggsave(P,filename = 'tRNA-100ROC.pdf',width = 3,height =4)

ys_result <- result
rs_result <- result
m_result <- result
t_result <- result
mi_result <- result
pi_result <- result
lnc_result <- result

get_top10 <- function(result){
  top10 <- data.frame(result$Best_Model$importance)
  top10$id <- rownames(top10)
  library(tidyverse)
  top10 <- top10 %>% arrange(desc(MeanDecreaseGini))
  top10 <- top10[1:10,]
  return(top10)
}

lnc_top10 <- get_top10(lnc_result)
mi_top10 <- get_top10(mi_result)
m_top10 <- get_top10(m_result)
pi_top10 <- get_top10(pi_result)
rs_top10 <- get_top10(rs_result)
ys_top10 <- get_top10(ys_result)
t_top10 <- get_top10(t_result)

write.csv(lnc_top10,'ROC-top10lncRNA.csv')
write.csv(mi_top10,'ROC-top10miRNA.csv')
write.csv(m_top10,'ROC-top10mRNA.csv')
write.csv(pi_top10,'ROC-top10piRNA.csv')
write.csv(rs_top10,'ROC-top10rsRNA.csv')
write.csv(ys_top10,'ROC-top10ysRNA.csv')
write.csv(t_top10,'ROC-top10tRNA.csv')

m_count <- data1[,m_top10$id]
mi_count <- data1[,mi_top10$id]
lnc_count <- data1[,lnc_top10$id]
pi_count <- data1[,pi_top10$id]
rs_count <- data1[,rs_top10$id]
ys_count <- data1[,ys_top10$id]
t_count <- data1[,t_top10$id]

ids <- as.data.frame(group$group)
rownames(ids) <- group$X
colnames(ids) <- 'group'
library(pheatmap)
library(ggplot2)
b <- t(pi_count)
b <- log2(b+1)
p <- pheatmap(b, 
              #annotation_row=dfGene, # （可选）指定行分组文件
              annotation_col=ids, # （可选）指定列分组文件
              #annotation_row = annotation_row,
              show_colnames = F, # 是否显示列名
              show_rownames = T,
              # 是否显示行名
              fontsize=9, # 字体大小
              color = colorRampPalette(c('navy','#ffffff','#d53e4f'))(50), # 指定热图的颜色
              annotation_legend=TRUE, # 是否显示图例
              breaks = seq(-2,2,length.out =50),
              #border_color='black',  # 边框颜色 NA表示没有
              scale = 'row', # 指定归一化的方式。"row"按行归一化，"column"按列归一化，"none"不处理
              cluster_rows = T, # 是否对行聚类
              cluster_cols = F ,# 是否对列聚类
              main = 'piRNA',
              border = F,
              #clustering_distance_rows = "euclidean", # 聚类距离方法
              #clustering_method = "ward.D2" # 聚类方法
)
p

ggsave(p,filename = 'F:\\R_script\\20230919 乳腺癌\\expr\\new1\\roc_piRNA_heatmap.pdf',width = 8,height = 6)


######整合前面筛选出所有RNA前十的表达量进行训练
all_count <- as.data.frame(cbind(lnc_count,mi_count,m_count,pi_count,rs_count,ys_count,t_count))
all_count$target <- target
result <- train_rf_100_times(all_count, 'target')
all_roc_value <- create_roc_value_dataframe(result$AUC_Matrix)

P <- ggplot(all_roc_value, aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  geom_point(aes(color = group), position = position_jitter(width = 0.2, height = 0), alpha = 0.5,size = 3) +
  scale_fill_manual(values = c("Train" = "lightblue", "Test" = "#DC143C")) +
  scale_color_manual(values = c("Train" = "blue", "Test" = "red")) +
  labs(title = "allRNA",
       x = NULL,
       y = "AUC Value") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.25),plot.title = element_text(hjust = 0.5)
  ) +theme_bw()
P
ggsave(P,filename = 'allRNA-100ROC.pdf',width = 3,height =4)

all_result <- result
all_top10 <- get_top10(all_result)
write.csv(all_top10,'ROC-top10allRNA.csv')

plot(ys_result$best_roc,col = '#7B68EE',lwd = 2)
plot(rs_result$best_roc,col = '#87CEEB',add = T,lwd = 2)
plot(t_result$best_roc,col = '#006400',add = T,lwd = 2)
plot(mi_result$best_roc,col = '#FF8C00',add = T,lwd = 2)
plot(lnc_result$best_roc,col = '#B22222',add = T,lwd = 2)
plot(m_result$best_roc,col = '#FF0000',add = T,lwd = 2)
plot(pi_result$best_roc,col = '#00008B',add = T,lwd = 2)
plot(all_result$best_roc,col = '#FFFF00',add = T,lwd = 2)
legend('bottomright', 
       legend = c('mRNA(AUC=0.910)', 'lncRNA(AUC=0.868)', 'miRNA(AUC=0.900)', 'piRNA(AUC=0.872)', 'rsRNA(AUC=0.930)', 'tRNA(AUC=0.962)', 'ysRNA(AUC=0.960)', 'allRNA(AUC=0.990)'), 
       col = c('#FF0000', '#B22222', '#FF8C00', '#00008B', '#87CEEB', '#006400', '#7B68EE', '#FFFF00'),
       lty = 1,
       lwd = 2,
       cex = 1,    # 调整字体大小
       pt.cex = 0.8, # 调整符号大小
       bty = 'n')    # 去掉图例边框
auc(ys_result$best_roc)
auc(rs_result$best_roc)
auc(t_result$best_roc)
auc(m_result$best_roc)
auc(mi_result$best_roc)
auc(pi_result$best_roc)
auc(lnc_result$best_roc)
auc(all_result$best_roc)

b <- all_count[,all_top10$id]
b <- t(b)
b <- log2(b+1)
p <- pheatmap(b, 
              #annotation_row=dfGene, # （可选）指定行分组文件
              annotation_col=ids, # （可选）指定列分组文件
              #annotation_row = annotation_row,
              show_colnames = F, # 是否显示列名
              show_rownames = T,
              # 是否显示行名
              fontsize=9, # 字体大小
              color = colorRampPalette(c('navy','#ffffff','#d53e4f'))(50), # 指定热图的颜色
              annotation_legend=TRUE, # 是否显示图例
              breaks = seq(-2,2,length.out =50),
              #border_color='black',  # 边框颜色 NA表示没有
              scale = 'row', # 指定归一化的方式。"row"按行归一化，"column"按列归一化，"none"不处理
              cluster_rows = T, # 是否对行聚类
              cluster_cols = F ,# 是否对列聚类
              main = 'allRNA',
              border = F,
              #clustering_distance_rows = "euclidean", # 聚类距离方法
              #clustering_method = "ward.D2" # 聚类方法
)
p

ggsave(p,filename = 'F:\\R_script\\20230919 乳腺癌\\expr\\new1\\roc_allRNA_heatmap.pdf',width = 8,height = 6)































