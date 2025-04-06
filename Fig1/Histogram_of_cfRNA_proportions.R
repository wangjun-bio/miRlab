
setwd('/Users/willow/Desktop/2025/cfRNA乳腺癌/DATA-cfRNA')
# 读取数据
data <- read.table("./fkrm.txt", header = FALSE, col.names = c("length", "sample_name"))

# 获取所有可能的RNA类型
all_RNA_types <- unique(c('missing', sapply(strsplit(data$sample_name[!grepl('trim$', data$sample_name)], '_'), function(x) gsub("rm", "", x[length(x)]))))

# 初始化结果数据框，包含所有可能的RNA类型列
result <- data.frame(sample_ID = character(), stringsAsFactors = FALSE)
for (type in all_RNA_types) {
  result[[type]] <- numeric()
}

previous_length <- 0

for (i in 1:nrow(data)) {
  current_length <- data$length[i]
  sample_name <- data$sample_name[i]
  
  # 提取更完整的样本ID
  sample_ID <- strsplit(sample_name, "_trim")[[1]][1]
  sample_ID <- gsub("_cutada", "", sample_ID)
  if (i > 1) {
    previous_row <- data[i - 1, ]
    previous_length <- as.numeric(previous_row[1])
    diff_length <- (previous_length - current_length) / 4

    if (!grepl('trim$', sample_name) || (length(unlist(strsplit(sample_name, 'trim'))) > 1 && grepl('trim$', sample_name))) {
      # 获取RNA类型名称
      if (!grepl('trim$', sample_name)) {
        RNA_type <- unlist(strsplit(sample_name, '_'))[length(unlist(strsplit(sample_name, '_')))]
        RNA_type <- gsub("rm", "", RNA_type)
      } else {
        RNA_type <-'missing'
      }
      
      if (!sample_ID %in% result$sample_ID) {
        new_row <- data.frame(sample_ID = sample_ID)
        for (col in colnames(result)) {
          if (col == "sample_ID") {
            next
          }
          if (col == RNA_type) {
            new_row[[col]] <- diff_length
          } else {
            new_row[[col]] <- NA
          }
        }
        result <- rbind(result, new_row)
      } else {
        result[result$sample_ID == sample_ID, RNA_type] <- diff_length
      }
    }
  }
}

# 保存为csv文件
write.csv(result, './processed_RNA_data.csv', row.names = FALSE)


# 加载必要的包
library(tidyverse)
processed_RNA_data <- read.csv('./processed_RNA_data.csv')
merged_data <- read.csv('./merged_ID.csv')
combined_data <- left_join(merged_data, processed_RNA_data, by = "sample_ID")

# 假设 combined_data 已经存在，排除 group 列为空的行
combined_data <- combined_data %>% filter(!is.na(group) & group != "")

# 找出需要进行长格式转换的列
rna_columns <- setdiff(colnames(combined_data), c("sample_ID", "group", "seq_ID", "stage"))

# 进行长格式转换
long_data <- combined_data %>%
  pivot_longer(cols = all_of(rna_columns), 
               names_to = "RNA_type", 
               values_to = "length")

# 根据 group 列进行分组
grouped_data <- split(long_data, long_data$group)

# 自定义 RNA 类型颜色，使用与前面示例一致的颜色方案
custom_colors <- c('#93c665CC','#f0f0f0CC','#8877b7CC', '#7eace0CC','#3673b1CC', '#849989DD','#d8744fAA')


# 定义绘图函数
plot_group <- function(group_df) {
  group_name <- unique(group_df$group)
  
  # 计算每个样本中各 RNA 类型的占比
  proportion_data <- group_df %>%
    group_by(sample_ID) %>%
    mutate(proportion = length / sum(length))
  
  # 绘制柱状占比图
  plot <- ggplot(proportion_data, aes(x = sample_ID, y = proportion, fill = RNA_type)) +
    geom_col(position = "stack") +  
    scale_fill_manual(values = custom_colors) +  # 使用自定义颜色
    labs(title = paste("RNA Type Proportions in Group", group_name),
         x = NULL,  # 不显示 x 轴标签
         y = "Proportion") +
    theme_minimal() + # 使用干净的主题
    theme(
      panel.grid = element_blank(), # 移除背景网格线
    )
  
  # 添加与前面示例一致的绘图风格
  plot <- plot +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(color = "#000000", linewidth = 0.05),
      strip.background = element_rect(color = "#000000", linewidth = 0.5),
      strip.text = element_text(size = 12, color = "#000000"),
      axis.text.y = element_text(color = "#000000", size = 12, face = "bold"),
      axis.title.y = element_text(color = "#000000", size = 15, face = "bold"),
      legend.background = element_rect(color = "#000000"),
      legend.text = element_text(face = "italic"),
      axis.text.x = element_blank(),  # 不显示 x 轴文本（样本名）
      axis.ticks.x = element_blank(), # 不显示 x 轴刻度线
      axis.title.x = element_blank(),  # 不显示 x 轴标题
      panel.border = element_blank(), 
      legend.position = "bottom",      # 将图例放置在图的底部
      legend.direction = "horizontal"  # 让图例呈水平排列
    )
  
  # 保存 PDF
  ggsave(filename = paste0("./figure/fig1/RNA_proportion_", group_name, ".pdf"), plot = plot, width = 8, height = 6)
  
  plot
}

# 分别计算 BEN 和 MAL 中每种 RNA 的占比
calculate_proportion_by_group <- function(data) {
  # 分别筛选 BEN 和 MAL 分组的数据
  ben_data <- data %>% filter(group == "BEN")
  mal_data <- data %>% filter(group == "MAL")
  
  # 定义一个内部函数来计算占比
  calculate_proportion <- function(group_data) {
    # 计算每种 RNA 的总长度
    total_length_by_rna <- group_data %>%
      group_by(RNA_type) %>%
      summarise(total_length = sum(length))
    
    # 计算所有 RNA 的总长度
    total_length_all <- sum(total_length_by_rna$total_length)
    
    # 计算每种 RNA 的占比
    proportion <- total_length_by_rna %>%
      mutate(proportion = total_length / total_length_all)
    
    return(proportion)
  }
  
  # 分别计算 BEN 和 MAL 的占比
  ben_proportion <- calculate_proportion(ben_data)
  mal_proportion <- calculate_proportion(mal_data)
  
  # 添加分组信息
  ben_proportion$group <- "BEN"
  mal_proportion$group <- "MAL"
  
  # 合并结果
  result <- rbind(ben_proportion, mal_proportion)
  
  return(result)
}

# 假设 grouped_data 已经定义
# 为每个分组绘制并保存柱状占比图
lapply(grouped_data, plot_group)

# 计算 BEN 和 MAL 中每种 RNA 的占比
proportions_by_group <- calculate_proportion_by_group(do.call(rbind, grouped_data))
print(proportions_by_group)
