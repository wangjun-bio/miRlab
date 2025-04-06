# 加载必要包
library(ggplot2)
library(dplyr)
library(tidyr)

# 设置工作目录
setwd('/Users/willow/Desktop/2025/cfRNA乳腺癌/DATA-cfRNA/')

# 读取分组信息
group <- read.csv('./merged_ID.csv')
# 去除 group 列中值为空字符串或者 NA 的行
group <- group[!(is.na(group$group) | group$group == ""), ]

# 获取 group 数据框中的两个分组名称
unique_groups <- unique(group$group)

# 假设只有两个分组
group1_name <- unique_groups[1]
group2_name <- unique_groups[2]

# 获取每个分组对应的 sample_ID
group1_samples <- group$sample_ID[group$group == group1_name]
group2_samples <- group$sample_ID[group$group == group2_name]

# 获取 pos 文件夹下所有的 csv 文件
file_list <- list.files(path = "./pos", pattern = "\\.csv$", full.names = TRUE)

# 自定义颜色
group2_color <- "#c47f7e"  # 可根据需求修改
group1_color <- "#7b8fad"  # 可根据需求修改

# 定义绘图函数
plot_frequency_M <- function(data_group1, data_group2, group1_name, group2_name, file_name) {
  
  # 转换为长格式 - 分组 1
  data_long1 <- data_group1 %>%
    pivot_longer(cols = -Length, names_to = "Sample", values_to = "Frequency") %>%
    mutate(Group = group1_name)
  
  # 转换为长格式 - 分组 2
  data_long2 <- data_group2 %>%
    pivot_longer(cols = -Length, names_to = "Sample", values_to = "Frequency") %>%
    mutate(Group = group2_name)
  
  # 合并两个分组的长格式数据
  combined_data <- rbind(data_long1, data_long2)
  
  # 计算每个 Length 和 Group 的统计量
  combined_stats <- combined_data %>%
    group_by(Length, Group) %>%
    summarise(
      Mean = mean(Frequency, na.rm = TRUE),
      Lower = quantile(Frequency, 0.05, na.rm = TRUE),
      Upper = quantile(Frequency, 0.95, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # 提取文件名中的 tRNA 亚型信息
  tRNA_subtype <- gsub(".*[/](.*)\\.csv", "\\1", file_name)
  
  title <- paste0(tRNA_subtype, ' Frequency Distribution with Confidence Band')
  
  # 绘图
  ggplot(combined_stats, aes(x = Length, color = Group, fill = Group)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, colour = NA) +
    geom_line(aes(y = Mean), linewidth = 0.8) +
    scale_x_continuous(limits = c(1, max(combined_stats$Length)), 
                       breaks = seq(1, max(combined_stats$Length), length.out = 8) %>% round()) +
    scale_color_manual(values = c(group1_color, group2_color)) +
    scale_fill_manual(values = c(group1_color, group2_color)) +
    labs(x = "DNA Length", 
         y = "Frequency",
         title = title) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),  # 去掉主要网格线
          panel.grid.minor = element_blank()   # 去掉次要网格线
          )
}

# 遍历每个 csv 文件
for (file in file_list) {
  data <- read.csv(file, check.names = FALSE)
  data1 <- melt(data)
  
  filename <- gsub(".*[/](.*)\\.csv", "\\1", file)
  filename <- paste0("./Figure/fig1/Table/", filename, '.csv')
  setwd('/Users/willow/Desktop/2025/cfRNA乳腺癌/')
  # 将合并后的数据保存为 CSV 文件
  write.csv(data1, file = filename, row.names = FALSE)
  if (grepl("RNY", file)) {
    # 找到除了 "Length" 之外的列名
    cols_to_modify <- setdiff(names(data), "Length")
    
    # 对这些列名在倒数第八个字符前添加下划线
    new_colnames <- sapply(cols_to_modify, function(col) {
      if (nchar(col) >= 8) {
        paste0(substr(col, 1, nchar(col) - 8), "_", substr(col, nchar(col) - 7, nchar(col)))
      } else {
        col
      }
    })
    # 更新列名
    names(data)[names(data) %in% cols_to_modify] <- new_colnames
  }
  # 从 data 数据框中提取对应分组的列（处理可能缺失的样本列）
  valid_group1_samples <- intersect(group1_samples, colnames(data))
  valid_group2_samples <- intersect(group2_samples, colnames(data))
  
  # 提取列并处理可能的缺失
  data_group1 <- data[, c("Length", valid_group1_samples)]
  data_group2 <- data[, c("Length", valid_group2_samples)]
  if (grepl("mt16S", file)) {
    mt16S_data_group1 <- data_group1
    mt16S_data_group2 <- data_group2
  }
  if (grepl("28S", file)) {
    R28S_data_group1 <- data_group1
    R28S_data_group2 <- data_group2
  }
  
  if (grepl("81S", file)) {
    R18S_data_group1 <- data_group1
    R18S_data_group2 <- data_group2
  }
  # 绘制图形
  p <- plot_frequency_M(data_group1, data_group2, group1_name, group2_name, file)
  print(p)
  filename <- gsub(".*[/](.*)\\.csv", "\\1", file)
  filename <- paste0(filename, '.pdf')
  setwd('/Users/willow/Desktop/2025/cfRNA乳腺癌/DATA-cfRNA/')
  ggsave(p,filename = filename, width = 10,height = 8)
  
}