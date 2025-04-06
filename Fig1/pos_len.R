setwd('/home/jhkuang/wj/BRC_WY_20230826/plasma/result_mapping/')

# 获取所有符合命名模式的文件
file_list <- list.files(pattern = "GDM.*_RNY.*_pos.txt")
head(file_list)
# 创建一个空列表来存储每个文件的频率数据
frequency_list <- list()
# 循环处理每个文件
for (file in file_list) {
  # 读取文件
  pos_table <- read.table(file)
  
  # 提取样本名和 RNA 亚型
  sample_name <- sub("RNY.*", "", file)
  print(file)
  rna_subtype <- sub(sample_name, "", sub("_pos.txt", "", file))
  sample_name <- sub("_", "", sample_name)
  # 统计每个数值出现的次数
  counts <- table(pos_table$V1)
  
  # 计算每个数值出现的频率
  frequencies <- prop.table(counts)
  
  # 将频率数据转换为数据框
  freq_df <- data.frame(Length = as.numeric(names(frequencies)), Frequency = as.numeric(frequencies))
  
  # 为频率列添加样本名作为列名
  colnames(freq_df)[2] <- sample_name
  
  # 将频率数据存储到列表中，以 RNA 亚型作为键
  if (!exists(rna_subtype, where = frequency_list)) {
    frequency_list[[rna_subtype]] <- freq_df
  } else {
    frequency_list[[rna_subtype]] <- merge(frequency_list[[rna_subtype]], freq_df, by = "Length", all = TRUE)
  }
}
setwd('/home/yangliu/')

# 保存每个 RNA 亚型的合并频率数据到单独的文件
for (subtype in names(frequency_list)) {
  subtype <- melt(subtype)
  output_file <- paste0(subtype, "_merged_frequencies.csv")
  write.csv(frequency_list[[subtype]], file = output_file, row.names = FALSE)
}
