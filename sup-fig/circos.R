if (!require("circlize")) {
  BiocManager::install("circlize")
  library(circlize)
}


# DCC
N_bed <-read.csv(file = "dcc_N_bed.csv", row.names = 1)
P_bed <-read.csv(file = "dcc_P_bed.csv", row.names = 1)
T_bed <-read.csv(file = "dcc_T_bed.csv", row.names = 1)

# 初始化 Circos 图，设置染色体 ideogram
circos.initializeWithIdeogram(species = "hg38")

# 绘制染色体密度图
circos.genomicDensity(N_bed, col = "#c12e20", track.height = 0.2)
circos.genomicDensity(P_bed, col = "#63a036", track.height = 0.2)
circos.genomicDensity(T_bed, col = "#d9942e", track.height = 0.2) 

# ciri2
N_bed <- read.csv(file = "ciri2_N_bed.csv", row.names = 1)
P_bed <- read.csv(file = "ciri2_P_bed.csv", row.names = 1)
T_bed <- read.csv(file = "ciri2_T_bed.csv", row.names = 1)

# 初始化 Circos 图，设置染色体 ideogram
circos.initializeWithIdeogram(species = "hg38")

# 绘制染色体密度图
circos.genomicDensity(N_bed, col = "#c12e20", track.height = 0.2)
circos.genomicDensity(P_bed, col = "#63a036", track.height = 0.2)
circos.genomicDensity(T_bed, col = "#d9942e", track.height = 0.2) 





N_bed <- read.csv(file = "55_N_bed.csv", row.names = 1)
T_bed <- read.csv(file = "55_T_bed.csv", row.names = 1)

# 初始化 Circos 图，设置染色体 ideogram
circos.initializeWithIdeogram(species = "hg38")

circos.genomicDensity(T_bed, col = "#d9942e", track.height = 0.3) 
cols <-  c("red", "blue")
bed_list <- list(N_bed, T_bed)
circos.genomicTrack(bed_list, 
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicPoints(region, value, pch = 16, cex = 1, col = cols[i], ...)
                    })

