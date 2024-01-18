if (!require("stringr")) {
  install.packages("stringr")
  library("stringr")
}
if (!require("dplyr")) {
  install.packages("dplyr")
  library("dplyr")
}
if (!require("ggplot2")) {
  BiocManager::install("ggplot2")
  library("ggplot2")
}
if (!require("WGCNA")) {
  install.packages("WGCNA")
  library("WGCNA")
}
if (!require("clusterProfiler")) {
  BiocManager::install("clusterProfiler")
  library("clusterProfiler")
}
if (!require("org.Rn.eg.db")) {
  BiocManager::install("org.Rn.eg.db")
  library("org.Rn.eg.db")
}

output_module_trait <- function(){
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  filename <- paste0(module, "_",traitName, "_Module_trait_relationships.pdf")
  ylab <-  paste0("Gene significance for ", traitName)
  # 打开PDF绘图设备，设置文件名
  pdf(filename, width = 15, height = 10)
  
  abs(geneModuleMembership[moduleGenes, column])
  abs(geneTraitSignificance[moduleGenes, 1])
  
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = ylab,
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module, pch = 20)
  abline(h=0.2,v=0.8,col="red",lwd=1.5)
  ### 输出结果文件
  # 关闭PDF绘图设备，保存文件
  dev.off()
}
output_geneInfo <- function(treated,treated_name){
  geneInfo0 = data.frame(substanceBXH = names(datExpr),
                         moduleColor = moduleColors,
                         geneTraitSignificance,
                         GSPvalue)
  
  modOrder = order(-abs(cor(MEs, treated, use = "p")))
  
  # Add module membership information in the chosen order
  for (mod in 1:ncol(geneModuleMembership)){
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                           MMPvalue[, modOrder[mod]]);
    
    names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                         paste("p.MM.", modNames[modOrder[mod]], sep=""))
  }
  geneInfo0 <- base::merge(geneInfo0, annot_matrix_ciri2, by.x = "substanceBXH", by.y = "chr")
  write.csv(geneInfo0, file = paste0(treated_name,"_geneInfo.csv"))
}







circ_and_hostGene <- read.csv("MA_geneInfo.csv")
expr <- read.csv("circRNAs.csv", row.names = "ID", header = TRUE, check.names = FALSE)
annot_matrix_ciri2 <- read.csv("ciri2_annotion_matrix.csv")

datExpr <- as.data.frame(t(expr))

gsg <- goodSamplesGenes(datExpr,verbose = 3)
gsg$allOK

datExpr <- datExpr[gsg$goodSamples,gsg$goodGenes]

sampleTree <- hclust(dist(datExpr),method = "average")

sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))

plot(sampleTree, main = "detect outliers", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)



h <- 800
abline(h = h,col = "red")


#去除离群值
clust <- cutreeStatic(sampleTree, cutHeight =h,minSize = 10)

keepSamples <- (clust == 1)
datExpr <- datExpr[keepSamples, ]
sampleTree <- hclust(dist(datExpr),method = "average")

plot(sampleTree, main = " sample Tree", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

datExpr[] <- lapply(datExpr, as.numeric)



group <- str_extract(rownames(datExpr), pattern = "[MNSH]{1}")
sample <- rownames(datExpr)
trait <- data.frame(sample = sample, group = group)

# 使用 mutate 进行条件判断并将 TRUE/FALSE 转为 1/0
trait <- trait %>%
  mutate(NC = as.integer(group == "N"),
         HA = as.integer(group == "H"),
         SA = as.integer(group == "S"),
         MA = as.integer(group == "M"))



datTraits <- trait[,-1]
datTraits <- as.data.frame(datTraits)
rownames(datTraits) <- trait[, 1]
datTraits <- datTraits[, -1]
collectGarbage()
## 转成颜色
traitsColors <- numbers2colors(datTraits)

plotDendroAndColors(sampleTree, traitsColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")

powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9
#Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))

text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
sft$fitIndices[,5]
## 构建网络 关联模块与表型

cor <- WGCNA::cor
power <- 8

net = blockwiseModules(datExpr, power = power,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "PAH_Rat_TOM", 
                       verbose = 3)


table(net$colors)
# open a graphics window
sizeGrWindow(12, 9)


# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");


dim(textMatrix) = dim(moduleTraitCor)

sizeGrWindow(10,6)
par(mar = c(6, 8.5, 3, 3));

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text =1.25,
               cex.lab.y = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

## 挑选相关模块基因
# NC
traitName <- "NC"
NC = as.data.frame(datTraits$NC);
names(NC) = traitName
NC
# names (colors) of the modules
modNames = substring(names(MEs), 3)
modNames

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));


names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, NC, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(NC), sep="");
names(GSPvalue) = paste("p.GS.", names(NC), sep="");

modules <- c( "greenyellow", "blue", "red", "pink", "purple")

for (module in modules) {
  ## 筛选与HA\SA\MCT相关模块基因
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  output_module_trait()
}

output_geneInfo(NC,traitName)


## 挑选相关模块基因
# HA
traitName <- "HA"
HA = as.data.frame(datTraits$HA);
names(HA) = traitName
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
MMPvalue

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, HA, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(HA), sep="");
names(GSPvalue) = paste("p.GS.", names(HA), sep="");

modules <- c( "magenta", "turquoise")

for (module in modules) {
  ## 筛选与HA\SA\MCT相关模块基因
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  output_module_trait()
}
output_geneInfo(HA,traitName)

## 挑选相关模块基因
# SA
traitName <- "SA"
SA = as.data.frame(datTraits$SA);
names(SA) = traitName
SA
# names (colors) of the modules
modNames = substring(names(MEs), 3)
modNames

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, SA, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(SA), sep="");
names(GSPvalue) = paste("p.GS.", names(SA), sep="");

modules <- c( "yellow")

for (module in modules) {
  ## 筛选与HA\SA\MCT相关模块基因
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  output_module_trait()
}
output_geneInfo(SA,traitName)

## 挑选相关模块基因
# MA
traitName <- "MA"
MA = as.data.frame(datTraits$MA);
names(MA) = traitName
# names (colors) of the modules
modNames = substring(names(MEs), 3)


geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));


names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, MA, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(MA), sep="");
names(GSPvalue) = paste("p.GS.", names(MA), sep="");

modules <- c("brown", "tan","lightcyan")

for (module in modules) {
  ## 筛选与HA\SA\MCT相关模块基因
  column = match(module, modNames);
  moduleGenes = moduleColors==modules;
  output_module_trait()
}

# Create the starting data frame

output_geneInfo(MA,traitName)


# Order modules by their significance for BLM

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
# geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.HA));
# geneInfo = geneInfo0[geneOrder, ]

circ_and_hostGene <- circ_and_hostGene[, -1]
m <- circ_and_hostGene[complete.cases(circ_and_hostGene) & complete.cases(circ_and_hostGene$host_gene_entrezgene_id), ]
selected_modules <- c("brown", "tan", "yellow", "magenta", "turquoise", "lightcyan")


m <- m %>%
  filter(moduleColor %in% selected_modules)

m <- m %>%
  mutate(group = case_when(
    moduleColor %in% c("brown", "tan","lightcyan") ~ "MA",
    moduleColor == "yellow" ~ "SA",
    moduleColor %in% c("magenta", "turquoise") ~ "HA",
    TRUE ~ NA_character_  # 添加默认条件，如果都不满足，返回NA
  ))
m$group <- factor(m$group, levels = c("HA", "SA", "MA"))

# GOenr = GOenrichmentAnalysis(circ_and_hostGene$moduleColor, entrezCodes =circ_and_hostGene$host_gene_entrezgene_id , organism = "rat", nBestP = 10);
###run go analysishttp://127.0.0.1:41421/graphics/plot_zoom_png?width=1279&height=988
OrgDb = "org.Rn.eg.db"
GO <- compareCluster(
  host_gene_entrezgene_id~group,
  data = m,
  fun = "enrichGO",
  OrgDb = OrgDb,
  ont = "ALL",  #One of "BP", "MF", and "CC"  or "ALL"
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

# KEGG
KEGG <- compareCluster(
  host_gene_entrezgene_id~group,
  data = m,
  fun = "enrichKEGG",
  organism="rno", 
  pvalueCutoff=1,
  qvalueCutoff = 1
)



### 绘制dotplot图
dotp <- dotplot(GO,
                showCategory=10,
                includeAll = FALSE, #将有overlap的结果也展示出来
                label_format=90,
                color = "pvalue")

dotp+scale_color_gradient(low="#55668E",high ="#862657")

dotp1 <- dotplot(KEGG,
                 showCategory=10,
                 includeAll = FALSE, #将有overlap的结果也展示出来
                 label_format=90,
                 color = "pvalue")


dotp1+scale_color_gradient(low="#55668E",high ="#862657")

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
TOM <- TOMsimilarityFromExpr(datExpr, power = power);
dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;

# Call the plot function
# sizeGrWindow(100000,100000)
# image <- TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

nSelect = 1000
set.seed(10);
select = sample(nGenes, size = nSelect);
# select <- nGenes
selectTOM = dissTOM[select, select];
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
sizeGrWindow(9,9)
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")


## Visualizing the network of eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes

# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, NC))
MET = orderMEs(cbind(MET, HA))
MET = orderMEs(cbind(MET, SA))
MET = orderMEs(cbind(MET, MA))

# Plot the relationships among the eigengenes and the trait
sizeGrWindow(15,15);
par(cex = 0.9)

plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)

# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)

plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)

plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)



## Exporting to Cytoscape
# Recalculate topological overlap if needed


# Select modules

modules = c("brown","yellow","turquoise");
# Select module probes
probes = names(datExpr)

inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = names(datExpr)[inModule];

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), "1.txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), "1.txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])


