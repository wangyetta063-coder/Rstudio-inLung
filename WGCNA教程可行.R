# script to perform WGCNA
#setwd("C:/Users/Administrator/Desktop")
install.packages("WGCNA")
library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
CorLevelPlot
library(gridExtra)
allowWGCNAThreads() 
#错误于 library(DESeq2): 不存在叫 ‘DESeq2’ 的程序包
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("DESeq2")
#错误于 library(GEOquery): 不存在叫 ‘GEOquery’ 的程序包
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("GEOquery")
#错误于 library(CorLevelPlot): 不存在叫 ‘CorLevelPlot’ 的程序包
install.packages("PerformanceAnalytics")
CorLevelPlot <- function(cor_mat, cex=0.8){
  library(ggplot2)
  library(reshape2)
  
  cor_m <- melt(cor_mat)
  ggplot(cor_m, aes(Var1, Var2, fill=value)) +
    geom_tile(color="white") +
    scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0) +
    geom_text(aes(label=round(value, 2)), size=3) +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=45, hjust=1))
}
#######以上是准备阶段已完成####准备开始WGCNA#######
# 1. Fetch Data ------------------------------------------------
#setwd("C:/Users/Administrator/Desktop")
data <- read.delim("C:/Users/Administrator/Desktop/GSE152418_p20047_Study1_RawCounts.txt",
                   header = TRUE,
                   check.names = FALSE)


getwd()   # 确认已经是 Desktop
# get metadata
geo_id <- "GSE152418"
gse <- getGEO(geo_id, GSEMatrix = TRUE)
phenoData <- pData(phenoData(gse[[1]]))
head(phenoData)
phenoData <- phenoData[,c(1,2,46:50)]###这个是我们对最后几列的数据重点观察

# prepare data
data[1:10,1:10]

data <- data %>% 
  gather(key = 'samples', value = 'counts', -ENSEMBLID) %>% 
  mutate(samples = gsub('\\.', '-', samples)) %>% 
  inner_join(., phenoData, by = c('samples' = 'title')) %>% 
  select(1,3,4) %>% 
  spread(key = 'geo_accession', value = 'counts') %>% 
  column_to_rownames(var = 'ENSEMBLID')

# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
data <- data[gsg$goodGenes == TRUE,]

#出图hhhhhhh   detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)

# 第三步pca - method 2

pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))


### NOTE: If there are batch effects observed, correct for them before moving ahead
# exclude outlier samples
samples.to.be.excluded <- c('GSM4615000', 'GSM4614993', 'GSM4614995')

data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]
#####第四步开始了%%%%%%%%%%%%%%%%%%%%%%%%%%

# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset

# exclude outlier samples
colData <- phenoData %>% 
  filter(!row.names(.) %in% samples.to.be.excluded)


# fixing column names in colData
names(colData)
names(colData) <- gsub(':ch1', '', names(colData))
names(colData) <- gsub('\\s', '_', names(colData))

# making the rownames and column names identical
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))


# create dds
dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~ 1) # not spcifying model



## remove all genes with counts < 15 in more than 75% of samples (31*0.75=23.25)
## suggested by WGCNA on RNAseq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >= 24,]
nrow(dds75) # 13284 genes


# perform variance stabilization
dds_norm <- vst(dds75)


# get normalized counts归一化计数
norm.counts <- assay(dds_norm) %>% 
  t()
####接下来进行质量监控并且排除异常值，然后对基因表达数据进行归一化
# 4. Network Construction构建网络，构建一个软阈值 ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)
####以上成功运行，选18作为软阈值############################
# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 18
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize###这里可以根据自己的数据进行选择
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000,
                          TOMType = "signed",
                          power = soft_power,###我们选择的软阈值
                          mergeCutHeight = 0.25,##切割高度
                          numericLabels = FALSE,##希望模型基因标签是颜色名称而不是数字
                          randomSeed = 1234,##设置随机种子提高可重复性
                          verbose = 3)##详细程度设置成3


cor <- temp_cor

# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs


# Print out a preview
head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)




# grey module = all genes that doesn't fall into other modules were assigned to the grey module

####到这里就进行第二个视频了#######
# 6A. Relate modules to traits --------------------------------------------------
# module trait associations



# create traits file - binarize categorical variables
traits <- colData %>% 
  mutate(disease_state_bin = ifelse(grepl('COVID', disease_state), 1, 0)) %>% 
  select(8)###这里确实There are only 8 columns.
# binarize categorical variables开始有点难###

colData$severity <- factor(colData$severity, levels = c("Healthy", "Convalescent", "ICU", "Moderate", "Severe"))##在呼叫数据中设置严重程度的级别

severity.out <- binarizeCategoricalColumns(colData$severity,
                                           includePairwise = FALSE,
                                           includeLevelVsAll = TRUE,##希望将一个级别与其他所有级别进行比较
                                           minCount = 1)####将最小值设为1


traits <- cbind(traits, severity.out)
# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)


module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)



# visualize module-trait association as a heatmap
########跳过#######################
heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')
###> 这里在console框输入后names(heatmap.data)，找到所有性状相关的联性数据，因此性状关联性数据从第19列到第23列
head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')




CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[18:22],####列名的向量
             y = names(heatmap.data)[1:17],####所有特征基因的名字
             col = c("blue1", "skyblue", "white", "pink", "red"))



module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'turquoise') %>% 
  rownames()

###############################分割线#########################

#✅ 📌 这是“完整版” CorLevelPlot（支持 x/y/col）
CorLevelPlot <- function(data, x, y, col = c("blue1", "skyblue", "white", "pink", "red")) {
  library(ggplot2)
  library(reshape2)
  
  # 选择行列
  sub_data_x <- data[, x, drop = FALSE]
  sub_data_y <- data[, y, drop = FALSE]
  
  # 相关矩阵
  cor_mat <- cor(sub_data_y, sub_data_x, use = "pairwise.complete.obs")
  
  # 转换为长格式
  cor_m <- melt(cor_mat)
  colnames(cor_m) <- c("Y", "X", "value")
  
  # 绘图
  ggplot(cor_m, aes(X, Y, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = col[1], mid = col[3], high = col[5], midpoint = 0) +
    geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank())
}

######￥￥￥￥￥￥￥￥￥￥￥￥￥￥￥￥￥￥


heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')




CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[18:22],
             y = names(heatmap.data)[1:17],
             col = c("blue1", "skyblue", "white", "pink", "red"))



module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'turquoise') %>% 
  rownames()
#####以上都对######
