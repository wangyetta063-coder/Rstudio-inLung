#使用harmony------2026/03/25-----------
# Finding markers in every cluster
# Finding conser//d markers 
# Finding markers DE between conditions

# setwd("~/Desktop/demo/single_cell_DEG/script")
set.seed(1234)

library(Seurat)
library(tidyverse)
library(harmony)
library(SeuratData)
library(ggplot2)
#官网下载harmony
if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("harmony")
# 1. 安装BiocManager
install.packages("BiocManager")
# 2. 安装harmony（正确方式）
BiocManager::install("harmony")
# 3. 加载
library(harmony)
#######然后开始了
library(Seurat)
library(SeuratData)
library(Seurat)
InstallData("ifnb")
# load dataset
LoadData("ifnb")
str(ifnb)
#以上都对
ifnb <- UpdateSeuratObject(ifnb)
# QC and filtering读取线粒体百分比
ifnb$mito.percent <- PercentageFeatureSet(ifnb, pattern = "^MT-")
View(ifnb@meta.data)
# explore QC可以看见绝大多数细胞没有线粒体百分比，但是没关系，接着下一步吧

# filter运行过滤程序，对细胞进行子集化，记录转录本计数超过800的细胞，保留所有基因数量大于200的细胞，线粒体百分比低于5%的细胞
ifnb
ifnb.filtered <- subset(ifnb, subset = nCount_RNA > 800 &
                           nFeature_RNA > 200 & 
                           mito.percent < 5)
ifnb.filtered <- NormalizeData(ifnb.filtered)
ifnb.filtered <- FindVariableFeatures(ifnb.filtered)
ifnb.filtered <- ScaleData(ifnb.filtered)
ifnb.filtered <- RunPCA(ifnb.filtered)
ElbowPlot(ifnb.filtered)
ifnb.filtered <- RunUMAP(ifnb.filtered, dims = 1:20, reduction = 'pca')
#想按条件分组看看STIM列包含了哪些信息。用STIM进行分组
before <- DimPlot(ifnb.filtered, reduction = 'umap', group.by = 'stim')
print(before)
# run Harmony -这步可以出图----------
ifnb.harmony <- ifnb.filtered %>%
   RunHarmony(group.by.vars = 'stim', plot_convergence = FALSE)

ifnb.harmony@reductions

ifnb.harmony.embed <- Embeddings(ifnb.harmony, "harmony")
ifnb.harmony.embed[1:10,1:10]
# Do UMAP and clustering using ** Harmony embeddings instead of PCA **
ifnb.harmony <- ifnb.harmony %>%
   RunUMAP(reduction = 'harmony', dims = 1:20) %>%
   FindNeighbors(reduction = "harmony", dims = 1:20) %>%
   FindClusters(resolution = 0.5)

# visualize 
after <- DimPlot(ifnb.harmony, reduction = 'umap', group.by = 'stim')

before|after
