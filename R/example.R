# Title     : 16 Days clustering pipeline 3.0
# Objective : Extract lncRNA sequences from the different populations at 3 different times and pool them together by population and timepoint
# Created by: DAVID
# Created on: 01/04/2020

library(tidyverse)
library(Seurat)
library(patchwork)

#Read the data
#-------------
day16.data <-Read10X(data.dir = "data/V3/FGF8plus_day16_V3/filtered_feature_bc_matrix")
day16 <- CreateSeuratObject(counts = day16.data, project="lncRNA", min.cells = 5, min.features = 100)
#min cells takes features present in at leas that many cells.
#min features includes cells with at least that many features

#Pre-Processing
#--------------
day16[["percent.mt"]]<-PercentageFeatureSet(day16, pattern = "^MT-")#calculate mt% to clean it later

#visualize features
VlnPlot(day16, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
#check nFeature plot, pretty wide


plot1 <- FeatureScatter(day16, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(day16, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2



day16 <- subset(day16, subset = nFeature_RNA > 1100 & nFeature_RNA < 4000 & percent.mt < 3)
#subset cells with at least 500 genes and max 4000 and less than 3% mt contamination


#Normalize the data
#------------------

day16 <- NormalizeData(day16,normalization.method = "LogNormalize", scale.factor = 10000)
#Is logarithmic normalization good enough(?)


#Feature selection
#-----------------
day16 <- FindVariableFeatures(day16, selection.method = "vst", nfeatures = 500)
#taking 3000 most variable features

top20 <- head(VariableFeatures(day16),20)
plot1 <- VariableFeaturePlot(day16)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot1 + plot2


#Scaling the data
#----------------
all.genes <- rownames(day16)
day16 <- ScaleData(day16, features = all.genes)


# Linear dimensional reduction
#-----------------------------
day16 <- RunPCA(day16, features = VariableFeatures(object = day16))

VizDimLoadings(day16, dims = 1:4, reduction = "pca")

DimPlot(day16, reduction = "pca")
DimHeatmap(day16, dims = 1, cells = 1000, balanced = TRUE)
DimHeatmap(day16, dims = 1:21, cells = 1000, balanced = TRUE)


#Dimensionality
#--------------

day16 <- JackStraw(day16,num.replicate = 100)
day16 <- ScoreJackStraw(day16, dims = 1:20)
JackStrawPlot(day16, dims = 1:20)
ElbowPlot(day16) #I'll take 3

#Clustering
#----------
day16 <- FindNeighbors(day16,dims = 1:3)
day16 <- FindClusters(day16, resolution = 0.05)

#Non-linear reduction
#--------------------
day16 <- RunUMAP(day16, dims = 1:3)
day16 <-RunTSNE(day16, dims = 1:3)
DimPlot(day16, reduction = "umap")#With 0.1 resolution looks better
DimPlot(day16, reduction = "tsne") #With 0.2 resolution looks better
#saveRDS(day16, file = "output/day16.rds")


#Find clusters
#-------------
day16.markers <- FindAllMarkers(day16, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2)
#cluster
day16.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
#Sort the clusters

#-------------
head(day16.markers,20)
VlnPlot(day16, features = c("TH"))
VlnPlot(day16, features = c("GFAP"))
top3<-day16.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)

FeaturePlot(day16, features = top3$gene, reduction = "tsne")
FeaturePlot(day16, features = top3$gene, reduction = "umap")
#--------------

#HEatmap
top15 <- day16.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)
DoHeatmap(day16, features = top15$gene) + NoLegend()
