#kidney rat
library(dplyr)
library(Seurat)
library(patchwork)
data1 <- Read10X(data.dir = "/Users/kanghongen/Desktop/博士课题/kidney/自发性高血压大鼠/SD/outs/filtered_feature_bc_matrix")
data2 <- Read10X(data.dir = "/Users/kanghongen/Desktop/博士课题/kidney/自发性高血压大鼠/SD1/outs/filtered_feature_bc_matrix")
data3 <- Read10X(data.dir = "/Users/kanghongen/Desktop/博士课题/kidney/自发性高血压大鼠/SHR/outs/filtered_feature_bc_matrix")
data4 <- Read10X(data.dir = "/Users/kanghongen/Desktop/博士课题/kidney/自发性高血压大鼠/SHR1/outs/filtered_feature_bc_matrix")
rat1 <- CreateSeuratObject(counts = data1, project = "SD")
rat2 <- CreateSeuratObject(counts = data2, project = "SD1")
rat3 <- CreateSeuratObject(counts = data3, project = "SHR")
rat4 <- CreateSeuratObject(counts = data4, project = "SHR1")

all <- merge(rat1,y=c(rat2,rat3,rat4),add.cell.ids = c("1","2","3","4"),project = "rat")
all
rownames(data1)[grep('^Mt-',rownames(data1))]
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
all[["percent.mt"]] <- PercentageFeatureSet(all, pattern = "^Mt-")
head(all@meta.data, 5)
# Visualize QC metrics as a violin plot
VlnPlot(all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

all <- subset(all, subset = nCount_RNA < 30000 & nFeature_RNA < 4000 & percent.mt < 50)

VlnPlot(all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
head(all@assays$RNA@counts)



all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 40000)

#Identification of highly variable features 
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(all), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(all)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
plot2

#Scaling the data
all.genes <- rownames(all)
all <- ScaleData(all, features = all.genes,)
all <- ScaleData(all, vars.to.regress = c("percent.mt","nCount_RNA"))

#Perform linear dimensional reduction
all <- RunPCA(all, features = VariableFeatures(object = all))
# Examine and visualize PCA results a few different ways
print(all[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(all, dims = 1:2, reduction = "pca")
DimPlot(all, reduction = "pca")
DimHeatmap(all, dims = 1:30, cells = 500, balanced = TRUE)

#Determine the ‘dimensionality’ of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
all <- JackStraw(all, num.replicate = 100)
all <- ScoreJackStraw(all, dims = 1:20)
JackStrawPlot(all, dims = 1:15)
ElbowPlot(all)


