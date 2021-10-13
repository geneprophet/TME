#kidney rat
library(Seurat)
library(dplyr)
library(Seurat)
library(patchwork)
data1 <- Read10X(data.dir = "D:/doctor/kidney/RAT/SD/outs/filtered_feature_bc_matrix")
data2 <- Read10X(data.dir = "D:/doctor/kidney/RAT/SD1/outs/filtered_feature_bc_matrix")
data3 <- Read10X(data.dir = "D:/doctor/kidney/RAT/SHR/outs/filtered_feature_bc_matrix")
data4 <- Read10X(data.dir = "D:/doctor/kidney/RAT/SHR1/outs/filtered_feature_bc_matrix")
rat1 <- CreateSeuratObject(counts = data1, project = "SD")
rat2 <- CreateSeuratObject(counts = data2, project = "SD1")
rat3 <- CreateSeuratObject(counts = data3, project = "SHR")
rat4 <- CreateSeuratObject(counts = data4, project = "SHR1")

all <- merge(rat1,y=c(rat2,rat3,rat4),add.cell.ids = c("1","2","3","4"),project = "rat")
all
