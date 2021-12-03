#utils
write.gmt <- function(gs,file){
  sink(file)
  lapply(names(gs), function(i){
    cat(paste(c(i,i,gs[[i]]),collapse = '\t'))
    cat('\n')
  })
  sink()
}
write.gmt(geneSet3,'new_signatures.txt')

library(clusterProfiler)
hallmarks = read.gmt('hallmarks.gmt')
length(unique(hallmarks$gene))
reactom = read.gmt('reactom.gmt')
length(unique(reactom$gene))
gbm_aracne = read.gmt('gbm_aracne.gmt')
length(unique(gbm_aracne$gene))

# ###########
# ## TCGAbiolinks
# library(TCGAbiolinks)
# library(DT)
# query = GDCquery(project = 'TCGA-SKCM',
#                 data.category = 'Transcriptome Profiling',
#                 data.type = "Gene Expression Quantification",
#                 #workflow.type = 'HTSeq - Counts',
#                 workflow.type = 'HTSeq - FPKM',
#                 sample.type = c("Primary Tumor","Metastatic")
#                 )
# datatable(
#     getResults(query), 
#     filter = 'top',
#     options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
#     rownames = FALSE )
# GDCdownload(query,directory = './GDC/')
# skcm.exp <- GDCprepare(query = query,directory = './GDC/', save = TRUE, save.filename = "skcmExp.rda")
# skcm.exp@assays@data$`HTSeq - FPKM`[1:3,1:3]
# skcm.exp@colData@rownames
# skcm.exp@rowRanges@elementMetadata$external_gene_name
# expressionMatrix <- skcm.exp@assays@data$`HTSeq - FPKM`
# expressionMatrix[1:3,1:3]
# colnames(expressionMatrix) = skcm.exp@colData@rownames
# rownames(expressionMatrix) = skcm.exp@rowRanges@elementMetadata$external_gene_name
# expressionMatrix[1:3,1:3]


#########
##fpkmToTpm
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}


##
library("rtracklayer")
gtf_data = import('gencode.v22.annotation.gtf') #gtf的路径
#这里使用import导入gtf文件， 生成一个GRangs对象
gtf_data = as.data.frame(gtf_data)
gtf_data <- gtf_data[which(gtf_data$type == "gene" & gtf_data$gene_type=="protein_coding"),]



######
library(Seurat)
library(loomR)
library(SeuratDisk)
load('G:/AD_seurat.Rdata')
N.loom <- as.loom(AD_seurat,filename="G:/AD.loom")
N.loom$close_all()
##
##conver from anndata to seurat
Convert("adata/adata_scvi_cm2.h5ad", dest = "h5seurat", overwrite = TRUE)
seurat <- LoadH5Seurat("adata/adata_scvi_cm2.h5seurat")
library(rhdf5)
structure <- h5ls("adata/adata_scvi_cm2.h5seurat")
cells <- h5read("adata/adata_scvi_cm2.h5seurat","cell.names")
features <- h5read("adata/adata_scvi_cm2.h5seurat","assays/RNA/counts")

###############

########
#use vega results to hclust and survival analysis
vega = read.csv('skcm/vega_results50.csv',header = F)
vega = t(vega)
colnames(vega) = colnames(HiSeqV2_PANCAN)[-1]
library('pheatmap')
# hierarchal cluster
#pheatmap(vega[-21,],cluster_rows =F,show_colnames = F, scale = 'column') -> res
pheatmap(vega,cluster_rows =F,show_colnames = F) -> res
cutree(res$tree_col,k=4) -> x1
###read the survival data
survivalData <- read_delim("TCGA/SKCM/survival_SKCM_survival.txt",
                           delim = "\t", escape_double = FALSE,
                           trim_ws = TRUE)
survivalData2 <- merge(survivalData,data.frame(sample=names(x1),cluster=x1),by='sample',all.x = T,sort = F)
survivalData2 <- survivalData2[which(!is.na(survivalData2$cluster)),]

library(survminer)
library(survival)
fit <- survfit(Surv(OS.time,OS) ~ cluster,data=survivalData2)
#fit <- survfit(Surv(PFI.time,PFI) ~ cluster,data=survivalData2)
#summary(fit)
ggsurvplot(fit, data=survivalData2,pval = TRUE,ggtheme = theme_minimal())
ggsurvplot(fit, data=survivalData2,pval = TRUE,ggtheme = theme_minimal(),xlim=c(0,2000),break.time.by=300)

##################
#Louvain R
library(igraph)

###################
#缩放表达值到0-1之后再作为VEGA输入进行Deep Learning
# expressionMatrix_scaled <- 
dim(HiSeqV2_PANCAN)

apply(HiSeqV2_PANCAN[,-1],2,function(x){return((x-min(x))/(max(x)-min(x)))}) -> HiSeqV2_Scaled
rownames(HiSeqV2_Scaled) = HiSeqV2_PANCAN$sample
write.table(HiSeqV2_Scaled,file = "skcm/HiSeqV2_Scaled",sep = '\t',
            quote = F,col.names = F,row.names = T)

#####################
#rank ssgsea scores
ssgsea <- read.csv('D:/projects/PycharmProjects/MFP/skcm/signature_scores_scaled.csv',header = T)
ssgsea[1:4,1:4]
ssgsea2 = ssgsea[,-1]
rownames(ssgsea2) = ssgsea$X
t(ssgsea2) -> ssgsea2
ssgsea2[1:4,1:4]
#pheatmap(ssgsea2,cluster_rows =F,show_colnames = F, scale = 'column') -> res
pheatmap(gsva.es,cluster_rows =F,show_colnames = F) -> res
cutree(res$tree_col,k=4) -> x1
###read the survival data
survivalData <- read_delim("TCGA/SKCM/survival_SKCM_survival.txt",
                           delim = "\t", escape_double = FALSE,
                           trim_ws = TRUE)
survivalData2 <- merge(survivalData,data.frame(sample=names(x1),cluster=x1),by='sample',all.x = T,sort = F)
survivalData2 <- survivalData2[which(!is.na(survivalData2$cluster)),]
library(survminer)
library(survival)
fit <- survfit(Surv(OS.time,OS) ~ cluster,data=survivalData2)
#fit <- survfit(Surv(PFI.time,PFI) ~ cluster,data=survivalData2)
summary(fit)
ggsurvplot(fit, data=survivalData2,pval = TRUE,pval.method = T,ggtheme = theme_minimal())
ggsurvplot(fit, data=survivalData2,pval = TRUE,ggtheme = theme_minimal(),xlim=c(0,2000),break.time.by=600)
res.cox <- coxph(Surv(OS.time,OS) ~ cluster,data=survivalData2)
summary(res.cox)

###roc
library(pROC)
data("aSAH")
head(aSAH)
roc(aSAH$outcome,aSAH$s100b,smooth=T,auc=T,ci=T) -> r
plot(r)
roc.list = roc(outcome ~ s100b + ndka + wfns,data=aSAH)
g.list = ggroc(roc.list)
g.list



##################
##检查分类与MFP肿瘤微环境分型的一致性
library(readr)
annotation_panmi <- read_delim("skcm/annotation-panmi.tsv", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)
survivalData3 <- merge(survivalData2,data.frame(sample=annotation_panmi$Sample,MFP=annotation_panmi$MFP),by.x = "_PATIENT",by.y='sample',all.x = T,sort = F)

skcm_final_mfp_clusters <- read_delim("skcm/skcm_final_mfp_clusters.tsv", 
                                      delim = "\t", escape_double = FALSE, 
                                      col_names = FALSE, col_types = cols(X2 = col_character()), 
                                      trim_ws = TRUE)
survivalData4 <- merge(survivalData2,data.frame(sample=skcm_final_mfp_clusters$X1,MFP=skcm_final_mfp_clusters$X2),by.x = "sample",by.y='sample',all.x = T,sort = F)

cluster1 <- survivalData3[survivalData3$cluster==1,]
table(cluster1$MFP)
cluster2 <- survivalData3[survivalData3$cluster==2,]
table(cluster2$MFP)
cluster3 <- survivalData3[survivalData3$cluster==3,]
table(cluster3$MFP)
cluster4 <- survivalData3[survivalData3$cluster==4,]
table(cluster4$MFP)
library(survminer)
library(survival)
fit2 <- survfit(Surv(OS.time,OS) ~ MFP,data=survivalData3)
ggsurvplot(fit2, data=survivalData3,pval = TRUE,ggtheme = theme_minimal())
ggsurvplot(fit2, data=survivalData3,pval = TRUE,ggtheme = theme_minimal(),xlim=c(0,1800),break.time.by=600)
ggsurvplot(fit, data=survivalData3,pval = TRUE,ggtheme = theme_minimal(),xlim=c(0,1800),break.time.by=600)
##only TCGA SKCM samples to classfication
fit4 <- survfit(Surv(OS.time,OS) ~ MFP,data=survivalData4)
ggsurvplot(fit4, data=survivalData4,pval = TRUE,ggtheme = theme_minimal())
ggsurvplot(fit4, data=survivalData4,pval = TRUE,ggtheme = theme_minimal(),xlim=c(0,1800),break.time.by=600)



