#utils
write.gmt <- function(gs,file){
  sink(file)
  lapply(names(gs), function(i){
    cat(paste(c(i,i,gs[[i]]),collapse = '\t'))
    cat('\n')
  })
  sink()
}
write.gmt(geneSet,'new_signatures.txt')

write.table(t(HiSeqV2_PANCAN[,-1]),file = 'matrix.mtx',sep = "\t",col.names = F,row.names = F)
write.table(names(HiSeqV2_PANCAN),file = 'barcodes.tsv',sep = '\n',col.names = F,row.names = F,quote = F)
write.table(HiSeqV2_PANCAN[,1],file = 'gene.tsv',sep = '\n',col.names = F,row.names = F,quote = F)


######
library(Seurat)
library(loomR)
library(SeuratDisk)
load('G:/AD_seurat.Rdata')
N.loom <- as.loom(AD_seurat,filename="G:/AD.loom")
N.loom$close_all()

########
#use vega results to hclust and survival analysis
vega = read.csv('skcm/vega_results.csv',header = F)
vega = t(vega)
colnames(vega) = colnames(HiSeqV2_PANCAN)[-1]
library('pheatmap')
# hierarchal cluster
#pheatmap(vega[-21,],cluster_rows =F,show_colnames = F, scale = 'column') -> res
pheatmap(vega,cluster_rows =F,show_colnames = F, scale = 'column') -> res
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
