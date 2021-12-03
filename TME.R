###get protein coding genes from ensembl
# library(biomaRt)
# listMarts(host="uswest.ensembl.org")
# ensembl_us_west = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org")
# head(listDatasets(ensembl_us_west))
# ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# head(listFilters(ensembl))
# head(listAttributes(ensembl))
# 
# pcgs <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','gene_biotype','chromosome_name','start_position','end_position'), filters =
#                 c('biotype','chromosome_name'), values =list(biotype="protein_coding",chromosome=c(as.character(seq(1:22)),"X","Y","MT")), mart = ensembl)
# #listFilters(ensembl)
# 

# signatures of 29 Fges
library(readr)
X29Fges <- read_delim("29Fges.txt", delim = "\t", 
                      escape_double = FALSE, trim_ws = TRUE)

all = unique(X29Fges$`Gene signature`)
geneSet = list()
for (i in all) {
  b = X29Fges$Gene[which(X29Fges$`Gene signature`==i)]
  #assign(i,b)
  geneSet[[i]] = b
}


########################################################

immunecelltypemarkers <- read_delim("immunecelltypemarkers",
                                    delim = "\t", escape_double = FALSE,
                                    trim_ws = TRUE)
geneSet2 = list()
for (i in immunecelltypemarkers$`cell type`) {
  trimws(unlist(strsplit(immunecelltypemarkers$markers[which(immunecelltypemarkers$`cell type`==i)],split = ","))) -> a
  geneSet2[[paste(i,"immune cell")]] = a
}
##impute the signatures
geneSet$`B cells`=union(geneSet$`B cells`,geneSet2$`Plasma B cells immune cell`)
geneSet$`T cells`= union(geneSet$`T cells`,geneSet2$`T cells naive immune cell`)
geneSet$Treg = union(geneSet$Treg,geneSet2$`T cells regulatory immune cell`)
geneSet$`NK cells` = union(geneSet$`NK cells`,geneSet2$`NK immune cell`)
geneSet[["INF-gamma"]] = c("ID01","CXCL10","CXCL9","HLA-DRA","STAT1","INFG")
#geneSet[["INF-gamma"]] = c("IDO1","PDCD1LG2","CD27","CD8A","LAG3","CD274","CXCR6","CMKLR1","NKG7","CCL5","PSMB10","IDO1","CXCL9","HLA-DQA1","CD276","STAT1","HLA-DRB1","HLA-E")
geneSet[["IMPRES_Activate"]] = c("OX40L","CD27","OX40L","CD40","CD28","CD137L","CD40","HVEM")
geneSet[["IMPRES_Inhibitor"]] = c("PD-1","CTLA4","CD86","CD80","PDL-1","VISTA","TIM-3","CD200","CD276")
geneSet[["myCAF"]] = c("MYL9","CALD1","MMP11","HOPX","BGN","IGFBP7","TPM2","CTHRC1","ACTA2","TAGLN","INHBA","COL10A1","TPM1","POSTN","GRP","CST1")
geneSet[["iCAF"]] = c("C3","DUSP1","FBLN1","LMNA","CLU","CCDC80","MYC","EFEMP1","HAS1","NR4A1","CFD","ANXA1","CXCL12","FGF7","KLF4","EMP1","GPRC5A","SRPX","MT2A","MEDAG","IGF1","MGST1","MCL1","CEBPD","S100A10","UAP1","TNXB","CEBPB","PNRC1","SOCS3","PTGDS","FOSB","NFKBIA","CXCL2","THBS1","CCL2","OGN","GSN","DPT","PLA2G2A","NAMPT","ITM2A","RGCC","JUND","NNMT","ZFP36","PIM1","CPE","GFPT2","SOD2","KDM6B","FSTL1","FBLN2","NR4A3","MFAP5","ABL2","SGK1","CILP","UGDH","FBLN5","ADAMTS1","ADH1B","WISP2","GPX3","S100A4","IL6","HAS2","PLAC9","IGFBP6","FBN1","BDKRB1","TPPP3","RASD1","MT1A","CXCL14","PI16","APOE","IL8","ARC","PTX3","TNFAIP6","MT1E","MT1X","CXCL1")
#geneSet$`Cancer-associated fibroblasts` = NULL

geneSet3 = c(geneSet,geneSet2[c(6:15,21:25)])

# # impute from the ESTIMATE  signature
# library(readr)
# SI_geneset <- read_delim("SI_geneset.gmt",
#                          delim = "\t", escape_double = FALSE,
#                          col_names = FALSE, trim_ws = TRUE)
# geneSet3$StromalSignature = as.character(SI_geneset[1,-c(1:2)])
# geneSet3$ImmuneSignature = as.character(SI_geneset[2,-c(1:2)])

# import SKCM expression matrix
###constract the expression matrix
# HiSeqV2_PANCAN <- read_delim("./skcm/HiSeqV2",
#                  delim = "\t", escape_double = FALSE,
#                  trim_ws = TRUE)
##gdc tcga
library(readr)
TCGA_SKCM_htseq_fpkm <- read_delim("skcm/TCGA-SKCM.htseq_fpkm.tsv",
                                   delim = "\t", escape_double = FALSE,
                                   trim_ws = TRUE)
##
library("rtracklayer")
gtf_data = import('gencode.v22.annotation.gtf') #gtf的路径
#这里使用import导入gtf文件， 生成一个GRangs对象
gtf_data = as.data.frame(gtf_data)
gtf_data <- gtf_data[which(gtf_data$type == "gene" & gtf_data$gene_type=="protein_coding"),]

TCGA_SKCM = merge(TCGA_SKCM_htseq_fpkm,data.frame(Ensembl_ID=gtf_data$gene_id,geneName=gtf_data$gene_name),by='Ensembl_ID',all.x = T,sort = F)
TCGA_SKCM[which(!is.na(TCGA_SKCM$geneName)),] -> TCGA_SKCM
TCGA_SKCM$Ensembl_ID = TCGA_SKCM$geneName
TCGA_SKCM <- subset(TCGA_SKCM,select= -geneName)
HiSeqV2_PANCAN <- TCGA_SKCM
##gdc tcga

expressionMatrix <- as.matrix(HiSeqV2_PANCAN[,-1])

# ###################
##log2(fpkm+1) to fpkm
2^expressionMatrix -1 -> expressionMatrix
##fpkm to tpm
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
fpkmToTpm(expressionMatrix) -> expressionMatrix
rownames(expressionMatrix) = HiSeqV2_PANCAN$Ensembl_ID
# ##################
# 
# rownames(expressionMatrix) = HiSeqV2_PANCAN$sample
expressionMatrix[1:5,1:5]
sapply(colnames(expressionMatrix), function(u){unlist(strsplit(u, split="-"))[4]}) -> type
expressionMatrix = expressionMatrix[, type=="06A" | type=="01A"]
dim(expressionMatrix)

# #############
# ##count to tpm
# expressionMatrix = 2^expressionMatrix - 1 
# expressionMatrix=GeoTcgaData::countToTpm_matrix(expressionMatrix)

# apply(expressionMatrix, 1, sum) -> check
# expressionMatrix = expressionMatrix[check!=0, ]
# dim(expressionMatrix)


# ###highly variable 10000 genes
# var(t(expressionMatrix)) -> v
# diag(v) -> v
# names(sort(v,decreasing = T)[1:10000]) -> hvg
# length(which(X29Fges$Gene %in% hvg))
# expressionMatrix[hvg,] -> expressionMatrix
# dim(expressionMatrix)


## ssgsea score
library("GSVA")
gsva.es <- gsva(expressionMatrix, geneSet,method="ssgsea", verbose=T)
# skcm_sclaed_score = read.csv(file = 'D:/projects/PycharmProjects/MFP/skcm/signature_scores_scaled.csv',row.names = 1)
# skcm_sclaed_score = t(skcm_sclaed_score)
# cor(gsva.es,skcm_sclaed_score) -> cc
# length(which(diag(cc)>0.9))
# plot(diag(cc))
dim(gsva.es)
gsva.es[1:5, 1:5]


############
##vega score
vega <- read.csv('D:/projects/PycharmProjects/vega/data/skcm/results50.csv',header = F)
dim(vega)
vega= t(vega)
colnames(vega) = colnames(gsva.es)
############

{
#################
# louvain clustering
##Generates a graph from the similarity_matrix (square dataframe). Each sample is a node, similarity - edge weight. 
##Only nodes with at least 1 edge with weight above the threshold will be present in the final graph
cor(skcm_sclaed_score)+1 -> cc
##上三角矩阵，并将对角线赋值为0
cc[!upper.tri(cc,diag = F)] <- 0
threshold = 1.51
which(cc > threshold, arr.ind = T) -> ii 
gd = c()
for (i in 1:nrow(ii)) {
  x=ii[i,]
  edge = c(rownames(cc)[x[1]],colnames(cc)[x[2]],round(cc[x[1],x[2]],digits = 2))
  gd = rbind(gd,edge)
}
dim(gd)
#Edges
edges <- as.data.frame(gd)
colnames(edges) <- c("from", "to", "weight")
#Create graph for the algorithms
library(igraph)
g <- graph_from_data_frame(edges, directed = FALSE)
# Louvain cluster
lc <- cluster_louvain(g)
modularity(lc)
membership(lc) -> clusterResult
communities(lc)
plot(lc, g)
table(clusterResult)
clusterResult ->x1
names(x1) = substr(names(x1),1,15)
}


library('pheatmap')
# hierarchal cluster
pheatmap(gsva.es,cluster_rows =F,show_colnames = F) -> res
cutree(res$tree_col,k=4) -> x1
names(x1) = substr(names(x1),1,15)

###read the survival data
survivalData <- read_delim("TCGA/SKCM/survival_SKCM_survival.txt",
                           delim = "\t", escape_double = FALSE,
                           trim_ws = TRUE)
# survivalData <- read_delim("skcm/TCGA-SKCM.survival.tsv",
#                            delim = "\t", escape_double = FALSE,
#                            trim_ws = TRUE)
survivalData2 <- merge(survivalData,data.frame(sample=names(x1),cluster=as.vector(x1)),by='sample',all.x = T,sort = F)
survivalData2 <- survivalData2[which(!is.na(survivalData2$cluster)),]

library(survminer)
library(survival)
fit <- survfit(Surv(OS.time,OS) ~ cluster,data=survivalData2)
#fit <- survfit(Surv(PFI.time,PFI) ~ cluster,data=survivalData2)
#summary(fit)
ggsurvplot(fit, data=survivalData2,pval = TRUE,ggtheme = theme_minimal())
ggsurvplot(fit, data=survivalData2,pval = TRUE,ggtheme = theme_minimal(),xlim=c(0,1800),break.time.by=300,risk.table = T)
ggsurvplot(fit, data=survivalData2,pval = TRUE,ggtheme = theme_minimal(),xlim=c(0,1800),break.time.by=300)
##HR of cluster 1,2,3,4
survivalData2$cluster[which(survivalData2$cluster!=4)] = rep(0,length(which(survivalData2$cluster!=1)))
res.cox <- coxph(Surv(OS.time,OS) ~ cluster,data=survivalData2)
HR_95CI <- function(x){ 
  x <- summary(x)
  HR <-signif(x$coef[2], digits=2);#exp(beta)
  HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
  HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
  # <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
  res<- c(HR,HR.confint.lower,HR.confint.upper)
  names(res)<-c("Hazard Ratio","95% Upper CI","95% lower CI")
  return(res)
}
HR_95CI(res.cox)
cox.zph(res.cox) -> temp
plot(temp)
###the log rank test is a popular test to test the null hypothesis of no difference in survival between two or more independent groups.
survdiff(Surv(OS.time,OS) ~ cluster,data=survivalData2)



##########################################
#pathwaycommons
geneSet = geneSet3
network <- read_delim("skcm/PathwayCommons12.All.hgnc.sif", 
                      delim = "\t", escape_double = FALSE, 
                      col_names = FALSE, trim_ws = TRUE)
dim(network)
names(network) <-c("PARTICIPANT_A","Type","PARTICIPANT_B")
network[1:3,1:3]
which(network$PARTICIPANT_A %in% rownames(expressionMatrix) & network$PARTICIPANT_B %in% rownames(expressionMatrix)) -> ii
network = network[ii, ]
dim(network)

# prepare edge sets
edgeSets = list()
for(k in 1:length(geneSet)){
  genes = geneSet[[k]]
  which(network$PARTICIPANT_A %in% genes & network$PARTICIPANT_B %in% genes ) -> ii
  paste0("E", ii) -> ee
  edgeSets[[ names(geneSet)[k] ]] = ee
}

crossEdges = c()
for(k1 in 1:(length(geneSet)-1) ){
  genes_1 = geneSet[[k1]]
  for(k2 in (k1+1):length(geneSet)){
    genes_2 = geneSet[[k2]]
    genes_3 = intersect(genes_1, genes_2)
    which(network$PARTICIPANT_A %in% genes_3 & network$PARTICIPANT_B %in% genes_3) -> ii
    crossEdges = c(crossEdges, ii)
  }
}
paste("E", crossEdges, sep="") -> ee
edgeSets[[ "crossEdges" ]] = ee

for(k1 in 1:(length(geneSet)-1) ){
  genes_1 = geneSet[[k1]]
  for(k2 in (k1+1):length(geneSet)){
    genes_2 = geneSet[[k2]]
    genes_3 = intersect(genes_1, genes_2)

    genes_1 = setdiff(genes_1, genes_3)
    genes_2 = setdiff(genes_2, genes_3)

    which(network$PARTICIPANT_A %in% genes_1 & network$PARTICIPANT_B %in% genes_2) -> ii_1
    which(network$PARTICIPANT_A %in% genes_2 & network$PARTICIPANT_B %in% genes_1) -> ii_2
    ii = union(ii_1, ii_2)
    if(length(ii) >= 5){
      paste("E", ii, sep="") -> ee
      edgeSets[[ paste("C_", k1, "_", k2, sep="") ]] = ee
    }
  }
}


#mean the expression of paired gene and perform ssgsea
### edge weights per sample
match(network$PARTICIPANT_A, rownames(expressionMatrix)) -> idx1
match(network$PARTICIPANT_B, rownames(expressionMatrix)) -> idx2

edgeWeight = c()
for(k in 1:ncol(expressionMatrix)){
  apply(cbind(expressionMatrix[idx1, k], expressionMatrix[idx2, k]), 1, mean) -> ee
  edgeWeight = cbind(edgeWeight, ee)
  cat(k,",",sep="")
}
dim(edgeWeight)
rownames(edgeWeight) = paste("E", 1:nrow(edgeWeight), sep="")
colnames(edgeWeight) = colnames(expressionMatrix)[1:ncol(edgeWeight)]
edgeWeight[1:5,1:5]
### edge set ssGSEA
gsva(edgeWeight, edgeSets, method="ssgsea") -> test2

# library(gplots)
# heatmap.2(test2, trace="none") -> fit
# cutree(as.hclust(fit$colDendrogram), 4) -> x2

pheatmap(test2,cluster_rows =F,show_colnames = F) -> res2
cutree(res2$tree_col,k=4) -> x2
names(x2) = substr(names(x2),1,15)

###louvion clustering
cor(test2)+1 -> cc2
##上三角矩阵，并将对角线赋值为0
cc2[!upper.tri(cc2,diag = F)] <- 0
threshold = 1.51
which(cc2 > threshold, arr.ind = T) -> ii 
gd = c()
for (i in 1:nrow(ii)) {
  x=ii[i,]
  edge = c(rownames(cc2)[x[1]],colnames(cc2)[x[2]],round(cc2[x[1],x[2]],digits = 2))
  gd = rbind(gd,edge)
}
dim(gd)
#Edges
edges <- as.data.frame(gd)
colnames(edges) <- c("from", "to", "weight")
#Create graph for the algorithms
library(igraph)
g <- graph_from_data_frame(edges, directed = FALSE)
# Louvain cluster
lc <- cluster_louvain(g)
modularity(lc)
membership(lc) -> clusterResult
plot(lc, g)
table(clusterResult)
clusterResult ->x1
names(x1) = substr(names(x1),1,15)
x2=x1


###survival analysis
survivalData3 <- merge(survivalData,data.frame(sample=names(x2),cluster=as.vector(x2)),by='sample',all.x = T,sort = F)

library(survminer)
library(survival)
fit2 <- survfit(Surv(OS.time,OS) ~ cluster,data=survivalData3)

#jpeg(file="myplot.jpeg",quality = 100)
ggsurvplot(fit2, data=survivalData3,pval = TRUE,ggtheme = theme_minimal())
ggsurvplot(fit2, data=survivalData3,pval = TRUE,ggtheme = theme_minimal(),xlim=c(0,2100),break.time.by=300)
save(list = ls(),file = 'skcm/edgeweight.Rdata')
#dev.off()

