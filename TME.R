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


for (path in list.dirs(path = './TCGA',recursive = F)) {
  expression_path = paste(path,list.files(path = path,pattern = "HiSeqV2"),sep = "/")
  survival_path =  paste(path,list.files(path = path,pattern = "survival"),sep = "/")
  print(path)
  
  #import expression matrix
  HiSeqV2 <- read_delim(expression_path,
                               delim = "\t", escape_double = FALSE,
                               trim_ws = TRUE)
  expressionMatrix <- as.matrix(HiSeqV2[,-1])
  rownames(expressionMatrix) = HiSeqV2$sample
  expressionMatrix[1:5,1:5]
  sapply(colnames(expressionMatrix), function(u){unlist(strsplit(u, split="-"))[4]}) -> type
  expressionMatrix = expressionMatrix[, type=="06" | type=="01"]
  dim(expressionMatrix)
  
  apply(expressionMatrix, 1, sum) -> check
  expressionMatrix = expressionMatrix[check!=0, ]
  dim(expressionMatrix)
  
  
  ## ssgsea score
  library("GSVA")
  gsva.es <- gsva(expressionMatrix, geneSet,method="ssgsea", verbose=FALSE)
  dim(gsva.es)
  gsva.es[1:5, 1:5]
  
  library('pheatmap')
  # hierarchal cluster 
  jpeg(file=paste(path,"pheatmap_29fges.jpeg",sep = "/"),quality = 100)
  pheatmap(gsva.es,cluster_rows =F,show_colnames = F) -> res
  dev.off()
  cutree(res$tree_col,k=4) -> x1
  
  ###read the survival data
  survivalData <- read_delim(survival_path, 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)
  survivalData2 <- merge(survivalData,data.frame(sample=names(x1),cluster=x1),by='sample',all.x = T,sort = F)
  survivalData2 <- survivalData2[which(!is.na(survivalData2$cluster)),]
  
  
  library(survminer) 
  library(survival) 
  fit <- survfit(Surv(OS.time,OS) ~ cluster,data=survivalData2)
  
  jpeg(file=paste(path,"survival_29fges.jpeg",sep = "/"),quality = 100)
  ggsurvplot(fit, data=survivalData2,pval = TRUE,ggtheme = theme_minimal())
  dev.off()
  
  ###########################
  #pathwaycommons 
  network <- read_delim("PathwayCommons12.All.hgnc.txt", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)
  dim(network)
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
  
  jpeg(file=paste(path,"pheatmap_edges.jpeg",sep = "/"),quality = 100)
  pheatmap(test2,cluster_rows =F,show_colnames = F) -> res2
  dev.off()
  
  cutree(res2$tree_col,k=4) -> x2
  
  ###survival analysis
  survivalData3 <- merge(survivalData,data.frame(sample=names(x2),cluster=x2),by='sample',all.x = T,sort = F)
  survivalData3 <- survivalData3[which(!is.na(survivalData3$cluster)),]
  
  library(survminer) # 加载包
  library(survival) # 加载包
  fit2 <- survfit(Surv(OS.time,OS) ~ cluster,data=survivalData3)
  
  jpeg(file=paste(path,"survival_edges.jpeg",sep = "/"),quality = 100)
  ggsurvplot(fit2, data=survivalData3,pval = TRUE,ggtheme = theme_minimal())
  dev.off()
  save(list = ls(),file = paste(path,"temp.Rdata",sep = "/"))
  
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
geneSet3 = c(geneSet,geneSet2[c(6:15,21:25)])
## impute from the ESTIMATE  signature
library(readr)
SI_geneset <- read_delim("SI_geneset.gmt", 
                         delim = "\t", escape_double = FALSE, 
                         col_names = FALSE, trim_ws = TRUE)
#geneSet3$StromalSignature = as.character(SI_geneset[1,-c(1:2)])
#geneSet3$ImmuneSignature = as.character(SI_geneset[2,-c(1:2)])
## import SKCM expression matrix
###constract the expression matrix


HiSeqV2_PANCAN <- read_delim("./skcm/HiSeqV2",
                 delim = "\t", escape_double = FALSE,
                 trim_ws = TRUE)
expressionMatrix <- as.matrix(HiSeqV2_PANCAN[,-1])
rownames(expressionMatrix) = HiSeqV2_PANCAN$sample
expressionMatrix[1:5,1:5]
sapply(colnames(expressionMatrix), function(u){unlist(strsplit(u, split="-"))[4]}) -> type
expressionMatrix = expressionMatrix[, type=="06" | type=="01"]
dim(expressionMatrix)

apply(expressionMatrix, 1, sum) -> check
expressionMatrix = expressionMatrix[check!=0, ]
dim(expressionMatrix)


## ssgsea score
library("GSVA")
gsva.es <- gsva(expressionMatrix, geneSet,method="ssgsea", verbose=FALSE)
dim(gsva.es)
gsva.es[1:5, 1:5]
#t(scale(t(gsva.es))) -> gsva.es2
#scale(gsva.es) -> gsva.es2

library('pheatmap')
# hierarchal cluster
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
#summary(fit)
ggsurvplot(fit, data=survivalData2,pval = TRUE,ggtheme = theme_minimal())
ggsurvplot(fit, data=survivalData2,pval = TRUE,ggtheme = theme_minimal(),xlim=c(0,2100),break.time.by=300)
##########################################
#pathwaycommons
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

###survival analysis
survivalData3 <- merge(survivalData,data.frame(sample=names(x2),cluster=x2),by='sample',all.x = T,sort = F)

library(survminer)
library(survival)
fit2 <- survfit(Surv(OS.time,OS) ~ cluster,data=survivalData3)

#jpeg(file="myplot.jpeg",quality = 100)
ggsurvplot(fit2, data=survivalData3,pval = TRUE,ggtheme = theme_minimal())
#dev.off()

