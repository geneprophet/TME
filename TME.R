###get protein coding genes from ensembl
library(biomaRt)
listMarts(host="uswest.ensembl.org")
ensembl_us_west = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org")
head(listDatasets(ensembl_us_west))
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
head(listFilters(ensembl))
head(listAttributes(ensembl))

pcgs <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','gene_biotype','chromosome_name','start_position','end_position'), filters =
                c('biotype','chromosome_name'), values =list(biotype="protein_coding",chromosome=c(as.character(seq(1:22)),"X","Y","MT")), mart = ensembl)
#listFilters(ensembl)


## import SKCM expression matrix
library(readr)
HiSeqV2_PANCAN <- read_delim("./skcm/HiSeqV2",
                  delim = "\t", escape_double = FALSE,
                  trim_ws = TRUE)

expressionMatrix <- as.matrix(HiSeqV2_PANCAN[,-1])
rownames(expressionMatrix) = HiSeqV2_PANCAN$sample
expressionMatrix[1:5,1:5]
sapply(colnames(expressionMatrix), function(u){strsplit(u, split="-")[[1]][4]}) -> type
expressionMatrix = expressionMatrix[, type=="06"]
dim(expressionMatrix)

apply(expressionMatrix, 1, sum) -> check
expressionMatrix = expressionMatrix[check!=0, ]
dim(expressionMatrix)

# signatures of 29 Fges
X29Fges <- read_delim("29Fges.txt", delim = "\t", 
                      escape_double = FALSE, trim_ws = TRUE)

all = unique(X29Fges$`Gene signature`)
geneSet = list()
for (i in all) {
  b = X29Fges$Gene[which(X29Fges$`Gene signature`==i)]
  #assign(i,b)
  geneSet[[i]] = b
}
library("GSVA")
gsva.es <- gsva(expressionMatrix, geneSet,method="ssgsea", verbose=FALSE)
dim(gsva.es)
gsva.es[1:5, 1:5]

heatmap.2(gsva.es, trace="none") -> fit
cutree(as.hclust(fit$colDendrogram), 4) -> x



library('pheatmap')
# hierarchal cluster 
# 
# mean(gsva.es[1,])
# sd(gsva.es[1,])
# kk <- function(x){
#   #print(x)
#   m = median(x)
#   s = sd(x)
#   return((x-m)/s)
# }
# scaled_ssgsea = apply(gsva.es,1,kk)


pheatmap(gsva.es,cluster_rows =F,show_colnames = F) -> res
cutree(res$tree_col,k=1)
plot(res$tree_col)
#pheatmap(t(scaled_ssgsea),cluster_rows =F,show_colnames = F)


##########################################
#pathwaycommons 
network <- read_delim("skcm/PathwayCommons12.All.hgnc.txt", 
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
library(gplots)
heatmap.2(test2, trace="none") -> fit
cutree(as.hclust(fit$colDendrogram), 4) -> x

match(names(x), anno[,1]) -> idx
cbind(x, anno[idx,]) -> x1
table(x1[,1], x1[,3])



#signatures_corr <- cor(gsva.es)
#pheatmap(signatures_corr,cluster_rows =F,show_colnames = F,show_rownames = F)

signatures_genes=unique(X29Fges$Gene)






###########################################################################################################


partner_a=vector(mode = "character",length = 0)
partner_b=vector(mode = "character",length = 0)
for (i in 1:(length(signatures_genes)-1)) {
  for (j in (i+1):length(signatures_genes)) {
    partner_a = c(partner_a,signatures_genes[i])
    partner_b = c(partner_b,signatures_genes[j])
  }
}
interactions = data.frame(partner_a,partner_b)

Cell = colnames(gsva.es)
cell_type = vector(mode = "character",length = 0)
for (i in 1:length(Cell)) {
  a=unlist(strsplit(Cell[i],split = "-"))
  if(a[4]=="01"){
    cell_type = c(cell_type,"Primary")
  }else if(a[4]=="06"){
    cell_type = c(cell_type,"Metastatic")
    
  }else if(a[4]=="07"){
    cell_type = c(cell_type,"Additional Metastatic")
    
  }else if(a[4]=="11"){
    cell_type = c(cell_type,"Solid Tissue Normal")
    
  }
  
}
meta = data.frame(Cell,cell_type)

#write.table(gsva.es,file = 'ssgsea.txt',sep = "\t")
write.table(interactions,file = "interactions.txt",sep = "\t",row.names = F)
write.table(meta,file = "meta.txt",sep = "\t",row.names = F)


