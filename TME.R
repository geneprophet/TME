#https://github.com/GSEA-MSigDB/ssGSEA-gpmodule
#https://github.com/broadinstitute/ssGSEA2.0
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GSVA")
library(GSVA)
p <- 10000 ## number of genes
n <- 30    ## number of samples
## simulate expression values from a standard Gaussian distribution
X <- matrix(rnorm(p*n), nrow=p,
            dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))
## sample gene set sizes
gs <- as.list(sample(10:100, size=100, replace=TRUE))
## sample gene sets
gs <- lapply(gs, function(n, p)
  paste0("g", sample(1:p, size=n, replace=FALSE)), p)
names(gs) <- paste0("gs", 1:length(gs))
gsva.es <- gsva(X, gs,method="ssgsea", verbose=FALSE)
dim(gsva.es)
gsva.es[1:5, 1:5]
####################################################################
## import SKCM expression matrix
library(readr)
HiSeqV2_PANCAN <- read_delim("HiSeqV2_PANCAN",
                  delim = "\t", escape_double = FALSE,
                  trim_ws = TRUE)

expressionMatrix <- as.matrix(HiSeqV2_PANCAN[,-1])
rownames(expressionMatrix) = HiSeqV2_PANCAN$Gene
expressionMatrix[1:5,1:5]

X29Fges <- read_delim("29Fges.txt", delim = "\t", 
                      escape_double = FALSE, trim_ws = TRUE)

all = unique(X29Fges$`Gene signature`)
geneSet = list()
for (i in all) {
  b = X29Fges$Gene[which(X29Fges$`Gene signature`==i)]
  #assign(i,b)
  geneSet[[i]] = b
}
gsva.es <- gsva(expressionMatrix, geneSet,method="ssgsea", verbose=FALSE)
dim(gsva.es)
gsva.es[1:5, 1:5]

expressionMatrix["FLT1",1:4]
write.csv(gsva.es,file = "ssgsea.txt")

library('pheatmap')

<<<<<<< HEAD
=======
# hierarchal cluster 
mean(gsva.es[1,])
sd(gsva.es[1,])
kk <- function(x){
  #print(x)
  m = median(x)
  s = sd(x)
  return((x-m)/s)
}
scaled_ssgsea = apply(gsva.es,1,kk)

pheatmap(gsva.es,cluster_rows =F,show_colnames = F)
pheatmap(t(scaled_ssgsea),cluster_rows =F,show_colnames = F)

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



#######
library(readr)
ccc_result <- read_delim("ccc_result.txt", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)
ccc = ccc_result[,-1]
rowname=paste(ccc$partner_a,ccc$partner_b,sep = '_')
ccc = ccc[,-c(1,2)]
rownames(ccc)=rowname
ccc2 = ccc[,c(1,6,11,16)]
rownames(ccc2)=rownames(ccc)

ccc3=apply(ccc2,1,function(x){if(min(x)<0.05) {return(x)} else{return()}})


>>>>>>> 3c2ba653b87c8620a952007b65c1106a95f566e4
