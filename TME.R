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
rownames(expressionMatrix) = HiSeqV2_PANCAN$sample
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


