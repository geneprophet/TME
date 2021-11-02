#
library(readr)
GSE91061 <- read_csv("D:/downloads/GSE91061_BMS038109Sample.hg19KnownGene.raw.csv")
library('org.Hs.eg.db')
columns(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
head(keys(org.Hs.eg.db, keytype = 'SYMBOL'))
symbol <- select(org.Hs.eg.db,keys = as.character(GSE91061$...1), columns = 'SYMBOL',keytype = 'ENTREZID')
GSE91061 <- cbind(gene=symbol$SYMBOL,GSE91061[,-1])

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

expressionMatrix <- as.matrix(GSE91061[,-1])
rownames(expressionMatrix) = GSE91061$gene
expressionMatrix[1:5,1:5]
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
pheatmap(gsva.es,cluster_rows =F,show_colnames = F) -> res
cutree(res$tree_col,k=2) -> x1


###read the survival data
library(readxl)
mmc2 <- read_excel("D:/downloads/mmc2.xlsx", 
                   skip = 2)
load('D:/downloads/Riaz_data.RData')
survivalData <- Riaz_data$Surv.info
sapply(names(x1),function(x){a=unlist(strsplit(x,split = "_"));if(a[2]=="On") {return(a[1])} else{return("Pre")}}) -> names(x1)
survivalData2 <- merge(survivalData,data.frame(PatientID=names(x1),cluster=x1),by='PatientID',all.x = T,sort = F)
survivalData2 <- survivalData2[which(!is.na(survivalData2$cluster)),]


library(survminer)
library(survival)
fit <- survfit(Surv(OS,OS_SOR) ~ cluster,data=survivalData2)
#summary(fit)
ggsurvplot(fit, data=survivalData2,pval = TRUE,ggtheme = theme_minimal())

