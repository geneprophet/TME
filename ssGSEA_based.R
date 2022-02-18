##ssGSEA based method
##

library(readr)
IIS_TIS_signature <- read_delim("marker/IIS_TIS_signature.txt", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE)
all = unique(IIS_TIS_signature$`Cell type`)
IIS_TIS_geneset = list()
for (i in all) {
  b = IIS_TIS_signature$Symbol[which(IIS_TIS_signature$`Cell type`==i)]
  #assign(i,b)
  #if(length(b)>1){
    IIS_TIS_geneset[[i]] = b
  #}
}

library("GSVA")
gsva.es <- gsva(data$TPM, IIS_TIS_geneset, method="ssgsea", verbose=T)

# TIS_score = 
TIS_sinature = c("CD8 T cells","T helper cells","Tcm cells","Tem cells","Th1 cells","Th2 cells","Th17 cells","Treg cells")
TIS_score = apply(gsva.es,2,function(x){return(sum(x[TIS_sinature]))})
IIS_score = colSums(gsva.es)

expMarker = merge(data$Samples,data.frame(Sample=names(TIS_score),IIS_score=IIS_score,TIS_score=TIS_score),by="Sample")
library(ggpubr)
ggboxplot(expMarker,x="Resp_NoResp",y="TIS_score",color = "Resp_NoResp",
          palette = "npg",add = "jitter") + stat_compare_means()
ggboxplot(expMarker,x="Resp_NoResp",y="IIS_score",color = "Resp_NoResp",
          palette = "npg",add = "jitter") + stat_compare_means()



##################
##innate anti-PD-1 resistance (IPRES) signature
# We row-normalized
# the GSVA scores of each gene set in the IPRES signature across the samples
# from the four cohorts
# The IPRES (enrichment) score was defined as the average
# Z score across all gene sets in the IPRES signature
library("qusage")
C2_chemical_genetic = read.gmt('marker/c2.cgp.v7.5.1.symbols.gmt')
C6_oncogenic = read.gmt('marker/c6.all.v7.5.1.symbols.gmt')
C7_immunologic = read.gmt('marker/c7.all.v7.5.1.symbols.gmt')
C_all = c(C2_chemical_genetic,C6_oncogenic,C7_immunologic)
IPRES_pathway = read_csv("marker/IPRES_signaures.txt", 
                                            col_names = FALSE)
length(which(IPRES_pathway$X1 %in% names(C_all)))
IPRES_signatures = list()
for (i in IPRES_pathway$X1) {
  print(i)
  IPRES_signatures[[i]] = C_all[[i]]
}
gsva.es <- gsva(data$TPM, IPRES_signatures, method="ssgsea", verbose=T)
gsva.es = scale(t(gsva.es))
IPRES_score = apply(gsva.es,1,mean)
expMarker = merge(data$Samples,data.frame(Sample=names(IPRES_score),IPRES_score=IPRES_score),by="Sample")
library(ggpubr)
ggboxplot(expMarker,x="Resp_NoResp",y="IPRES_score",color = "Resp_NoResp",
          palette = "npg",add = "jitter") + stat_compare_means()


##############
##xCell: 64 cell type signatures
# devtools::install_github('dviraran/xCell')
library(xCell)
#expression matrix with genes in rows and samples in columns
xCell = xCellAnalysis(data$TPM)




#############
##immune microenvironment score (IMS) from gastric cancer
IMS_signature <- read_delim("marker/IMS_signature.txt", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)
IMS_signature_geneset = list()
all = unique(IMS_signature$`Immune cell type`)
for (i in all) {
  b = IMS_signature$Gene[which(IMS_signature$`Immune cell type` == i)]
  IMS_signature_geneset[[i]] = b
}
IMS_signature_meta <- read_delim("marker/IMS_signature_meta.txt", 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)
IMS_signature_HR_more = IMS_signature_meta$`Immune Cell`[which(IMS_signature_meta$HR>1)]
IMS_signature_HR_less = IMS_signature_meta$`Immune Cell`[which(IMS_signature_meta$HR<1)]

gsva.es <- gsva(data$TPM, IMS_signature_geneset, method="ssgsea", verbose=T)

IMS_score <- apply(gsva.es, 2, function(x){
  NES1 = sum(x[which(is.element(names(x),IMS_signature_HR_more))])
  NES2 = sum(x[which(is.element(names(x),IMS_signature_HR_less))])
  return(NES2-NES1)
})
expMarker = merge(data$Samples,data.frame(Sample=names(IMS_score),IMS_score=IMS_score),by="Sample")
library(ggpubr)
ggboxplot(expMarker,x="Resp_NoResp",y="IMS_score",color = "Resp_NoResp",
          palette = "npg",add = "jitter") + stat_compare_means()





