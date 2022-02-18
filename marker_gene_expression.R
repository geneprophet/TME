#marker gene expression can predict ICB response
##CD274(PD-L1), PDCD1(PD-1),PDCD1LG2(PD-L2)

data = Gide_data
expression = data$TPM["CD274",]
#names(expression)
expMarker = merge(data$Samples,data.frame(Sample=names(expression),PDL1=expression),by="Sample")
# t.test(PDL1 ~ Resp_NoResp, data = expMarker)
# wilcox.test(PDL1 ~ Resp_NoResp, data = expMarker)
# boxplot(PDL1 ~ Resp_NoResp, data = expMarker)
library(ggpubr)
ggboxplot(expMarker,x="Resp_NoResp",y="PDL1",color = "Resp_NoResp",
          palette = "npg",add = "jitter") + stat_compare_means(method.args = list(alternative = "greater"))



expression = data$TPM["PDCD1",]
#names(expression)
expMarker = merge(data$Samples,data.frame(Sample=names(expression),PD1=expression),by="Sample")
t.test(PD1 ~ Resp_NoResp, data = expMarker)
wilcox.test(PD1 ~ Resp_NoResp, data = expMarker)
boxplot(PD1 ~ Resp_NoResp, data = expMarker)
ggboxplot(expMarker,x="Resp_NoResp",y="PD1",color = "Resp_NoResp",
          palette = "npg",add = "jitter") + stat_compare_means()



expression = data$TPM["PDCD1LG2",]
#names(expression)
expMarker = merge(data$Samples,data.frame(Sample=names(expression),PDL2=expression),by="Sample")
t.test(PDL2 ~ Resp_NoResp, data = expMarker)
wilcox.test(PDL2 ~ Resp_NoResp, data = expMarker)
boxplot(PDL2 ~ Resp_NoResp, data = expMarker)
ggboxplot(expMarker,x="Resp_NoResp",y="PDL2",color = "Resp_NoResp",
          palette = "npg",add = "jitter") + stat_compare_means()



############
###cytolytic activity (CYT) scores
CYT_gene = c("GZMA","PRF1")
CYT_gene = intersect(CYT_gene,rownames(data$TPM))
expression = data$TPM[CYT_gene,]
mean_expression = apply(expression, 2, mean)
expMarker = merge(data$Samples,data.frame(Sample=names(mean_expression),CYT_score=as.numeric(mean_expression)),by="Sample")
ggboxplot(expMarker,x="Resp_NoResp",y='CYT_score',color = "Resp_NoResp",
          palette = "npg",add = "jitter") + stat_compare_means()



##############
##MHC-II
expression = data$TPM["HLA-DRA",]
expMarker = merge(data$Samples,data.frame(Sample=names(expression),MHCII=as.numeric(expression)),by="Sample")
ggboxplot(expMarker,x="Resp_NoResp",y='MHCII',color = "Resp_NoResp",
          palette = "npg",add = "jitter") + stat_compare_means()



# ##EMT scores
# EMT_signature <- read_csv("marker/EMT_signature.csv")



##INF-gamma signature
#INF-gamma_6 = c("ID01","CXCL10","CXCL9","HLA-DRA","STAT1","INFG")
INF_gamma_18 = c("TIGIT","PDCD1LG2","CD27","CD8A","LAG3","CD274","CXCR6","CMKLR1","NKG7","CCL5","PSMB10","IDO1","CXCL9","HLA-DQA1","CD276","STAT1","HLA-DRB1","HLA-E")
expression = data$TPM[INF_gamma_18,]
#names(expression)
expression = as.data.frame(t(expression))
expression["INF_gama"] = apply(expression, 1, function(x){return(mean(as.numeric(x)))})
expression["Sample"] = rownames(expression)
expMarker = merge(data$Samples,expression[c("Sample","INF_gama")],by="Sample")
ggboxplot(expMarker,x="Resp_NoResp",y="INF_gama",color = "Resp_NoResp",
          palette = "npg",add = "jitter") + stat_compare_means()


##Immunophenoscore (IPS)
source("marker/Immunophenogram/IPS.R")
DF = calculateIPS(data$TPM)
expMarker = merge(data$Samples,data.frame(Sample=DF$SAMPLE,IPS=DF$IPS),by="Sample")
ggboxplot(expMarker,x="Resp_NoResp",y="IPS",color = "Resp_NoResp",
          palette = "npg",add = "jitter") + stat_compare_means()


##Immune Signature (IS) 105 genes

##immuno-predictive score (IMPRES) coding in MATLAB
#15 pairs of gene,0-15, High scores predict good response
Gene1 = c("PDCD1","CD27","CTLA4","CD40","CD86","CD28","CD80","CD274","CD86","CD40","CD86","CD40","CD28","CD40","TNFRSF14")
Gene2 = c("TNFSF4","PDCD1","TNFSF4","CD28","TNFSF4","CD86","TNFSF9","VSIR","HAVCR2","PDCD1","CD200","CD80","CD276","CD274","CD86")
IMPRES_pairs = data.frame(Gene1,Gene2)
FS = apply(IMPRES_pairs, 1, function(x){
  if(x[1] %in% rownames(data$TPM) & x[2] %in% rownames(data$TPM)){
    expression1 = data$TPM[x[1],]
    expression2 = data$TPM[x[2],]
    return(as.numeric(expression1<expression2))
  }else{
    return(rep(0,ncol(data$TPM)))
  }
})
rownames(FS) = colnames(data$TPM)
IMPRES_score = apply(FS,1,sum)
expMarker = merge(data$Samples,data.frame(Sample=names(IMPRES_score),IMPRES_score=as.numeric(IMPRES_score)),by="Sample")
ggboxplot(expMarker,x="Resp_NoResp",y="IMPRES_score",color = "Resp_NoResp",
          palette = "npg",add = "jitter") + stat_compare_means()




######################
##TILs score,Total TILs score was calculated as the average of all cell scores with >0.6 correlations with CD45
library(readr)
TILs_cell_type_markers <- read_csv("marker/TILs cell type markers.csv")
TILs_cell_type_markers <- TILs_cell_type_markers[which(is.element(TILs_cell_type_markers$Gene,rownames(data$TPM))),]
celltypes = unique(TILs_cell_type_markers$`Cell Type`)
celltype_score = matrix(NA,nrow = ncol(data$TPM),ncol = length(celltypes))
dimnames(celltype_score) = list(colnames(data$TPM),celltypes)
for (celltype in celltypes) {
  #print(celltype)
  marker_genes = TILs_cell_type_markers$Gene[which(TILs_cell_type_markers$`Cell Type` == celltype)]
  if(length(marker_genes)==1){
    celltype_score[,celltype] = mean(log(data$TPM[marker_genes,]))
  }else{
    celltype_score[,celltype] = colMeans(log(data$TPM[marker_genes,]))
  }
}
# calc allTotTils as elsewhere:
high.cor.w.cd45 = cor(celltype_score)[,"CD45"]>0.6
alltottils = rowMeans(celltype_score[,high.cor.w.cd45])
celltype_score = cbind(alltottils,celltype_score)
####to be continue






##################
##MGAE-A : (MAGEA3, CSAG3, CSAG2,MAGEA2, MAGEA2B, CSAG1, MAGEA12, MAGEA6)
##CTLA-4 Blockade
MGAEA_genes = c("MAGEA3", "CSAG3", "CSAG2","MAGEA2", "MAGEA2B", "CSAG1", "MAGEA12", "MAGEA6")
expression = data$TPM[MGAEA_genes[1],]
#names(expression)
expMarker = merge(data$Samples,data.frame(Sample=names(expression),MGAEA=expression),by="Sample")
ggboxplot(expMarker,x="Resp_NoResp",y='MGAEA',color = "Resp_NoResp",
          palette = "npg",add = "jitter") + stat_compare_means(method = "t.test")


#####
##MHC proteins: SOX10, HLA-A, HLA-B, HLA-C, B2M
MHC_genes = c("SOX10", "HLA-A", "HLA-B", "HLA-C", "B2M")
expression = data$TPM[MHC_genes[2],]
#names(expression)
expMarker = merge(data$Samples,data.frame(Sample=names(expression),MHC=expression),by="Sample")
ggboxplot(expMarker,x="Resp_NoResp",y='MHC',color = "Resp_NoResp",
          palette = "npg",add = "jitter") + stat_compare_means(method = "t.test")



#######
#overall expression
discretize<-function(v,n.cat){
  q1<-quantile(v,seq(from = (1/n.cat),to = 1,by = (1/n.cat)))
  u<-matrix(nrow = length(v))
  for(i in 2:n.cat){
    u[(v>=q1[i-1])&(v<q1[i])]<-i
  }
  return(u)
}

get.semi.random.OE <- function(r,genes.dist.q,b.sign,num.rounds = 1000,full.flag = F){
  # Previous name: get.random.sig.scores
  sign.q<-as.matrix(table(genes.dist.q[b.sign]))
  q<-rownames(sign.q)
  idx.all<-c()
  B<-matrix(data = F,nrow = length(genes.dist.q),ncol = num.rounds)
  Q<-matrix(data = 0,nrow = length(genes.dist.q),ncol = num.rounds)
  for (i in 1:nrow(sign.q)){
    num.genes<-sign.q[i]
    if(num.genes>0){
      idx<-which(is.element(genes.dist.q,q[i]))
      for (j in 1:num.rounds){
        idxj<-sample(idx,num.genes) 
        Q[i,j]<-sum(B[idxj,j]==T)
        B[idxj,j]<-T
      }  
    }
  }
  rand.scores<-apply(B,2,function(x) colMeans(r$zscores[x,]))
  if(full.flag){return(rand.scores)}
  rand.scores<-rowMeans(rand.scores)
  return(rand.scores)
}

get.OE.bulk <- function(r,gene.sign = NULL,num.rounds = 1000,full.flag = F){
  set.seed(1234)
  r$genes.mean<-rowMeans(r$tpm)
  r$zscores<-sweep(r$tpm,1,r$genes.mean,FUN = '-')
  r$genes.dist<-r$genes.mean
  r$genes.dist.q<-discretize(r$genes.dist,n.cat = 50)
  r$sig.scores<-matrix(data = 0,nrow = ncol(r$tpm),ncol = length(gene.sign))
  sig.names<-names(gene.sign)
  colnames(r$sig.scores)<-sig.names
  r$sig.scores.raw<-r$sig.scores
  rand.flag<-is.null(r$rand.scores)|!all(is.element(names(gene.sign),colnames(r$rand.scores)))
  if(rand.flag){
    print("Computing also random scores.")
    r$rand.scores<-r$sig.scores
  }
  for (i in sig.names){
    b.sign<-is.element(r$genes,gene.sign[[i]])
    if(sum(b.sign)<2){next()}
    if(rand.flag){
      rand.scores<-get.semi.random.OE(r,r$genes.dist.q,b.sign,num.rounds = num.rounds)
    }else{
      rand.scores<-r$rand.scores[,i]
    }
    raw.scores<-colMeans(r$zscores[b.sign,])
    final.scores<-raw.scores-rand.scores
    r$sig.scores[,i]<-final.scores
    r$sig.scores.raw[,i]<-raw.scores
    r$rand.scores[,i]<-rand.scores
  }
  if(full.flag){return(r)}
  sig.scores<-r$sig.scores
  return(sig.scores)
}

data$tpm = data$TPM
data$genes = rownames(data$TPM)

geneSet = list()
load("D:/projects/R Projects/ImmuneResistance-master/Results/Signatures/cell.type.sig.full.RData")
OE = get.OE.bulk(data, gene.sign = cell.sig)
OE = as.data.frame(OE)
OE["Sample"] = colnames(data$TPM)
expMarker = merge(data$Samples,OE,by="Sample")
ggboxplot(expMarker,x="Resp_NoResp",y='stroma',color = "Resp_NoResp",
          palette = "npg",add = "jitter") + stat_compare_means()




############ EMT/Stroma_core signature: 8 genes were included in the EdgeSeq expression panel 
##Patients with high CD8 infiltration and low EMT/Stroma core
##gene expression had the highest response rates and longest PFS and OS
EMT_Stroma_core_signature = c("FLNA","EMP3","CALD1","FN1","FOXC2","LOX","FBN1","TNC")
expression = data$TPM[EMT_Stroma_core_signature,]
#names(expression)
expression = as.data.frame(t(expression))
expression["EMT_Stroma_core_signature"] = apply(expression, 1, function(x){return(mean(as.numeric(x)))})
expression["Sample"] = rownames(expression)
expMarker = merge(data$Samples,expression[c("Sample","EMT_Stroma_core_signature")],by="Sample")
ggboxplot(expMarker,x="Resp_NoResp",y="EMT_Stroma_core_signature",color = "Resp_NoResp",
          palette = "npg",add = "jitter") + stat_compare_means()



#############
##risk score, 11 IRGs (immune-related genes)
IRGs = c("LEPR","PRLHR","NR2F2","PRL","NRP1","ANGPTL5","IGF1","TNFRSF10B","TNFRSF10A","PLAU","IFI30")
coeff = c(0.32196,-0.64921,-0.32677,0.23573,0.39005,0.38166,-0.03522,0.02975,0.39830,0.14607,-0.68625)
expression = data$TPM[IRGs,]
riskScore = apply(expression, 2, function(x){return(sum(x*coeff))})
expMarker = merge(data$Samples,data.frame(Sample=names(riskScore),RiskScore=riskScore),by="Sample")
ggboxplot(expMarker,x="Resp_NoResp",y="RiskScore",color = "Resp_NoResp",
          palette = "npg",add = "jitter") + stat_compare_means(method.args = list(alternative = "less"))




############
#Tertiary lymphoid structures (TLS) gene signature score
TLS_gene = c("CD79B","CD1D","CCR6","LAT","SKAP1","CETP","EIF1AY","RBP5","PTGDS")
TLS_gene = intersect(TLS_gene,rownames(data$TPM))
expression = data$TPM[TLS_gene,]
mean_expression = apply(expression, 2, mean)
expMarker = merge(data$Samples,data.frame(Sample=names(mean_expression),TLS_score=as.numeric(mean_expression)),by="Sample")
ggboxplot(expMarker,x="Resp_NoResp",y='TLS_score',color = "Resp_NoResp",
          palette = "npg",add = "jitter") + stat_compare_means()



##########
##CXCL9
expression = data$TPM["CXCL9",]
expMarker = merge(data$Samples,data.frame(Sample=names(expression),CXCL9=as.numeric(expression)),by="Sample")
ggboxplot(expMarker,x="Resp_NoResp",y='CXCL9',color = "Resp_NoResp",
          palette = "npg",add = "jitter") + stat_compare_means()




##########
##HRH1
expression = data$TPM["HRH1",]
expMarker = merge(data$Samples,data.frame(Sample=names(expression),HRH1=as.numeric(expression)),by="Sample")
ggboxplot(expMarker,x="Resp_NoResp",y='HRH1',color = "Resp_NoResp",
          palette = "npg",add = "jitter") + stat_compare_means()








