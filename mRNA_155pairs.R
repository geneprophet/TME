library(readr)
mRNA_count_155pair <- read_csv("mRNA_count_155pair.csv")
#View(mRNA_count_155pair)

genes = unlist(lapply(strsplit(mRNA_count_155pair$EnsemblGene_GeneSymbol,split = '_'),function(x) {return(x[2])}))
mRNA_count_155pair$EnsemblGene_GeneSymbol = genes
write.csv(mRNA_count_155pair[,1:20],file = "mRNA_count_155pair_symbol.csv",row.names = F)
###get protein coding genes from ensembl
library(biomaRt)
listMarts(host="uswest.ensembl.org")
ensembl_us_west = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org")
head(listDatasets(ensembl_us_west))
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
head(listFilters(ensembl))
head(listAttributes(ensembl))
pcgs <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','gene_biotype','chromosome_name','start_position','end_position'), filters =
                      'biotype', values ="protein_coding", mart = ensembl)

#pcgs[which(pcgs$hgnc_symbol %in% signatures_genes),]
mRNA_count_155pair = mRNA_count_155pair[which(mRNA_count_155pair$gene %in% pcgs$hgnc_symbol),]
mRNA_count_155pair[1:5,1:5]
mRNA_count_155pair$EnsemblGene_GeneSymbol=mRNA_count_155pair$gene
mRNA_count_155pair = mRNA_count_155pair[,-312]
row_name = mRNA_count_155pair$EnsemblGene_GeneSymbol
mRNA_count_155pair = mRNA_count_155pair[,-1]
mRNA_count_155pair = log2(mRNA_count_155pair+1)
mRNA_count_155pair$genes = row_name
mRNA_count_155pair = mRNA_count_155pair[,c(311,1:310)]

write.table(mRNA_count_155pair,file = "mRNA_count_155pair.tsv",row.names = F,sep = '\t')

######################################

library(survminer) # 加载包
library(survival) # 加载包
data(lung)
View(lung)
attach(lung) # 绑定数据集
Surv(time,status) # 创建生存对象
fit <- survfit(Surv(time,status) ~ sex,  # 创建生存对象 
               data = lung) # 数据集来源
fit # 查看拟合曲线信息
summary(fit)
ggsurvplot(fit, data = lung)
#ggsurvplot(fit, data = lung,conf.int = T, pval=T,surv.median.line = "hv",add.all=T)
library(readxl)
clnical <- read_excel("clnical.xlsx")
library(readr)
final_clusters <- read_delim("D:/projects/PycharmProjects/MFP/final_clusters.tsv", 
                             delim = "\t", escape_double = FALSE, 
                             col_names = FALSE, trim_ws = TRUE)
final_clusters = final_clusters[endsWith(final_clusters$X1,suffix="T"),]
final_clusters$X1 = substr(final_clusters$X1,7,9)
library(stringr)
clnical$sampleid = str_pad(as.character(clnical$`sample ID           （BDESCC2-)`),3,side = "left","0")
rownames(clnical) = clnical$sampleid
survival155 = clnical[final_clusters$X1,]
allsurvival = cbind(final_clusters,survival155)
allsurvival$TME = allsurvival$X2
allsurvival$`overall survival        time` = as.numeric(allsurvival$`overall survival        time`)
allsurvival$`overall survival              （1, dead; 0,alive)`[which(allsurvival$`overall survival              （1, dead; 0,alive)` == "1")] = 2
allsurvival$`overall survival              （1, dead; 0,alive)`[which(allsurvival$`overall survival              （1, dead; 0,alive)` == "loss")] = 1
allsurvival$`overall survival              （1, dead; 0,alive)` = as.numeric(allsurvival$`overall survival              （1, dead; 0,alive)`)
attach(allsurvival) # 绑定数据集
Surv(`overall survival        time`,`overall survival              （1, dead; 0,alive)`) # 创建生存对象
fit <- survfit(Surv(`overall survival        time`,`overall survival              （1, dead; 0,alive)`) ~ TME,  # 创建生存对象 
               data = allsurvival) # 数据集来源
fit # 查看拟合曲线信息
summary(fit)
ggsurvplot(fit, data = allsurvival)
