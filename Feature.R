library(readr)
TCGA_SKCM_varscan2_snv <- read_delim("D:/博士课题/TME/TCGA.SKCM/skcm/TCGA-SKCM.varscan2_snv.tsv", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)
#how many samples
length(unique(TCGA_SKCM_varscan2_snv$Sample_ID))

rownames(expressionMatrix)[which(rownames(expressionMatrix) %in%  signatures_genes)]

TCGA_EB_A3Y6 = TCGA_SKCM_varscan2_snv[which(TCGA_SKCM_varscan2_snv$Sample_ID=="TCGA-EB-A3Y6-01A"),]

mutation = vector(mode = "numeric",length = 0)
for (i in signatures_genes) {
  a = TCGA_EB_A3Y6[which(TCGA_EB_A3Y6$gene==i)[1],]
  if(is.na(a)){
    mutation=c(mutation,0)
  }else if(a$effect=="synonymous_variant"){
    mutation=c(mutation,0)
  }else if(a$effect=="missense_variant"){
    mutation=c(mutation,1)
  }else if(a$effect=="stop_gained"){
    mutation=c(mutation,2)
  }else {
    mutation=c(mutation,0)
  }
}

expression = vector(mode = "numeric",length = 0)
sample_e = HiSeqV2_PANCAN$`TCGA-EB-A3Y6-01`
for (i in signatures_genes) {
  if(i %in% HiSeqV2_PANCAN$Gene){
    expression = c(expression,sample_e[which(HiSeqV2_PANCAN$Gene==i)])
  }else{
    expression = c(expression,0)
  }
}

features = data.frame(gene=signatures_genes,mutation=mutation,expression=expression)

write.csv(features,file = 'features.txt',row.names = F,col.names = F)


all_sample_exp = SKCM_htseq_counts[,-1]
write.csv(all_sample_exp,file = 'all_sample_exp.txt',row.names = F,col.names = F)

#################################################
#Expression counts
library(readr)
TCGA_SKCM_htseq_counts <- read_delim("D:/博士课题/TME/TCGA.SKCM/skcm/TCGA-SKCM.htseq_counts.tsv", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)

library('org.Hs.eg.db')
keytypes(org.Hs.eg.db)
genes = unlist(lapply(strsplit(TCGA_SKCM_htseq_counts$Ensembl_ID,split = '.',fixed = T),function(x){return(x[1])}))
TCGA_SKCM_htseq_counts$Ensembl_ID = genes

genes2 = select(org.Hs.eg.db, keys=signatures_genes, columns="ENSEMBL", keytype="SYMBOL")
gene3 = genes2$ENSEMBL[which(!duplicated(genes2$SYMBOL))]
SKCM_htseq_counts = TCGA_SKCM_htseq_counts[which(TCGA_SKCM_htseq_counts$Ensembl_ID %in% gene3),]
rownames(SKCM_htseq_counts) = SKCM_htseq_counts$Ensembl_ID
SKCM_htseq_counts=SKCM_htseq_counts[gene3,]
SKCM_htseq_counts$Ensembl_ID=signatures_genes
write.table(SKCM_htseq_counts,file = 'expression_count.txt',row.names = F,sep = '\t')

