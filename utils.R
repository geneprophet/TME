#utils
write.gmt <- function(gs,file){
  sink(file)
  lapply(names(gs), function(i){
    cat(paste(c(i,'tmp',gs[[i]]),collapse = '\t'))
    cat('\n')
  })
  sink()
}
write.gmt(geneSet,'sink_example.txt')

write.table(t(HiSeqV2_PANCAN[,-1]),file = 'matrix.mtx',sep = "\t",col.names = F,row.names = F)
write.table(names(HiSeqV2_PANCAN),file = 'barcodes.tsv',sep = '\n',col.names = F,row.names = F,quote = F)
write.table(HiSeqV2_PANCAN[,1],file = 'gene.tsv',sep = '\n',col.names = F,row.names = F,quote = F)


######
library(Seurat)
library(loomR)
library(SeuratDisk)
load('G:/AD_seurat.Rdata')
N.loom <- as.loom(AD_seurat,filename="G:/AD.loom")
N.loom$close_all()


