### data.1
tmp = readLines("D:/BIG/work/test/29.txt")
genesets = list()
for(k in 1:length(tmp)){
	a = strsplit(tmp[k],split="\t")[[1]]
	genesets[[ a[1] ]] = a[-1]
}

### data.2
network = read.table("D:/BIG/work/19-scGWAS/PathwayCommons12.All.hgnc.exPCDHA.tsv")

### data.3
dat = read.table("D:/BIG/work/test/HiSeqV2", as.is=T, header=T)
expr = as.matrix(dat[,-1])
rownames(expr) = dat[,1]

sapply(colnames(expr), function(u){strsplit(u, split="\\.")[[1]][4]}) -> type
expr = expr[, type=="06"]

apply(expr, 1, sum) -> check
expr = expr[!is.na(check), ]
apply(expr, 1, sum) -> check
expr = expr[check!=0, ]

sapply(colnames(expr), function(u){strsplit(u, split="\\.")[[1]] -> v; paste(v[1:3], collapse="-")}) -> x1
colnames(expr) = x1

### data.4
anno = read.table("D:/BIG/work/test/TCGA-SKCM_final_cluster.txt", as.is=T)

### data.5
surv = read.table("D:/BIG/work/test/survival_SKCM_survival.txt", header=T, sep="\t")

########################################################################################

### gene set ssGSEA
gsva(expr, genesets,method="ssgsea") -> test
heatmap.2(test, trace="none") -> fit
cutree(as.hclust(fit$colDendrogram), 4) -> x

match(names(x), anno[,1]) -> idx
cbind(x, anno[idx,]) -> x1
table(x1[,1], x1[,3])

### survival
library("survival")
library("survminer")

surv.data = surv[match(colnames(expr), surv[,2]),]
Y1 = Surv(surv.data[,"OS.time"], surv.data[,"OS"])

dat = data.frame(cbind(surv.data,x) )
fit = survfit( Surv(surv.data[,"OS.time"], surv.data[,"OS"]) ~ x, data = dat)
coxph(Y1 ~ xvector)

ggsurvplot(fit, data=surv.data , risk.table = TRUE,pval = TRUE,ggtheme = theme_minimal())

########################################################################################

### edge set ssGSEA

which(network[,1] %in% rownames(expr) & network[,2] %in% rownames(expr)) -> ii
network = network[ii, ]

### prepare edge sets
edgesets = list()
for(k in 1:length(genesets)){
	genes = genesets[[k]]
	which(network[,1] %in% genes & network[,2] %in% genes ) -> ii
	paste("E", ii, sep="") -> ee
	edgesets[[ names(genesets)[k] ]] = ee
}

crossedges = c()
for(k1 in 1:(length(genesets)-1) ){
	genes_1 = genesets[[k1]]
	for(k2 in (k1+1):length(genesets)){
		genes_2 = genesets[[k2]]
		genes_3 = intersect(genes_1, genes_2)
		which(network[,1] %in% genes_3 & network[,2] %in% genes_3) -> ii
		crossedges = c(crossedges, ii)
	}
}
paste("E", crossedges, sep="") -> ee
edgesets[[ "crossedges" ]] = ee

for(k1 in 1:(length(genesets)-1) ){
	genes_1 = genesets[[k1]]
	for(k2 in (k1+1):length(genesets)){
		genes_2 = genesets[[k2]]
		genes_3 = intersect(genes_1, genes_2)
		
		genes_1 = setdiff(genes_1, genes_3)
		genes_2 = setdiff(genes_2, genes_3)
		
		which(network[,1] %in% genes_1 & network[,2] %in% genes_2) -> ii_1
		which(network[,1] %in% genes_2 & network[,2] %in% genes_1) -> ii_2
		ii = union(ii_1, ii_2)
		if(length(ii) >= 5){
			paste("E", ii, sep="") -> ee
			edgesets[[ paste("C_", k1, "_", k2, sep="") ]] = ee
		}
	}
}


### edge weights per sample
match(network[,1], rownames(expr)) -> idx1
match(network[,2], rownames(expr)) -> idx2

### code below, row 93 -- 98, takes time...

edge_weight = c()
for(k in 1:ncol(expr)){
	apply(cbind(expr[idx1, k], expr[idx2, k]), 1, mean) -> ee
	edge_weight = cbind(edge_weight, ee)
	cat(k,",",sep="")
}
rownames(edge_weight) = paste("E", 1:nrow(edge_weight), sep="")
colnames(edge_weight) = colnames(expr)[1:ncol(edge_weight)]

### edge set ssGSEA
gsva(edge_weight, edgesets, method="ssgsea") -> test2

heatmap.2(test2, trace="none") -> fit
cutree(as.hclust(fit$colDendrogram), 4) -> x

match(names(x), anno[,1]) -> idx
cbind(x, anno[idx,]) -> x1
table(x1[,1], x1[,3])


### survival analysis

surv.data = surv[match(gsub("\\.","-",colnames(expr)), surv[,1]),]

new3 = cbind(drug.response, surv.data)
Y1 = Surv(surv.data[,"OS.time"], surv.data[,"OS"])

dat = data.frame(cbind(surv.data,x) )
fit = survfit( Surv(surv.data[,"OS.time"], surv.data[,"OS"]) ~ x, data = dat)
coxph(Y1 ~ xvector)

ggsurvplot(fit, data=surv.data , risk.table = TRUE,pval = TRUE,ggtheme = theme_minimal())

