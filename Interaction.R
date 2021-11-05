library(readr)
PathwayCommons11_All_hgnc <- read_delim("skcm/PathwayCommons12.All.hgnc.sif", 
                                        delim = "\t", escape_double = FALSE, 
                                        col_names = FALSE, trim_ws = TRUE)
# PathwayCommons11_All_hgnc <- read_delim("pathwaycommons/PathwayCommons11.All.hgnc.txt", 
#                                         delim = "\t", escape_double = FALSE, 
#                                         col_names = T, trim_ws = TRUE)
PathwayCommons11_All_hgnc[1:3,1:3]
signatures_genes = unique(X29Fges$Gene)
signatures_genes[1:4]

#TODO:extract singnature interactions
length(which(PathwayCommons11_All_hgnc$X1 %in% signatures_genes))
length(which(PathwayCommons11_All_hgnc$X3 %in% signatures_genes))

interaction = PathwayCommons11_All_hgnc[which(PathwayCommons11_All_hgnc$X1 %in% signatures_genes & PathwayCommons11_All_hgnc$X3 %in% signatures_genes),]

interaction2 =  interaction[which(!duplicated.data.frame(interaction)),]

table(interaction2$X2)
library(igraph)

nodes = data.frame(name=signatures_genes)
edges = data.frame(from=interaction2$X1,to=interaction2$X3,type=interaction2$X2)

#undirected graph
interaction_graph = graph_from_data_frame(edges,directed = F,vertices = nodes)
plot(interaction_graph)
cluster_louvain(interaction_graph) -> cluster_results
l <- layout_randomly(interaction_graph,dim=2) 
plot(interaction_graph, layout=l)

adj = as_adjacency_matrix(interaction_graph,names = T)
adj = as.matrix(adj)
dim(adj)
adj[1:5,1:5]
adj[adj>1] = 1
adj[1:5,1:5]

signatures_genes[1:5]

write.csv(adj,file = "adj.txt",row.names = T)
write.csv(adj,file = "adj_vgae.txt",row.names = F)





