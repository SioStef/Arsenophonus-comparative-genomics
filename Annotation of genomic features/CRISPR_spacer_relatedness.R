library(igraph)
library(reshape2)

# prepare data
M <- read.table(file="relatedness_matrix.tsv", header = T, sep = "\t")
M <- M[-c(1)]
rownames(M)<-colnames(M)

M = as.matrix(M)

links<- setNames(melt(M), c('from', 'to', 'weight'))
nodes <-read.csv("nodes.tsv", sep = "\t", header=T, as.is=T)

net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)
net<-simplify(net,remove.multiple = T, remove.loops = T, edge.attr.comb ="max")
net<-delete.edges(net, E(net)[weight=0])

# colour nodes based on Cas type
colrs <- c("darkcyan", "tomato", "gold")
V(net)$color <- colrs[V(net)$type]

# scale node size based on number of spacers
#V(net)$size <- V(net)$spacers

# scale edge size based on weight
E(net)$width <- E(net)$weight*5

# plot
l <- layout_nicely(net)
plot(net, edge.arrow.size=.1, layout=l, vertex.label.dist=2, vertex.size=10)

legend(x=-1.5, y=-0.7, c("typeI-F","typeI-Fv", "typeI-E"), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)
