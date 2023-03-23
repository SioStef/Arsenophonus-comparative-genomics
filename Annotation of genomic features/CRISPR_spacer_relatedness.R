library(igraph)
library(reshape2)

# read node data
nodes <-read.csv("nodes.tsv", sep = "\t", header=T, as.is=T)

# prepare link data from blast output
b<-read.table(file="blast_out.tsv", header = T, sep = "\t")
b2<-subset(b, idn >= 98 & cov >= 97, select=c(Q,S))
b3<-table(b2)
b3<-matrix(b3, ncol=ncol(b3), dimnames=dimnames(b3))

M<-melt(b3)

my_fun <- function(x,y,z) {                    # Compute relatedness score
  z/min(nodes$spacers[nodes$strain==x],nodes$spacers[nodes$strain==y])
}

M$value<-mapply(my_fun, M$Q, M$S, M$value)

M<-acast(M,Q~S)

links<- setNames(M, c('from', 'to', 'weight'))

# network
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)
net<-simplify(net,remove.multiple = T, remove.loops = T, edge.attr.comb ="max")
net<-delete.edges(net, E(net)[weight==0])

# node colors based on Cas type
colrs <- c("darkcyan", "tomato", "gold")
V(net)$color <- colrs[V(net)$type]

# edge size based on weight
E(net)$width <- E(net)$weight*5

# plot network
l <- layout_nicely(net)
plot(net, edge.arrow.size=.1, layout=l, vertex.label.dist=2, vertex.size=10)
#plot(net, edge.arrow.size=.1, layout=l, vertex.label.dist=2)
#plot(net, edge.arrow.size=.1, layout=l, vertex.label.dist=2, vertex.size=35, rescale = FALSE, ylim=c(-3.5,4.5),xlim=c(2,2))

# add legend
legend(x=-1.5, y=-0.7, c("typeI-F","typeI-Fv", "typeI-E"), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)
