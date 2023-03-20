library(ggplot2)
library(ggridges)
library(ggtree)
library(aplot)
library(seqinr)
library(reshape2)
library(treeio)



files<-dir(pattern =".fasta")
anot<- read.csv(file="tmp_table.txt", sep = "\t", header = T)

# create an empty data frame and populate it with GC and aa-lengths per CDS per genome
df<-NULL
for (q in 1:length(files)){
  f <- seqinr::read.fasta(file = files[q])
  g<-NA
  l<-NA
  s<-rep(gsub('-CDS.*fasta','', files[q]),length(f))
  r<-rep(grepl("-p",files[q]),length(f))
  for (i in 1:length(f)){
    g[i]<-GC(f[[i]])
    l[i]<-getLength(f[[i]])/3
  }
  dft<-data.frame(g,l,s,r)
  df = rbind(df, dft)
}

# ridgeline plot of CDS GC content distribution per replicon 
p2 <- ggplot(df, aes(x = g, y = s, fill = r, alpha = 0.7)) +
  geom_density_ridges() +
  theme_ridges() +
  scale_y_discrete(name="", labels=c()) +
  scale_x_continuous(limits=c(0.15,0.65), name="GC") +
  scale_fill_manual(name="replicon", values = c("#FF6E00", "#009999"))

# ridgeline plot of CDS length (aa) distribution per replicon
p3 <- ggplot(df, aes(x = l, y = s, fill = r, alpha = 0.7)) +
  geom_density_ridges(aes(vline_color = r), quantile_lines=TRUE, quantile_fun=function(x,...)median(x)) +
  scale_vline_color_hue(l = 20) +
  theme_ridges() +
  scale_y_discrete(name="", labels=c()) +
  scale_x_continuous(limits=c(0,1000), name="gene length (aa)") +
  scale_fill_manual(name="replicon", values = c("#FF6E00", "#009999")) +
  theme(legend.position = "none")

# barplots of genome size and fraction of interrupted genes per Arsenophonus strain
p4<-ggplot(anot, aes(x = size/1000000, y = strain)) +
    geom_col() +
    scale_y_discrete(name="", labels=c()) +
    scale_x_continuous(name="genome size (Mb)") +
    theme(legend.position="none", axis.title = element_text(color="black", size=16),
        axis.text = element_text(size = 12))

p5<-ggplot(anot, aes(x = pseudo, y = strain)) +
    geom_col() +
    scale_y_discrete(name="", labels=c()) +
    scale_x_continuous(name="fraction of interrupted genes") +
    theme(legend.position="none", axis.title = element_text(color="black", size=16),
        axis.text = element_text(size = 12))

# Heat map of plasmid and phage element abundance
anot2<-melt(anot, id.vars="strain", measure.vars=c("plasmids","phages"))
p6<-ggplot(anot2, aes(x=variable, y=strain)) + 
  geom_tile(aes(fill=value),color = "white", linewidth = 1.5) + scale_fill_viridis_c(option = "B", na.value = 'white') + 
  geom_text(aes(label = value), color = "grey", size = 4) +
  scale_y_discrete(name="", labels=c()) +
  theme_minimal(base_size=16)


###status and crispr
#anot3<-melt(anot, id.vars="strain", measure.vars=c("status","crispr"))
#p7<-ggplot(anot3, aes(x=variable, y=strain)) + 
#  geom_tile(aes(fill=value),color = "white", linewidth = 1.5) + scale_fill_viridis_d(option = "B", na.value = 'white') + 
#  scale_y_discrete(name="", labels=c()) +
#  theme_minimal()

# Core phyloTree with support values
tree <- read.newick(file = "tree.nwk", node.label = "label")
root <- rootnode(tree)  
p1<-ggtree(tree) + 
  geom_tiplab(align=TRUE, color='black') + xlim(0,0.4) + 
  geom_point2(aes(subset=!isTip & node != root,fill=cut(as.numeric(label), c(0, 70, 95, 100))), 
              shape=21, size=2.5) + 
  theme_tree(legend.position=c(0.1, 0.9)) + 
  scale_fill_manual(values=c("black", "grey", "white"), guide='legend', 
                    name='bootstrap support', 
                    breaks=c('(95,100]', '(70,95]', '(0,70]'), 
                    labels=expression("95-100", "70-95", "< 70")) +
  geom_treescale(x=0, y=10)


# make combined plot
p4 %>% insert_left(p1, width = 4) %>% insert_right(p5) %>% insert_right(p3) %>% insert_right(p2) %>% insert_right(p6,width = 0.4)

# we can rename the tree tip labels using the manes from the anot table 
#library(dplyr)
#n<-anot %>% select(strain, name)
#tree$tip.label<-n[[2]][match(tree$tip.label, n[[1]])]
#plot(tree)
