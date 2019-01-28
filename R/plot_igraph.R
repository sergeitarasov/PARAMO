ontology<-ONT
# processing ontology to incorporate character statements
#ontology<-onto_process(ONT, H$char, do.annot = F)
ontology$name_characters<-H$char
names(ontology$name_characters)<-H$char_id
# embedding manual annotations
ontology$annot_characters<-ONT$terms_selected_id

H$char
ONT$terms_selected_id

# exporting
cyto<-export_cytoscape(ontology, annotations = ontology$annot_characters, is_a = c("is_a"),
                       part_of = c("BFO:0000050"))
write.csv(cyto, file="HAO_cyto.csv")


########
cyto<-as_tibble(cyto)
net <- graph.data.frame(select(cyto, id_from, id_to), c(cyto$id_from, cyto$id_to) %>% unique(), directed=FALSE)
plot(net, vertex.label=NA, vertex.size=0.2, edge.width=.5, edge.curved=.4,
     vertex.color="orange", vertex.frame.color="black"
     )
str(V(net))

V(net)$color<-"black"
V(net)$color[1:394]<-"orange"
plot(net, vertex.label=NA, vertex.size=0.2, edge.width=.5, edge.curved=.4,
     vertex.color="orange", vertex.frame.color="black"
)
V(net)$frame.color<-V(net)$color
V(net)$frame.color<-"black"
V(net)$frame.color[1:394]<-"white"

V(net)$size<-.5
V(net)$size[1:394]<-4
plot(net, vertex.label=NA, edge.width=1, edge.color=rgb(0, 0, 0, alpha=.1), edge.curved=0, layout=layout.kamada.kawai(net) )
plot(net, vertex.label=NA, vertex.size=2, edge.width=.1, edge.curved=.6, layout=layout.kamada.kawai(net) )

plot(net, vertex.label=NA, vertex.size=.1, edge.width=.5, edge.curved=.4, layout=l)
l<-layout.forceatlas2(net, iterations=100, plotstep=500)
layout.spring
V(net)$color<-rep("orange", length(V(net)$name))[V(net)$name %>% grepl("CHAR",.)]

##
cg1<-net
# construct the nodes and edges data for gexf conversion
nodes <- data.frame(cbind(V(cg1), as.character(V(cg1))))
edges <- t(Vectorize(get.edge, vectorize.args='id')(cg1, 1:ecount(cg1)))

# do the conversion
write.gexf(nodes, edges, output = "gephi")   

write.gexf(nodes, edges)
###
layouts <- grep("^layout\\.", ls("package:igraph"), value=TRUE)
# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama", layouts)]
par(mfrow=c(3,3))

for (layout in layouts) {
  print(layout)
  l <- do.call(layout, list(net))
  plot(net, vertex.label=NA, vertex.size=.1, edge.width=.5) }
dev.off()
##
rexp(1000, 2) %>% log %>% hist
