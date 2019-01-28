setwd("~/my-papers-2017/HAO_new")
# read dependency table
chrs.uber<-read.csv("Dep_table_Jan2018.csv", header=T,  stringsAsFactors = F, na.strings = "")

# checking if dependecies are simple and correctly represented
chrs.uber[,c(1,3)]
paste(chrs.uber[,1], chrs.uber[,3]) %>% unique() ->tmp1
  strsplit(tmp1, " ") %>% lapply(., function(x) x[1]) %>% unlist() ->tmp
tmp1 %>% length()
tmp %>% unique() %>% length()
tmp %>% length()
tmp[duplicated(tmp)]

# quick check by grap
# construct matrix out of tmp1
strsplit(tmp1, " ") %>% lapply(., function(x) x[1]) %>% unlist() %>% paste("CHAR:",., sep="") -> col2
strsplit(tmp1, " ") %>% lapply(., function(x) x[2]) %>% unlist()%>% paste("CHAR:",., sep="") -> col1
chrs.uber=cbind(col1, col2)

# convert the table to matrix
#chrs.uber<-as.matrix(chrs.uber[,c(1,2)])

# make igraph object
g=graph_from_edgelist(chrs.uber)
# plot all dependecies
plot(g, vertex.color="green", vertex.size=2,
     vertex.label.dist=0.2, vertex.label.cex=0.5, vertex.label=NULL, edge.arrow.size=.2 )
components(g, "weak")
isomorphism_class(g)

# now let's select some components to get matrices for them
# there are a few ways to do it:

###########
# Assigning dependecies into groups
#######
# remove chars with two layers of dependecy 
# chars 199, 200, 201
#tb<-tb[!tb[,1] %in% c("199", "200", "201"),]
#
tb[,7]=read.csv("Dep_table_19Jan2018_no-199-200-201.csv", header=T,  stringsAsFactors = F, na.strings = "")[,2]
tb[,8]=read.csv("Dep_table_19Jan2018_no-199-200-201.csv", header=T,  stringsAsFactors = F, na.strings = "")[,4]

tb[,7]=read.csv("Dep_table_19Jan2018.csv", header=T,  stringsAsFactors = F, na.strings = "")[,2]
tb[,8]=read.csv("Dep_table_19Jan2018.csv", header=T,  stringsAsFactors = F, na.strings = "")[,4]

tb<-read.csv("Dep_table_19Jan2018_no-199-200-201.csv", header=T,  stringsAsFactors = F, na.strings = "")
tb<-read.csv("Dep_table_19Jan2018.csv", header=T,  stringsAsFactors = F, na.strings = "")
#tb<-read.csv("Dep_table_Jan2018.csv", header=T,  stringsAsFactors = F, na.strings = "")


paste(tb[,1], tb[,3], sep="<") -> dp.gr
dp.gr%>% unique()%>% match(dp.gr, .) ->dp.gr.num
tb=cbind(tb, dp.gr, dp.gr.num)
tb[,1]=paste0("CHAR:", tb[,1])
tb[,3]=paste0("CHAR:", tb[,3])

# check consistency of dependecies
coding_states_mt(char_matrix)->coding_states_matrix
ln=unique(dp.gr.num)

# check if controlling chars have only one state in depen.tb
lapply(ln, function(x) tb[tb[,6]==x,][,4]%>% unique() %>% length()) %>% unlist() ->contr.n.state
tb[tb[,6]==which(contr.n.state>1),] #show it

i=1
tb.check=c()
for (i in seq_along(ln)){
  tb[tb[,6]==i,]->tb.focal
  
  # check if states of controlled char in dep. table and in matrix are the same 
  set.controlled=setequal(coding_states_matrix[[tb.focal[1,1]]]%>%unname(), tb.focal[,2])
  
  # check if states of controlling char in dep. table are in matrix 
  con.in.mt=is.element(tb.focal[,4] %>% unique(),
  coding_states_matrix[[tb.focal[1,3]]]%>%unname() )
  
  tb.check=rbind(tb.check,
  cbind(ln[i], set.controlled, con.in.mt)
  )
  
}

tb.check[,2]
coding_states_matrix["CHAR:15"]
coding_states_matrix["CHAR:17"]

coding_states_matrix["CHAR:272"]

############################################################
# # requirements: 1. dependency table, 2. read matrix
# 
# 
# ed=list2edges(coding_states_matrix)
# ed[,2]=paste0("val:", ed[,2])
# ed=cbind(ed, paste0(ed[,1], ed[,2]))
# 
# # from matrix
# # chrs present in dependency table adn their full states from matrix
# ed.red=ed[ed[,1]%in%chrs[,1],]
# 
# # constructing nodes
# # chrs from dep tb
# chrs=c(tb[,1], tb[,3]) %>% unique()
# chrs=cbind(chrs, chrs, "present", NA)
# colnames(chrs)<-c("id", "label", "in.matrix", "in.dep.tb")
# abs.chr=chrs[,1][!chrs[,1]%in%ed[,1]] # charatcers absent in matrix
# chrs[!chrs[,1]%in%ed[,1],][,3]<-"absent"
# 
# # states
# sts=c( # from dpendency table
#   paste0("val:", tb[,2])%>% paste0(tb[,1], .),
#   paste0("val:", tb[,4])%>% paste0(tb[,3], .)
# ) %>% cbind(c(tb[,2], tb[,4]), .) 
# sts=sts[!duplicated(sts[,2]),]
# sts[,1]=paste0("val:", sts[,1])
# 
# # which states are absent in matrix
# #sts[,2][!sts[,2]%in%ed[,4]]
# sts[,2][!sts[,2]%in%ed[,3]]
# 
# # from matrix
# # chrs present in dependency table adn their full states from matrix
# ed.red=ed[ed[,1]%in%chrs[,1],]
# 
# all.st=rbind(ed.red[,c(4,2)], sts[,c(2,1)]) # all relevant states
# all.st[duplicated(all.st[,1]),] #all states unique
# all.st=cbind(all.st, "present", "present")
# colnames(all.st)<-c("id", "label", "in.matrix", "in.dep.tb")
# 
# # states present in matrix but absent in dep tb
# ed.red[!ed.red[,4]%in%sts[,2],][,4] %>% match(., all.st[,1])->v
# all.st[v,][,4]<-"absent"
#   
# # states present in dp table but absent in matrix
# sts[,2][!sts[,2]%in%ed.red[,4]] %>% match(., all.st[,1])->v
# all.st[v,][,3]<-"absent"
# 
# #combine chrs and states
# nodes=rbind(chrs, all.st)

#####################################
# constructing edges
###########################
#requirements: 1. dependency table, 2. read matrix

#makeng edge table and tailored dep table from dependecy table
tb[,2]=paste0("val:", tb[,2])
tb[,4]=paste0("val:", tb[,4])
dep.t=cbind(tb[,1:2], paste0(tb[,1], tb[,2]),
      tb[,3:4],
paste0(tb[,3], tb[,4])
)

colnames(dep.t)<-rep("a", 6)
dep.t.ed=rbind(dep.t[,1:3], dep.t[,4:6])
# selecting unique
dep.t.ed=dep.t.ed[!duplicated(dep.t.ed[,3]),]
colnames(dep.t.ed)<-rep("a", 3)
################
# getting edge table from matrix data
ed=list2edges(coding_states_matrix)
ed[,2]=paste0("val:", ed[,2])
ed=cbind(ed, paste0(ed[,1], ed[,2]))

# reducing edges to thos present in dependecy
ed.red=ed[ed[,1]%in%dep.t.ed[,1],]
colnames(ed.red)<-rep("a", 3)
# getting tb of all chrs links
all.chrs=rbind(ed.red, dep.t.ed, make.row.names = F)
all.chrs=all.chrs[!duplicated(all.chrs[,3]),]

# making edges
edges=rbind(
  data.frame(from=all.chrs[,1],
           to=all.chrs[,3],
           type="char"),
data.frame(from=dep.t[,6],
           to=dep.t[,3],
           type="state")
)
# making nodes
all.chrs.n=rbind(ed.red, dep.t.ed, make.row.names = F)
# which chrs from dep absent in matrix
dep.t.ed.u=dep.t.ed[!duplicated(dep.t.ed[,1]),]
ed.red.u=ed.red[!duplicated(ed.red[,1]),]
chr.absent=dep.t.ed.u[!(dep.t.ed.u[,1] %in% ed.red.u[,1]),]

# which states from dep absent in matrix
dep.t.ed.u.s=dep.t.ed[!duplicated(dep.t.ed[,3]),]
ed.red.u.s=ed.red[!duplicated(ed.red[,3]),]
st.abs.in.mt=dep.t.ed.u.s[!(dep.t.ed.u.s[,3] %in% ed.red.u.s[,3]),]
# which states from matrix absent in dep table
st.abs.in.dep=ed.red.u.s[!(ed.red.u.s[,3] %in% dep.t.ed.u.s[,3]),]

#id for characters
chr.u=c(dep.t.ed[,1], ed.red[,1]) %>% unique()
nodes=rbind(
       data.frame(id=chr.u, label=chr.u
),
# id for states
data.frame(id=all.chrs[,3], label=all.chrs[,2])
) 

types.n=c()
types.n[nodes[,1]%in% chr.u]<-"lightblue"

types=c()
types[nodes[,1] %in% chr.absent[,1]]<-"chr.a.in.mt" #"chr.a.in.mt"
types[nodes[,1] %in% st.abs.in.mt[,3]]<-"st.abs.in.mt"  #"st.abs.in.mt"
types[nodes[,1] %in% st.abs.in.dep[,3]]<-"st.abs.in.dep" #"st.abs.in.dep"

nodes[,3]<-types.n
colnames(nodes)<-c("id", "label",   "color")

library("plyr", lib.loc="~/.local/R/site-library")
edges1=edges
edges1[,3]=revalue(edges[,3], c("char"=NA, "state"="to"))
edges1[,4]=revalue(edges[,3], c("char"=TRUE, "state"=FALSE))
edges1[,4]=as.logical(edges1[,4])
edges1[,5]=as.numeric(as.character(revalue(edges[,3], c("char"=4, "state"=2))))
colnames(edges1)<-c("from", "to",   "arrows", "dashes", "width")
visNetwork(nodes, edges1, height = "500px", width = "100%")

visNetwork(nodes, edges1, height = "5000px", width = "5000px") %>%
  visIgraphLayout(layout = "layout_with_fr") %>%
  visNodes(size = 20, font="20px") 
 #visExport() ->x

layout = "layout_with_fr"
str(edges1)


######
# igrap
# removing chars 399 and 400
iedge<-edges1
iedge<-iedge[-grep("CHAR:399", iedge[,1]),]
iedge<-iedge[-grep("CHAR:399", iedge[,2]),]

iedge<-iedge[-grep("CHAR:400", iedge[,1]),]
iedge<-iedge[-grep("CHAR:400", iedge[,2]),]

inodes<-nodes
inodes<-inodes[-grep("CHAR:400", inodes[,1]),]
inodes<-inodes[-grep("CHAR:399", inodes[,1]),]

net <- graph_from_data_frame(d=iedges, vertices = inodes, directed=TRUE)

pdf(file = paste0("DepGraph_", layout, ".pdf"), width=7, height=7)
plot(net, vertex.label=NA, vertex.size=6, edge.width=2, edge.arrow.size=.4, vertex.frame.color="white", layout=layout.fruchterman.reingold,
     rescale=T )
dev.off()

# with vertex labels
pdf(file = "Dependecy_graph_names.pdf", width=17, height=17)
plot(net, vertex.label.cex=.8, vertex.label.color="black",  vertex.size=4, edge.width=2, edge.arrow.size=.4, vertex.frame.color="white", layout=layout.fruchterman.reingold,
     rescale=T )
dev.off()


tkplot(net, vertex.label=NA, vertex.size=5, edge.width=2, edge.arrow.size=.4, vertex.frame.color="white", layout=layout.fruchterman.reingold)
net

### layots
###
layouts <- grep("^layout\\.", ls("package:igraph"), value=TRUE)
# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama", layouts)]
par(mfrow=c(3,3))

for (layout in layouts) {
  print(layout)
  l <- do.call(layout, list(net))
  #plot(net, vertex.label=NA, vertex.size=.1, edge.width=.5) 
  png(filename = paste0("DepGraph_", layout, ".png"))
  plot(net, vertex.label=NA, vertex.size=5, edge.width=2, edge.arrow.size=.2, vertex.frame.color="white", layout=l)
  dev.off()
  }
dev.off()
##


###########
# color
V(net)$color
V(net)$color[1:83]<-"orange"
V(net)$color[84:268]<-"grey"
V(net)$color[84:268]<-"#525252"
#

# e color
E(net)$dashes[1:185]
E(net)$color<-"grey"
E(net)$color<-"#252525"
E(net)$color<-"#E3191B"
E(net)$color[1:185]<-"orange"

# arrows
E(net)$arrows[186:305]
E(net)$arrow.mode<-0
E(net)$arrow.mode[186:305]<-2

###########################################
# check in de tb: if states in dependent chars are unique
edges2list(cbind(tb[,1], tb[,2]))->tmp
lapply(tmp, function(x) length(unique(x))==length(x))%>%list2edges()

##

paste(char_matrix[["CHAR:227"]],
char_matrix[["CHAR:193"]], sep="") %>% unique()

paste(char_matrix[["CHAR:78"]],
      char_matrix[["CHAR:79"]], sep="") %>% unique()

paste(char_matrix[["CHAR:273"]],
      char_matrix[["CHAR:276"]], sep="") %>% unique() #good

paste(char_matrix[["CHAR:180"]],
      char_matrix[["CHAR:183"]], sep="") %>% unique()

paste(char_matrix[["CHAR:238"]],
      char_matrix[["CHAR:239"]], sep="") %>% unique()

paste(char_matrix[["CHAR:14"]],
      char_matrix[["CHAR:15"]], sep="") %>% unique()

paste(char_matrix[["CHAR:166"]],
      char_matrix[["CHAR:167"]], sep="") %>% unique()

paste(char_matrix[["CHAR:25"]],
      char_matrix[["CHAR:23"]], sep="") %>% unique()

paste(char_matrix[["CHAR:295"]],
      char_matrix[["CHAR:296"]], sep="") %>% unique()

paste(char_matrix[["CHAR:403"]],
      char_matrix[["CHAR:37"]], sep="") %>% unique()

# # make report
# # check if controlling chars have only one state in depen.tb
# lapply(ln, function(x) tb[tb[,6]==x,][,4]%>% unique() %>% length()) %>% unlist() ->contr.n.state
# tb[tb[,6]==which(contr.n.state>1),] #show it
# 
# 
# i=3
# #tb1=tb
# #tb=tb1
# unique_comb<-NA
# susp_dependencies<-NA
# susp_inappl<-NA
# report<-NA
# 
# tb=cbind(tb, unique_comb, susp_dependencies, susp_inappl, report)
# for (i in 1:max(tb$dp.gr.num)){
#   i.mt=tb[which(tb$dp.gr.num==i),]
#   
#   # unique coding combinations from matrix
#   coding=paste( char_matrix[[i.mt[1,3]]], 
#          char_matrix[[i.mt[1,1]]], sep=":") %>% unique()
#   tb[tb$dp.gr.num==i,][1,9]<-paste(coding, collapse = ", ") 
#   
#   
#   coding.l=strsplit(coding, ":")
#   
#   # filter out non-numeric
#   coding[
#   lapply(coding.l, function(x) all(!(x==("-")|x==("?"))) ) %>% unlist()
#   ]->coding.un
#   
#   # suspicious coding
#   coding.dep=paste(i.mt[,8], i.mt[,7], sep=":") %>% unique() 
#   coding.sus=coding.un[!(coding.un %in% coding.dep)]
#   coding.sus=paste(coding.sus, collapse=", ")
#   tb[tb$dp.gr.num==i,][1,10]<-coding.sus
#   
#   # identify weird characters with "-", "?:-", "controlling charactre:-"
#   coding[ lapply(coding.l, function(x) any(x[1]=="-")) %>% unlist()  ]->coding.inap
#   if (any(coding %in% "?:-")) coding.inap=c(coding.inap, "?:-")
#   
#   # if states which are not in dependency used with inapplicable
#   #coding[ lapply(coding.l, function(x) any(x[2]=="-")) %>% unlist()  ]->coding.inap.dep
#   unique(char_matrix[[i.mt[1,3]]])%>%as.character()->ques
#   ques=ques[!(ques=="-" | ques=="?")]
#   paste(ques[!(ques %in% i.mt[,8]%>% unique())], ":-", sep="")->ques 
#   ques%in%coding->coding.inap.dep
#   if (any(coding.inap.dep)) coding.inap=c(coding.inap, paste(ques[coding.inap.dep], collapse=", ") )
#   # if (any(coding %in% paste(i.mt[,7], ":-", sep = ""))){
#   #   ques=coding[coding %in% paste(i.mt[,7], ":-", sep = "")]
#   #   coding.inap=c(coding.inap, ques)}
#   #check if dependent and controlled have overlapping states
#   #if (any(i.mt[,7]%in%i.mt[,8]))  coding.inap=c(coding.inap, paste("overlap", "(", i.mt[,7][i.mt[,7]%in%i.mt[,8]], ")", sep="", collapse=", "))
#   
#   tb[tb$dp.gr.num==i,][1,11]<-paste(coding.inap, collapse = ", ") 
#   
#   #controlls report
#   tb[tb$dp.gr.num==i,][1,12]<-paste(i.mt[1,3], "{", paste(unique(i.mt[,8]), collapse=", "), "} > ", i.mt[1,1], "{", paste(unique(i.mt[,7]), collapse=", "), "}",
#         sep="")
# }
# 
# write.csv(file="dependecy_report.csv", tb)


##############################33
# New version

############
charC=char_matrix[[i.mt[1,3]]]
charD=char_matrix[[i.mt[1,1]]]

comb_char_matrix<-function(charC, charD){ # charC contorlling, charD dependent
  # unique coding combinations from matrix
  coding=paste( charC, 
                charD, sep=":") %>% unique()
  coding.l=strsplit(coding, ":")
  
  lapply(coding.l, function(x) x[1]) %>% unlist()->vec
  lapply(coding.l, function(x) x[2]) %>% unlist()->names(vec)

  return(list(coding, vec))
  
}
comb_char_matrix(charC, charD)
#########################


i=49
#tb1=tb
#tb=tb1
type1<-NA
type2<-NA
type2a<-NA
type3<-NA
type3a<-NA
type4<-NA
type4a<-NA
type5<-NA
report<-NA
unique_comb<-NA

i=49
tb=cbind(tb, type1, type2, type2a, type3, type4, type4a, type5, report, unique_comb, type3a)
for (i in 1:max(tb$dp.gr.num)){
#for (i in 1:35){
  i.mt=tb[which(tb$dp.gr.num==i),]
  
  # # unique coding combinations from matrix
  # coding=paste( char_matrix[[i.mt[1,3]]], 
  #               char_matrix[[i.mt[1,1]]], sep=":") %>% unique()
  
  vec=comb_char_matrix(char_matrix[[i.mt[1,3]]], char_matrix[[i.mt[1,1]]])
  comb2=paste(char_matrix[[i.mt[1,3]]], char_matrix[[i.mt[1,1]]], sep=":")
  # Type 1
  if (any(vec[[1]]=="-:-")) {tb[tb$dp.gr.num==i,][1,9]="-:-"
### !! recoding
id2inap<-comb2=="-:-"
char_matrix[[i.mt[1,3]]][id2inap]<-"?"
char_matrix[[i.mt[1,1]]][id2inap]<-"?"
####  
  }
  # Type 2
  if (any(vec[[1]]=="-:?")) tb[tb$dp.gr.num==i,][1,10]="-:?"
  
  
  # Type 2a
  if (any(vec[[1]]=="?:-")) {tb[tb$dp.gr.num==i,][1,11]="?:-"
  
  ### !! recoding
  id2inap<-comb2=="?:-"
  char_matrix[[i.mt[1,3]]][id2inap]<-"?"
  char_matrix[[i.mt[1,1]]][id2inap]<-"?"
  ####  

  }
  
  # Type 3
  contr.state=i.mt[,8] %>% unique()
  unique(vec[[2]])[
    !(unique(vec[[2]])==contr.state | unique(vec[[2]])=="-" | unique(vec[[2]])=="?")]->ncontr.state
  
  vec[[2]][names(vec[[2]])=="-"] %in% ncontr.state ->res
  if (any(res)){
    tb[tb$dp.gr.num==i,][1,12]<-paste(vec[[1]][names(vec[[2]])=="-"][res], collapse = ", ")
  }
  
  # Type 3a
  contr.state=i.mt[,8] %>% unique()
  vec[[1]] == paste0(contr.state, ":-") ->res
  if (any(res)){
    tb[tb$dp.gr.num==i,][1,"type3a"]<-vec[[1]][res]
    
    ### !! recoding
    id2inap<-comb2==vec[[1]][res]
    char_matrix[[i.mt[1,1]]][id2inap]<-"?"
    ####  
  }
  
  
  #Type 4
  vec[[2]][!(names(vec[[2]])=="-" | names(vec[[2]])=="?")]->t4
  #if (any(vec[[2]]=="-")) tb[tb$dp.gr.num==i,][1,13]="-:x"
  if (any(t4=="-")) tb[tb$dp.gr.num==i,][1,13]<-paste(paste(t4[t4=="-"], ":", names(t4[t4=="-"]), sep=""), collapse = ", " )
  
  #Type 4a
  vec[[2]][!(names(vec[[2]])=="-" | names(vec[[2]])=="?")]->t4a
  if (any(t4a=="?")) tb[tb$dp.gr.num==i,][1,14]<-paste(paste(t4a[t4a=="?"], ":", names(t4a[t4a=="?"]), sep=""), collapse = ", " )
  
  #Type 5
  paste(i.mt[,8], i.mt[,7], sep=":")->t5
  # filter out - and ? and check
  vec[[1]][
    !(grepl("-", vec[[1]]) | grepl("\\?", vec[[1]])) ] -> in.mt
  if (!all(in.mt %in% t5)){
    tb[tb$dp.gr.num==i,][1,15]<-paste(in.mt[!(in.mt %in% t5)], collapse = ", ")
  }
  
  #controlls report
  tb[tb$dp.gr.num==i,][1,]$report<-paste(i.mt[1,3], "{", paste(unique(i.mt[,8]), collapse=", "), "} > ", i.mt[1,1], "{", paste(unique(i.mt[,7]), collapse=", "), "}",
                                    sep="")
  # unique combs
  tb[tb$dp.gr.num==i,][1,]$unique_comb=paste(vec[[1]], collapse = ", ")
  
} #end loop
#############################################  
  
write.csv(file="dependecy_report_19Jan-final.csv", tb)

write.csv(file="matrix_19Jan-final.csv", char_matrix)
