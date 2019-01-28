library("magrittr")
library("igraph")

##############################################################################################
# FUNCTIONS
##############################################################################################

#' @title Combining two matrices
#' @description Combining two matrices. The parametric schem of matrice is defined by nattural
#' numbers; different numbers = different rate parameters
#' @param M1 matrix; if dependency true thenM1 controls M2
#' @param M2 matrix; if dependency true then: M2 depends on those states of M1 specified in dependent.state
#' @param dependent.state state(s) of M1 that switches on matrix M2 
#' @param name.sep separator for state names
#' @param diag.as hpopulate main diagonal with
#' @return Matrix
#' @examples
#' M1<-matrix(c(-1,1,  2,-2),2,2,byrow=TRUE)
#' rownames(M1)<-colnames(M1)<-c("0","1")
#' M2<-matrix(c(-3,3,  4,-4),2,2,byrow=TRUE)
#' rownames(M2)<-colnames(M2)<-c("0","1")
#' comb2matrices(M1, M2, dependent.state=NULL)
#' comb2matrices(M1, M2, dependent.state=2)
# if dependency true then: M2 depends on M1 states specified in dependent.state

M1=init_matrices(g1, char_matrix, diag.as=0)[[1]]
M2=init_matrices(g1, char_matrix, diag.as=0)[[2]]

comb2matrices(M1, M2, dependent.state=2)

comb2matrices<-function(M1,M2, dependent.state=NULL, name.sep="", diag.as=""){
  
  if (!is.null(dependent.state)){
    matrix.diag<-rep(0, ncol(M1))
    matrix.diag[dependent.state]<-1
    matrix.diag<-diag(matrix.diag)
  }
  
  if (is.null(dependent.state)){
    matrix.diag<-diag(nrow(M1))
  }
  
  M_kr=(M1%x%diag(nrow(M2)))+ (matrix.diag%x%M2)
  
  
  #getting colnames
  
  col=paste(colnames(kronecker(M1, diag(nrow(M2)), make.dimnames = T)),
            colnames(kronecker(diag(nrow(M1)), M2, make.dimnames = T)), sep="")
  col=gsub("::", name.sep, col, fixed=T)
  
  # merging two names
  rownames(M_kr)<-colnames(M_kr)<-col
  if (diag.as!="") diag(M_kr)<-diag.as
  
  return(M_kr)
}



#' @title Combining multiple matrices
#' @description Liely as independently evolving
#' @param list.matrices Lsit of matrices
#' @return Matrix
#' @examples
combNmatrices<-function(list.matrices,  ...){
  comb.matrix<-list.matrices[[1]]
  
  for (i in 1:(length(list.matrices)-1)){
    comb.matrix=comb2matrices(comb.matrix, list.matrices[[i+1]], dependent.state=NULL, diag.as = 0)
  }
  
  return(comb.matrix)
}

#' @title Initialize binary matrices given graph
#' @description Call matrices are populated with different parameters
#' @param graph igraph object
#' @return List of matrices
#' @examples
#' init_binary_matrix(g)
init_binary_matrix<-function(graph){
  matrix.list=list()
  n.matrix=vcount(graph)
  vertices=V(graph)$name
  
  param.scheme<-matrix(seq(1:(2* n.matrix )), ncol = 2, byrow=TRUE)
  for (i in 1:n.matrix){
    matrix.list[[vertices[i]]]<-matrix(c(-param.scheme[i,1],param.scheme[i,1], 
                                         param.scheme[i,2],-param.scheme[i,2]),2,2,byrow=TRUE)
    rownames(matrix.list[[vertices[i]]])<-colnames(matrix.list[[vertices[i]]])<-c("0","1")
  }
  
  return(matrix.list)
}

############### intialize any Q-matirx for a set of chars
char.matrix<-char_matrix
graph=subgraph.list[[2]]
graph=g1
graph=g[1:3]
class(graph)
graph=c("CHAR:7", "CHAR:8")

init_matrices(graph, char_matrix.rem, diag.as=0)

init_matrices<-function(graph, char_matrix, diag.as=0){
  
  matrix.list=list()
  
  if (class(graph)=="character"){
    n.matrix=length(graph)
    vertices=graph
  }
  
  if (class(graph)=="igraph"){
  
  n.matrix=vcount(graph)
  vertices=V(graph)$name
  }
  
  
  tmp=apply(char_matrix[,vertices], 2, list)
  lapply(tmp, function(x) unique(x[[1]])) %>%
  #apply(char_matrix[,vertices], 2, unique) %>% 
    lapply(., function(x) as.numeric(x) %>% na.omit() %>% as.character() %>% as.numeric(x) %>% sort()) -> l.state
  
  # total number different rate cats
  #lapply(l.state, function(x) (length(x)^2)-length(x)) %>% unlist() %>% sum() %>% seq(1, ., 1)->total.state
  
  #unlist(l.state) %>% length %>% seq(1, ., 1)->total.state
  
  l.param<-l.state
  
  k=1
  
  for (i in 1:length(l.state)){
    l.param[[i]]<-seq(k, k-1+length(l.param[[i]])^2-length(l.param[[i]]), 1)
    #k=l.param[[i]][length(l.param[[i]])]
    k=k+length(l.param[[i]])
  }
  
  # for (i in 1:length(l.state)){
  #   l.param[[i]]<-total.state[k:(k+length(l.param[[i]])-1)]
  #   k=k+length(l.param[[i]])
  # }

  matrix.list<-l.state
  for (i in 1:length(l.state)){
    matrix.list[[i]]<-init_char_matrix(char.state=l.state[[i]], rate.param=l.param[[i]], diag.as=diag.as)
  }
  return(matrix.list)
}
####################

## get matrix for a character
char.state<-c(1,2,3)
rate.param<-c(1:6) # rate parameters = (n.states^2)-n
init_char_matrix(char.state, rate.param, diag.as=0)

init_char_matrix<-function(char.state, rate.param, diag.as=NA){
  n.state<-length(char.state)
  Q=matrix(ncol = n.state, nrow=n.state,byrow=TRUE)
  Q[xor(lower.tri(Q, diag = FALSE), upper.tri(Q, diag = FALSE))]<-rate.param
  Q<-t(Q)
  diag(Q)<-diag.as
  rownames(Q)<-colnames(Q)<-as.character(char.state)
return(Q)
}

#' @title Get all dependency matrices given a dependecy graph
#' @description Construct dependency matrices and their correponding attributes
#' @param graph igraph object of ontology terms
#' @return List of matrices and their attributes
#' @examples
#' get_graph_matrix(g)

########## Structure of the List
$binary.matrices # intial binary matrices assigned to each node of graph

$comb.matrices$matrix # combined matrix for each node
$comb.matrices$state.string # vector of states [1] "00" "01" "10" "11"
$comb.matrices$state.ident # specifies the order of ontology terms in each state [1] "UBERON:0007829" "UBERON:2000663"
$comb.matrices$state.observable # ids and names of "observable" states. In red-blue tail notation refers to blue and red states
$comb.matrices$state.hidden # ids and names of "hidden" states. In red-blue tail notation refers to "blue absent"
# and "red absent"

$nodes.sorted # topologically sorted nodes
$vertex.hier # hierrachy of the nodes
###########
dep.tb<-tb
g<-g1

M1=init_matrices(g1, char_matrix, diag.as=0)[[1]]
M2=init_matrices(g1, char_matrix, diag.as=0)[[2]]
M1>M2
comb2matrices(M1, M2, dependent.state=2)
get_graph_matrix_any(graph, char_matrix, tb)

i=5
get_graph_matrix_any(graph, char_matrix, dep.tb, new.polymorph="&")

get_graph_matrix_any<-function(graph, char_matrix, dep.tb, new.polymorph=" ", Lnew.polymorph=NULL, Rnew.polymorph=NULL){

  g=graph  
 
  
  # dependent chars object
  complex.char<-list()
  complex.char$binary.matrices<- init_matrices(g, char_matrix, diag.as = 0)
  complex.char$comb.matrices<-list()
  
  # traverse graph
  topo=topo_sort(g, mode = c("out")) %>% names()
  complex.char$nodes.sorted=topo
  vertex.hier=ego(g, order=1, nodes = topo, mode = c("in"), mindist = 1)
  names(vertex.hier)<-topo
  complex.char$vertex.hier=vertex.hier
  
  i=1
  for (i in seq_along(topo)){
    focal.v=complex.char$vertex.hier[[i]]
    
    if (length(focal.v)==0){
      complex.char$comb.matrices[[topo[i]]]$matrix=complex.char$binary.matrices[[topo[i]]]
      complex.char$comb.matrices[[topo[i]]]$state.string=complex.char$binary.matrices[[topo[i]]] %>% row.names()
      complex.char$comb.matrices[[topo[i]]]$state.ident=topo[i]
      #complex.char$comb.matrices[[topo[i]]]$dependency.true=2
      #complex.char$comb.matrices[[topo[i]]]$state.observable=integer(0)
      #complex.char$comb.matrices[[topo[i]]]$state.hidden=integer(0)
      
      complex.char$comb.matrices[[topo[i]]]$matrix.new=complex.char$binary.matrices[[topo[i]]]
      
      # create new char coding
      entire.char=paste(char_matrix[[topo[i]]])
      new.mtchar<-entire.char
      names(new.mtchar)<-rownames(char_matrix)
      ### sort species to make the 1st one havinh non-polymorphic state
      #sort=grepl(new.polymorph, new.mtchar)
      #new.mtchar<-c(new.mtchar[!sort], new.mtchar[sort])
      ##
      complex.char$comb.matrices[[topo[i]]]$new.coding<-new.mtchar
      
    }
    
    if (length(focal.v)>0){
      contr.st.val<-dep.tb[which(dep.tb[,1]==topo[i])[1], 8]
      dep.st.val<-dep.tb[which(dep.tb[,1]==topo[i]), 7]
      
      if (length(focal.v)==1){ # if length =1 the dependency is chain like
        MC=complex.char$comb.matrices[[names(focal.v)]]$matrix
        M=complex.char$binary.matrices[[topo[i]]]
        #which(colnames(MC)==contr.st.val)
        dps=which(colnames(MC)==contr.st.val) # dependent state id
        cmb=comb2matrices(MC, M, dependent.state=dps, diag.as = 0)
      }
      
      
      # if (length(focal.v)>1){
      #   # sequentially combine multiple matrices as independently coevolving
      #   list.matrices=lapply(names(focal.v), function(x) complex.char$comb.matrices[[x]]$matrix)
      #   #names(list.matrices)=names(focal.v)
      #   #list.matrices=list.matrices[1:2]
      #   comb.mt=combNmatrices(list.matrices)
      #   
      #   # combine matrices  from above with the focal node matrix;
      #   # the state bearing dependency is where all entities=1, i.e. the last state
      #   M=complex.char$binary.matrices[[topo[i]]]
      #   dps=ncol(comb.mt)
      #   cmb=comb2matrices(comb.mt, M, dependent.state=dps, diag.as = 0)
      # }
      
   
      # adding attributes
      complex.char$comb.matrices[[topo[i]]]$matrix=cmb
      complex.char$comb.matrices[[topo[i]]]$state.string<-r.name<-cmb %>% row.names()
      
      # classify Qmatrix states in d and n
      name.splt=strsplit(row.names(cmb), "")
      Qmt.state.class<-list()
      #lapply(name.splt, function(x) x[dps]==contr.st.val) %>%unlist %>% as.character() -> Qmt.state.class$cat
      lapply(name.splt, function(x) x[1]==contr.st.val) %>%unlist %>% as.character() -> Qmt.state.class$cat
      
      Qmt.state.class$cat<-revalue(Qmt.state.class$cat, c("TRUE"="d", "FALSE"="n"))
      names(Qmt.state.class$cat)<- row.names(cmb)
      # creating new states
      #Qmt.state.class$new.state<-c(1:length(Qmt.state.class$cat))
      Qmt.state.class$new.state<-c(0:(length(Qmt.state.class$cat)-1))
      names(Qmt.state.class$new.state)<-row.names(cmb)
      complex.char$comb.matrices[[topo[i]]]$Qmt.state.classes<-Qmt.state.class
      
      # getting char matrix classes
      mt.state.class<-classify_state(char_matrix[[names(focal.v)]], char_matrix[[topo[i]]], contr.st.val, dep.st.val)
      complex.char$comb.matrices[[topo[i]]]$mt.state.class<-mt.state.class
      
      # creating remapping map
      new.char<-c()
      # recoding class n
      new.char[mt.state.class$cats=="n"]<-
        #paste(Qmt.state.class$new.state[Qmt.state.class$cat=="n"], collapse=new.polymorph)
        paste(Qmt.state.class$new.state[Qmt.state.class$cat=="n"], collapse=new.polymorph) %>%
        paste(Lnew.polymorph, ., Rnew.polymorph, sep="")
      
      # recoding class d (not dp, dq)
      new.st<-match(mt.state.class$comb.merged[mt.state.class$cats=="d"],
            names(Qmt.state.class$cat))
      #new.char[mt.state.class$cats=="d"]<-new.st
      new.char[mt.state.class$cats=="d"]<-Qmt.state.class$new.state[new.st]
      
      # recoding class dp, !!! this recoding work only for dp
      new.st<-match(mt.state.class$polymorph,
                    names(Qmt.state.class$cat))
      #new.char[mt.state.class$cats=="dp"]<-paste(new.st, collapse=new.polymorph)
      new.char[mt.state.class$cats=="dp"]<-
        paste(new.st, collapse=new.polymorph) %>% paste(Lnew.polymorph, ., Rnew.polymorph, sep="")
      
      # recoding class dq
       #new.char[mt.state.class$cats=="dq"]<-
       #paste(Qmt.state.class$new.state[Qmt.state.class$cat=="d"], collapse=new.polymorph)
      new.char[mt.state.class$cats=="dq"]<-
      paste(Qmt.state.class$new.state[Qmt.state.class$cat=="d"], sep="", collapse=new.polymorph) %>%
        paste(Lnew.polymorph, ., Rnew.polymorph, sep="")
      
      # recoding class q
      new.char[mt.state.class$cats=="q"]<-"?"
      
      names(new.char)<-mt.state.class$comb.sep
      complex.char$comb.matrices[[topo[i]]]$char.mapping<-new.char
      
      # create new char coding
      entire.char=paste(char_matrix[[names(focal.v)]], char_matrix[[topo[i]]], sep = ":")
      new.mtchar<-mapvalues(entire.char, from = names(new.char), to = new.char)
      names(new.mtchar)<-rownames(char_matrix)
      ### sort species to make the 1st one havinh non-polymorphic state
      sort=grepl(new.polymorph, new.mtchar)
      new.mtchar<-c(new.mtchar[!sort], new.mtchar[sort])
      ##
      complex.char$comb.matrices[[topo[i]]]$new.coding<-new.mtchar
      
      # rename and sort Q state names
      mt.new.name=mapvalues(colnames(cmb), from=names(Qmt.state.class$new.state), to=Qmt.state.class$new.state)
      cmb.new=cmb
      colnames(cmb.new)<-rownames(cmb.new)<-mt.new.name
      #colnames(cmb.new)<-rownames(cmb.new)<-c(3,4,1,2,5,6)
      colnames(cmb.new) %>% order()->ord
      cmb.new=cmb.new[ord,ord]
      complex.char$comb.matrices[[topo[i]]]$matrix.new=cmb.new
      
      ##########
      # st.iden=lapply(names(focal.v), function(x) complex.char$comb.matrices[[x]]$state.ident) %>% unlist()
      # complex.char$comb.matrices[[topo[i]]]$state.ident<-st.iden<-c(st.iden, topo[i])
      # 
      # # get observable and "hidden" states of the focal node
      # ln=length(st.iden)-1 # observable a/p are only those which have all states from other chrs=1
      # obs=which(substr(r.name, 1, ln)==paste(rep(1, ln), collapse=""))
      # names(obs)=r.name[obs]
      # complex.char$comb.matrices[[topo[i]]]$state.observable=obs
      # 
      # hid=(1:length(r.name))[-obs]
      # names(hid)<-r.name[hid]
      # complex.char$comb.matrices[[topo[i]]]$state.hidden=hid
      
      
      
    } #end if (length(focal.v)>0)
  } #end all
  
  return(complex.char)
} # end function


#################
# END FUNCTIONS
################################################################################################################################



#########################
# 
# TUTORIAL
#
#########################
# remove chars with two layers of dependecy 
# chars 199, 200, 201
#tb.red<-tb[!tb[,1] %in% c("CHAR:199", "CHAR:200", "CHAR:201"),]
# remove chars 101, 103, 102 from dependecy
#tb.red<-tb[!tb[,1] %in% c("CHAR:101", "CHAR:103", "CHAR:102"),]


gr.id=match(1:max(tb.red$dp.gr.num), tb.red$dp.gr.num)
g.edges<-tb[gr.id,c(3,1)] %>% as.matrix()

gr.id=match(1:max(tb$dp.gr.num), tb$dp.gr.num)
g.edges<-tb[gr.id,c(3,1)] %>% as.matrix()

# make igraph object
#c(g.edges[,1], g.edges[,2]) %>% unique()
# remove chars 101, 103, 102 from dependecy
#g.edges=g.edges[!g.edges[,1] %in% c("CHAR:101"),]

g=graph_from_edgelist(g.edges)
# plot all dependecies
plot(g, vertex.color="green", vertex.size=1,
     vertex.label.dist=0.5, vertex.label.cex=0.5, vertex.label=NULL, edge.arrow.size=.5 )

# recode all dependencies
#char.recode=get_graph_matrix_any(g, char_matrix, dep.tb, new.polymorph="&")
#save(char.recode, file="char.recode.RData")

str(char.recode)
# now let's select some components to get matrices for them
# there are a few ways to do it:

# 1. working with connected components
con.comp=components(g, "weak")
# let's get a component with three nodes
#com.id=which(con.comp$csize==5)
com.id=1
comp=con.comp$membership[con.comp$membership==com.id] %>% names()
# crating a new graph for this component
g1=subgraph(g, comp)
g2=subgraph(g, comp)
# plotting it
plot(g1)
char.recode=get_graph_matrix_any(g1, char_matrix.rem, dep.tb, new.polymorph=" ",  Lnew.polymorph=NULL, Rnew.polymorph=NULL)


subgraph.list<-list()
# all subgraphs of graph
for (i in 1:con.comp$no){
  com.id=i
  comp=con.comp$membership[con.comp$membership==com.id] %>% names()
  # crating a new graph for this component
  subgraph.list[[i]]=subgraph(g, comp)
}

plot(subgraph.list[[8]])
#for (i in 1:con.comp$no){
  for (i in 1:29){
  get_graph_matrix_any(subgraph.list[[i]], char_matrix.rem, dep.tb, new.polymorph="&")
  }


# getting all matrices for g1
char.recode=get_graph_matrix_any(g, char_matrix, dep.tb, new.polymorph="&")
# get a table of recodings
ln=length(char.recode$nodes.sorted )
i=2
mt.recoded<-matrix(ncol = 4, nrow = ln)
for (i in 1:ln){
  ch=char.recode$nodes.sorted[i]
  mt.recoded[i, 1]=ch
  tmp=names(char.recode$comb.matrices[[ch]]$Qmt.state.classes$new.state)
  if (!is.null(tmp)){
  mt.recoded[i, 2]=tmp %>% paste(., collapse=", ")
  mt.recoded[i, 3]=char.recode$comb.matrices[[ch]]$mt.state.class$comb.sep %>% paste(., collapse=", ")
  mt.recoded[i, 4]=char.recode$comb.matrices[[ch]]$char.mapping %>% paste(., collapse=", ")
  }
}

write.csv("char_recode.csv", mt.recoded)

####
mt=get_graph_matrix(g1)
str(mt)
# matrix for a term
mt$comb.matrices$`UBERON:2002076`$matrix



# 2. working with subgraphs
# let's have a look on the subgraph of 16 nodes
com.id=which(con.comp$csize==16)
comp=con.comp$membership[con.comp$membership==com.id] %>% names()
# craeting a new graph for this component
g2=subgraph(g, comp)
plot(g2)
# now I make a subgraph for "UBERON:2000663" that includes only the nearest neighbors
g3=make_ego_graph(g, order=1, nodes = "UBERON:2000663", mode = c("all"), mindist = 0)
g3=g3[[1]]
plot(g3)

# getting all matrices for g3
mt=get_graph_matrix(g3)

# get a list of all combined matrices
comb.mt=lapply(mt$nodes.sorted, function(x) mt$comb.matrices[[x]]$matrix)
names(comb.mt)=mt$nodes.sorted
comb.mt


#####################
#
# JUNK CODE
#
#####################

# faked toy exmple
chrs.uber<-matrix(c(
  "UBERON:2201587", "UBERON:2102027",
  "UBERON:2201587", "UBERON:4300092",
  "UBERON:2202028", "UBERON:2201587",
  "UBERON:4200103", "UBERON:2201587",
  "UBERON:4200105", "UBERON:2201587"
),
ncol=2, byrow = T
)



#sbc=subcomponent(g, "UBERON:2001537", mode="in")

# #sorting connected components and saving them
# ub.ids=cbind(names(con.comp$membership), as.numeric(unname(con.comp$membership)))
# ub.ids=ub.ids[order(ub.ids[,2]),]
# colnames(ub.ids)<-c("id", "component_category")
# 
# write.csv(ub.ids, file="chrs_as_connected_components.csv")
##########

# # sorting graph to navigate
# node.sorted=topo_sort(g, mode = c("out"))
# node.sorted=names(node.sorted)
# chrs.uber=chrs.uber[factor(chrs.uber[,1], levels=node.sorted) %>% order(),]
# # now graph form orderd table
# #g=graph_from_edgelist(chrs.uber)

# read dependecy table
chrs.depen<-read.csv("dependencies.txt", header=T,  stringsAsFactors = F, na.strings = "")

# get a dependent pair of chars
char.id<-2
seq.of.chars<-which(chrs.depen[,1]==char.id)
seq.of.chars<-chrs.depen[seq.of.chars,3]%>% unique()
# first chrs in the vecto depends on the following chars
seq.of.chars<-c(char.id, seq.of.chars)

#intialize N matrices where N=number of binary chars involved in a given dependecy scheme,
# all matrices are populated with different parameters
# this can be better done using array()
matrix.list=list()
#
# each row is the id of rate parameters which will be plugged in rate matrices
param.scheme<-matrix(seq(1:(2*length(seq.of.chars))), ncol = 2, byrow=TRUE)
for (i in seq_along(seq.of.chars)){
  matrix.list[[i]]<-matrix(c(-param.scheme[i,1],param.scheme[i,1], 
                             param.scheme[i,2],-param.scheme[i,2]),2,2,byrow=TRUE)
  rownames(matrix.list[[i]])<-colnames(matrix.list[[i]])<-c("0","1")
}

# combine matrices together given their dependecies
# this matrix can be used for inference with hidden models
combined.matrix=comb2matrices(matrix.list[[2]], matrix.list[[1]], dependent.state=2, diag.as=0)

# check lumpability in respect to given partitioning scheme for states reflecting meaningfull
# absence, so we can use non-hidden models
# I will work on this function to more to get it working properly
#part_scheme=list(c(1, 2), c(3), c(4))
#is_strg_lumpable(combined.matrix, part_scheme)
# if matrix is lumpable then we construct the aggregated chain, let's imagine the matrix
# is lumpable
is.lumpable=T

if (is.lumpable==T) {
  #identify which states are not meaningfull given absence
  # I will work on it to extend for arbitrary number of combined matrices
  # by far I just create a matrix how it should look like
  aggregated.matrix<-matrix(c(0,1,2,  3,0,4, 5, 6, 0),3,3,byrow=TRUE)
  rownames(aggregated.matrix)<-colnames(aggregated.matrix)<-c("0","10", "11")
  
}

### Outputs
combined.matrix # combined matrix
aggregated.matri # aggregated matrix
seq.of.chars # pair of chars where #1 depends on #2 (swith off dependency), also useful to 
#know character order in state names
