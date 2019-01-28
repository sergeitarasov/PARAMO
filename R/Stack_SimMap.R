library(plyr)
library("tibble", lib.loc="~/.local/R/site-library")
library("stringdist") # Hamming distance between strings stringdist(c("35"), c("10"), method = "hamming")

#############################################
#
# FUNCTIONS
#
#####################

sim1 = read.simmap(file="CHAR-7.sstm", format="phylip")
sim2 = read.simmap(file="CHAR-30.sstm", format="phylip")

#### stack two discrete stm's lists; x,y are the list of state names (i.e. maps)
stack2<-function(x,y){
  mapply(function(x,y) 
  {paste(x,y, sep="") },
  x=x, y=y )
}

st1=lapply(stack.L[[4]]$maps, function(x) names(x))
st2=lapply(stack.L[[5]]$maps, function(x) names(x))
cbind(stack.L[[4]]$edge.length, stack.L[[5]]$edge.length)

st1=lapply(sim1$maps, function(x) names(x))
st2=lapply(sim2$maps, function(x) names(x))
stack2(st1, st2)
#######################
# stack N simmap objects
# sim.list is a list of Simmap object (trees with mapping)
stm.list<-stack.L[4:5]
template.tree<-stack.L[[1]]

lapply(stack.L[[4]]$maps, function(x) length(x)) %>% unlist ->a
lapply(stack.L[[5]]$maps, function(x) length(x)) %>% unlist ->b
a==b
cbind(a,b)
stack.L[[4]]$maps[[8]]
stack.L[[5]]$maps[[8]]

stack_stm<-function(stm.list){
  M<-lapply(stm.list, function(x) x$maps)
  M<-lapply(M, function(x) lapply(x, function(y) names(y)))
  M<-Reduce(stack2, M)
  
  M.out<-mapply(function(x,y) 
  {setNames(x, y) },
  x=stm.list[[1]]$maps, y=M )
  
  out<-stm.list[[1]]
  out$maps<-M.out
  return(out)
  
}
length(M)
# read first 10 char maps into list sim
c=char.rep[1:10,1]
sim=lapply(c, function(x) read.simmap(file=paste0(x, ".sstm"), format="phylip"))

# stack
z=stack_stm(sim)
z$maps
# states
lapply(z$maps, names) %>% unlist %>% unique->states
length(states)
# plot
plot(z, setNames(getPalette(length(states)), states))
plot.phylo(z, show.node.label=T)
#########################3
#
# calculate number of chages on discretized stack of simmaps over branches

tr<-stack_dist(z)
tr$maps
plot(tr)
tree<-stack.L[[311]]
stack.L[[311]]$maps[[150]]
stack.L[[311]]$maps[[4]]
i=150
sim.d[[1]][[311]]$maps[[150]]
sim.d[[2]][[311]]$maps[[150]]
sim[[2]][[311]]$maps[[150]]

stack_dist<-function(tree, init.val=0){

maps.out<-tree$maps
for (i in 1:nrow(tree$edge)){
  E<-tree$edge[i,]
  A<-anc_edge(E, tree$edge)
  
  if (A =="root"){ # at root the first bin recives init.val - changes cannot be detected there
    ME<-names(tree$maps[[i]])
      
      if (length(ME)>1)
        names(maps.out[[i]])<- stringdist(ME[2:length(ME)], ME[1:(length(ME)-1)], method = "hamming") %>%
            c(init.val, .)
      if (length(ME)==1)
        names(maps.out[[i]])<-init.val
  
  } else {
    MEA<- names(tree$maps[[A]])
    ME<-c(MEA[length(MEA)], names(tree$maps[[i]]) ) 
    names(maps.out[[i]])<-stringdist(ME[2:length(ME)], ME[1:(length(ME)-1)], method = "hamming")
    
    #stringdist(ME[2:length(ME)], ME[1:(length(ME)-1)], method = "hamming")
    #nchar(MEA)
    #nchar(ME)
  }
  
}
tree$maps<-maps.out
return(tree)
}
#########################3
#
# calculate mean number of chages on a sample of discretized stack of simmaps over branches
# the sample is given as a list 
library("reshape")

tree<-stack_dist(z)
tree$maps
plot(tree)
treeL<-list(tree, tree, tree)
tr<-stack_mean(treeL, norm=1)
tr$maps

stack_mean<-function(treeL, norm="mean"){
if (norm=="mean")
  norm<-length(treeL)  
  
maps.out<-treeL[[1]]$maps
M<-lapply(treeL, function(x) x$maps)

for (i in 1:length(maps.out)){
  BR<-lapply(M, function(x) x[i] %>% unlist %>% names %>% as.numeric)
  names(maps.out[[i]])<- Reduce("+", BR)/norm      
}

tree<-treeL[[1]]
tree$maps<-maps.out
return(tree)
}
###
####################################
# get ancetral edge of a focal eadge
edge<-z$edge
E<-edge[2,]
anc_edge(E, edge)

anc_edge<-function(E, edge){
  A<-which(edge[,2]==E[1])
  if (length(A)==0)
    return("root")
  else return(A)
  }
##########################################################
### merge the same discretized char categories over branch

merge_branch_cat<-function(br){
  i=2
  while (i<=length(br)){
    if( (names(br[i])) == names( br[i-1] )) {
      br[i-1]<-br[i-1]+br[i]
      br<-br[-i]
    } else{
      i=i+1
    }
  }
  return(br)
}

br=z$maps[[153]]
br
merge_branch_cat(br)

z$maps<-lapply(z$maps, merge_branch_cat)
plot(z, setNames(getPalette(length(states)), states))

z$maps %>% unlist
br=z$maps[[172]]
a=array(rbind(br, NA, NA), dim=c(3,4))
colnames(a)<-names(br)
###

####################################################################################
####################################################################################

# plotting

plotSimmap(sim[[1]],setNames(c("blue","red"), c(0,1)),ftype="off",lwd=4, node.numbers =T)
plotSimmap(sim[[2]],setNames(c("blue","red"), c(0,1)),ftype="off",lwd=4, node.numbers =T)
plot(sim1)
plot(sim2)
str(sim1)
sim1$maps
sim2$maps
###

# read in undesritezed trees
setwd("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/Run1_all/stm_R500")
setwd("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/Run1_all/stm_R100")
#c=char.rep[1:1,1]
sim=lapply(c, function(x) read.simmap(file=paste0(x, ".stmR"), format="phylip"))

# descritize trees
sim.d<-lapply(sim, function(x) discr_Simmap_all(x, 300))
sim.d[[1]][[1]]$maps

#######################
# Stacking archivi
#
# read descritized trees from archive
#setwd("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/Run1_all/stm_R500_discr")
setwd("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/Run1_all/stm_R500_discr_new")

#rm(F)
F1<-F.obj
# mouthparts
c<-sub(":", "-",
       get_descendants_chars(ONT, annotations="manual", F1$Level_1_ids[4] )  )
ntrees<-50
tr<-vector("list", ntrees)

i=2
for (i in 1:ntrees){
  
  fl<-paste0(c, "_", i, ".rds")  
  j=1
  stack.L<-vector("list", length(fl))
  for (j in 1:length(fl)){
    
  #con<-unz("CHAR-7.zip", filename="CHAR-7_1.rds")
    print(paste0("Reading ", paste0(c[j], ".zip"), " and ", fl[j]))
    con<-unz(paste0(c[j], ".zip"), filename=fl[j])
    con2 <- gzcon(con)
    stack.L[[j]] <- readRDS(con2)
    close(con)
  }
 
  tr[[i]]<- stack_stm(stack.L)
}

plot(tr[[1]])

stack_stm(stack.L[c(1,6)])
stack.L[4]

closeAllConnections()

countSimmap(tree=mat2, states=NULL, message=TRUE)

showConnections (all=T)
#####################################################################
###################
# get N random smples of trees and conver them in N stacks
n<-101
stack.L<-vector("list", n)
for (i in 1:n){
# get a random sample of trees
S<-lapply(sim.d, function(x) x[[i]])
#plot(S[[3]])

# stack the sample
stack.L[[i]]<-stack_stm(S)
#z$maps
}

####################

# calculate hamming for each stack
tr<-lapply(stack.L, function(x) stack_dist(x))
tr[[1]]$maps
length(tr)
# check Inf
lapply(tr, function(x) x$maps[[150]][1]) %>% unlist->In
which(names(In)=="Inf")
length(In)

# claculate the mean of hamming
tr1<-stack_mean(tr, norm="mean")
tr1$maps
tr1$maps[[150]]
plot(tr1)
#stack.L[[1]]$maps

tr<-tr[[2]]
str(tr)
tr$mapped.edge %>% head()
lapply(tr$maps, names) %>% unlist %>% unique->states
lapply(tr$maps, function(x) which(names(x)=="Inf"))->t
lapply(t, function(x) length(x)!=0) %>% unlist %>% which(.==T)

states
length(states)

states.c<-states %>% as.numeric() #%>% round(2)
states.c[states.c==Inf]<-0
max(states.c)
seq=seq(min(states.c),  max(states.c), length.out=100)
W<-sapply(states.c, function(x) which(abs(seq-x)==min(abs(seq-x)) ) ) %>% unlist
bluecols2[100:200][W]

# plot
plot(tr, setNames(getPalette(length(states)), states))
plot(tr, setNames(bluecols2[W], states))
plot(tr, setNames(bluecols2[c(20:69, 151:200)][W], states))

####
# Palettes
bluecols <- brewer.pal(11, 'Spectral')
bluecols <- brewer.pal(9, 'YlOrRd')
bluecols <- brewer.pal(9, 'Blues')
newcol <- colorRampPalette(bluecols)
ncols <- 200
bluecols2 <- newcol(ncols)#apply the function to get 100 colours
pie(rep(1, ncols), col = bluecols2, border = NA, labels = c(1:100))
#
