#rm(alphadata)
# working dir
#setwd("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/Run1_all/stm_R500_discr")
# which(c=="CHAR-377")
# length(c)
c=char.rep[,1]

setwd("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/Run1_all/stm_R500_discr_new")
# dir to read
dir= ("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/Run1_all/stm_R500/")

i=39

##################################################
# char 377 did not perform well; work it out separately
i=

setwd("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/Run1_all/stm")
chr377<-read_Simmap_Rev(file="CHAR-377.stm", 
                start=700, end=1500, 
                save = NULL) %>% read.simmap(text=., format="phylip")

tr<-chr377
errors<-c()

#### check if disretization is ok
i=1
for (i in 1:length(tr)) {
  
  print(paste0("Reading ", i) )
 
  tr[[i]]$maps %>% unlist->MM
  names(MM)->MMn
  
  c( 
    !(MM %>% is.na() %>% any),
    (MM > 0) %>% all,
    (MM!=Inf) %>% all,
    
    !(MMn %>% is.na() %>% any),
    (MMn!=Inf) %>% all
  ) ->out
  
  out<-!out
  
  if (any(out)){
    errors<-rbind(errors, c(i, !out) )
  }
}
##
# select only good trees
errors
tr<-tr[-errors[,1]]
tr<-tr[201:700]
sim<-tr
#c(tr, tr[sample(c(1:490), 10)])
#save
setwd("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/Run1_all/chr377")
write.simmap(tr, file="CHAR-377.stmR", append=FALSE, map.order=NULL, quiet=FALSE)
#############################
######################################### # Descritize: looping across all chars
#for (i in 1:length(c)){ 
  for (i in 378:length(c)){ 
  
  # read in undesritezed trees
  print(paste0("Reading ", c[i]))
  sim=read.simmap(file=paste0(dir, c[i], ".stmR"), format="phylip")
  
        # descritize trees by looping over sample and saving as rds
        #j=1
  

        for (j in 1:length(sim)){
        tryCatch({
          
          print(paste0("Descritizing tree ", j))
          
          ## errors with na
          
          ##
          
          ##### make trees equal with template
          sim.d<-make_tree_eq(tree.tmp.final, sim[[j]], round=5)
          ###
          
          #sim.d<-discr_Simmap_all(sim[[j]], 1000)
          sim.d<-discr_Simmap_all(sim.d, 1000)
          
          saveRDS(sim.d, file =  paste0(c[i], "_", j, ".rds") )
          
        }, error=function(e){
                  cat("ERROR :",conditionMessage(e), "\n")
                  #errors<-rbind(errors, c(ii,jj))
                          }  )
          
        } 
        
  # putting rds files into archive
  files<-paste0(c[i], "_", c(1:length(sim)), ".rds")
  zip(paste0(c[i], ".zip"), files=files)
  file.remove(files)

}
###################
# showConnections (all=T)
# closeAllConnections()
#########################

####################################################
# check if plotting of descrized and non descritezed trees is the same

# par(mfrow=c(1,2))
# plot(sim[[500]])
# plot(sim.d)
# sim.d$maps
# sim[[500]]$maps

setwd("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/Run1_all/stm_R500_discr_new")
con<-unz("CHAR-377.zip", filename = "CHAR-377_1.rds")
con2 <- gzcon(con)
tr.disc <- readRDS(con2)
close(con)

# use tr object from above (prediscretized)
par(mfrow=c(1,2))
plot(tr[[1]])
plot(tr.disc)
plot(sim.d)
tr.disc$maps

###
#tree.tmp<-read.tree("Hymenoptera_br_resolved.tre")

setwd("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/Run1_all/stm_R500")
tree.tmp1<-read.simmap(file="CHAR-363.stmR", format="phylip")
tree.tmp1[[1]]$maps
tree.tmp<-tree.tmp1[[1]]
# round edges and maps
tree.tmp$edge.length<-round(tree.tmp$edge.length, 5)
# make final tree tmp
tree.tmp.final<-make_tree_eq(tree.tmp, tree.tmp, round=5)
plot(tree.tmp.final)
all.equal(lapply(tree.tmp.final$maps, function(x) sum(x)) %>% unlist , tree.tmp.final$edge.length )
lapply(tree.tmp.final$maps, function(x) sum(x)) %>% unlist == tree.tmp.final$edge.length 

# target tree
target.tr1<-read.simmap(file="CHAR-364.stmR", format="phylip")
target.tr<-target.tr1[[2]]
plot(target.tr)

tree.test<-make_tree_eq(tree.tmp.final, target.tr, round=5)
plot(tree.test)
all.equal(lapply(tree.test$maps, function(x) sum(x)) %>% unlist , tree.tmp$edge.length )
lapply(tree.test$maps, function(x) sum(x)) %>% unlist == tree.tmp$edge.length 

# Test
# make trees equall, descritize
d1<-discr_Simmap(tree.tmp.final, 1000)
d2<-discr_Simmap(tree.test, 1000)
dd<-list(d1,d2)
stack_stm(dd)
# WORKS

dd1<-list(discr_Simmap(tree.tmp1[[1]], 1000), discr_Simmap(target.tr1[[2]], 1000) )
stack_stm(dd1)
# this doesnt

# make branch length and maps equal, use tree.tmp with rounded (3) edge.length
make_tree_eq<-function(tree.tmp, target.tr, round=3){
  target.tr$edge.length<-tree.tmp$edge.length
  
  target.tr$maps<-lapply(target.tr$maps, function(x) round(x,3) )
  Maps.targ<-lapply(target.tr$maps, function(x) sum(x)) %>%unlist()
  
  #k<-tree.tmp$edge.length-Maps.targ
  k<-tree.tmp$edge.length/Maps.targ
  k<-as.list(k)
  
  maps.out<-mapply(function(x,y) 
  {x*y },
  x=target.tr$maps, y=k )
  
  target.tr$maps<-lapply(maps.out, function(x) round(x, round))
  return(target.tr)
}
########
target.tr<-read.simmap(file="CHAR-364.stmR", format="phylip")
tt<-make_tree_eq_ALL(tree.tmp.final, target.tr, round=3)
stack_stm(tt)

make_tree_eq_ALL<-function(tree.tmp.final, target.tr, round=3){
  out<-lapply(target.tr, function(x)  make_tree_eq(tree.tmp.final, x, round=round) )
  return(out)
}

##############################################




