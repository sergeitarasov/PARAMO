library("phytools")
c=paste0("C", AN$CHAR_ID2)

setwd("~/Documents/Recon-Anc_Anat/Supplementary_materials/STEP_4/RevBayes/output")
# dir to write files
dir= ("~/Documents/Recon-Anc_Anat/Supplementary_materials/STEP_4/RevBayes/Discr_maps/")

# read a sample of 100 maps from .stm files and save them in the poper format .stmR
i=1
for (i in 1:length(c))
{
  tree<-read_Simmap_Rev(paste0(c[i], ".stm"),
                        start=400, end=500,
                        save = NULL) %>% read.simmap(text=., format="phylip")
  
  
  write.simmap(tree, file=paste0(dir, c[i], ".stmR"))
}
##########

tree.tmp.final<-read.tree("Hymenoptera_br_resolved.tre")
for (i in 1:length(c))
  { 
  
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
      
      saveRDS(sim.d, file =  paste0(dir,c[i], "_", j, ".rds") )
      
    }, error=function(e){
      cat("ERROR :",conditionMessage(e), "\n")
      #errors<-rbind(errors, c(ii,jj))
    }  )
    
  } 
  
  # putting rds files into archive
  files<-paste0(dir, c[i], "_", c(1:length(sim)), ".rds")
  zip(paste0(c[i], ".zip"), files=files)
  file.remove(files)
  
}
###################
# showConnections (all=T)
# closeAllConnections()
#########################

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

###########################
############################################
# Reading unsummarized Stoch Map files from ReVBayes

#' @param file file
#' @param start start from tree
#' @param end end with tree
#' @param save save to file. if NULL reads in R
#' 
#file="CHAR-1.stm"
#start=1
#end=2
#save="/home/tarasov/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/Run1_all/stm_R/CHAR-1.stmR"

#text <- scan(file=file, sep = "\n", what = "character")

#read_Simmap_Rev(file, start=1400, end=1500)

#sim1 = read.simmap(file=save, format="phylip")
#sim1 = read.simmap(text=trees, format="phylip")
#plot(sim1)

# sm2<-read_Simmap_Rev(paste0("CHAR-363", ".stm"), 
#                      start=1001, end=1500, 
#                      save = NULL) %>% read.simmap(text=., format="phylip")
file="CHAR-363.stm"
start=1
end=1
cat(trees, file="chr363.test", sep="\n")

read_Simmap_Rev<-function(file, start=1, end=1, save=NULL){
  
  skip=start+2
  max2read=end-start+1
  
  text <- scan(file=file, sep = "\n", what = "character", skip=skip, nlines=max2read)
  
  trees<-c()
  for (i in 1:length(text)){
    
    #trees[i]<-strsplit(text[i], "\\}\t\\(")[[1]][2]
    
    ss=regexpr("\\}\t\\(",  text[i])[1]
    trees[i]<-substring(text[i], first=ss+2)
  }
  
  if (is.null(save)){
    return(trees)
  }else{
    
    cat(trees, file=save, sep="\n")
    print(paste0("Tree(s) are saved to ", save))
  }
}
#########


# for one tree
discr_Simmap<-function(tree, res){
  
  steps <- 0:res/res * max(phytools:::nodeHeights(tree))
  H <- phytools:::nodeHeights(tree)
  maps.n <- vector(mode = "list", length = nrow(tree$edge))
  
  # i=170
  for (i in 1:nrow(tree$edge)) {
    YY <- cbind(c(H[i, 1], steps[intersect(which(steps > H[i, 1]), which(steps < H[i, 2]))]), 
                c(steps[intersect(which(steps > H[i, 1]), which(steps < H[i, 2]))], H[i, 2])) -  H[i, 1]
    
    
    TR<-cumsum(tree$maps[[i]])
    # TR[length(TR)]<-YY[nrow(YY), 2] # this to make the length equal as it sometiems does not hold
    # YY[,1]-YY[,2]
    ######
    #sprintf("%.54f", c(TR[length(TR)], YY[nrow(YY), 2]) )
    # TR[1]==TR[2]
    # all.equal(TR[1], TR[2])
    # all.equal(TR)
    # duplicated(TR)
    #length(int.out)
    #all.equal(TR[length(TR)], YY[nrow(YY), 2])
    #TR[length(TR)]==YY[nrow(YY), 2]
    ######
    
    #TR[length(TR)]==YY[nrow(YY), 2]
    int.out=findInterval(YY[,2], c(0,TR), left.open=T, rightmost.closed = FALSE, all.inside = TRUE)
    #int.out=findInterval(YY[,2], c(0,TR), left.open=T, rightmost.closed = FALSE)
    #findInterval(seq(0.1, 4, .1), c(0, 0.5, 0.7, 1.5, 1.6, 4 ), left.open=T, rightmost.closed = F)
    maps.n[[i]]<-setNames(YY[,2]-YY[,1], names(tree$maps[[i]])[int.out])
  }
  tree$maps<-maps.n
  return(tree)        
}

discr_Simmap_all<-function(tree, res){
  
  if (class(tree)[1]=="simmap") {
    tree<-discr_Simmap(tree, res)
  }
  
  if (class(tree)[1]=="multiSimmap") {
    
    for (j in 1:length(tree)){
      tree[[j]]<-discr_Simmap(tree[[j]], res)
    }
  }
  return(tree)
}



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

#########################