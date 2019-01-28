library(phytools)

# Discretize Simmap Tree by classifying states into bins

# tree<-read_Simmap_Rev(file, start=1500, end=1500) %>% read.simmap(text=., format="phylip")
# plot(tree)
# dtree=discr_Simmap(tree,1000)
# plot(dtree)
#dtree$maps

# for multiphylo
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
##
# discr_Simmap_all(sim[[j]], 1000)
tree<-sim.d
tree<-sim[[1]]
res<-1000
sim.d$maps[[172]]
sim.d$maps[[49]][[2]]==0
sim.d$maps %>% unlist %>% names
is.na(tree$maps %>% unlist %>% names) %>% any
lapply(tree$maps, function(x) which(is.na(x%>%names)==TRUE) )
tree$maps[[163]]
tree$edge.length[[169]]

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
##
#
# Read in undesritized sample (multiphylo) of trees and save to discritized format

tree<-read.simmap(file = , format="phylip")
# plot(tree)
# dtree=discr_Simmap(tree,1000)
# plot(dtree)
#dtree$maps

dir="/home/tarasov/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/Run1_all/stm_R500_discr/"

i=1
for (i in 1:length(char.rep[,1]) ){
  tree=read.simmap(file=paste0(char.rep[i,1], ".stmR"), format="phylip")
  tree<-discr_Simmap_all(tree, 10)
  #tree<-discr_Simmap_all(tree[[1]], 10)

  write.simmap(tree, file=paste0(dir, char.rep[i,1], ".stmR") )
}

class(tree[[1]])
tree[[100]]$maps
tree$maps
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









