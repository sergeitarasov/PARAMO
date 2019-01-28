devtools::install_github("beast-dev/RBeast")
setwd("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/PolytomyResolver")
library("ape", lib.loc="~/.local/R/site-library")
# Polytomy resolver
str(hym.drop)
source("PolytomyResolver.R") 
PolytomyResolver(hym.dropA, error=.5, file.out="HymPolytomy05")

# read the best tree
tmp<-read.nexus("Hymenoptera_best_tree.nex")
tmp<-c()
class(tmp)
plot(tmp[[293]])
write.tree(tmp[[293]], file="Hymenoptera_best_tree.tre")

