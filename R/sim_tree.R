## first load packages & source code
library(phytools)
library(geiger)
source("fitPagel.R")
.check.pkg<-phytools:::.check.pkg
## now let's simulate some uncorrelated data
tree<-pbtree(n=5,scale=3)
plot(tree)
tree$tip.label<-row.names(taxa)
plot(tree)
write.tree(file="test.tree", tree)
write.nexus(file="test.tree.nex", tree)

Q<-matrix(c(-1,.1,.2, .01, -1, 0.05,  0.16, 0.007, -1),3,3, byrow=T)


Q<-matrix(c(-1,1,1,-1),2,2)
rownames(Q)<-colnames(Q)<-c(1:2)
tt1<-sim.history(tree,Q, nsim=3)
plot(tt1)
z=density(tt1, bw=2)
z=densityMap(tt1)
plot(z)
str(z)
z$tree$maps


taxa1=cbind(row.names(taxa)[1:5], runif(5, 0, 3) %>% round())
taxa1=cbind(row.names(taxa)[1:5], runif(5, 0, 3) %>% round())
runif(87, 1, 4) %>% round()->taxa1


taxa1=cbind(names(tt2$states), setNames(tt2$states, NULL))
write.table(taxa1, file="char-test.char", quote = F, sep = " ", col.names = F, row.names = F)

str(tt1)
new.col=c(0.24314420/2, 0.24314420/2)
names(new.col)<-c("a", "b")
tt1$maps[[2]]<-new.col
tt1$edge.length
tt1$edge
tt1$mapped.edge[2,2]<-0.24314420/2
tt1$mapped.edge[2,1]<-0.24314420/2
tt1$states
tt1$node.states
write.simmap(tt1, file="simmap.writte.nex", map.order="left-to-right")
read.simmap(file="simmap.writte.nex", format="phylip")

tt2<-sim.history(tree,Q)
plot(tt2)
str(tt2)

make.simmap(tree, x, model="SYM", nsim=1, ...)
par(mfrow=c(1,2))
plotSimmap(tt1,setNames(c("blue","red"),letters[1:2]),ftype="off",lwd=1, node.numbers =T)
str(tt1)

plotSimmap(tt2,setNames(c("blue","red"),letters[1:2]),ftype="off",lwd=1,direction="leftwards")

Q<-matrix(c(-2,1,1,1,-2,1,1,1,-2),3,3)
rownames(Q)<-colnames(Q)<-letters[1:3]
tree<-pbtree(n=26,tip.label=LETTERS,scale=1)
tree

x<-sim.history(tree,Q)$states
x<-x[tree$tip.label] 
x[1]<-c("c/a")
tr1<-make.simmap(tree,x,model="ER")
plotSimmap(tr1,lwd=3)
