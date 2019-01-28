# read in undesritezed trees
#setwd("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/Run1_all/stm_R500")
setwd("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/Run1_all/stm_R500")
#c=char.rep[1:1,1]
library("ggridges", lib.loc="~/.local/R/site-library")

#ST<-countSimmap(tree=sim[[1]], states=NULL, message=TRUE)
hist(ST$Tr[,1])
read.simmap(file="CHAR-199.stmR")
##################
#color<-c()

ST<-vector("list", nrow(F))
names(ST)<-F$Level_1
char.rep[,2]

i=11
for (i in 1:nrow(F)){

# get files to read
  c<-sub(":", "-",
         get_descendants_chars(ONT, annotations="manual", F$Level_1_ids[i] )  )
# c<-sub(":", "-",
# get_descendants_chars(ONT, annotations="manual", get_onto_id("fore wing", ONT) ) )

# read trees
sim=lapply(c, function(x) read.simmap(file=paste0(x, ".stmR"), format="phylip"))

# count number of changes
TM<-lapply(sim, function(x) countSimmap(tree=x, states=NULL, message=TRUE) )

# get as distribution
ST[[i]]<- lapply(TM, function(x) x$Tr[,1]) %>% unlist

}

# mean
lapply(ST, function(x) sum(x)/(500))->exp.ch
exp.ch<-list2edges(exp.ch)
F$mean.changes<-as.numeric(exp.ch[,2])
#F<-add_column(F, mean.changes=as.numeric(exp.ch[,2]))
F<-add_column(F,mean.changes.per.state=F$mean.changes/al$Level_1_nstates)
plot(al$Level_1_nstates, F$mean.changes )
plot(al$Level_1_nchars, F$mean.changes )
plot(al$Level_1_nstates, F$mean.changes.per.state)

# sd
lapply(ST, function(x) sd(x))->sd.ch
sd.ch<-list2edges(sd.ch)
F<-add_column(F, sd.changes=as.numeric(sd.ch[,2]))
plot(al$Level_1_nstates, F$sd.changes )
plot(al$Level_1_nchars, F$sd.changes )

plot(F$mean.changes , F$sd.changes )

write.csv(F, file="mean.changes.csv")
######
TM[[12]]$Tr[,1]


ST$`fore leg` %>% length
13*500
which(ST$`fore leg`==max(ST$`fore leg`))
ST$`fore leg`[5500:6000]
ch<-read.simmap(file=paste0(c[12], ".stmR"), format="phylip")
plot(ch[[1]])
ch[[100]]$maps
ch[[1]]
countSimmap(ch, states=NULL, message=TRUE)
lapply(sim, function(x) countSimmap(tree=x, states=NULL, message=TRUE)$Tr[,1] %>% mean )


#########
d <- data.frame(x = rep(1:5, 3), y = c(rep(0, 5), rep(1, 5), rep(2, 5)),
                height = c(0, 1, 3, 4, 0, 1, 2, 3, 5, 4, 0, 5, 4, 4, 1))
ggplot(d, aes(x, y, height = height, group = y)) + geom_ridgeline(fill = "lightblue")

ggplot(d, aes(x, y, height = height, group = y)) + 
  geom_density_ridges(stat = "identity", scale = 1)

ggplot(iris, aes(x = Sepal.Length, y = Species)) + geom_density_ridges()

###
ST.p<-setNames(unlist(ST[1:10], use.names=F),rep(names(ST[1:10]), lengths(ST[1:10])))
ST.p<-setNames(unlist(ST, use.names=F),rep(names(ST), lengths(ST)))
ST.p<-enframe(ST.p)

ggplot(ST.p, aes(x = value, y = name)) + geom_density_ridges()
hist(ST[[13]])
ST[[2]] %>% max
