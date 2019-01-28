library(plyr)
library("tibble", lib.loc="~/.local/R/site-library")
# Hamming distance between strings
library("stringdist")
stringdist(c("35"), c("10"), method = "hamming")

tree<-pbtree(n=5,scale=2)

Q<-matrix(c(-1,1,1,-1),2,2)
rownames(Q)<-colnames(Q)<-c(0,1)
sim<-sim.history(tree,Q, nsim = 1)
plot(sim)
str(tree)

plotSimmap(sim[[1]],setNames(c("blue","red"), c(0,1)),ftype="off",lwd=4, node.numbers =T)
plotSimmap(sim[[2]],setNames(c("blue","red"), c(0,1)),ftype="off",lwd=4, node.numbers =T)

sim[[1]]

sim1 = read.simmap(file="CHAR-7test.sstm", format="phylip")
sim2 = read.simmap(file="CHAR-30test.sstm", format="phylip")

sim1 = read.simmap(file="CHAR-7.sstm", format="phylip")
sim2 = read.simmap(file="CHAR-30.sstm", format="phylip")

plot(sim1)
plot(sim2)

str(sim1)
sim1$maps
sim2$maps
st1=lapply(sim1$maps, function(x) names(x))
st2=lapply(sim2$maps, function(x) names(x))

cmb=mapply(function(x,y) 
  {paste(x,y, sep="") },
 x=st1, y=st2 )

#### stack two discrete stm's lists; x,y are the list of state names
stack2<-function(x,y){
mapply(function(x,y) 
{paste(x,y, sep="") },
x=x, y=y )
}
stack2(st1, st2)
#######################
st3=Reduce(stack2, list(st1, st2))

#########
stm.list<-list(sim1, sim2)
z=stack_stm(stm.list)
z$maps
plot(z)
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
############

# read first 10 char maps into list
c=char.rep[1:10,1]
sim=lapply(c, function(x) read.simmap(file=paste0(x, ".sstm"), format="phylip"))

z=stack_stm(sim)
z$maps
lapply(z$maps, names) %>% unlist %>% unique->states
length(states)
plot(z, setNames(getPalette(length(states)), states))
######

### merge the same discretized char categories over branch
br=z$maps[[153]]
br
merge_branch_cat(br)

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
################
z$maps<-lapply(z$maps, merge_branch_cat)
plot(z, setNames(getPalette(length(states)), states))

z$maps %>% unlist
br=z$maps[[172]]
a=array(rbind(br, NA, NA), dim=c(3,4))
colnames(a)<-names(br)

ld=ldply(z$maps, rbind)
kmeans(ld,2)

# unlist list and back
H<-c(1: length(z$maps))
W<-lapply(z$maps, length) %>% unlist

v=unlist(z$maps, use.names = F)

LL=z$maps


for (i in 1:length(LL)){
  LL[[i]]
}

W<-lapply(LL, length) %>% unlist
H<-c(1: length(LL))
CL<-mapply(rep, H, W) %>% unlist

LL.U<-unlist(LL)

out<-list()
for (i in 1:max(CL)){
  out[[i]]<-LL.U[CL==i]
}

# classify states for likelihood
LL<-sim$maps 
LL.n=lapply(LL, function(x) c(rep("D", length(x)-1), "E") )
# edges id which are tips
t.e=match(1:length(sim$tip.label), sim$edge[,2])

for (i in t.e){ 
  LL.n[[i]] %>% length->en
  LL.n[[i]][en]<-"T"
  }

D<-data.frame(time=unlist(LL),
              state=names(unlist(LL)),
              Lcat=unlist(LL.n)
)

D %>%
  group_by(state, Lcat) %>%
  summarize(mean_vals = mean(time))

filter(D, state==0 & Lcat=="E")

iris <- as_tibble(iris) # so it prints a little nicer
select(iris, starts_with("Petal"))
select(iris, ends_with("Width"))

exp(-2*30)*dexp(1e-10, 2)
dexp(30, 2)
###########

##########################
#
# Classify mapped states into E, D, T; and return as data frame

#' @param LL maps from Simmap

LL<-sim$maps 
class_maps(sim$maps)

class_maps<-function(LL){
  LL.n=lapply(LL, function(x) c(rep("D", length(x)-1), "E") )
  
  # edges id which are tips
  t.e=match(1:length(sim$tip.label), sim$edge[,2])
  
  for (i in t.e){ 
    LL.n[[i]] %>% length->en
    LL.n[[i]][en]<-"T"
  }
  
  # making datat frame
  D<-data.frame(time=unlist(LL),
                state=names(unlist(LL)),
                Lcat=unlist(LL.n)
  )
  
  return(D)
}
##

####################
#
# ML of Exponential distribution
r=rexp(100, 2)
dexp(r, 2) %>% prod() %>% log()
L_exp(2, r)

L_exp<-function(lambda, data, log=F){
  L<-log(lambda^length(data))-(lambda*length(data)*mean(data))
  if (!log) { L<-exp(L)}
  return(L)
}
##

#####################
#
# ML of Possion seeing zero changes
r=c(1,2,3)
L_pois0(2, r)

L_pois0<-function(lambda, data, log=F){
  L<-exp(-lambda*sum(data))
  if (log) {log(L)}
  return(L)
}
##
#####################################
#
# Initilize state data frame S given mumber of categoriues to estimate

D<-class_maps(sim$maps)
as.tibble(D)
S<-distinct(D, state) %>% as.tibble()
#S<-add_column(S, prob=NA)
#S<-c()

ncat<-3
S=init_S(D, 1)

init_S<-function(D, ncat=1){
  S<-distinct(D, state) %>% as.tibble()

  m=matrix(0, ncol=ncat*2+1, nrow=nrow(S)) %>% as.tibble()
  colnames(m)<-c( paste0("prob", 1:ncat), paste0("post", 1:ncat), "sum.of.comps"  )
  S<-bind_cols(S, m)
return(S)

}

###################################################
#
# Claculate states probailities

#' @param S state data frame
#' @param D category data frame with times
D<-class_maps(sim$maps)
S<-distinct(D, state)
S<-add_column(S, prob=NA)
S<-c()

S=init_S(D, 2)
state_probs(S, D, lambda=c(1,2), alpha=c(1,2))

state_probs<-function(S, D, lambda, alpha){

for (i in 1:nrow(S)){
    for (j in 1:((ncol(S)-1)/2)){
  col<-paste0("prob", j)
  s.val<-S$state[i]
  prod(
  filter(D, state==s.val & Lcat=="D")$time %>% L_exp(lambda = lambda[j], .),
  filter(D, state==s.val & Lcat=="E")$time %>% L_pois0(lambda = lambda[j], .), # removed ^2
  filter(D, state==s.val & Lcat=="T")$time %>% L_pois0(lambda = lambda[j], .)
  )*alpha[j] -> S[[col]][i]
    
  }
}
  return(S)
}
##

###################################################
#
# Claculate posterioirs

S=init_S(D, ncat=2)
S<-state_probs(S, D, lambda=c(1,2), alpha=c(1,2))
S<-state_post(S)

state_post<-function(S){
  
  S$sum.of.comps<-apply(S[grep("prob", colnames(S))]
                                   , 1, sum)
  
  S[grep("post", colnames(S))]<-
    apply(S[grep("prob", colnames(S))], 2, function(x) x/S$sum.of.comps)
  
  sum.of.comps.ln.sum<- sum(S$sum.of.comps) %>% log() %>% sum
  
  return(list(S=S, loglik=sum.of.comps.ln.sum))
}
###

#######################
S<-S$S
m_step(S=S$S, D)

m_step <- function(S, D) {
  
  comp.n<-apply(S[grep("post", colnames(S))], 2, sum)
  
  #comp1.n <- sum(posterior.df[, 1])
  #comp2.n <- sum(posterior.df[, 2])

  D %>% group_by(state) %>% summarize(sum = sum(time)) %>% 
    .$sum * S[grep("post", colnames(S))] ->PostD
  
  lambda.n<-comp.n/apply(PostD, 2, sum)
  
  #comp1.lambda <- comp1.n / sum(posterior.df[, 1] * x)
  #comp2.lambda <- comp2.n / sum(posterior.df[, 2] * x)
  
  alpha.n=comp.n/nrow(S)
  # comp1.alpha <- comp1.n / length(x)
  # comp2.alpha <- comp2.n / length(x)
  
  list("lambda" = lambda.n,
       "alpha" = alpha.n)
}

####

for (i in 1:50) {
  if (i == 1) {
    # Initialization
    
    #e.step <- e_step(data, data.summary.df[["lambda"]], data.summary.df[["alpha"]])
    # e-step
    S=init_S(D, ncat=2)
    S<-state_probs(S, D, lambda=c(.1,.2), alpha=c(.5,.5))
    S<-state_post(S)
    ##
    
    # m-step
    #m.step <- m_step(data, e.step[["posterior.df"]])
    m.step <-m_step(S=S$S, D)
    
    #cur.loglik <- e.step[["loglik"]]
    #loglik.vector <- e.step[["loglik"]]
    
    cur.loglik <- S[["loglik"]]
    loglik.vector <- S[["loglik"]]
  
    } else {
    # Repeat E and M steps till convergence
    
    #e.step <- e_step(data, m.step[["lambda"]], m.step[["alpha"]])
      S<-state_probs(S$S, D, lambda=m.step[["lambda"]], alpha=m.step[["alpha"]])
      S<-state_post(S)
    
    #m.step <- m_step(data, e.step[["posterior.df"]])
      m.step <-m_step(S=S$S, D)
    
    loglik.vector <- c(loglik.vector, S[["loglik"]])
    
    loglik.diff <- abs((cur.loglik - S[["loglik"]]))
    if (loglik.diff < 1e-6) {
      break
    } else {
      cur.loglik <- S[["loglik"]]
    }
  }
}
loglik.vector
