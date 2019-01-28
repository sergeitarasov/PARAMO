Q<-matrix(c(0,3,  3,0),2,2,byrow=TRUE)
rownames(Q)<-colnames(Q)<-c("0","1")
diag(Q)<-rowSums(Q)*-1
Q

tree<-pbtree(n=500,scale=100)
tt<-sim.history(tree,Q)
plot(tt)

recon <- rayDISC(tt, cbind(names(tt$states), tt$states), model="ARD",
                 node.states="marginal")
recon <- rayDISC(tt, st, model="ARD",
                 node.states="marginal")
str(recon)
recon$solution[1,2]

st=cbind(names(tt$states), tt$states)
st[,2]<-"0"

ln=c()
i=1
for (i in 1:nrow(st)){
  chrs=st
  chrs[i,2]="1"
  recon <- rayDISC(tt, chrs, model="ARD",
                   node.states="marginal")
  ln=rbind(ln, c(recon$solution[1,2], recon$solution[2,1]))
}

dev.off()
hist(ln[,1])
min(ln[,1])
1/sum(tt$edge.length)

1000/sum(tt$edge.length)

hist(tt$edge.length, breaks=30)

#1. ntaxa 100 scale 100
> 1/sum(tt$edge.length)
[1] 0.0004738088

> min(ln[,1])
[1] 0.0004745806

#2. ntaxa 100, scale=200
> min(ln[,1])
[1] 0.0001907497
> 1/sum(tt$edge.length)
[1] 0.000189283

#3. ntaxa 50, scale 300
> min(ln[,1])
[1] 0.0002437459
> 1/sum(tt$edge.length)
[1] 0.0002385987



######
gmodel<-function(param, low.bound, mean.val){
  (pgamma(low.bound, shape=mean.val*param, rate=param)-0.05)^2
}

emodel<-function(lambda, low.bound, quantl){
  (pexp(low.bound, lambda)- quantl)^2
}


optimize(gmodel, interval=c(0, 4), maximum=F, low.bound=0.0005, mean.val=0.5 )
0.5*.8

optimize(emodel, interval=c(0, 100), maximum=F, low.bound=0.0156, quantl=0.5)

pgamma(0.0001, shape=0.5*.703, scale=.703)
0.5*.8
0.35/.7

curve(dnorm(x, 2, 1), 0,5)
curve(dexp(x, 3.2), 0,5)

dgamma(0.2, shape=0.4, rate=0.6)
mean(rgamma(10000, .4,.8))
#######
dexp(0.1, 0.0157)
pexp(0.0156, 6.7)
curve(emodel(x, 0.0156, 0.5), 0, 100)
##########
#EXP
> optimize(emodel, interval=c(0, 100), maximum=F, low.bound=0.0156, quantl=0.05)
$minimum
[1] 3.288038

> optimize(emodel, interval=c(0, 100), maximum=F, low.bound=0.0156, quantl=0.1)
$minimum
[1] 6.753877

> optimize(emodel, interval=c(0, 100), maximum=F, low.bound=0.0156, quantl=0.3)
$minimum
[1] 22.86378

> optimize(emodel, interval=c(0, 100), maximum=F, low.bound=0.0156, quantl=0.5)
$minimum
[1] 44.43253



