library("klaR", lib.loc="~/.local/R/site-library")
library("EnsCat", lib.loc="~/.local/R/site-library")
### a 5-dimensional toy-example:

## generate data set with two groups of data:
set.seed(1)
x <- rbind(matrix(rbinom(250, 2, 0.25), ncol = 5),
           matrix(rbinom(250, 2, 0.75), ncol = 5))
colnames(x) <- c("a", "b", "c", "d", "e")
x
## run algorithm on x:
(cl <- kmodes(x, 2))

## and visualize with some jitter:
plot(jitter(x), col = cl$cluster)
points(cl$modes, col = 1:5, pch = 8)

### The running is time consuming
### Run hamming distance
data(alphadata)
str((alphadata))
dis0<-hammingD(alphadata)
### Save as distance format
REDIST<-as.dist(dis0)
### Run a hierarchical clustering using average linkage
hc0 <- hclust(REDIST,method = "average")
### plot the dendrogram
plot(hc0)

# test
rbind(
  sp1=c(0,0,0,0),
  sp2=c(0,0, 1, 3),
  sp3=c(0,2, 1, 3)
) %>% data.frame() ->t
dis0<-hammingD(t)
REDIST<-as.dist(dis0)
hc0 <- hclust(REDIST,method = "average")
plot(hc0)

