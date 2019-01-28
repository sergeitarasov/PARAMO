
############################################################
# Making character report
#
#

runs.done<-fl.id
runs.done=sort(runs.done)
#######
# Character table
char.rep<-rbind(
cbind(c(anal3.revision, indep.n), "nondep"),
cbind(contr.chars, "contr"),
cbind(dep.chars1, "dep"),
cbind(add.chars, "add")
)

#write.csv(char.rep, "character_table.csv")
unique(char.rep[,1]) %>% length()
char.rep[,1][duplicated(char.rep[,1])]
# reading char table  of completed chars
char.rep=read.csv("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/dataset/Data_final/character_table.csv", stringsAsFactors =  FALSE, header=T)
CH<-as_tibble(char.rep)
#char.rep=cbind(file=gsub(":", "-", char.rep$char.id), char.rep)

cat(file="file_list.txt",
paste0("\"", char.rep$file, "\"", collapse=", ") %>% paste0("[", ., "]")
)


#############################################################
# Read RevBayes stochastic mapping trees and svae 100 of 
# them to R friendly format

source("Intervals.R")
setwd("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/Run1_all/stm")

dir="/home/tarasov/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/Run1_all/stm_R500/"

# Read rev files and save them to file

i=1
#for (i in 1:3){
for (i in 171:length(char.rep[,1])){
  read_Simmap_Rev(paste0(char.rep[i,1], ".stm"), 
                  start=1001, end=1500, 
                        save = paste0(dir, char.rep[i,1], ".stmR"))
}
#

##################################################################

# 
# Read all mcmc outpur and save the last 1000 samples to mcmc.per.char obj
#
library(coda) 
library("dplyr", lib.loc="~/.local/R/site-library")
library("tibble", lib.loc="~/.local/R/site-library")

thin=500

mcmc.per.char<-vector(mode="list", length = length(char.rep[,1]))
names(mcmc.per.char)<-char.rep[,2]

i=1
for (i in 1:length(char.rep[,1])){

    file=paste0(char.rep[i,1], ".log")
    chain1=read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    chain1<-as_tibble(chain1) %>% select(-c(1:4))
    chain1=slice(chain1, 5002:n())
    mcmc.per.char[[char.rep[i,2]]]$run1<-chain1
}

lapply(mcmc.per.char, length) %>% unlist->tmp
which(tmp==0)
#save(file="mcmc.per.char.RData", mcmc.per.char)

unlist(mcmc.per.char)->tmp
hist(tmp, breaks=500,  freq = T, xlim = c(0,0.1))
max(tmp)
#chain1=mcmc(chain1, thin=thin)
#HPDinterval(chain1)

#### getting numer of states per character
setwd("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/Run1_all")
load(file="mcmc.per.char.RData")
mcmc.per.char$`CHAR:7`$run1 %>% ncol

lapply(mcmc.per.char, function(x) x$run1 %>% ncol)
################################################################
plot(hym.res)



Sharkey_2011_annot[10]
a=list2edges(Sharkey_2011_annot)

char.recode$comb.matrices %>% length
char.recode.cor$comb.matrices %>% length
char.recode$comb.matrices["CHAR:376"]

char.recode.cor$comb.matrices[["CHAR:201"]]
char.recode$comb.matrices[["CHAR:201"]]

load("char.recode.cor.all.RData")
load("char.recode.cor.RData")
load("char.recode.RData")
load("char.recode-v2.RData")
#################
#
# MERGING STATES IN DEPENDENT CHARACTERS ON SIMMAPS
#
#

################
get_states2merge<-function(XX){
N=sapply(names(XX), function(x) substr(x, 2,2))
L<-vector("list",length(unique(N)))
names(L)<-unique(N)

for (i in seq_along(unique(N))){
L[[i]]<-XX[which(N==N[i])]
}
return(L)
}

XX<-char.recode.cor$comb.matrices[["CHAR:376"]]$Qmt.state.classes$new.state
char.recode.cor$comb.matrices[["CHAR:24"]]$Qmt.state.classes$new.state
Rec=get_states2merge(XX)

# getting number of states per character
lapply(char.recode.cor$comb.matrices, function(x) x$Qmt.state.classes$new.state %>% get_states2merge(.) %>% length)->X1
X1<-unlist(X1)
CH<-add_column(CH, X1=NA)
CH$X1<-X1[match(CH$char.id, (X1 %>% names()))]


char.recode.cor.all$init.reoder %>% length()
lapply(char.recode.cor.all$init.reoder, function(x) x %>% nrow)->X1
X1<-unlist(X1)
CH<-add_column(CH, X2=NA)
CH$X2<-X1[match(CH$char.id, (X1 %>% names()))]

char.recode.cor$binary.matrices %>% length()
lapply(char.recode.cor$binary.matrices, function(x) x %>% nrow)->X1
X1<-unlist(X1)
CH<-add_column(CH, X3=NA)
CH$X3<-X1[match(CH$char.id, (X1 %>% names()))]

M<-read.csv("matrix_1Feb-final_sp-selected.csv", stringsAsFactors =  FALSE, header=T, row.names=1)
apply(M, 2, unique)->z

# remove ?
lapply(z, function(x) paste(x[x!="?"], collapse=" ") ) %>% unlist ->z1
z1
X1<-z1
CH<-add_column(CH, X5=NA)
CH$X5<-X1[match(CH$char.id, (gsub("\\.", ":", X1 %>% names()) ) )]


lapply(z, function(x) x[x!="?"] %>% length) %>% unlist ->z2
X1<-z2
CH<-add_column(CH, X4=NA)
CH$X4<-X1[match(CH$char.id, (gsub("\\.", ":", X1 %>% names()) ) )]

write.csv(CH, file="char_rep_withStatesN.csv", na = "", row.names = FALSE)
##################
#############
# merge multiple states in multiphylo
tree<-read.simmap(file="CHAR-376.stmR", format="phylip")
XX<-char.recode.cor$comb.matrices[["CHAR:376"]]$Qmt.state.classes$new.state
Rec=get_states2merge(XX)

treeM=merge_mult_states_All(tree, Rec)
tree[[1]]$maps
treeM[[1]]$maps

i=4
merge_mult_states_All<-function(tree, Rec){
for (i in seq_along(tree)){
  tree[[i]]<-merge_mult_states(tree[[i]], Rec)
}
  return(tree)
}
################
#############
# merge multiple states in one tree
i=172
merge_mult_states<-function(tree, Rec){
# remove states whicg are not present in Simmap
S=tree$mapped.edge %>% colnames()
Rec=lapply(Rec, function(x) x[x%in% S])

# merge states
#i=1
for (i in seq_along(Rec)){
  tree=mergeMappedStates(tree, Rec[[i]], paste0("m", names(Rec[i])) )
}

# rename states back to leave only numbers
for (i in seq_along(Rec)){
  tree=mergeMappedStates(tree, paste0("m", names(Rec[i])), names(Rec[i]) )
}
return(tree)
}
############################
library("dplyr", lib.loc="~/.local/R/site-library")
head(char.rep)
DP<-char.rep %>% filter(dep=="dep")
DP[,2][!DP[,2] %in% (char.recode.cor$comb.matrices %>% names)]
(char.recode.cor$comb.matrices %>% names)[! (char.recode.cor$comb.matrices %>% names) %in% DP[,2]]

# removng two-level dep to make it manually later
DP<-DP %>% filter(file!="CHAR-201" & file!="CHAR-199" & file!="CHAR-200")


i=6
dir="/home/tarasov/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/Run1_all/stm_R500D/"

for (i in 7:length(DP[,1])){
  tree=read.simmap(file=paste0(DP[i,1], ".stmR"), format="phylip")
  XX<-char.recode.cor$comb.matrices[[DP[i,2]]]$Qmt.state.classes$new.state
  Rec=get_states2merge(XX)
  
  treeM=merge_mult_states_All(tree, Rec)
  write.simmap(treeM, file=paste0(dir, DP[i,1], ".stmR"))
}

########################
# working on chars file!="CHAR-201" & file!="CHAR-199" & file!="CHAR-200")
#
char.recode.cor$comb.matrices[["CHAR:201"]]

XX<-char.recode.cor$comb.matrices[["CHAR:201"]]$matrix %>% colnames()
names(XX)<-char.recode.cor$comb.matrices[["CHAR:201"]]$matrix.new %>% colnames()
XX

tree=read.simmap(file=paste0("CHAR-201", ".stmR"), format="phylip")
tree[[1]]
treeM[[1]]
plot(treeM[[1]])

# CHAR 199
x0=names(XX)[which(substr(XX, 1,1)=='0')]
x1=names(XX)[which(substr(XX, 1,1)=='1')]
Rec=list('0'=x0, '1'=x1)
treeM=merge_mult_states_All(tree, Rec)
write.simmap(treeM, file=paste0(dir, "CHAR-199", ".stmR"))

# CHAR 200
x0=names(XX)[which(substr(XX, 2,2)=='0')]
x1=names(XX)[which(substr(XX, 2,2)=='1')]
Rec=list('0'=x0, '1'=x1)
treeM=merge_mult_states_All(tree, Rec)
write.simmap(treeM, file=paste0(dir, "CHAR-200", ".stmR"))

# CHAR 201
x0=names(XX)[which(substr(XX, 3,3)=='0')]
x1=names(XX)[which(substr(XX, 3,3)=='1')]
Rec=list('0'=x0, '1'=x1)
treeM=merge_mult_states_All(tree, Rec)
write.simmap(treeM, file=paste0(dir, "CHAR-201", ".stmR"))

################################################
#
# READ ANNOTAIONS
#
setwd("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/dataset/Data_final")
H<-read.csv("finest_HAO.csv", header=T,  stringsAsFactors = F, na.strings = "")
H<-read.csv("fines_HAO_shortlist.csv", header=T,  stringsAsFactors = F, na.strings = "")
H<-read.csv("finest_HAO_shortlist3.csv", header=T,  stringsAsFactors = F, na.strings = "")
###
H[,2]
H<-cbind(H, get_onto_name(H[,2], HAO) )
H<-cbind(H, NA )

H[!is.na(H[,3]), 7]<-get_onto_name(H[,3], HAO)
H[!is.na(H[,3]), 7] %>% length()
get_onto_name(H[,3], HAO)

which(is.na(H[,3])==F) %>% length
for ( i in 1:nrow(H)){
 if (!is.na(H[i,3]))
  H[i,7]<- get_onto_name(H[i,3], HAO)
}

sapply(H[,2], function(x) get_onto_name(x, HAO) %>% length)->z
sapply(H[,3], function(x) get_onto_name(x, HAO) %>% length)->z
which(z==0)

get_onto_name(c("HAO:0001140"), HAO)
#write.csv(H, file="finest_HAO_shortlist3a.csv")
##
C<-read.csv("Terms4Graphics.csv", header=T,  stringsAsFactors = F, na.strings = "")
# C<-cbind(C, C[,1]==get_onto_name(C[,2], HAO), get_onto_name(C[,2], HAO) )
# focal terms for graphics and their layer id
F<-set_names(C[,3], C[,2])
str(F)
F<-F[!is.na(F)]
##
H<-H[,-4]
H[,1]<-paste0("CHAR:", H[,1])
H<-set_names(table2list(H), H[,1])

ONT<-HAO
ONT$terms_selected_id<-H

stat<-sapply(names(F), function(x)
get_descendants_chars(ONT, annotations="manual", terms=x) %>% length )

col<-cbind(F, stat)
sum(col[,2])
col<-as.data.frame(col)
col<-as_tibble(col)
col<-tibble::add_column(col, hao=names(F), .before =T)
#cbind(
#get_onto_name(names(stat), HAO), stat)

##
get_descendants_chars(ONT, annotations="manual", terms="HAO:0000351") #fore wing
get_descendants_chars(ONT, annotations="manual", terms="HAO:0000400") #hind wing
get_descendants_chars(ONT, annotations="manual", terms="HAO:0001089") # wing

get_descendants_chars(ONT, annotations="manual", terms="HAO:0000020") # abdomenal segment 2
adb=get_descendants_chars(ONT, annotations="manual", terms="HAO:0000015") # abdomen
get_descendants_chars(ONT, annotations="manual", terms="HAO:0000626") # metasoma
get_onto_name(get_descendants(ONT, roots="HAO:0000626"),
              ONT)

mp=get_descendants_chars(ONT, annotations="manual", terms="HAO:0000604") # metapectal propodeal complex

adb %in% mp

chars_per_term(ONT, annotations = "manual")

# read in focal terms
TR<-read.csv("Terms4Graphics.csv", header=T,  stringsAsFactors = F, na.strings = "")
# query number of chars per each term
x<-sapply(TR[,2], function(x) get_descendants_chars(ONT, annotations="manual", terms=x) %>% length)
TR<-cbind(TR, x)

# sum chars over graphical terms and check if all summed chars are present in original annot
x<-sapply(TR[,2], function(x) get_descendants_chars(ONT, annotations="manual", terms=x))
x<-unlist(x) %>% unique()
names(H)[!names(H) %in% x]

# Problems
# 1. Metasoma