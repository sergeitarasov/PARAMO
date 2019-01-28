# non dep chars
c=c("CHAR:103", "CHAR:102")
mt.c=init_matrices(c, char_matrix=char_matrix.rem, diag.as = 0)

#mm=read.csv(file = "matrix_1Feb-final_sp-selected.csv")
#char_matrix.rem[["CHAR:102"]]<-mm[["CHAR.102"]]
char_matrix.rem[["CHAR:102"]]<-revalue(char_matrix.rem[["CHAR:102"]], c("1"="0", "2"="1"))

char_id="CHAR:102"

char.recode.cor.all$comb.matrices[[char_id]]$matrix.new<-mt.c[[char_id]]
##################################33
# chars 199 200, 201

c=c("CHAR:199", "CHAR:200", "CHAR:201")
mt.c=init_matrices(c, char_matrix=char_matrix.rem, diag.as = 0)

m2=comb2matrices(mt.c[[1]], mt.c[[2]], dependent.state=2, diag.as = 0)
m3=comb2matrices(m2, mt.c[[3]], dependent.state=4, diag.as = 0)

matrix=m3

matrix.new=m3
colnames(matrix.new)<-rownames(matrix.new)<-seq(0, 7, 1)

# rcoding
paste(
char_matrix.rem[["CHAR:199"]],
char_matrix.rem[["CHAR:200"]],
char_matrix.rem[["CHAR:201"]],
sep=":")->ch
unique(ch)
ch.cor=revalue(ch, c("1:1:1"=7, "1:1:0"=6, "0:?:?"="0&1&2&3", "?:?:?"="?", "1:0:-"="4&5"))
names(ch.cor)<-rownames(char_matrix.rem)
ch.rev=revalue(ch, c("1:1:1"=7, "1:1:0"=6, "0:?:?"="0 1 2 3", "?:?:?"="?", "1:0:-"="4 5"))
names(ch.rev)<-rownames(char_matrix.rem)

char_id="CHAR:201"


char.recode.cor$comb.matrices[[char_id]]$matrix.new=matrix.new
char.recode$comb.matrices[[char_id]]$matrix.new=matrix.new

char.recode.cor$comb.matrices[[char_id]]$matrix=matrix
char.recode$comb.matrices[[char_id]]$matrix=matrix

char.recode.cor$comb.matrices[[char_id]]$new.coding=ch.cor
char.recode$comb.matrices[[char_id]]$new.coding=ch.rev

char.recode.cor$comb.matrices[[char_id]]$mt.state.class$comb.sep=unique(ch)
char.recode.cor$comb.matrices[[char_id]]$mt.state.class$comb.sep=unique(ch)

$mt.state.class$comb.sep
$matrix
$new.coding
$matrix.new
###########
char_id="CHAR:201"
mt.c[[char_id]]
char.recode.cor$binary.matrices[[char_id]]<-mt.c[[char_id]]

char.recode$binary.matrices<-char.recode.cor$binary.matrices

#########################################################################333
#chars plymorphs 17, 399

c=c("CHAR:17") # depends on 14-1
char_id="CHAR:24"
char.recode.cor.all$comb.matrices[[char_id]]
#char.recode.cor.all$binary.matrices[[char_id]]
#char.recode$binary.matrices<-char.recode.cor$binary.matrices

char.recode.cor$comb.matrices[[char_id]]$matrix.new.recon
char.recode$comb.matrices[[char_id]]$matrix.new

char.recode.cor$comb.matrices[[char_id]]$matrix
char.recode$comb.matrices[[char_id]]$matrix

char.recode.cor$comb.matrices[[char_id]]$new.coding
char.recode$comb.matrices[[char_id]]$new.coding

char.recode.cor$comb.matrices[[char_id]]$mt.state.class$comb.sep
char.recode.cor$comb.matrices[[char_id]]$mt.state.class$comb.sep

$mt.state.class$comb.sep
$matrix
$new.coding
$matrix.new
###########

paste(
  char_matrix.rem[["CHAR:400"]],
  char_matrix.rem[["CHAR:399"]],
  sep=":")->ch
unique(ch)

$char.mapping
1:1     1:2   1:1/2     0:?     ?:?     1:0     1:? 
"4"     "5"   "5&6"->"4&5" "0&1&2"     "?"     "3" "3&4&5" 

$char.mapping
1:0   1:1   0:-   ?:? 1:0/1   1:? 
"2"   "3" "0&1"   "?" "3&4"->"2&3" 

ch.cor=revalue(char.recode.cor$comb.matrices[[char_id]]$new.coding, c("3&4"="2&3"))
char.recode.cor$comb.matrices[[char_id]]$new.coding<-ch.cor
# names(ch.cor)<-rownames(char_matrix.rem)
ch.rev=revalue(char.recode$comb.matrices[[char_id]]$new.coding,  c("3 4"="2 3"))
char.recode$comb.matrices[[char_id]]$new.coding<-ch.rev
#######################################################################
# polymorphic chars independent
poly=c("CHAR:352", "CHAR:361", "CHAR:362", "CHAR:369", "CHAR:384", "CHAR:388")

#char_id="CHAR:399  "

i=2
for (i in seq_along(poly)){
  char_id=poly[i]
  print(char_id)
  char_matrix.rem[[char_id]]<-gsub("/", "&", char_matrix.rem[[char_id]])
}

char.recode.cor.all$comb.matrices

mt.c=init_matrices(poly, char_matrix=char_matrix.rem, diag.as = 0)

for (i in seq_along(poly)){
  char_id=poly[i]
  print(char_id)
char.recode.cor.all$comb.matrices[[char_id]]$matrix.new<-mt.c[[char_id]]
}
