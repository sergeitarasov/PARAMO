#Get ML estimates
char.recode.all=get_graph_matrix_any(g, char_matrix.rem, dep.tb, new.polymorph=" ",  Lnew.polymorph=NULL, Rnew.polymorph=NULL)

# getting a vector of all chars in dep table + char193 that has to be removed
c(g.edges[,1], g.edges[,2]) %>% c(., "CHAR:193") %>% unique() ->tb.chars

indep.chars=colnames(char_matrix.rem)[!colnames(char_matrix.rem) %in% tb.chars]

# remove chars with two layers of dependecy 
# chars 199, 200, 201
#tb.red<-tb[!tb[,1] %in% c("CHAR:199", "CHAR:200", "CHAR:201"),]
#
gr.id=match(1:max(tb.red$dp.gr.num), tb.red$dp.gr.num)
g.edges<-tb[gr.id,c(3,1)] %>% as.matrix()

gr.id=match(1:max(tb$dp.gr.num), tb$dp.gr.num)
g.edges<-tb[gr.id,c(3,1)] %>% as.matrix()

exp=matrix(c("c1", "c2", "c3", NA, NA, NA), 3,2)
exp=matrix(c("c1", "c2", "c3"), 3,1)
g1=graph_from_edgelist(exp)
plot(g1)

# make igraph object
g=graph_from_edgelist(g.edges)
# plot all dependecies
plot(g, vertex.color="green", vertex.size=1,
     vertex.label.dist=0.5, vertex.label.cex=0.5, vertex.label=NULL, edge.arrow.size=.5 )

char_id="CHAR:14"
char.recode.cor$comb.matrices[[char_id]]$matrix.new
char.recode.cor$comb.matrices[[char_id]]$matrix.new.recon$solution
####################################################################################
# ML estimates
hym.res=read.tree(file="Hymenoptera_br_resolved.tre")

char.recode.cor.all<-list()
char.recode.cor.all$init=init_matrices(indep.chars, char_matrix.rem, diag.as = 0)
tmp=lapply(char.recode.cor.all$init, function(x) nrow(x)) %>% unlist()
tmp[tmp>2] %>% length
char.recode.cor.all$init[["CHAR:46"]]

init_matrices(char.revision, char_matrix=char_matrix.rem, diag.as = 0)

#char_matrix.rem[,indep.chars]

char_id="CHAR:7"
i=1

anal3.revision
char_id="CHAR:8"
contr.chars

#for (i in seq_along(anal3.revision)){
#for (i in seq_along(indep.chars)){
  #for (i in seq_along(contr.chars)){
for (i in seq_along(poly)){
  tryCatch({
  
  #char_id=indep.chars[i]
    char_id=poly[i]
  
  print(paste0("Working on ", char_id))
  
  
  #mt=char.recode.cor.all$init[[char_id]]
  
  #mt=char.recode.cor.all$init[[char_id]]
  mt=char.recode.cor.all$comb.matrices[[char_id]]$matrix.new
  
  #mt=char.recode$comb.matrices[[char_id]]$matrix.new
  mt[mt==0]<-NA
  ####
  #rownames(mt)<-colnames(mt)<-c(0,1)
  ####
  
  taxa=char_matrix.rem[char_id]
  
  taxa=as.data.frame(cbind(rownames(taxa), taxa),stringsAsFactors =F)
  
  recon=c()
  
  #recon <- rayDISC(hym.res, taxa, rate.mat=mt, node.states="marginal", model="ARD", root.p=NULL)
  recon <- rayDISC(hym.res, taxa, rate.mat=NULL, node.states="marginal", model="ER", root.p="maddfitz")
  recon <- rayDISC(hym.res, taxa, rate.mat=NULL, node.states="marginal", model="ER", root.p=c(0,0,0,0,1))
  
  char.recode.cor.all$comb.matrices[[char_id]]$matrix.new.recon<-recon
  #str(recon)
  #recon$solution
  #plotRECON(hym.res,recon$states,title="rayDISC Example")
  
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
#save(char.recode.cor.all, file = "char.recode.cor.all.RData")

plotRECON(hym.res, char.recode.cor.all$comb.matrices[[char_id]]$matrix.new.recon$states,title="rayDISC Example")
recon1=recon
str(recon)
char.recode.cor.all$comb.matrices[[char_id]]$matrix.new.recon$solution

char.recode.cor.all$comb.matrices

rates<-lapply(char.recode.cor.all$comb.matrices, function(x) x$matrix.new.recon$solution[!is.na(x$matrix.new.recon$solution) &
                                                                                       (x$matrix.new.recon$solution)>0 ] %>%
                unique()
)

unlist(rates)->rates
hist(rates, breaks=150)

char.recode.cor.all$comb.matrices[[char_id]]$matrix.new.recon
lapply(char.recode.cor.all$comb.matrices, function(x) x$matrix.new.recon$solution %>% nrow()) %>% unlist ->tmp

length(tmp)
abs=indep.chars[!indep.chars %in% names(tmp)]
######
# reorder matrix rows and cols
i=2
indep.chars.cor=char.recode.cor.all$comb.matrices %>% names()

for ( i in seq_along(indep.chars.cor)){
char_id=indep.chars.cor[i]
init=char.recode.cor.all$init[[char_id]]
rc=char.recode.cor.all$comb.matrices[[char_id]]$matrix.new.recon$solution

char.recode.cor.all$init.reoder[[char_id]]<-init[match(rownames(init), rownames(rc)), match(rownames(init), rownames(rc))]

}
#####

#########
# char where enumeration does not startt from 0
lapply(char.recode.cor.all$comb.matrices, function(x) rownames(x$matrix.new.recon$solution)[2]) %>% unlist ->tmp
lapply(char.recode.cor.all$comb.matrices, function(x) rownames(x$matrix.new.recon$solution)[5]) %>% unlist ->tmp
which(tmp!="4")

nm=names(char.recode.cor.all$comb.matrices)
sapply(nm, function(x) char_matrix.rem[[x]] %>% unique %>% as.character() %>% grepl("0", .) %>% any() ) %>% all()

lapply(char.recode.cor.all$comb.matrices, function(x) rownames(x$matrix.new.recon$solution)[1])->tmp
tmp
non0=which(tmp!="0") %>% names
char_id=non0[3]
char.recode.cor.all$init.reoder[[char_id]]
char.recode.cor.all$comb.matrices[[char_id]]

char_matrix.rem.tmp=char_matrix.rem
for (i in non0){
  
  char_matrix.rem[[i]]=revalue(char_matrix.rem[[i]], c("1"="0", "2"="1"))
}

for (i in non0){
  colnames(char.recode.cor.all$comb.matrices[[i]]$matrix.new.recon$solution)<-
    rownames(char.recode.cor.all$comb.matrices[[i]]$matrix.new.recon$solution)<-c("0", "1")
  
  colnames(char.recode.cor.all$init.reoder[[i]])<-
    rownames(char.recode.cor.all$init.reoder[[i]])<-c("0", "1")
  
  #char.recode.cor.all$init.reoder[[i]]
}

#write.csv(char_matrix.rem, file = "matrix_1Feb-final_sp-selected.csv")
###################

###
#char.recode$comb.matrices$`CHAR:239`
char_id="CHAR:7"
char.recode$comb.matrices[[char_id]]$matrix.new
str(char.recode$comb.matrices[[char_id]])
lapply(char.recode$comb.matrices, function(x) x$matrix.new %>% nrow()) %>% unlist ->tmp
which(tmp>=4)

templates=c("ASR_templ")
templates.fl=paste0(templates, ".Rev")
#simulations$anal_names=paste(simulations$run.ids, "-", templates.fl, sep="")
# !!! simulations$anal_names["CHAR:101"] is # 8
#simulations$anal_names=simulations$anal_names[-8]

indep.chars.cor
char.recode.cor.all$init.reoder[[char_id]]
char.recode.cor.all$comb.matrices[[char_id]]
str(char.recode.cor.all)
char.recode.cor.all$comb.matrices[1]
# add "fake" sublist to use in the follwoing function
for (i in indep.chars.cor){
  char.recode.cor.all$comb.matrices[[i]]$matrix.new<-char.recode.cor.all$init.reoder[[i]]
}
char.recode.cor.all$comb.matrices[2]

###################################
# Revsions of Analys3: recoding chars that did not work + chars with matrices dim>2
#
##
anal3.mt=init_matrices(anal3.revision, char_matrix=char_matrix.rem, diag.as = 0)

i=1
for (i in seq_along(anal3.mt)){
  char_id=names(anal3.mt[i])
char.recode.cor.all$comb.matrices[[char_id]]$matrix.new<-anal3.mt[[i]]
}

#####################
anal3.revision[anal3.revision %in% indep.chars.cor]
indep.n=indep.chars.cor[!indep.chars.cor %in% anal3.revision]
lapply(char.recode.cor.all$comb.matrices, function(x) x$matrix.new %>% nrow()) %>% unlist()
sapply(indep.n, function(x) char.recode.cor.all$comb.matrices[[x]]$matrix.new %>% nrow()==2) %>% all()
#######


##########
char_id=anal3.revision[1]
i=1
#for (i in 9:80){
#for (i in seq_along(char.revision)){

#for (i in seq_along(indep.chars.cor)){
  # write character file
#for (i in seq_along(anal3.revision)){
  #for (i in seq_along(indep.n)){
    for (i in seq_along(poly)){
  
  #char_id=indep.chars.cor[i]
    #char_id=anal3.revision[i]
    
  char_id=poly[i]
  
  char_id.n=gsub(":", "-", char_id)
  
  #taxa=cbind(char.recode$comb.matrices[[char_id]]$new.coding)
  taxa=char_matrix.rem[char_id]
  write.table(taxa, file=paste("data/", char_id.n, ".char", sep=""), quote = F, sep = " ", col.names = F)
  
  # write Rev file
  fl.in  <- readLines(paste0("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/script-final/", templates.fl))
  fl.in  <- gsub(pattern = "@analysis_name@", replace = paste(char_id.n, sep=""),
                 x = fl.in)
  #fl.in  <- gsub(pattern = "@tree_2_read@", replace = simulations$names.trees[i], x = fl.in)
  fl.in <- gsub(pattern = "@chrs_2_read@", replace = paste("data/", char_id.n, ".char", sep=""), x = fl.in)
  
  fl.in <- gsub(pattern = "@NUM_STATES@", replace = char.recode.cor.all$comb.matrices[[char_id]]$matrix.new %>% nrow(),
                x = fl.in)
  
  fl.in <- gsub(pattern = "@rate-prior@", replace = Rev_rate_prior(char_id, rateUnif=NULL, rateExp=NULL, Ml.prior=T, 
                                                                  char.recode.cor=char.recode.cor.all, char.recode=char.recode.cor.all), x = fl.in)
  
  #### Unifform prior
  #fl.in <- gsub(pattern = "@rate-prior@", replace = Rev_rate_prior(char_id, rateUnif=0.1, rateExp=NULL, Ml.prior=F, 
    #                                                               char.recode.cor=char.recode.cor.all, char.recode=char.recode.cor.all), x = fl.in)
  ####
  
  fl.in <- gsub(pattern = "@rate-moves@", replace = Rev_rate_moves(char_id, rateUnif=F, Ml.prior=T, 
                                                                  char.recode.cor=char.recode.cor.all, char.recode=char.recode.cor.all), x = fl.in)
  
  # Uniform moves
  #fl.in <- gsub(pattern = "@rate-moves@", replace = Rev_rate_moves(char_id, rateUnif=T, Ml.prior=F, 
   #                                                                char.recode.cor=char.recode.cor.all, char.recode=char.recode.cor.all), x = fl.in)
  
  
  fl.in <- gsub(pattern = "@rates2matrix@", replace = Rev_rates2Q(char_id, Ml.prior=T, 
                                                                 char.recode.cor=char.recode.cor.all, char.recode=char.recode.cor.all), x = fl.in)
  
 # Uniform
  #fl.in <- gsub(pattern = "@rates2matrix@", replace = Rev_rates2Q(char_id, Ml.prior=F, 
   #                                            char.recode.cor=char.recode.cor.all, char.recode=char.recode.cor.all), x = fl.in)

  
  
  
   fl.in <- gsub(pattern = "@root_freq@", replace = paste0("simplex(",
                                                          paste(rep(1, char.recode.cor.all$comb.matrices[[char_id]]$matrix.new %>% nrow()), collapse = ", "), ")") , 
                x = fl.in)
  
  
  
  #writeLines(fl.in , con=simulations$anal_names[i])
  cat(file=paste(char_id.n, ".Rev", sep=""), sep="\n", fl.in)
  
  
} # end loop
##############


# Bash scripts

### 
# lapply(char.recode.cor.all$comb.matrices, function(x) x$matrix.new %>% nrow()) %>% unlist ->tmp
# which(tmp<6) %>% names(.) ->chrsL6st
# diff=chrsL6st[!chrsL6st %in% simulations$anal_names[1:29]]
# diff=diff[-1]

sim.nodep=indep.n
mt.nodep=matrix(sim.nodep, ncol = 6, nrow = 43, byrow = F)
mt.nodep[ matrix(duplicated(c(mt.nodep)), ncol = 6, nrow = 43)  ] <-NA
mt.nodep=gsub(":", "-", mt.nodep)
#####

sim.nodep=poly
mt.nodep=matrix(sim.nodep, ncol = 1, nrow = 6, byrow = F)
mt.nodep[ matrix(duplicated(c(mt.nodep)), ncol = 1, nrow = 6)  ] <-NA
mt.nodep=gsub(":", "-", mt.nodep)
#####

sim.nodep=anal3.revision
mt.nodep=matrix(sim.nodep, ncol = 2, nrow =28, byrow = F)
mt.nodep[ matrix(duplicated(c(mt.nodep)), ncol = 2, nrow = 28)  ] <-NA
mt.nodep=gsub(":", "-", mt.nodep)
#####

# make bash script
bash.name=c("Abies_run2.sh")
cat("#!/bin/bash\n", file=bash.name)

i=1
for (i in 1:nrow(mt.nodep)){
  #for (i in seq_along(simulations$anal_names)){
  str.in<-paste("./rb ", mt.nodep[i,], ".Rev > ", mt.nodep[i,], ".screen", sep="") %>% paste(., collapse=" && ") %>%
    paste0(., " &")
  
  cat(str.in, file=bash.name, append=T)
  cat("\n", file=bash.name, append=T)
}


#####################
