
setwd("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/script-final")
# write Rev scripts using templates
# make up file names for rev

##
#char.recode.cor=get_graph_matrix_any(g, char_matrix.rem, dep.tb, new.polymorph="&")
#char.recode=get_graph_matrix_any(g, char_matrix.rem, dep.tb, new.polymorph=" ",  Lnew.polymorph=NULL, Rnew.polymorph=NULL)
######
sapply(dep.chars, function(x) char.recode.cor$comb.matrices[[x]]$matrix.new.recon$solution %>% nrow()) %>% unlist ->dd
sapply("CHAR:15", function(x) char.recode.cor$comb.matrices[[x]]$matrix.new.recon$solution %>% nrow()) %>% unlist
unique(tb[,1])[!unique(tb[,1]) %in% dep.chars]
dep.chars[!dep.chars %in% names(dd)]

#char.recode$comb.matrices$`CHAR:239`
contr.chars
char_id="CHAR:17"
char.recode$comb.matrices[[char_id]]$matrix.new
char.recode$comb.matrices[[char_id]]$matrix
char.recode$binary.matrices[[char_id]]
comb.matrices[[char_id]]
char.recode.cor$comb.matrices[[char_id]]$matrix.new.recon$solution

str(char.recode$comb.matrices[[char_id]])
lapply(char.recode$comb.matrices, function(x) x$matrix.new %>% nrow()) %>% unlist ->tmp
which(tmp>=4)

templates=c("ASR_templ")
templates.fl=paste0(templates, ".Rev")
#simulations$anal_names=paste(simulations$run.ids, "-", templates.fl, sep="")
# !!! simulations$anal_names["CHAR:101"] is # 8
#simulations$anal_names=simulations$anal_names[-8]
i=2
#for (i in 9:80){
for (i in seq_along(dep.chars)){
#for (i in seq_along(contr.chars)){
  # write character file
  #char_id=contr.chars[i]
  tryCatch({
  
  #char_id=dep.chars[i]
  char_id="CHAR:17"
  
  char_id.n=gsub(":", "-", char_id)
  
    taxa=cbind(char.recode$comb.matrices[[char_id]]$new.coding)
    write.table(taxa, file=paste("data/", char_id.n, ".char", sep=""), quote = F, sep = " ", col.names = F)

  # write Rev file
  fl.in  <- readLines(paste0("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/script-final/", templates.fl))
  fl.in  <- gsub(pattern = "@analysis_name@", replace = paste(char_id.n, sep=""),
                 x = fl.in)
  #fl.in  <- gsub(pattern = "@tree_2_read@", replace = simulations$names.trees[i], x = fl.in)
  fl.in <- gsub(pattern = "@chrs_2_read@", replace = paste("data/", char_id.n, ".char", sep=""), x = fl.in)

  fl.in <- gsub(pattern = "@NUM_STATES@", replace = char.recode$comb.matrices[[char_id]]$matrix.new %>% nrow(),
                x = fl.in)
  
  fl.in <- gsub(pattern = "@rate-prior@", replace = Rev_rate_prior(char_id, rateUnif=NULL, rateExp=NULL, Ml.prior=T, 
                                                         char.recode.cor=char.recode.cor, char.recode=char.recode), x = fl.in)
  
  fl.in <- gsub(pattern = "@rate-moves@", replace = Rev_rate_moves(char_id, rateUnif=F, Ml.prior=T, 
                                                                   char.recode.cor=char.recode.cor, char.recode=char.recode), x = fl.in)

  
  fl.in <- gsub(pattern = "@rates2matrix@", replace = Rev_rates2Q(char_id, Ml.prior=T, 
                                                                  char.recode.cor=char.recode.cor, char.recode=char.recode), x = fl.in)
  
  fl.in <- gsub(pattern = "@root_freq@", replace = paste0("simplex(",
                  paste(rep(1, char.recode$comb.matrices[[char_id]]$matrix.new %>% nrow()), collapse = ", "), ")") , 
                x = fl.in)
  
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})  
  
  #writeLines(fl.in , con=simulations$anal_names[i])
  cat(file=paste(char_id.n, ".Rev", sep=""), sep="\n", fl.in)

  
} # end loop
##############

# make bash script
bash.name=c("Temp_1-1.sh")

cat("#!/bin/bash\n", file=bash.name)

i=1

for (i in 11:20){
#for (i in seq_along(simulations$anal_names)){
  char_id=simulations$anal_names[i]
  char_id.n=gsub(":", "-", char_id)
  
  str.in=paste("./rb ", paste0(char_id.n, ".Rev"), " > ", char_id.n, ".screen", " &", sep="")
  cat(str.in, file=bash.name, append=T)
  cat("\n", file=bash.name, append=T)
}
#########################
# New code for 2018

### Abies_st4_35runs.sh; 35 runs from dependecy table; chars with states < 6
lapply(char.recode$comb.matrices, function(x) x$matrix.new %>% nrow()) %>% unlist ->tmp
which(tmp<6) %>% names(.) ->chrsL6st
diff=chrsL6st[!chrsL6st %in% simulations$anal_names[1:29]]
diff=diff[-1]
sim.nodep=diff
mt.nodep=matrix(sim.nodep, ncol = 1, nrow = 35, byrow = F)
mt.nodep[ matrix(duplicated(c(mt.nodep)), ncol = 1, nrow = 35)  ] <-NA
mt.nodep=gsub(":", "-", mt.nodep)
#####



### Alcies_states6_15anal.sh; Alcies run for slow datase which have 6 and more states; altogether 15 runs; all these from 
#dependency table
lapply(char.recode$comb.matrices, function(x) x$matrix.new %>% nrow()) %>% unlist ->tmp
which(tmp>=6) %>% names(.) ->chrsL6st
sim.nodep=chrsL6st
mt.nodep=matrix(sim.nodep, ncol = 2, nrow = 8, byrow = F)
mt.nodep[ matrix(duplicated(c(mt.nodep)), ncol = 2, nrow = 8)  ] <-NA
mt.nodep=gsub(":", "-", mt.nodep)
#####

#### this is Abies 1-29 run with no dependecies
sim.nodep=simulations$anal_names[1:29]
mt.nodep=matrix(sim.nodep, ncol = 3, nrow = 10, byrow = F)
mt.nodep[ matrix(duplicated(c(mt.nodep)), ncol = 3, nrow = 10)  ] <-NA
mt.nodep=gsub(":", "-", mt.nodep)
############3

#### this is Rerun of the majoprity chars from Abies 1-29; this chars are controlling
#sim.nodep=contr.chars
sim.nodep=rerun.analys_2
mt.nodep=matrix(sim.nodep, ncol = 1, nrow = 10, byrow = F)
mt.nodep[ matrix(duplicated(c(mt.nodep)), ncol = 1, nrow = 10)  ] <-NA
mt.nodep=gsub(":", "-", mt.nodep)
############3

#### this is Rerun of the majoprity chars from Abies 1-29; this chars are controlling
#sim.nodep=dep.chars
sim.nodep=contr.chars
#sim.nodep=rerun.analys_2
mt.nodep=matrix(sim.nodep, ncol = 1, nrow = 29, byrow = F)
mt.nodep[ matrix(duplicated(c(mt.nodep)), ncol = 1, nrow = 29)  ] <-NA
mt.nodep=gsub(":", "-", mt.nodep)
############3


# make bash script
bash.name=c("Abies_contr_chars.sh")
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

lapply(char.recode$comb.matrices, function(x) rownames(x$matrix.new)[1]) %>% unlist ->tmp
tmp=names(tmp[tmp!=0])
rerun=tmp
contr.chars=char.recode$vertex.hier[1:29] %>% names()

contr.chars %in% tmp
tmp %in% contr.chars

contr.chars=contr.chars[contr.chars %in% simulations$anal_names]
#non0=which(tmp!="0") %>% names

char.recode$comb.matrices[[char_id]]$matrix.new
######
# reorder matrix controlling chars rows and cols



i=19
indep.chars.cor=char.recode.cor.all$comb.matrices %>% names()

for ( i in seq_along(contr.chars)){
  char_id=contr.chars[i]
  init=char.recode$comb.matrices[[char_id]]$matrix.new
  rc=char.recode.cor$comb.matrices[[char_id]]$matrix.new.recon$solution
  
  char.recode$comb.matrices[[char_id]]$matrix.new<-init[match(rownames(rc), rownames(init)), match(rownames(rc), rownames(init))]
  
  
}
#####
char.recode$comb.matrices[["CHAR:272"]]$matrix.new
char.recode.cor$comb.matrices[["CHAR:272"]]$matrix.new.recon$solution
#save(char.recode, file="char.recode-v2.RData")

# Making up char report
contr.chars

char.rep=rbind(
cbind(simulations$anal_names[1:29], "controll-char"),
cbind(simulations$anal_names[30:79], "depen-char"),
cbind(indep.chars.cor, "non-dep.tb")
)

which(char.rep[,1]=="CHAR:193")
char.rep[c(47,48),2]<-"non-dep.tb"
#write.csv(char.rep, file="Chars_used_in_RevBayes.csv")
