# Check RevBayes results



matches=c("completed", "Error:\tProblem", "|")
names(matches)<-c("Done", "Error", "Processing")
file="CHAR-24.screen"
Rev_check_screen(file="CHAR-500.screen", matches)
file.exists("CHAR-500.screen")
  
Rev_check_screen<-function(file, matches){
#  tryCatch({

  if (file.exists(file)){
print(paste0("Reading ", file))
     
text <- scan(file, sep = "\n", what = "character")
tail=length(text)
line=text[tail-1]

line=strsplit(line, " ")[[1]]
line=line[!line==""]
out=names(matches)[matches %in% line]
#sapply(line, function(x) grep(x, matches)) %>% unlist->result
#out=names(matches)[result[1]]
return(out)

#}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

} else return("No file found")

}


as.file=sapply(char.rep[,1], function(x) gsub(":", "-", x)) %>% setNames(., NULL)
char.rep=cbind(as.file, char.rep)
char.rep %>% nrow

Rev.result=sapply(char.rep[,1], function(x) Rev_check_screen(file=paste0(x, ".screen"), matches))
Rev.result[Rev.result=="Error"] %>% length
Rev.result[Rev.result=="Processing"] %>% length
Rev.result[Rev.result=="Done"] %>% length
tmp=cbind(char.rep, Rev.result)

unlist(Rev.result[Rev.result=="Processing"] %>% names()) %>% gsub("-", ":", .)->a
length(a)
anal3.revision [!anal3.revision %in% a]

tmp[tmp[,4]=="Error",]->err
lapply(char.recode.cor.all$comb.matrices, function(x) nrow(x$matrix.new)) %>% unlist->nr
c(names(nr[nr>2]), err[,2]) %>% unique() ->anal3.revision
err[,2]
#rerun.analys_2=tmp[tmp[,4]=="Error",2]

tmp[tmp[,4]=="Processing",]
char_id="CHAR:8"
char.recode.cor$comb.matrices[[char_id]]$matrix.new 
char.recode.cor$comb.matrices[[char_id]]$matrix.new.recon$solution
#char.recode$comb.matrices[[char_id]]$new.coding
char.recode.cor$comb.matrices[[char_id]]$matrix
############
x1=melt(recon$solution)
x2=melt(char.recode.cor$comb.matrices[[char_id]]$matrix.new.recon$solution)
x3=melt(char.recode.cor$comb.matrices[[char_id]]$matrix.new)
comp=cbind(x3, x2[,3], x1[,3])
comp
recon


###################################
#
# Checking convergence
#
##################################
library(coda) 

mcmc_check_batch<-function(analyses, dir1, dir2){
  out<-list()
  for (i in seq_along(analyses)){
    print(paste0("Working on ", analyses[i]))
    analyses[i] %>%  gsub(":", "-", .) %>%  paste0(., ".log") %>% mcmc_check(., dir1, dir2) ->R
    out[[analyses[i]]]<- R
  }
  return(out)
}
########################


mcmc_check<-function(file, dir1, dir2, thin=500){
  chain1=read.table(paste0(dir1,"/", file), header=TRUE, sep="\t", stringsAsFactors=FALSE)
  chain1=cbind(chain1[,5:ncol(chain1)])
  chain1=mcmc(chain1, thin=thin)
  
  chain2=read.table(paste0(dir2,"/", file), header=TRUE, sep="\t", stringsAsFactors=FALSE)
  chain2=cbind(chain2[,5:ncol(chain2)])
  chain2=mcmc(chain2, thin=thin)
  
  combinedchains<-NULL
  tryCatch({
    combinedchains = mcmc.list(chain1, chain2)
    return(gelman.diag(combinedchains))
  }, error=function(e){return(list(mpsrf=101, psrf=matrix(c(101, 101), nrow=1, ncol=2 ) ))})
  #}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}
################################

#######
# anal3.revision
# all converged
####
file="CHAR-306.log"
dir1="/home/tarasov/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/analys_3_revision/Empirical_prior/output"
dir2="/home/tarasov/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/analys_3_revision/Empirical_prior_run2/output"
mcmc_check(file, dir1, dir2)

anal3.revision

ss=mcmc_check_batch(anal3.revision, dir1, dir2)

#check mpsrf - multivariate pot red factor
lapply(ss, function(x) x$mpsrf) %>% unlist %>% min
lapply(ss, function(x) x$mpsrf<1.1) %>% unlist->mp
names(mp)[!mp]

#check psrf - multivariate pot red factor
lapply(ss, function(x) all(x$psrf[,1]<1.1)) %>% unlist->pp
length(pp)
names(pp)[!pp]

#[1] "CHAR:58"  "CHAR:177" "CHAR:306"
#####################################


#######
# contr.chars
# all converged
####

dir1="/home/tarasov/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/analys_5_DepChars/contr_chars/output"
dir2="/home/tarasov/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/analys_5_DepChars/contr_chars_run2/output"

file="CHAR-14.log"
mcmc_check(file, dir1, dir2)

contr.chars

ss=mcmc_check_batch(contr.chars, dir1, dir2)

#check mpsrf - multivariate pot red factor
lapply(ss, function(x) x$mpsrf) %>% unlist %>% min
lapply(ss, function(x) x$mpsrf<1.1) %>% unlist->mp
names(mp)[!mp]

#check psrf - multivariate pot red factor
lapply(ss, function(x) all(x$psrf[,1]<1.1)) %>% unlist->pp
names(pp)[!pp]



#####################################


#######
# dep.chars
# Error: CHAR-19 CHAR-34 CHAR-24 

Rev.result=sapply(char.rep[,1], function(x) Rev_check_screen(file=paste0(x, ".screen"), matches))
Rev.result[Rev.result=="Error"] %>% length
Rev.result[Rev.result=="Processing"] %>% length
Rev.result[Rev.result=="Done"] %>% length

dir1="/home/tarasov/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/analys_5_DepChars/dep_chars/output"
dir2="/home/tarasov/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/analys_5_DepChars/dep_chars_run2/output"

#dep.chars1<-dep.chars
which(dep.chars1=="CHAR:17")
which(dep.chars1=="CHAR:399")
dep.chars1<-dep.chars1[-c(2,47)]
ss=mcmc_check_batch(dep.chars1, dir1, dir2)

#check mpsrf - multivariate pot red factor
lapply(ss, function(x) x$mpsrf) %>% unlist %>% min
lapply(ss, function(x) x$mpsrf<1.1) %>% unlist->mp
names(mp)[!mp]

#check psrf - multivariate pot red factor
lapply(ss, function(x) all(x$psrf[,1]<1.1)) %>% unlist->pp
names(pp)[!pp]

#########
#######
# add.chars
# Done, all converged


add.chars=c(poly, "CHAR:17", "CHAR:102", "CHAR:103", "CHAR:201", "CHAR:399")
Rev.result=sapply(add.chars %>%  gsub(":", "-", .), function(x) Rev_check_screen(file=paste0(x, ".screen"), matches))

dir1="/home/tarasov/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/additional_chars/add_chars/output"
dir2="//home/tarasov/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/additional_chars/add_chars_run2/output"


ss=mcmc_check_batch(add.chars, dir1, dir2)

#check mpsrf - multivariate pot red factor
lapply(ss, function(x) x$mpsrf) %>% unlist %>% min
lapply(ss, function(x) x$mpsrf<1.1) %>% unlist->mp
names(mp)[!mp]

#check psrf - multivariate pot red factor
lapply(ss, function(x) all(x$psrf[,1]<1.1)) %>% unlist->pp
names(pp)[!pp]

#######
# analr3_revision nondep chars 2 states
# 
# all converged

dir1="/home/tarasov/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/analys_3_revision/nondep-chars_2states/output"
dir2="/home/tarasov/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/analys_3_revision/nondep-chars_2states_run2/output"

indep.n
ss=mcmc_check_batch(indep.n, dir1, dir2)

#check mpsrf - multivariate pot red factor
lapply(ss, function(x) x$mpsrf) %>% unlist %>% min
lapply(ss, function(x) x$mpsrf<1.1) %>% unlist->mp
names(mp)[!mp]

#check psrf - multivariate pot red factor
lapply(ss, function(x) all(x$psrf[,1]<1.1)) %>% unlist->pp
names(pp)[!pp]

# all converged
#################

fl=list.files()
strsplit(fl, "-") %>% lapply(., function(x) x[2]) %>% unlist %>% strsplit(., "\\.") %>% lapply(., function(x) x[1]) %>% unlist->fl.id
fl.id=as.numeric(fl.id)
c(1:401)[!c(1:401) %in% fl.id]
