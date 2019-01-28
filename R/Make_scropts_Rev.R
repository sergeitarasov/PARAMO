library(reshape)

char.recode=get_graph_matrix_any(g1, char_matrix.rem, dep.tb, new.polymorph=" ",  Lnew.polymorph=NULL, Rnew.polymorph=NULL)
# save data
taxa=cbind(char.recode$comb.matrices$`CHAR:239`$new.coding)
write.table(taxa, file="Hym-char239-1.nex", quote = F, sep = " ", col.names = F)

char.recode$comb.matrices$`CHAR:239`$new.coding


char.recode.cor.all$comb.matrices[["CHAR:43"]]$matrix.new.recon$solution
char.recode.cor.all$comb.matrices[["CHAR:43"]]$matrix.new

char.recode$comb.matrices$`CHAR:239`$matrix.new
mt<-char.recode.cor.all$comb.matrices[[char_id]]$matrix.new
mt.r=char.recode.cor.all$comb.matrices[[char_id]]$matrix.new.recon$solution
##############3
Rev_rates2Q<-function(char_id, Ml.prior=F, char.recode.cor=char.recode.cor, char.recode=char.recode){
mt<-char.recode$comb.matrices[[char_id]]$matrix.new
mt.r=char.recode.cor$comb.matrices[[char_id]]$matrix.new.recon$solution

#if (!rateUnif){

if (Ml.prior){
  mt.r[is.na(mt.r)]<-0
  bb.mt.r<-melt(mt.r)
  bb.mt.r[,1]=bb.mt.r[,1]+1
  bb.mt.r[,2]=bb.mt.r[,2]+1
  
  bb.mt<-melt(mt)
  bb.mt[,1]=bb.mt[,1]+1
  bb.mt[,2]=bb.mt[,2]+1
  
  bb.mt= bb.mt[bb.mt[,3]>0 & bb.mt.r[,3]>0,]
  
  return(
    apply(bb.mt, 1, function(x) paste("rates[", x[1], "][", x[2], "]:=r", x[3], sep="")) %>%
      paste(., collapse = "\n")
  )

} else {

bb.mt<-melt(mt)
bb.mt[,1]=bb.mt[,1]+1
bb.mt[,2]=bb.mt[,2]+1

bb.mt=bb.mt[bb.mt[,3]>0,]
return(
apply(bb.mt, 1, function(x) paste("rates[", x[1], "][", x[2], "]:=r", x[3], sep="")) %>%
paste(., collapse = "\n")
)
}


}

###########################################



##################################
##############3
log(2)/0.0002437692
1/0.0002437692
char.recode$comb.matrices[[char_id]]$matrix.new %>% nrow()
Rev_rate_prior<-function(char_id, rateUnif=NULL, rateExp=NULL, Ml.prior=F, char.recode.cor=char.recode.cor, char.recode=char.recode){
  mt<-char.recode$comb.matrices[[char_id]]$matrix.new
  mt.r=char.recode.cor$comb.matrices[[char_id]]$matrix.new.recon$solution
  
  if (Ml.prior){
    mt.r[is.na(mt.r)]<-0
    rt<-mt[(c(mt)>0) & (c(mt.r)>0)] %>% unique()
    #rt.mt<-c(mt)[c(mt)>0] %>% unique()
    
    return(
    sapply(rt, function(x) paste("r", x, " ~ dnExp(", round(log(2)/mt.r[mt==x][1], 2), ")", sep=""))%>%
      paste(., collapse = "\n")
    )
    
  } #else {
  
  
  if (!is.null(rateExp)){
  rt<-c(mt)[c(mt)>0] %>% unique()
  return(
  sapply(rt, function(x) paste("r", x, " ~ dnExp(", rateExp, ")", sep=""))%>%
    paste(., collapse = "\n")
  )
  }
  
  if (!is.null(rateUnif)){
    rt<-c(mt)[c(mt)>0] %>% unique()
    sapply(rt, function(x) paste("r", x, " ~ dnUniform(", 0, ", ", rateUnif, ")", sep=""))%>%
      paste(., collapse = "\n")
  }

}
##################################

##############3
Rev_rate_moves<-function(char_id, rateUnif=F, Ml.prior=F, char.recode.cor=char.recode.cor, char.recode=char.recode){
  mt<-char.recode$comb.matrices[[char_id]]$matrix.new
  #rt<-c(mt)[c(mt)>0] %>% unique()
  mt.r=char.recode.cor$comb.matrices[[char_id]]$matrix.new.recon$solution
  
  if (!rateUnif){
  
  if (Ml.prior){
    mt.r[is.na(mt.r)]<-0
    rt<-mt[(c(mt)>0) & (c(mt.r)>0)] %>% unique()
    #rt.mt<-c(mt)[c(mt)>0] %>% unique()
    
    return(
    sapply(rt, function(x) paste("moves[++mvi] = mvScale(",
                                 "r", x, ", lambda=1, tune=true, weight=2)",
                                 sep=""))%>%
      paste(., collapse = "\n")
    )
    
  } else {
    
    rt<-c(mt)[c(mt)>0] %>% unique()
    
    return(
  
  sapply(rt, function(x) paste("moves[++mvi] = mvScale(",
                               "r", x, ", lambda=1, tune=true, weight=2)",
                               sep=""))%>%
    paste(., collapse = "\n")
    )
  }
    
  }
  
  if (rateUnif){
    rt<-c(mt)[c(mt)>0] %>% unique()
    
    sapply(rt, function(x) paste("moves[++mvi] = mvSlide(",
                                 "r", x, ", delta=1, tune=true, weight=2)",
                                 sep=""))%>%
      paste(., collapse = "\n")
     
  }
  
}
##################################
char.recode.cor$comb.matrices[[char_id]]$matrix.new.recon
char_id="CHAR:17"
x=c()
Rev_MLprior<-function(char_id, char.recode.cor=char.recode.cor){
  mt.r=char.recode.cor$comb.matrices[[char_id]]$matrix.new.recon$solution
  mt=char.recode.cor$comb.matrices[[char_id]]$matrix.new
  
  # both intial rates and recon rates have to be higher than zero
  mt.r[is.na(mt.r)]<-0
  rt<-mt[(c(mt)>0) & (c(mt.r)>0)] %>% unique()
  #rt.mt<-c(mt)[c(mt)>0] %>% unique()
  
  sapply(rt, function(x) paste("r", x, " ~ dnExp(", round(1/mt.r[mt==x][1], 2), ")", sep=""))%>%
    paste(., collapse = "\n")
}


###########################
#Get ML estimates
char.recode=get_graph_matrix_any(g, char_matrix.rem, dep.tb, new.polymorph=" ",  Lnew.polymorph=NULL, Rnew.polymorph=NULL)

# get resolved tree
hym.res=read.tree(file="Hymenoptera_br_resolved.tre")

char.recode.cor=get_graph_matrix_any(g, char_matrix.rem, dep.tb, new.polymorph="&",  Lnew.polymorph=NULL, Rnew.polymorph=NULL)

# get names of chars
simulations<-list()
char.recode.cor$comb.matrices %>% names() -> simulations$anal_names

#### separating dependent chars
dep.chars=simulations$anal_names[!simulations$anal_names %in% contr.chars]


# ML estimates
char.revision=c("CHAR:34", "CHAR:35", "CHAR:260")
char.recode.cor$binary.matrices[[char_id]]
i=1
"CHAR:17"  "CHAR:399"
char_id="CHAR:17"
contr.chars
  
for (i in seq_along(dep.chars)){
#for (i in seq_along(simulations$anal_names)){
 # for (i in seq_along(contr.chars)){
  tryCatch({
  
  #char_id="CHAR:201"
  #char_id=contr.chars[i]
  
print(paste0("Working on ", char_id))
  
mt=char.recode.cor$comb.matrices[[char_id]]$matrix.new
#char.recode.cor$comb.matrices[[char_id]]$matrix.new<-mt
char.recode$comb.matrices[[char_id]]$matrix.new<-mt
#mt=char.recode.cor$binary.matrices[[char_id]]
#mt=char.recode.cor.all$comb.matrices[[char_id]]$matrix.new

#rownames(mt)<-colnames(mt)<-c(0,1)
#char.recode.cor$comb.matrices[[char_id]]$matrix.new<-mt
####
mt[mt>20]<-23
mt[mt==19]<-NA

mt[3,1]<-NA
mt[2,3]<-NA
###
mt[mt==0]<-NA

taxa=char.recode.cor$comb.matrices[[char_id]]$new.coding
taxa=as.data.frame(cbind(names(taxa), taxa),stringsAsFactors =F)

##############
#taxa=char_matrix.rem[[char_id]] %>% as.character() 
#names(taxa)=rownames(char_matrix.rem)

#taxa=revalue(taxa, c("-"="?", "1"="0", "2"="1"))
#char.recode.cor$comb.matrices[[char_id]]$new.coding<-taxa

#taxa=as.data.frame(cbind(names(taxa), taxa),stringsAsFactors =F)
#taxa[,2] %>% unique()
###########
recon=c()

#recon <- rayDISC(hym.res, taxa, rate.mat=mt, node.states="marginal", model="SYM", root.p=NULL, verbose=T)
#recon <- rayDISC(hym.res, taxa, rate.mat=NULL, node.states="marginal", model="ARD", root.p=NULL,diagn=T, verbose=T)
recon <- rayDISC(hym.res, taxa, rate.mat=mt, node.states="marginal", model="ARD", root.p=NULL)

char.recode.cor$comb.matrices[[char_id]]$matrix.new.recon<-recon
#str(recon)
#recon$solution
#plotRECON(hym.res,recon$states,title="rayDISC Example")

  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}
########################
x=c()
str(recon)
recon$solution.se
char.recode.cor$comb.matrices[["CHAR:17"]]$matrix.new.recon$solution
recon1=recon
char.recode.cor$comb.matrices[[char_id]]$matrix.new<-mt
#setwd("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/dataset/Data_final")
#save(char.recode.cor, file="char.recode.cor.RData")
char.recode$comb.matrices[[char_id]]$matrix.new %>% nrow()

tmp=lapply(char.recode.cor$comb.matrices, function(x) x$matrix.new %>% nrow() %>% unique() )
tmp=unlist(tmp)
which(tmp==1)

rates<-lapply(char.recode.cor$comb.matrices, function(x) x$matrix.new.recon$solution[!is.na(x$matrix.new.recon$solution) &
                                                                                (x$matrix.new.recon$solution)>0 ] %>%
                unique()
              )

char.recode.cor$comb.matrices$`CHAR:376`$matrix.new.recon$solution

rates.v=unlist(rates)
hist(rates.v, breaks=150)
mean(rates.v)
var(rates.v)
sd(rates.v)
median(rates.v)
min(rates.v)
max(rates.v)

str(recon)
recon$solution.se
char.recode.cor$comb.matrices$`CHAR:110`$matrix.new.recon$solution
