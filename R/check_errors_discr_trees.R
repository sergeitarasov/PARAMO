char.rep=read.csv("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/dataset/Data_final/character_table.csv", stringsAsFactors =  FALSE, header=T)


c=char.rep[,1]
length(c)

setwd("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/Run1_all/stm_R500_discr_new")
#c<-sub(":", "-",
#       get_descendants_chars(ONT, annotations="manual", F1$Level_1_ids[4] )  )

# dir to read
#dir= ("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/RevBayes/Run1_all/stm_R500/")

i=1
errors<-c()
ntrees<-500
#tr<-vector("list", ntrees)
rm(tr)

# tr[[1]]$maps %>% unlist %>% is.numeric()
# c(1, 0.1, NA) %>% is.na() %>% any
# (c(1, 0.1, Inf)==Inf) %>% any
# 
# tr[[1]]$maps[[1]] %>% unlist %>% names %>% as.numeric()

i=2
#for (i in 1:length(c)){
  for (i in 40:length(c)){
  
  fl<-paste0(c[i], "_", c(1:ntrees), ".rds")  
  #j=1
  for (j in 1:ntrees) {
    
    print(paste0("Reading ", c[i], "   ", fl[j]))
    
    con<-unz(paste0(c[i], ".zip"), filename = fl[j])
    con2 <- gzcon(con)
    tr <- readRDS(con2)
    close(con)
    
    tr$maps %>% unlist->MM
    names(MM)->MMn
    
    c( 
     !(MM %>% is.na() %>% any),
     (MM > 0) %>% all,
     (MM!=Inf) %>% all,
     
     !(MMn %>% is.na() %>% any),
     (MMn!=Inf) %>% all
    ) ->out
    
    out<-!out
    
    if (any(out)){
      errors<-rbind(errors, c(i,j, !out) )
    }
  }
}
 

errors
#
######## only 443 trees in char377
[1] "Reading CHAR-377   CHAR-377_19.rds"
Error in readRDS(con2) : cannot open the connection
In addition: Warning message:
  In readRDS(con2) :
  cannot locate file 'CHAR-377_19.rds' in zip file 'CHAR-377.zip'
#############