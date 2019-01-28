# Creating data files to run with RevBayes
str(MT)

# creating chracter files
#setwd("~/Documents/Recon-Anc_Anat/Supplementary_materials/STEP_4/RevBayes/data")
for (i in 1:ncol(MT))
{
  C.rev<-MT[,i]
  C.rev<-gsub("&", " ", C.rev)
  
  out<-cbind(rownames(MT), C.rev)
  write.table(file=paste0(colnames(MT[i]), ".char"), out, quote=F, sep=" ", 
              row.names=F, col.names=F)
}



# write Rev file for two-state characters
setwd("~/Documents/Recon-Anc_Anat/Supplementary_materials/STEP_4/RevBayes/")

fl.in  <- readLines("PARAMO2_templ.Rev")

for (i in 1:ncol(MT))
{
  fl.in  <- readLines("PARAMO2_templ.Rev")
  fl.in  <- gsub(pattern = "@analysis_name@", replace = paste0(colnames(MT[i])),
                 x = fl.in)
  fl.in <- gsub(pattern = "@chrs_2_read@", replace = paste0("data/", colnames(MT[i]), ".char"), x = fl.in)
  
  cat(file=paste0(colnames(MT[i]), ".Rev"), sep="\n", fl.in)

}

# write Rev file for dependent foru-state character C3-2
setwd("~/Documents/Recon-Anc_Anat/Supplementary_materials")
dd <- getwd()
knitr::opts_knit$set(root.dir= paste(dd,'/../../')) 

# I use precooked set of functions for constracting SMM
source("STEP_4/SMM_functions.R")

###################################
# same SMMs as for the tail color problem
###################################
char.state<-c("a", "p")
rate.param<-c(1, 1)
TL<-init_char_matrix(char.state, rate.param, diag.as=0)
char.state<-c("r", "b")
rate.param<-c(1, 1)
COL<-init_char_matrix(char.state, rate.param, diag.as=0)

#SMM-ind
TC.ind<-comb2matrices(TL, COL, controlling.state=NULL, name.sep="", diag.as="")
TC.ind
in.rev<-Mk_Rev(TC.ind)
cat(in.rev) # COPY the output and insert in Rev template PARAMO2_templ.Rev
#cat(in.rev, file="STEP_4/input_Rev.txt") # or save this outout to a file and then copy to Rev template



