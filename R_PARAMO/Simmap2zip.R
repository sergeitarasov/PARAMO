library("phytools")
c=paste0("C", AN$CHAR_ID2)

setwd("~/Documents/Recon-Anc_Anat/Supplementary_materials/STEP_4/RevBayes/output")
# dir to write files
dir= ("~/Documents/Recon-Anc_Anat/Supplementary_materials/STEP_4/RevBayes/Discr_maps/")

# read a sample of 100 maps from .stm files and save them in the poper format .stmR
i=1
for (i in 1:length(c))
{
  tree<-read_Simmap_Rev(paste0(c[i], ".stm"),
                        start=400, end=500,
                        save = NULL) %>% read.simmap(text=., format="phylip")
  
  
  write.simmap(tree, file=paste0(dir, c[i], ".stmR"))
}
##########

# Read stmR and save each map as a separate rds file; in turn all rds file for a chracter are stored 
# in zip archive

tree.tmp.final<-read.tree("Hymenoptera_br_resolved.tre")
setwd("~/Documents/Recon-Anc_Anat/Supplementary_materials/STEP_4/RevBayes/Discr_maps")


for (i in 1:length(c))
  { 
  
  # read in undesritezed trees
  print(paste0("Reading ", c[i]))
  sim=read.simmap(file=paste0(dir, c[i], ".stmR"), format="phylip")
  
  # descritize trees by looping over sample and saving as rds
  #j=1
  
  
  for (j in 1:length(sim)){
    tryCatch({
      
      print(paste0("Descritizing tree ", j))
      
      ## errors with na
      
      ##
      
      ##### make trees equal with template
      sim.d<-make_tree_eq(tree.tmp.final, sim[[j]], round=5)
      ###
      
      #sim.d<-discr_Simmap_all(sim[[j]], 1000)
      sim.d<-discr_Simmap_all(sim.d, 1000)
      
      saveRDS(sim.d, file =  paste0(dir,c[i], "_", j, ".rds") )
      
    }, error=function(e){
      cat("ERROR :",conditionMessage(e), "\n")
      #errors<-rbind(errors, c(ii,jj))
    }  )
    
  } 
  
  # putting rds files into archive
  files<-paste0( c[i], "_", c(1:length(sim)), ".rds")
  zip(paste0(c[i], ".zip"), files=files)
  file.remove(files)
  
}
###################
# showConnections (all=T)
# closeAllConnections()
#########################

