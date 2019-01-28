library("tibble")
library("dplyr")

# objects
# H annotation of chars with HAO
# F list of focal terms
# GL graphic links


# H obj
#setwd("/Volumes/easystore/NIMBIOS_comp/my-papers-2017/Onto-Phylo/onto_phylo/data/working/dataset/Data_final")
setwd("~/my-papers-2017/Onto-Phylo/onto_phylo/data/working/dataset/Data_final")
#H<-read.csv("finest_HAO_shortlist3a.csv", header=T,  stringsAsFactors = FALSE, na.strings = c("", NA))
H<-read.csv("finest_HAO.csv", header=T,  stringsAsFactors = FALSE, na.strings = c("", NA))
str(H)
H<-as_tibble(H)

# H<-add_column(H, n1=get_onto_id(H$correct, HAO, names = F))
# saving th H objects
# write.csv(H, "finest_HAO.csv")

# getting ids for names
# checking if we can identify all terms
get_onto_id(H$correct, HAO, names = F) %>% length
get_onto_id(H$X2, HAO, names = F) %>% .[!is.na(.)] %>% length
H$X2 %>% .[!is.na(.)]  %>% length
get_onto_id(H$X3, HAO, names = F)

# merging
H<-add_column(H,
           n1=get_onto_id(H$X1, HAO, names = F),
           n2=get_onto_id(H$X2, HAO, names = F),
           n3=get_onto_id(H$X3, HAO, names = F)
           )



#H$char_id<-paste0("CHAR:", H$char_id)
#H[,c(2, 10:12)]
#H<-set_names(table2list(H[,c(2, 10:12)]), H$char_id)

# making list
H<-set_names(table2list(H[,c(1, 4)]), H$char_id)
######
########

####### GL obj
C<-read.csv("Terms4Graphics.csv", header=T,  stringsAsFactors = FALSE, na.strings = c("", NA))
C<-as_tibble(C)
GL<-vector("list", nrow(C))
names(GL)<-C$ID

for (i in 1:nrow(C)){
  GL[[i]]$name<-C$Term[i]
  GL[[i]]$layer<-strsplit(C$pic_id[i], ", ")[[1]]
}
########

# F obj, focal terms for graphics and their layer id
# level 1, elementary focal terms
al<-read.csv(file="alluvial_plot_data.csv", header=T,  stringsAsFactors = FALSE, na.strings = c("", NA) )
al<-as_tibble(al)

# I renamed F to F.obj objects
F.obj<-filter(C, !is.na(pic_id))[,c(1,2)]
F.obj<-select(al, Level_1, Level_1_ids)
##

############
#
# Getting statistics
#

# creating ontology with annotations
ONT=get_OBO("STEP_3/HAO.obo", extract_tags="everything", propagate_relationships = c("BFO:0000050", "is_a"))
#get_onto_name("HAO:x000576", ONT)

#ONT<-HAO
ONT$terms_selected_id<-H

# getting statistics: number of chars per each term
stat<-sapply(F.obj$ID, function(x)
  get_descendants_chars(ONT, annotations="manual", terms=x) %>% length )
F.obj<-add_column(F, stat=stat)

# getting statistics: number of chars per each term
stat<-sapply(al$Level_1_ids, function(x)
  get_descendants_chars(ONT, annotations="manual", terms=x) %>% length )

stat1<-sapply(al$Level_2_ids, function(x)
  get_descendants_chars(ONT, annotations="manual", terms=x) %>% length )

al<-add_column(al, Level_1_nchars=stat, Level_2_nchars=stat1, Level_3_nchars=394)

#write.csv(al, file="alluvial_plot_data.csv", row.names =FALSE )
#F<-add_column(F, stat)

## getting statistics: number of states per each term
stat<-sapply(al$Level_1_ids, function(x)
  get_descendants_chars(ONT, annotations="manual", terms=x))
lapply(stat, function(x) CH[CH$char.id %in% x, "Nstates" ] %>% sum) %>% unlist->stat

stat1<-sapply(al$Level_2_ids, function(x)
  get_descendants_chars(ONT, annotations="manual", terms=x))
lapply(stat1, function(x) CH[CH$char.id %in% x, "Nstates" ] %>% sum) %>% unlist->stat1

al<-add_column(al, Level_1_nstates=stat, Level_2_nstates=stat1, Level_3_nstates=CH$Nstates %>% sum)

###

# getting color for statistics
hm.palette <- colorRampPalette(brewer.pal(9, 'Spectral') %>% rev, space='Lab')
color<-hm.palette(max(F$stat)+1)
F<-add_column(F, color=color[F$stat+1])

# random color
hm.palette <- colorRampPalette(brewer.pal(9, 'Spectral') %>% rev, space='Lab')
hm.palette <- colorRampPalette(brewer.pal(9, 'Accent') %>% rev, space='Lab')
color<-hm.palette(nrow(F))
sample(color, 16)
F<-add_column(F, color=color)
F$color<-sample(color, 16)


############
#
# Getting layer ids to draw
#
library("ontologyIndex", lib.loc="~/.local/R/site-library")


##

Layer<-sapply(F$Level_1_ids, function(x) GL[[x]]$layer %>% as.numeric)


Layer<-layer2terms(F$Level_1_ids, GL, ONT)
layer2terms(al$Level_2_ids, GL, ONT)

# get layes ids for each terms
layer2terms<-function(terms, GL, ONT){
out<-vector("list", length(terms))
names(out)<-terms    
  #i=11
  for (i in 1:length(terms)){
X<-names(GL)[names(GL) %in% get_descendants(ONT, terms[i])]
sapply(X, function(x) GL[[x]]$layer %>% as.numeric) %>% unlist %>% c(.)->X
out[[i]]<-X[!is.na(X)]
  }
return(out)
}
####
# random color
hm.palette <- colorRampPalette(brewer.pal(9, 'Spectral') %>% rev, space='Lab')
hm.palette <- colorRampPalette(brewer.pal(9, 'Accent') %>% rev, space='Lab')
hm.palette <- colorRampPalette(brewer.pal(10, 'Paired') %>% rev, space='Lab')
hm.palette <- colorRampPalette(brewer.pal(9, 'Blues'), space='Lab')
color<-hm.palette(nrow(F))
sample(color, 16)
F<-add_column(F, color=color)
F$color<-sample(color, 16)
F$color<-color

color<-hm.palette(14)
color<-color[F$mean.changes.per.state %>% round(0)]
color[is.na(color)]<-hm.palette(14)[14]
F$color<-color
# scheme 2
color<-hm.palette(150)
color<-color[F$mean.changes.per.state %>% round(0)%>%.^2]
color[is.na(color)]<-hm.palette(150)[150]
F$color<-color
# scheme 2

############
###
# DRAW
F
petal@paths[c(1, 2)]$path@rgb
i=2

draw_pic<-function(F, Layer, petal){
  i=1
for (i in 1:nrow(F)){
  LA<-Layer[[F[[i,2]] ]]
  if (length(LA)>0){
  #j=1
  for (j in 1:length(LA)){
    petal@paths[LA[j]]$path@rgb<-F$color[i]
  }
  }
}
  return(petal)
}
# color bar
par(mar=c(bottom=3, left=3, top=19, right=35))
color.bar(color, 0, max(F$stat), ticks=c(0, 37, 75), title="N characters")
# color bar
par(mar=c(bottom=3, left=3, top=19, right=35))
color.bar(hm.palette(14), 0, max(F$mean.changes.per.state), ticks=c(1, 86), title="")
##
F$color<-sample(color, 16)
#F$color<-color
pic<-draw_pic(F, Layer, petal)
grid.picture(pic)
graphics.off()

# F2
F2<-select(al, Level_2, Level_2_ids)
color<-hm.palette(nrow(F2))
F2<-add_column(F2, color=color)
Layer2<-layer2terms(terms=unique(al$Level_2_ids), GL, ONT)
F2$color<-sample(color, 16)
pic<-draw_pic(F2, Layer2, petal)
grid.picture(pic)
graphics.off()
#####################

F
H<-read.csv("finest_HAO_shortlist3a.csv", header=T,  stringsAsFactors = FALSE, na.strings = c("", NA))
H<-as_tibble(H)
H$X1 %>% table()

filter(H, !is.na(X2)) -> X
paste(X$X1, X$X2) %>% unique
c(X$X1[1], X$X2[1]) %>% order
mapply(function(x, y){c(x, y) %>% sort(.) %>% paste(., collapse=" ") }, x=X$X1, y=X$X2, SIMPLIFY = FALSE, USE.NAMES = FALSE)->A 
A %>% unlist %>% unique ->A

write.csv(cbind(A), "tmp-test.csv")

# check descndnats
sapply(al$Level_2_ids, function(x)
  get_descendants(ONT, x))->z
lapply(z, function(x) x[x %in%  al$Level_1_ids] %>% get_onto_name(., ONT))

# check number of chars
sapply(al$Level_1_ids, function(x)
  get_descendants_chars(ONT, annotations = "manual", x) %>% length)->z
names(z)<-get_onto_name(names(z), ONT)
z
sapply(al$Level_2_ids, function(x)
  get_descendants_chars(ONT, annotations = "manual", x) %>% length)->z
names(z)<-get_onto_name(names(z), ONT)
z

sapply(al$Level_3_ids, function(x)
  get_descendants_chars(ONT, annotations = "manual", x) %>% length)

get_descendants_chars(ONT, annotations = "manual", x)
get_ancestors(ONT, get_onto_id("mesosoma", ONT)) %>% get_onto_name(., ONT)

H$annot_name %>% table
