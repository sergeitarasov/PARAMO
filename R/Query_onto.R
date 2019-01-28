library("ontologyIndex", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")

HAO<-get_OBO("STEP_3/HAO.obo", extract_tags="everything", propagate_relationships = c("BFO:0000050", "is_a"))

char_id<-paste0("CHAR:", AN$CHAR_ID2)
annot<-set_names(table2list(AN[,c(2,5)]), char_id)

ONT=get_OBO("STEP_3/HAO.obo", extract_tags="everything", propagate_relationships = c("BFO:0000050", "is_a"))

ONT$terms_selected_id<-annot


level2<-set_names(c("HAO:0000397", "HAO:0001089", "HAO:0000494"), c("head", "wings", "legs") )
level3<-set_names(c("HAO:0000012"), c("whole_organism") )
get_descendants_chars(ONT, annotations="manual", terms="HAO:0000012")

stat<-lapply(level2, function(x)
  get_descendants_chars(ONT, annotations="manual", terms=x)  )

stat<-lapply(level3, function(x)
  get_descendants_chars(ONT, annotations="manual", terms=x)  )


# getting statistics: number of chars per each term
stat<-sapply(F.obj$ID, function(x)
  get_descendants_chars(ONT, annotations="manual", terms=x) %>% length )
F.obj<-add_column(F, stat=stat)


#