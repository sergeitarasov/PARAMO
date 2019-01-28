# Scripts to run batch analyses of RevBayes

setwd("/home/tarasov/Documents/RevBayes/hmm_dpp/alces_run/mny_runs")

# make up file names for rev
anal_names=paste("run_", c(1:20), sep="")

for (i in seq_along(anal_names)){
fl.in  <- readLines("hmm_dpp_templ.Rev")
fl.out  <- gsub(pattern = "@analysis_name@", replace = anal_names[i], x = fl.in)
writeLines(fl.out , con=paste0(anal_names[i], ".Rev"))
}


# make bash script


cat("#!/bin/bash\n", file="analyses20.sh")

i=1
for (i in seq_along(anal_names)){
  str.in=paste("./rb ", anal_names[i], ".Rev", " > ", "screen_", i, ".log", " &", sep="")
  cat(str.in, file="analyses20.sh", append=T)
  cat("\n", file="analyses20.sh", append=T)
}


~/Documents/RevBayes/hmm_dpp/alces_run$ ./rb hmm_dpp_alce.Rev

"output_dpp_hmm/hidden_dpp-alces.log"

#gsub(pattern = "@analysis_name@", replace = "run_1", x = "analysis_name = $analysis_name$")

?readLines
paste(getwd(), "01")
