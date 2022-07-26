# this is script is here just to rename the rmarkdown output file, not sure how to do it in a nicer way
setwd('Documents/git_repos/ril_epistasis/RILepi/rqtl/')
rmarkdown::render("README.Rmd", output_file = "README.html")



