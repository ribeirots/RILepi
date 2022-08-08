# script to run scantwo
library(qtl)
setwd('Documents/git_repos/ril_epistasis/RILepi/rqtl/')
etoh = read.cross("csv", , "FR_etoh_rqtl.csv", sep=",", alleles=c("Z","F"), 
                  genotypes = c("FF","FZ","ZZ"))

etoh_genoprob2 <- calc.genoprob(etoh,step=0, err=0.001, )

etoh_s2 <- scantwo(etoh_genoprob2, verbose = FALSE, pheno.col=2)

