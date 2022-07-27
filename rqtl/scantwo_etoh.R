# script to run scantwo
library(qtl)
etoh = read.cross("csv", , "FR_etoh_rqtl.csv", sep=",", alleles=c("Z","F"), 
                  genotypes = c("FF","FZ","ZZ"))

etoh_genoprob2 <- calc.genoprob(etoh,step=0, err=0.001, )
