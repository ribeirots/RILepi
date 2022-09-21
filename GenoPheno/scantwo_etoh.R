# arg1: input file

args = commandArgs(trailingOnly=TRUE)


library(qtl)
finame=args[1]
fname = paste(finame,'.csv',sep="")

ril_data <- read.cross(format='csv', dir=".", file=fname, genotypes = c('ZZ','FZ','FF'), alleles=c('Z','F'))
ril_data_genoprob <- calc.genoprob(ril_data,step=0)

# Scanone


############ SCANTWO
print("starting calc.geno.")
ril_data_genoprob2 <- calc.genoprob(ril_data,step=0, err=0.001)

print("starting scantwo")
data_s2 <- scantwo(ril_data_genoprob2, verbose = FALSE, pheno.col=2)

#perm_etoh_s2 <- scantwo(ril_data_genoprob2, n.perm=10)


#write.table(etoh_s2,file="etoh_mr_scantwo_3r.csv",row.names = T,sep=",")
#write.table(summary(perm_etoh_s2),file="etoh_scantwo_perm10_summary.csv",row.names = T,sep=",")


#jpeg(filename=paste(finame,"_step0.jpg",sep=""),width = 1000, height = 648)
#par(mar=c(5,5,4,2))
#plot(etoh_s2, chr=c(1,2,3,4,5))
#plot(etoh_s2, chr=c(1,2,3,4,5))
#dev.off()

save.image(file = paste(finame,"_avoid.RData",sep=""))
