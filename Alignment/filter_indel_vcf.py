# This script reads a bed file and remove regions from such file from a VCF file.
# It removed the regions/bases by changing the genotypes to GT ./.
# tribeiro@wisc.edu

# arg 1: indel bed file
# arg 2: input vcf
import re, sys

# each contig - chrm arm - will have a list of the positions to be masked according to a bed file
indel_contig = {}
with open(sys.argv[1]) as indels:
    for row in indels:
        row = re.split('\t',row[:-1])
        if row[0] in indel_contig.keys():
            indel_contig[row[0]] += list(range(int(row[1]),int(row[2])))
        else:
            indel_contig[row[0]] = list(range(int(row[1]),int(row[2])))

# read the vcf file and output a new one masking the positions in the dictionary above
output = open(sys.argv[2][:-4]+'_indelfree.vcf','w')
with open(sys.argv[2]) as vcffile:
    for row in vcffile:
        if row[0] != '#': # if it is not a header
            row = re.split('\t',row[:-1])
            pos = int(row[1])
            if pos in indel_contig[row[0]]: # the dictionary key is the contig, checks if the position is in the region to be masked. If yes, replaces GT with ./.
                row[8] = 'GT'
                row[9] = r'./.'
            elif row[5] == r'.':
                row[8] = 'GT'
                row[9] = r'./.'
            elif float(row[5]) < 32: # also masks out low quality sites
                row[8] = 'GT'
                row[9] = r'./.'
            outting = '\t'.join(row)+'\n'
            output.write(outting)
output.close()
