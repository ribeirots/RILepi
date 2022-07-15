# script that used INDEL VCF to generate bed file with indel locations
# tribeiro@wisc.edu
# arg1: input INDELS.vcf

import sys, re

flank = 3 # number of sites surrounding the indel that will be removed

output = open(sys.argv[1][:-4]+'_indelfilter.vcf','w')
outting = []
poscheck = []
chrm = 'N'
with open(sys.argv[1]) as vcfindel:
    for row in vcfindel:
        if row[0] != '#': # skips header
            row = re.split('\t', row[:-1])
            if chrm != row[0]:
                poscheck = []
            chrm = row[0] # chromosome ID
            startpos = int(row[1])-(flank+1) # starting position of the bed region is the position of the indel in the VCF minus flanking
            reflen = len(row[3]) # get the length of the reference allele
            refalt = len(row[4]) # get the length of the alt allele
            endpos = int(row[1])+max(reflen,refalt)+(flank-1) # end positions of the bed region is the position of the indel plus whichever is larger between the reference and alt allele plus the flanking region
            for indrange in range(startpos, endpos): # print the vcf rows
                if indrange not in poscheck:
                    poscheck.append(indrange)
                    filtpos = indrange
                    outting = [chrm, str(filtpos), r'.', 'N', r'.', r'.', r'.', r'.', 'GT', r'./.'] # output
                    output.write('\t'.join(outting)+'\n')
        else:
            output.write(row)
output.close()
