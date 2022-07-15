# script to split the large fasta file with all the 16 contigs into 5 files for the 5 chrm arms of interest
# tribeiro@wisc.edu

# arg1: input
# arg2: output prefix (e.g. EF43N)

import re, sys

outputid = -1

# ordered in the genome files
chrmarms = ['Chr2L','ChrX','Chr3L','Chr2R','Chr3R']
chrmid = ['3', '4', '5', '7', '8']

with open(sys.argv[1]) as mainfasta:
    for r in mainfasta:
        if r[0] == '>': # start of the contig
            if outputid >= 0: # if a new contig is found and the id is already >= 0 it means there is a file to be finalized and output.
                output.write(outting+'\n')
                output.close()
                outputid = -1
            if r[1:-1] in chrmid: # check if the contig was found. r[] from the 1st character after > until the last character before linebreak - it is, the contig ID
                outting = ''
                for findid in range(0,5):
                    if r[1:-1] == chrmid[findid]:
                        outputid = findid
                output = open(sys.argv[2]+'_'+chrmarms[outputid]+'.fasta','w')
            else:
                outputid = -1 # if contig is found but is not one of the right chrm arms, change this variable that is also used as key
        elif outputid >= 0: 
            outting += r[:-1]


if outputid >= 0: # if the id is already >= 0 by the end of the loop it means there is a file to be finalized and output.
    output.write(outting+'\n')
    output.close()
