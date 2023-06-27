# script to compare RIL against two parental genomes
# Args 1: parental dict suffix (pkl file)
# 2: ril vcf

import re, sys, pickle
from infordepth_functions import *
from pvcf_ril import *
#from readPairs import *

# turn screenLoadKey to 1 to print loading% on screen
screenLoadKey = 0

# VCF
vcfF = open(sys.argv[2], 'r') # focalvcf

#chromosomes = ['4'] # chrm X
chromosomes = ['3','4','5','7','8'] # chrms 2L, X, 3L, 2R, ad 3R from Dmel.

# returns 1 dict per chrm arm passed in the 3rd arg
#xdict = parentalVCF(vcfLowP, vcfHighP, chromosomes)[0]
#print(len(xdict.keys()))
pxp_dicts = []
for i in range(0,5):
    with open(sys.argv[1]+str(i)+'.pkl','rb') as readInput:
        pxp_dicts.append(pickle.load(readInput))

dictId = sys.argv[1]
rilId = sys.argv[2]

dictId = re.split('/',dictId)[-1]
rilId = re.split('/',rilId)[-1]

dictId = re.split('_', dictId)[0]
rilId = re.split('_', rilId)[0]

mismatches = onexpp(vcfF, pxp_dicts, chromosomes)
mismatches = [rilId, sys.argv[1]] + mismatches

outfile =  'qc_'+dictId+'_'+rilId+'.csv'
output = open(outfile,'w')
outting = '\t'.join(list(map(str,mismatches)))+'\n'
output.write(outting)
output.close()
print(mismatches)
