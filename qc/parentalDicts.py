# script to compare RIL against two parental genomes
# Args 1 and 2: parental Low and High sites.vcf files
# arg 3: outputPrefix
import re, sys, pickle
from infordepth_functions import *
from pvcf_ril import *
#from readPairs import *

# turn screenLoadKey to 1 to print loading% on screen
screenLoadKey = 0

# Parental VCFs
vcfP1 = open(sys.argv[1], 'r') # Parent1
vcfP2 = open(sys.argv[2], 'r') # Parent2

outfilename = re.split(r'/',sys.argv[3])[-1]
#chromosomes = ['4'] # chrm X
chromosomes = ['3','4','5','7','8'] # chrms 2L, X, 3L, 2R, ad 3R from Dmel.

# returns 1 dict per chrm arm passed in the 3rd arg
#xdict = parentalVCF(vcfLowP, vcfHighP, chromosomes)[0]
#print(len(xdict.keys()))
pxp_dicts = PF_VCF(vcfP1, vcfP2, chromosomes)
print(len(pxp_dicts[1].keys()))

for i in range(0,len(pxp_dicts)):
    with open (sys.argv[3]+str(i)+'.pkl', 'wb') as output:
        pickle.dump(pxp_dicts[i], output)

