# script to cound info depth per pool
# Args 1 and 2 (-pl and -ph): parental Low and High sites.vcf files
# Args 3 and 4 (-pli and -phi): parental Low and High INDELS.vcf files
import re, sys, argparse
from infordepth_functions import *
from parentalVCF_dict import *
from readPairs import *

parser = argparse.ArgumentParser()

parser.add_argument("-pl", "--parentLow", help="vcf file from the parent with the low phenotype (full path)")
parser.add_argument("-ph", "--parentHigh", help="vcf file from the parent with the high phenotype (full path)")
parser.add_argument("-pli", "--parentLowIndel", help="INDELS vcf file from the parent with the low phenotype (full path)")
parser.add_argument("-phi", "--parentHighIndel", help="INDELS vcf file from the parent with the high phenotype (full path)")

args = parser.parse_args()

# turn screenLoadKey to 1 to print loading% on screen
screenLoadKey = 0

# Parental VCFs
vcfLowP = open(args.parentLow, 'r') # Low-Parent
vcfHighP = open(args.parentHigh, 'r') # High-Parent

indelLowP = open(args.parentLowIndel, 'r') # Low-Parent INDEL
indelHighP = open(args.parentHighIndel, 'r') # High-Parent INDEL

#chromosomes = ['4'] # chrm X
chromosomes = ['4','3','7','5','8'] # chrms 3:2L, 4:X, 5:3L, 7:2R, ad 8:3R from Dmel.
Qual = 20 # minimum quality to accept a variant in the indel vcf file

# returns 1 dict per chrm arm passed in the 3rd arg
#xdict = parentalVCF(vcfLowP, vcfHighP, chromosomes)[0]
#print(len(xdict.keys()))

chrDicts = parentalVCF(vcfLowP, vcfHighP, chromosomes)
xdict, l2dict, r2dict, l3dict, r3dict = chrDicts
#print('Parental VCFs:')
#print(len(l2dict.keys()), len(xdict.keys()), len(l3dict.keys()), len(r2dict.keys()), len(r3dict.keys()))

# for each chrm, gets the VCF dictionary and the sam files to create output with info depths
# still not designed to include all the chrms
# results in the creation of 2 dictionaries: readDictL and readDictH with all the read pairs from each offspring pool for each chrm
#for chrm in ['X']:
#outputL = open('lowOutput1read.out','w')
#outputH = open('highOutput1read.out','w')
for chrmID in chromosomes:
    if chrmID == '4':
        chrm = 'X'
        dictsnps = xdict
        output = open('X.panel','w')
        print(chrm)
    elif chrmID == '3':
        chrm = '2L'
        dictsnps = l2dict
        output = open('2L.panel','w')
        print(chrm)
    elif chrmID == '7':
        chrm = '2R'
        dictsnps = r2dict
        output = open('2R.panel','w')
        print(chrm)
    elif chrmID == '5':
        chrm = '3L'
        dictsnps = l3dict
        output = open('3L.panel','w')
        print(chrm)
    elif chrmID == '8':
        chrm = '3R'
        dictsnps = r3dict
        output = open('3R.panel','w')
        print(chrm)

    # this function uses low and high parent indel vcfs to remove indel sites from SNP dict (also remove sites around indels) 
    dict_to_use = indelVCFCleaning(indelLowP, indelHighP, dictsnps, chrmID, Qual, 3) # lowParentINDEL, highParentINDEL, SNP dict, chrmID, Quality, surroudning sites to rm
    dictsnps = 0

    # biallelic alleles only
    for k in dict_to_use.keys():
        genoL = dict_to_use[k][0]
        genoH = dict_to_use[k][1]
        alleles = 0
        A = 0
        T = 0
        C = 0
        G = 0
        if 'A' in genoL or 'A' in genoH:
            alleles += 1
            A = 1
        if 'T' in genoL or 'T' in genoH:
            alleles += 1
            T = 1
        if 'C' in genoL or 'C' in genoH:
            alleles += 1
            C = 1
        if 'G' in genoL or 'G' in genoH:
            alleles += 1
            G = 1

        if alleles == 2: # if 1: not a SNP, if 0: missing data, if > 2: not bi-allelic
            if A and T:
                allele1 = 'A'
                allele2 = 'T'
            elif A and C:
                allele1 = 'A'
                allele2 = 'C'
            elif A and G:
                allele1 = 'A'
                allele2 = 'G'
            elif T and C:
                allele1 = 'T'
                allele2 = 'C'
            elif T and G:
                allele1 = 'T'
                allele2 = 'G'
            elif C and G:
                allele1 = 'C'
                allele2 = 'G'
            else:
                print('Error: SNP not found. Quitting.')
                quit()

            if allele1 in genoL:
                genoL0 = 1
            else:
                genoL0 = 0
            if allele2 in genoL:
                genoL1 = 1
            else:
                genoL1 = 0
            if allele1 in genoH:
                genoH0 = 1
            else:
                genoH0 = 0
            if allele2 in genoH:
                genoH1 = 1
            else:
                genoH1 = 0

            SNP = [chrm, k, allele1, allele2, genoL0, genoL1, genoH0, genoH1]
            outting = list(map(str,SNP))
            outting = '\t'.join(outting)+'\n'
            output.write(outting)
    output.close()
