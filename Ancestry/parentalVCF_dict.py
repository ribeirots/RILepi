# This function takes as input two VCF files (and a list of empty dicts - one per chrm) and generates dictionaries for each chromosome dict gave to it
import re, sys
from infordepth_functions import *

# PL_vcf (path to VCF file from Low Parent
# PH_vcf (path to VCF file from High Parent
# chrmID (list with the IDs of the chromosomes to be used as present in the VCF file (e.g. ['4'] for X, or ['3','4','5','7','8'] for 5chrm arms in Dmel (2L, X, 3L, 2R, 3R)
def parentalVCF(PL_vcf, PH_vcf, chrmID):
    Q = 20 # quality threshold
    parentDicts = [] # will hold one dict per chrmID

    for i in chrmID:
        parentDicts.append({})

    for PL_site, PH_site in zip(PL_vcf, PH_vcf):
        # check they have the same number of headers
        if (PL_site[0] != '#' and PH_site[0] == '#') or (PL_site[0] == '#' and PH_site[0] != '#'):
            print('vcf headers had different number of lines.\nQuitting script.')
            quit()
        elif PL_site[0] != '#':
            PL_site = re.split('\t',PL_site[:-1])
            PH_site = re.split('\t',PH_site[:-1])
            if PL_site[0] in chrmID and PL_site[5] != '.' and PH_site[5] != '.':
                if PL_site[:2] != PH_site[:2]:
                    print('vcf files dont match. \tQuitting.')
                    quit()
                if float(PL_site[5]) >= Q and float(PH_site[5]) >= Q:
                    if len(PL_site[4]) == 1 and len(PH_site[4]) == 1 and (PL_site[4] != '.' or PH_site[4] != '.'):
    #                    print('SNP1: '+ str(c[1]))
                        # looping through this section in case there is more than one sample in the vcf for the same parental genotype (multiple sequencing)
                        genotypeL = gtParentalVCF(PL_site) # returns genotype or -1 if missing
                        if genotypeL != -1:
                            genotypeH = gtParentalVCF(PH_site) # returns genotype or -1 if missing

                            if genotypeH != -1 and genotypeL != genotypeH and not (len(genotypeL) == 2 and len(genotypeH) == 2):
                                if (len(genotypeL) == 1 and len(genotypeH) == 2 and genotypeL in genotypeH) or (len(genotypeL) == 2 and len(genotypeH)== 1 and genotypeH in genotypeL) or (len(genotypeL)==1 and len(genotypeH)==1):
                                    for i in range(len(chrmID)):
                                        if PL_site[0] == chrmID[i]:
                                            parentDicts[i][PL_site[1]] = [genotypeL, genotypeH]
    print('Parental SNP dictionaries created.')
    return parentDicts



# function to find the list of sites to filter because of 1 indel
# it returns a list of the sites to be removed from SNP dictionary because of indels
# variant = list object with the vcf information for a indel
# surround = how many bases around the indel should be filtered out
def indelFilter(variant, surround):
    pos = variant[1]
    start = int(pos) - surround # the first position to filtered out, the start of the indel minus the chosen surrounding area
    # next lines decided whether ref or alt are longer, the longer one is considered the entire length of the indel and will be added the chosen surronding area to be removed
    ref = len(variant[3])
    alt = len(variant[4])
    if ref > alt:
        end = int(pos)+ref+surround
    else:
        end = int(pos)+alt+surround
    filterRange = []
    # makes a list with all the positions to be filtered out. Output them as str to be compared with the dictionary keys of positions
    for i in range(start, end):
        filterRange.append(str(i))
    return filterRange


# return parental SNP dictionary without indels and regions neighboring indels
def indelVCFCleaning(indelPL, indelPH, snpP, c, Q, neighb):
    chrmID = c
    Qual = Q

    for dkt in [indelPL, indelPH]: # remove all indel sites from indelPL and then all the sites from indelPH
        for indel in dkt:
            if indel[0] != '#': # skips rows that start with #, those are headers/additional information
                indel = re.split('\t',indel[:-1])
                if indel[0] == chrmID and float(indel[5]) >= Qual: # checks if this row contains information about the current chromosome arm and if the variant passes the minimum quality check
                    # the function below obtains a list of the positions to be filtered out (indel + surrounding)
                    indelSites = indelFilter(indel, neighb) # neighb is the number of bases around the indel that will also be filtered out
                    for indelSite in indelSites:
                        if indelSite in snpP.keys(): # if the site to be filtered out is in the
                            del snpP[indelSite]

    return snpP

# returns a list of positions to be removed for a given chrm art
# input: indel VCF, chrom, Quality, surrounding regions size
def indel2rm(indelVCF, c, Q, neighb):
    chrmID = c
    Qual = Q
    indelRM = []

    for indel in indelVCF:
        if indel[0] != '#': # skips rows that start with #, those are headers/additional information
            indel = re.split('\t',indel[:-1])
            if indel[0] == chrmID and float(indel[5]) >= Qual: # checks if this row contains information about the current chromosome arm and if the variant passes the minimum quality check
                # the function below obtains a list of the positions to be filtered out (indel + surrounding)
                indelSites = indelFilter(indel, neighb) # neighb is the number of bases around the indel that will also be filtered out
                for indelSite in indelSites:
                    indelRM.append(indelSite)
    
    return indelRM

