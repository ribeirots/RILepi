# This function takes as input two VCF files (and a list of empty dicts - one per chrm) and generates dictionaries for each chromosome dict gave to it
import re, sys, pickle
from infordepth_functions import *

# this function takes in the dictionary with the SNPs between the 2 potential parental strainst being investigated and remove the sites in which the Base from OffStrain matches the OutParent
# OffStrain = the strain that is not believed to have been used to create the RILs
# OutParent = the Parental strain not being investigated (e.g. Investigating whether ZI418N is the right parent in ZI418NxFR320N, FR320N is the OutParent, and a line such as ZI251N or ZI413N can be used as OffStrain)
# outP_vcf (path to VCF file from Parent not being checked)
# pp_dicts (dictionary with mismatches among the parental genomes)
def outParentvcfFilter(outP_vcf, pp_dicts):
    Q = 20 # quality threshold
    sitesAnalyzed = 0
    countP0 = 0
    countP1 = 0
    print10 = 0
    p1list = []
    p2list = []

    for F_site in outP_vcf:
        readSite = 1
        if F_site[0] != '#':
            F_site = re.split('\t',F_site[:-1])
            F_site[0] = str(F_site[0])
            if str(F_site[0]) == '3':
                dictKey = 0
            elif str(F_site[0]) == '4':
                dictKey = 1
            elif str(F_site[0]) == '5':
                dictKey = 2
            elif str(F_site[0]) == '7':
                dictKey = 3
            elif str(F_site[0]) == '8':
                dictKey = 4
            else:
                readSite = 0

            if readSite:
                if str(F_site[1]) in pp_dicts[dictKey].keys() and F_site[5] != '.': # site is mismatch between parentals and isn't missing data
                    if float(F_site[5]) >= Q: # matches quality
                        if len(F_site[4]) < 4: # len == 1 when 1 alternate base, == 3 when alternate is heterozygous and neither base match the reference
                            # looping through this section in case there is more than one sample in the vcf for the same parental genotype (multiple sequencing)
                            genotypeF = gtParentalVCF(F_site) # returns genotype or -1 if missing
                            if genotypeF != -1:
                                gtP1 = pp_dicts[dictKey][str(F_site[1])][0]
                                gtP2 = pp_dicts[dictKey][str(F_site[1])][1]
                                if len(genotypeF) == 1:
                                    if genotypeF in gtP1 and genotypeF in gtP2:
                                        print("Error generating parental mismatch!")
                                        print('Quitting!')
                                        quit()
                                    elif genotypeF in gtP1:
                                        countP1 += 1
                                        if print10:
                                            if len(p1list) < 10:
                                                p1list.append([F_site[0],F_site[1],genotypeF, gtP1, gtP2])
                                            elif len(p2list) >= 10:
                                                print(p1list)
                                                print(p2list)
                                                quit()
                                    elif genotypeF in gtP2:
                                        del pp_dicts[dictKey][F_site[1]]
                                        if print10:
                                            if len(p2list) < 10:
                                                p2list.append([F_site[0],F_site[1],genotypeF, gtP1, gtP2])
                                            elif len(p1list) >= 10:
                                                print(p1list)
                                                print(p2list)
                                                quit()
                                    else:
                                        countP0 += 1
                                elif len(genotypeF) == 2:
                                    if genotypeF[0] in gtP1 or genotypeF[1] in gtP1:
                                        if genotypeF[0] not in gtP2 and genotypeF[1] not in gtP2:
                                            countP1 += 1
                                            if print10:
                                                if len(p1list) < 10:
                                                    p1list.append([F_site[0],F_site[1],genotypeF, gtP1, gtP2])
                                                elif len(p2list) >= 10:
                                                    print(p1list)
                                                    print(p2list)
                                                    quit()
                                    elif genotypeF[0] in gtP2 or genotypeF[1] in gtP2:
                                        del pp_dicts[dictKey][F_site[1]]
                                        if print10:
                                            if len(p2list) < 10:
                                                p2list.append([F_site[0],F_site[1],genotypeF, gtP1, gtP2])
                                            elif len(p1list) >= 10:
                                                print(p1list)
                                                print(p2list)
                                                quit()
                                    else:
                                        countP0 += 1
                                else:
                                    print('Error generating offspring genotype!')
                                    print('Quitting.')
                                    quit()

    return pp_dicts


# compare two SNP list of dictionaries
def rilxpp(dict1list, dict2list):
    countPa = 0 # in dict1 but not a SNP in dict2 (could be due to missing data or not being a SNP in dict2
    countP = 0 # dict1 matching both, at least partially
    countP0 = 0 # dict1 matching neither
    countP1 = 0 # dict1 match only P1
    countP2 = 0 # dict1 match only P2
    for i in range(0,len(dict1list)):
        for k in dict1list[i].keys():
            if k in dict2list[i].keys() and k != '-1':
                gtRIL = dict1list[i][k][1] # genotype RIL
                gtP1 = dict2list[i][k][0] # gt parent 1
                gtP2 = dict2list[i][k][1] # gt parent 2
                if len(gtRIL) == 1:
                    if gtRIL in gtP1 and gtRIL not in gtP2:
                        countP1 += 1
                    elif gtRIL not in gtP1 and gtRIL in gtP2:
                        countP2 += 1
                    elif gtRIL not in gtP1 and gtRIL not in gtP2:
                        countP0 += 1
                    else:
                        countP += 1
                elif len(gtRIL) == 2:
                    if gtRIL[0] in gtP1 or gtRIL[1] in gtP1:
                        if gtRIL[0] not in gtP2 and gtRIL[1] not in gtP2:
                            countP1 += 1
                        else:
                            countP += 1
                    elif gtRIL[0] in gtP2 or gtRIL[1] in gtP2:
                        countP2 += 1
                    else:
                        countP0 += 1
                else:
                    print('Weird ril genotype. \nQuitting.')
                    quit()
            else:
                countPa += 1 # SNP between RIL and offParent not a SNP between the 2 parents being investigated
    return [countPa, countP0, countP, countP1, countP2]


# compare one file against a SNP list of dictionaries with genotypes from SNPs among two other vcf files
def onexpp(P_vcf, dict2list, chrmID):
    Q = 20 # quality threshold

    countP = 0 # focalvcf matching both, at least partially
    countP0 = 0 # focalvcf matching neither
    countP1 = 0 # focalvcf match only P1
    countP2 = 0 # focalvcf match only P2

    for P_site in P_vcf:
        if P_site[0] != '#':
            P_site = re.split('\t',P_site[:-1])
            P_site[0] = str(P_site[0])
            if P_site[0] in chrmID and P_site[5] != '.':
                i = chrmID.index(P_site[0]) # get if of the dictionary for the correct chrm arm 
                if P_site[1] in dict2list[i].keys():
                    k = P_site[1]
                    if float(P_site[5]) >= Q:
                        if len(P_site[4]) < 4: # len == 1 when 1 alternate base, == 3 when alternate is heterozygous and neither base match the reference
                            # looping through this section in case there is more than one sample in the vcf for the same parental genotype (multiple sequencing)
                            genotypeP = gtParentalVCF(P_site) # returns genotype or -1 if missing
                            if genotypeP != -1:
                                gtRIL = genotypeP
                                gtP1 = dict2list[i][k][0] # gt parent 1
                                gtP2 = dict2list[i][k][1] # gt parent 2
                                if len(gtRIL) == 1:
                                    if gtRIL in gtP1 and gtRIL not in gtP2:
                                        countP1 += 1
                                    elif gtRIL not in gtP1 and gtRIL in gtP2:
                                        countP2 += 1
                                    elif gtRIL not in gtP1 and gtRIL not in gtP2:
                                        countP0 += 1
                                    else:
                                        countP += 1
                                elif len(gtRIL) == 2:
                                    if gtRIL[0] in gtP1 or gtRIL[1] in gtP1:
                                        if gtRIL[0] not in gtP2 and gtRIL[1] not in gtP2:
                                            countP1 += 1
                                        else:
                                            countP += 1
                                    elif gtRIL[0] in gtP2 or gtRIL[1] in gtP2:
                                        countP2 += 1
                                    else:
                                        countP0 += 1
                                else:
                                    print('Weird ril genotype. \nQuitting.')
                                    quit()
    return [countP0, countP, countP1, countP2]

# F_vcf (path to VCF file from Parent
# pp_dicts (dictionary with mismatches among the parental genomes)
def ppCheck(F_vcf, pp_dicts):
    Q = 20 # quality threshold
    sitesAnalyzed = 0
    countP0 = 0
    countP1 = 0
    countP2 = 0
    print10 = 0
    p1list = []
    p2list = []

    for F_site in F_vcf:
        readSite = 1
        if F_site[0] != '#':
            F_site = re.split('\t',F_site[:-1])
            F_site[0] = str(F_site[0])
            if str(F_site[0]) == '3':
                dict_to_use = pp_dicts[0]
            elif str(F_site[0]) == '4':
                dict_to_use = pp_dicts[1]
            elif str(F_site[0]) == '5':
                dict_to_use = pp_dicts[2]
            elif str(F_site[0]) == '7':
                dict_to_use = pp_dicts[3]
            elif str(F_site[0]) == '8':
                dict_to_use = pp_dicts[4]
            else:
                readSite = 0

            if readSite:
                if str(F_site[1]) in dict_to_use.keys() and F_site[5] != '.': # site is mismatch between parentals and isn't missing data
                    if float(F_site[5]) >= Q: # matches quality
                        if len(F_site[4]) < 4: # len == 1 when 1 alternate base, == 3 when alternate is heterozygous and neither base match the reference
                            # looping through this section in case there is more than one sample in the vcf for the same parental genotype (multiple sequencing)
                            genotypeF = gtParentalVCF(F_site) # returns genotype or -1 if missing
                            if genotypeF != -1:
                                gtP1 = dict_to_use[str(F_site[1])][0]
                                gtP2 = dict_to_use[str(F_site[1])][1]
                                if len(genotypeF) == 1:
                                    if genotypeF in gtP1 and genotypeF in gtP2:
                                        print("Error generating parental mismatch!")
                                        print('Quitting!')
                                        quit()
                                    elif genotypeF in gtP1:
                                        countP1 += 1
                                        if print10:
                                            if len(p1list) < 10:
                                                p1list.append([F_site[0],F_site[1],genotypeF, gtP1, gtP2])
                                            elif len(p2list) >= 10:
                                                print(p1list)
                                                print(p2list)
                                                quit()
                                    elif genotypeF in gtP2:
                                        countP2 += 1
                                        if print10:
                                            if len(p2list) < 10:
                                                p2list.append([F_site[0],F_site[1],genotypeF, gtP1, gtP2])
                                            elif len(p1list) >= 10:
                                                print(p1list)
                                                print(p2list)
                                                quit()
                                    else:
                                        countP0 += 1
                                elif len(genotypeF) == 2:
                                    if genotypeF[0] in gtP1 or genotypeF[1] in gtP1:
                                        if genotypeF[0] not in gtP2 and genotypeF[1] not in gtP2:
                                            countP1 += 1
                                            if print10:
                                                if len(p1list) < 10:
                                                    p1list.append([F_site[0],F_site[1],genotypeF, gtP1, gtP2])
                                                elif len(p2list) >= 10:
                                                    print(p1list)
                                                    print(p2list)
                                                    quit()
                                    elif genotypeF[0] in gtP2 or genotypeF[1] in gtP2:
                                        countP2 += 1
                                        if print10:
                                            if len(p2list) < 10:
                                                p2list.append([F_site[0],F_site[1],genotypeF, gtP1, gtP2])
                                            elif len(p1list) >= 10:
                                                print(p1list)
                                                print(p2list)
                                                quit()
                                    else:
                                        countP0 += 1
                                else:
                                    print('Error generating offspring genotype!')
                                    print('Quitting.')
                                    quit()

    return [countP0, countP1, countP2]

# P_vcf (path to VCF file from Parent
# F_vcf (path to VCF file from Offspring
# chrmID (list with the IDs of the chromosomes to be used as present in the VCF file (e.g. ['4'] for X, or ['3','4','5','7','8'] for 5chrm arms in Dmel (2L, X, 3L, 2R, 3R)
def PF_VCF(P_vcf, F_vcf, chrmID):
    Q = 20 # quality threshold
    parentDicts = [] # will hold one dict per chrmID
    sitesAnalyzed = 0

    for i in chrmID:
        parentDicts.append({})

    for P_site, F_site in zip(P_vcf, F_vcf):
        # check they have the same number of headers
        if (P_site[0] != '#' and F_site[0] == '#') or (P_site[0] == '#' and F_site[0] != '#'):
            print('vcf headers had different number of lines.\nQuitting script.')
            quit()
        elif P_site[0] != '#':
            P_site = re.split('\t',P_site[:-1])
            P_site[0] = str(P_site[0])
            F_site = re.split('\t',F_site[:-1])
            F_site[0] = str(F_site[0])
            if P_site[0] in chrmID and P_site[5] != '.' and F_site[5] != '.':
                if P_site[:2] != F_site[:2]:
                    print('vcf files dont match. \tQuitting.')
                    quit()
                if float(P_site[5]) >= Q and float(F_site[5]) >= Q:
                    if len(P_site[4]) < 4 and len(F_site[4]) < 4: # len == 1 when 1 alternate base, == 3 when alternate is heterozygous and neither base match the reference
    #                    print('SNP1: '+ str(c[1]))
                        # looping through this section in case there is more than one sample in the vcf for the same parental genotype (multiple sequencing)
                        genotypeP = gtParentalVCF(P_site) # returns genotype or -1 if missing
                        if genotypeP != -1:
                            genotypeF = gtParentalVCF(F_site) # returns genotype or -1 if missing, i.e. sites with genotype information for both vcf files
                            if genotypeF != -1:
                                sitesAnalyzed += 1 # count total number of sites analyzed
                                if genotypeP != genotypeF:
                                    if (len(genotypeP) == 1 and len(genotypeF) == 2 and genotypeP not in genotypeF) or (len(genotypeP) == 2 and len(genotypeF)== 1 and genotypeF not in genotypeP) or (len(genotypeP)==1 and len(genotypeF)==1):
                                        for i in range(len(chrmID)):
                                            if P_site[0] == chrmID[i]:
                                                parentDicts[i][P_site[1]] = [genotypeP, genotypeF]
                                    elif len(genotypeP) == 2 and len(genotypeF) == 2:
                                        if genotypeP[0] not in genotypeF and genotypeP[1] not in genotypeF:
                                            for i in range(len(chrmID)):
                                                if P_site[0] == chrmID[i]:
                                                    parentDicts[i][P_site[1]] = [genotypeP, genotypeF]
#    print('Parental SNP dictionaries created.')
    parentDicts[0]['-1'] = sitesAnalyzed # total
    return parentDicts
