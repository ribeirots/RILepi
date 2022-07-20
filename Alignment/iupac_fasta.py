# script to check output fasta file for parental genomes to ensure they were correctly generated
# tribeiro@wisc.edu

# argv1 = sample ID (fasta files and _sites.vcf need to be present in the same folder)
# argv2 = chrm arm, e.g. X, 2L, 2R, 3L, 3R

import re, sys, random


# Chrm size based on bam file header. chrmid (Chrmosome ID) is the number of the contig relative to each chrm arm.
if sys.argv[2] == 'X':
    armsize = 22422827
    chrmid = '4'
elif sys.argv[2] == '2L':
    armsize = 23011544
    chrmid = '3'
elif sys.argv[2] == '2R':
    armsize = 21146708
    chrmid = '7'
elif sys.argv[2] == '3L':
    armsize = 24543557
    chrmid = '5'
elif sys.argv[2] == '3R':
    armsize = 27905053
    chrmid = '8'
else:
    print('No chrm arm input as expected (X, 2L, 2R, 3L, 3R). Quitting script.')
    quit()

def iupacode(altref):
    if "A" in altref:
        if "C" in altref:
            return "M"
        elif "G" in altref:
            return "R"
        elif "T" in altref:
            return "W"
        else:
            print('IUPAC code error1')
            quit()
    elif "C" in altref:
        if "G" in altref:
            return "S"
        elif "T" in altref:
            return "Y"
        else:
            print('IUPAC code error2')
            print(altref)
            print(r)
            quit()
    elif "G" in altref and "T" and altref:
        return "K"
    else:
        print('IUPAC code error3')
        quit()

# finding other sites to check (missing data, ht iupac code, lowQ filter, alternate allele (A, in this code))
iupac_dict = {}
iupac_dict['M'] = []
iupac_dict['R'] = []
iupac_dict['W'] = []
iupac_dict['S'] = []
iupac_dict['Y'] = []
iupac_dict['K'] = []

with open(sys.argv[1]+'_sites.vcf') as sitevcf:
    for r in sitevcf:
        if r[0] != '#': # skips headers
            r = re.split('\t',r[:-1])
            if int(r[0]) > int(chrmid): # all the sites for the chrmid were already checked
                break
            elif r[0] == chrmid: # right chrm arm
                if r[-1] != './.': # site is not missing data
                    if float(r[5]) >= 32: # site is high quality
                        if r[4] != '.': # found a site with an alternate allele
                            if r[-1][:3] != '1/1': # not homozygous for the alternate allele
                                if r[-1][:3] == '1/0' or r[-1][:3] == '0/1': # heterozygous ALT/REF
                                    basepair = []
                                    basepair.append(r[3])
                                    basepair.append(r[4])
                                    iupac_dict[iupacode(basepair)].append(int(r[1])-1) # the function returns the IUPAC code
print('ambiguity dict done.')
#checking the fasta file for all the conditions listed in the header
with open(sys.argv[1]+'_Chr'+sys.argv[2]+'.fasta') as fastaf: # open the fasta file
    for r in fastaf:
        r = [c for c in r[:-1]]
        if len(r) != armsize:
            print('Fasta file length = '+str(len(r))+'. While the chrm arm expected length is: '+str(armsize))
            print('quitting script.')
            quit()
        print("chrm arm length is correct.")
        
        for k in iupac_dict.keys():
            for site in iupac_dict[k]:
                r[site] = k
        output = open(sys.argv[1]+'_Chr'+sys.argv[2]+'_diploid.fasta','w')
        output.write(''.join(r)+'\n')

