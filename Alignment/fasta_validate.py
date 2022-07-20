# script to check output fasta file for parental genomes to ensure they were correctly generated
# the script will check (1) their size, (2) if they were masked for indels, missing data, and low-quality data, (3) IUPAC code for ht sites, and (4) alternate allele (A) used instead of reference allele.
# tribeiro@wisc.edu

# argv1 = sample ID (fasta files, _INDELS.vcf and _sites.vcf need to be present in the same folder)
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

# find list of indels:
indels_to_test = []
with open(sys.argv[1]+'_INDELS.vcf') as indels:
    for r in indels:
        if r[0] != '#': # skips headers
            r = re.split('\t', r[:-1])
            if r[0] == chrmid and len(indels_to_test) <= 10:
                indels_to_test.append(int(r[1])) # append one position of the indel
            elif r[0] == chrmid: #it will only test 10 indels, if there are more than 10 in the file it will decide whether to store a new one and where in the test vector
                if random.randrange(0,100) < 2: # 2% change to keep new indel
                    repl_index = random.randrange(0,10) # the indel to be replaced is random
                    indels_to_test[repl_index] = int(r[1])
            elif int(r[0]) > int(chrmid):
                break # breaks the loop if the chr arm in the vcf file is larger than the chrm ID, because by then all the sites have been sampled already

if indels_to_test[9] > 5:
    indels_to_test[9] = indels_to_test[9] - 2 # I just picked one random large position to reduce by 2 to test if sites downstream of the indel were masked as well.

# finding other sites to check (missing data, ht iupac code, lowQ filter, alternate allele (A, in this code))
missing_to_test = []
lowq_to_test = []
a_to_test = []
iupac_to_test = []
with open(sys.argv[1]+'_sites.vcf') as sitevcf:
    for r in sitevcf:
        if r[0] != '#': # skips headers
            r = re.split('\t',r[:-1])
            if int(r[0]) > int(chrmid): # all the sites for the chrmid were already checked
                break
            elif r[0] == chrmid and len(missing_to_test) <= 10 and r[-1] == './.': # is in the right chrm arm, vector is small enough, site is missing data
                missing_to_test.append(int(r[1])-1) # subtracting 1 because vcf is not 0 based. I didn't do it for indels because surrounding sites are also expected to be filtered out.
            elif r[0] == chrmid and r[-1] == './.': # is in the right chrm arm, vector is NOT small enough, site is missing data
                if random.randrange(0,1000) < 1: # 0.1% chance to replace old missing site
                    repl_index = random.randrange(0,10)
                    missing_to_test[repl_index] = int(r[1])-1
            elif r[0] == chrmid and len(lowq_to_test) <= 10 and float(r[5]) < 32: # is in the right chrm arm, vector is small enough, site is low quality
                lowq_to_test.append(int(r[1])-1)
            elif r[0] == chrmid and float(r[5]) < 32: # is in the right chrm arm, vector is NOT small enough, site is low quality
                if random.randrange(0,1000) < 1: # 0.1% chance to replace old missing site
                    repl_index = random.randrange(0,10)
                    lowq_to_test[repl_index] = int(r[1])-1
            elif r[0] == chrmid and r[4] != '.': # found a site with an alternate allele
                if r[-1][:3] == '1/1': # homozygous for the alternate allele
                    if r[4] == 'A': # randomly chose to only check alternate A, for simplicity
                        if len(a_to_test) <= 100: # this is a larger vector because some sites will be indels and wont be tested
                            a_to_test.append(int(r[1])-1)
                        else:
                            if random.randrange(0,1000) < 1:
                                repl_index = random.randrange(0,100)
                                a_to_test[repl_index] = int(r[1]) - 1
                elif r[-1][:3] == '1/0' or r[-1][:3] == '0/1': # heterozygous
                    if len(iupac_to_test) <= 10:
                        iupac_to_test.append(int(r[1])-1)
                    else:
                        if random.randrange(0,1000) < 1:
                            repl_index = random.randrange(0,10)
                            iupac_to_test[repl_index] = int(r[1]) - 1

#checking the fasta file for all the conditions listed in the header
with open(sys.argv[1]+'_Chr'+sys.argv[2]+'_diploid.fasta') as fastaf: # open the fasta file
    for r in fastaf:
        r = [c for c in r[:-1]]
        if len(r) != armsize:
            print('Fasta file length = '+str(len(r))+'. While the chrm arm expected length is: '+str(armsize))
            print('quitting script.')
            quit()

        # check indels
#        for site in indels_to_test:
#            if r[site] != 'N':
#                print('This site is an indel. It should have been filtered out and replaced with N. Site: '+str(site)+'. Allele: '+r[site])
#                quit()
        for site in indels_to_test:
            if r[site] != 'N':
                print('This site should have been an indel and filtered out. It should be N. Site:'+str(site)+'. Allele: '+r[site])
                quit()

        for site in missing_to_test:
            if r[site] != 'N':
                print('This site has missing data. It should have been filtered out and replaced with N. Site: '+str(site)+'. Allele: '+r[site])
                quit()

        for site in lowq_to_test:
            if r[site] != 'N':
                print('This site has low quality. It should have been filtered out and replaced with N. Site: '+str(site)+'. Allele: '+r[site])
                quit()

        for site in a_to_test:
            if r[site] != 'A' and r[site] != 'N': # the site could be N if it was filtered out for being indel. This validate script doesnt check all the indels.
                print('This site has alternate allele. It should have been A. Site: '+str(site)+'. Allele: '+r[site])
                quit()

        for site in iupac_to_test:
            if r[site] == 'A' or r[site] == 'C' or r[site] == 'G' or r[site] == 'T':
                print('This site is heterozygous. It should have been had IUPAC code. Site: '+str(site)+'. Allele: '+r[site])
                quit()
        print('No issues with fasta file: '+sys.argv[1]+'_Chr'+sys.argv[2])
