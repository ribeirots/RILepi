# functions needed to calculate informative depth
import re, sys

# Will return the pool vote (Low-Parent, High-Parent, uNknown) and the vote weight (1 - certain pool identity, 0.5 - heterozygous)
# baseF = offspring base that will be compared against the parental bases.
# basePL and basePH: bases from Low-Parent and High-Parent respectively, could be a single base (homozygous) or a list with two bases (heterozygous)
def baseVote(baseF, basePL, basePH):
    if len(basePL) == 1 and len(basePH) == 1:
        if baseF == basePL:
            return ['L',1]
        elif baseF == basePH:
            return ['H',1]
        else:
            return ['N',0]
    elif len(basePL) == 1:
        if baseF == basePL and baseF not in basePH:
            return ['L',1]
        elif baseF == basePL and baseF in basePH:
            return ['L',0.5]
        elif baseF in basePH:
            return ['H', 1]
        else:
            return ['N',0]
    elif len(basePH) == 1:
        if baseF == basePH and baseF not in basePL:
            return ['H',1]
        elif baseF == basePH and baseF in basePL:
            return ['H',0.5]
        elif baseF in basePL:
            return ['L', 1]
        else:
            return ['N',0]
    elif len(basePL) == 2 and len(basePH) == 2:
        if baseF in basePL and baseF not in basePH:
            return ['L',1]
        elif baseF in basePH and baseF not in basePL:
            return ['H',1]
        else:
            return ['N',0]
    else:
        return ['N',0]

# Will return a summary of the votes from a given read, including position of the SNPs within the read and their votes
# parses read sequence in accordance with its CIGAR string information
# seq = DNA sequence
# start = starting position of the sequence
# cigar = CIGAR string
# chrm_dict = dictionary with all the SNPs in a chromosome, the key is the position in the chromosome and the content are the genotypes from the Low-Parent and the High-Parent (in this order)
def ancestry_information(inputRead, chrm_dict):
    seq = inputRead[9]
    start = int(inputRead[3])
    cigar = inputRead[5]

    bases = re.split('',seq)[1:-1]
    cigList = re.split('(\d+)',cigar)[1:]
    cigInt = cigList[::2]
    cigInt = list(map(int,cigInt))
    cigCode = cigList[1::2]
    readIndex = 0
    refIndex = 0
    Msize = 0
    SNPsize = 0
    readVotes = []
    for cigI, cigC in zip(cigInt, cigCode):
#        print(cigI)
#        print(cigC)
        if cigC == 'M':
            for bpIndex in range(0,int(cigI)):
#                walk normally on both looping I
                Msize += 1
                b = bases[readIndex]
                ref = int(start)+refIndex
                if str(ref) in chrm_dict.keys():
                    SNPsize += 1
                    pair = chrm_dict[str(ref)]
#                    print(int(start))
#                    print(refIndex)
#                    print(b)
#                    print(readIndex)
#                    print(bases)
#                    print(pair)
#                    print(baseVote(b,pair[0],pair[1]))
                    readVotes.append([int(start)+refIndex]+baseVote(b,pair[0],pair[1]))
#                    print(readVotes)
#                    print(b)
#                    print(pair)
                readIndex += 1
                refIndex += 1
        elif cigC == 'I' or cigC == 'S':
#                walk read, not ref
            for bpIndex in range(0,int(cigI)):
                readIndex += 1
        elif cigC == 'D':
#                walk ref, not read
            for bpIndex in range(0,int(cigI)):
                refIndex += 1
        elif cigC != 'H': # if it is H, just loops back to the start, it is a hard clip and the sites have been removed from this read already
            print("CIGAR code: "+cigC+" is not recognized. Quitting script.")
            print(cigInt)
            print(cigCode)
            print(seq)
            quit()
    return [[int(start),readVotes], Msize, SNPsize] # [read=[pos, [snp,snpsnp]] , Msize = properly aligned bases] # previous version had "+refIndex for the position, I don't know why.

# Will return the final vote of the read pair, based on the summary of how the SNPs in each read in the pair voted. Assigning the pair to Low-Parent, High-Parent, or Unknown.
# Only used if both reads are fully in one window
# r1 and r2 are the summaries for each read in the pair
def voteFullPair(r1,r2):
    vote1_L = 0
    vote1_H = 0
    vote05_L = 0
    vote05_H = 0
    if len(r1[1]) > 0:
        for snp in r1[1]:
            if snp[1] == 'H' and snp[2] == 1:
                vote1_H += 1
            elif snp[1] == 'H' and snp[2] == 0.5:
                vote05_H += 1
            elif snp[1] == 'L' and snp[2] == 1:
                vote1_L += 1
            elif snp[1] == 'L' and snp[2] == 0.5:
                vote05_L += 1
    if len(r2[1]) > 0:
        for snp in r2[1]:
            if snp[1] == 'H' and snp[2] == 1:
                vote1_H += 1
            elif snp[1] == 'H' and snp[2] == 0.5:
                vote05_H += 1
            elif snp[1] == 'L' and snp[2] == 1:
                vote1_L += 1
            elif snp[1] == 'L' and snp[2] == 0.5:
                vote05_L += 1
    if vote1_H > 0 or vote1_L > 0:
        if vote1_H > vote1_L:
            return 'H'
        elif vote1_H < vote1_L:
            return 'L'
        else:
            return 'N'
    elif vote05_H > 0 and vote05_L == 0: # these will count as half informative in the final count
        return 'h'
    elif vote05_H == 0 and vote05_L > 0:
        return 'l'
    else:
        return 'N'

# Will return the final vote of the read pair, based on the summary of how the SNPs in each read in the pair voted. Assigning the pair to Low-Parent, High-Parent, or Unknown.
# Only used if both reads are NOT in one window only. It is, they spam two or more windows.
# r1 and r2 are the summaries for each read in the pair
# wStart and wEnd, the boundaries of the window that will be considered for the read pair final vote
def votePartialWindow(r1,r2,wStart,wEnd):
    vote1_L = 0
    vote1_H = 0
    vote05_L = 0
    vote05_H = 0
    if len(r1[1]) > 0:
        for snp in r1[1]:
            if snp[0] >- wStart and snp[0] <= wEnd:
                if snp[1] == 'H' and snp[2] == 1:
                    vote1_H += 1
                elif snp[1] == 'H' and snp[2] == 0.5:
                    vote05_H += 1
                elif snp[1] == 'L' and snp[2] == 1:
                    vote1_L += 1
                elif snp[1] == 'L' and snp[2] == 0.5:
                    vote05_L += 1
    if len(r2[1]) > 0:
        for snp in r2[1]:
            if snp[0] >- wStart and snp[0] <= wEnd:
                if snp[1] == 'H' and snp[2] == 1:
                    vote1_H += 1
                elif snp[1] == 'H' and snp[2] == 0.5:
                    vote05_H += 1
                elif snp[1] == 'L' and snp[2] == 1:
                    vote1_L += 1
                elif snp[1] == 'L' and snp[2] == 0.5:
                    vote05_L += 1
    if vote1_H > 0 or vote1_L > 0:
        if vote1_H > vote1_L:
            return 'H'
        elif vote1_H < vote1_L:
            return 'L'
        else:
            return 'N'
    elif vote05_H > 0 and vote05_L == 0:
        return 'h'
    elif vote05_H == 0 and vote05_L > 0:
        return 'l'
    else:
        return 'N'

# returns genotype for a given site in a vcf file.
# can return just one single nucleotide code if homozygous, or a list with two nucleotides if heterozygoues
# will identify genotype with best quality if multiple are present
def gtParentalVCF(vcfSite): # a row of the vcf file
    genotypes = vcfSite[9:] # e.g. list of 0/0:10 geno:count
    gtCode = vcfSite[8] # series of codes separated by :, including DP which specifies the column with count in the information on vcfSite[9] shown above
    gtCode = re.split(':',gtCode)
    refAlt = vcfSite[3:5] # [ref, alt] alleles
    genotype = -1

    if len(vcfSite[4]) == 3: # this means 2 alleles that don't match the reference allele, in this case it will those as the genotype for this site
        altHt = re.split(',',vcfSite[4])
        if altHt[0] in ['A','T','C','G'] and altHt[1] in ['A','T','C','G']:
            refAlt = altHt
        else:
            print('Check this genotype:  '+altHt)
            print(vcfSite)
            quit()
            return -1
    
    gqual = -1
    genoCode = -1
    run_multiple = 1 # will turn to 0 when there were multiple alignments in the one being checked isnt better than the ones before
    for gt in genotypes: # it can be multiple when there were multiple bam files
        gt = re.split(':',gt)
        if len(genotypes) > 1:
            if gt[0] != './.' and int(gt[gtCode.index('DP')]) > gqual: # highest genotype quality
                gqual = int(gt[gtCode.index('DP')])
                genoCode = gt[0]
            else:
                run_multiple = 0
        elif gt[0] != './.':
            genoCode = gt[0]

        if genoCode != -1 and run_multiple == 1:
            if genoCode == '0/0':
                genotype = refAlt[0]
            elif genoCode == '1/1':
                genotype = refAlt[1]
            else:
                if 'A' in refAlt and 'T' in refAlt:
                    genotype = ['A','T']
                elif 'A' in refAlt and 'C' in refAlt:
                    genotype = ['A','C']
                elif 'A' in refAlt and 'G' in refAlt:
                    genotype = ['A','G']
                elif 'T' in refAlt and 'C' in refAlt:
                    genotype = ['T','C']
                elif 'T' in refAlt and 'G' in refAlt:
                    genotype = ['T','G']
                elif 'C' in refAlt and 'G' in refAlt:
                    genotype = ['C','G']
                else:
                    genotype = -1
    return genotype

# function to assign the vote to a window
def winVote(window, SNPvote):
    if SNPvote == 'H':
        window[-1] += 1
    elif SNPvote == 'L':
        window[-2] += 1
    elif SNPvote == 'h': # lower case = half informative
        window[-1] += 0.5
    elif SNPvote == 'l':
        window[-2] += 0.5
    elif SNPvote == 'NoSNP':
        window[-4] += 1 
    else: # has SNP but it is not informative to distinguish ancestry
        window[-3] += 1

    return window
