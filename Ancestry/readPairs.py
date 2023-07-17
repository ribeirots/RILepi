# script to cound info depth per pool
# Args 1 and 2: parental Low and High sites.vcf files
# Args 3 and 4: offspring pool Low and High sam file prefix
# Arg 5: window
import re, sys
from infordepth_functions import *
from parentalVCF_dict import *

# for each chrm, gets the VCF dictionary and the sam files to create output with info depths
# still not designed to include all the chrms
# results in the creation of 2 dictionaries: readDict and read_dictH with all the read pairs from each offspring pool for each chrm

# low or highSam (path to sam files for each pool, sam file for the chrm to be used and not from the whole alignment)
# call the func one pool at a time, the returned SNPdict is pool and chrm specific
def SNPdict(samOffspring, dictToUse, screenLoad=0): # screenLoad determines whether to display loading % or not
    totalMsize = 0
    totalSNPsize = 0
    Msize = 0
    SNPsize = 0
    readDict = {}
    readTotalCount = 0
    lowQualCount = 0
    read1count = 0
    read2count = 0
    read3Pcount = 0
    noSNP = 0
    atLeast1SNP = 0
    blackKeyList = [] # after a read ID is remove from the set of keys for containing more than 2 reads (likely a chimera/supplementary alignment) the ID is saved here to ensure it is not used again
    with open(samOffspring) as samf:
        for r in samf:
            readTotalCount += 1
            if screenLoad:
                z = round(readTotalCount / 5000000, 2)
                print(z, end='\r')
            r = re.split('\t',r[:-1])
            if ( abs(int(r[7]) - int(r[3])) <= 10000) and (float(r[4]) >= 20): # r[4] = Quality
                if r[0] not in readDict.keys():
                    if r[0] not in blackKeyList:
                        read1count += 1
                        readDict[r[0]] = [[],[]] # [[read1], [read2]]
                        readDict[r[0]][0], Msize, SNPsize = ancestry_information(r, dictToUse)
                        totalMsize += Msize
                        totalSNPsize += SNPsize
                else:
                    # when this is false (goes into the following else) there are more than 2 reads with the same ID in the sam file, it is likely a case of supplementary alignment
                    if len(readDict[r[0]][1]) == 0:
                        read2count += 1
                        readDict[r[0]][1], Msize, SNPsize = ancestry_information(r, dictToUse)
                        totalMsize += Msize
                        totalSNPsize += SNPsize
                        if len(readDict[r[0]][0][1]) == 0 and len(readDict[r[0]][1][1]) == 0:
                            noSNP += 1
                            if 1 == 0:
                                print(noSNP)
                                print(readDict[r[0]])
                                print(r)
                                if noSNP == 3:
                                    quit()
                        elif len(readDict[r[0]][0][1]) != 0 or len(readDict[r[0]][1][1]) != 0:
                            atLeast1SNP += 1
                    else:
                        read3Pcount += 1
                        blackKeyList.append(r[0])
                        del readDict[r[0]]
            else:
                # this if is true if the 2nd read in the pair is the one violating the distance or lowQual criteria
                if r[0] in readDict.keys():
                    del readDict[r[0]]
                lowQualCount += 1
    print('##### Pool of SNPs analized:')
    print('Total Read Count: '+str(readTotalCount))
    print('Read1 Count: '+str(read1count))
    print('Read2 Count: '+str(read2count))
    print('Supplementary Alignment encountered: '+str(read3Pcount))
    print('Low Quality (filtered out): '+str(lowQualCount))
    print('Read pairs with at least one SNP: '+str(atLeast1SNP))
    print('Read pairs without any SNP: '+str(noSNP))
    print('Final number of read pairs to be analized for informative depth: '+str(len(readDict.keys())))
    meanMsize = totalMsize/(read1count + read2count)
    meanSNPsize = totalSNPsize/(read1count + read2count)
    print('Average read size (considering correctly mapped bases only, removing indels and H/S clips): '+str(meanMsize))
    print('Average number of SNPs per read: '+str(meanSNPsize))
    print('#################################')
    
    return readDict

#use read_dicts to count votes per windows
def windowCountVote(windows, readDictPool, screenLoad=0):
    print('Counting votes.\n')
    missingPairL = 0
    noSNP_L2 = 0
    vFull1 = 0
    vFull2 = 0
    v1Read = 0
    vOutOfB = 0
    vWinPair = 0
    multiWinNoVote = 0
    noWin = 0
    kCount = 0
    
    for k in readDictPool.keys():
        kCount += 1
        if screenLoad:
            screenLoading = kCount/len(readDictPool.keys())
            print(screenLoading, end='\r')
        r1Start = int(readDictPool[k][0][0])
        wFoundCheck = 0
        for i in range(0, len(windows)):
            if wFoundCheck == 1:
                break

            w1Start = int(windows[i][0])
            w1End = int(windows[i][1])
            
            if i == (len(windows)-1):
                w2Start = 0
                w2End = 0
            else:
                w2Start = int(windows[i+1][0])
                w2End = int(windows[i+1][1])

            if r1Start >= w1Start and r1Start <= w1End:
                wFoundCheck = 1
                read1 = readDictPool[k][0]
                read2 = readDictPool[k][1]
                if len(read1) > 0 and len(read2) > 0: # Has both reads from paired end reads
                    if len(read1[1]) > 0 and len(read2[1]) > 0:
                        ##### what if one of the reads dont have SNPs? need to account for that in the code
                        if read2[1][-1][0] <= w1End: # -1 = lastSNP: [Pos,Vote,Weight]
                            # fully in window 1: read result is added to window 1
                            vote = voteFullPair(read1,read2)
                            vFull1 += 1
                            windows[i] = winVote(windows[i], vote)
                        elif read1[1][0][0] >= w2Start and read2[1][-1][0] <= w2End: # read1[1] = 1st SNP read1
                            # all SNPs in window 2: read result is added to window 2
                            vote = voteFullPair(read1,read2)
                            vFull2 += 1
                            windows[i+1] = winVote(windows[i+1], vote)
                        elif read2[1][-1][0] > w2End:
                            # do nothing - read two out of bounds, beyond window 2
                            vOutOfB += 1
                            vote = 'x'
                        else:
                            # will count how many snps are voting for each window and the vote weight
                            vWinPair += 1
                            w11 = 0
                            w105 = 0
                            w21 = 0
                            w205 = 0
                            for read in [read1, read2]:
                                for snp in read[1]:
                                    if snp[0] <= w1End:
                                        if snp[2] == 1:
                                            w11 += 1
                                        elif snp[2] == 0.5:
                                            w105 += 1
                                    else:
                                        if snp[2] == 1:
                                            w21 += 1
                                        elif snp[2] == 0.5:
                                            w205 += 1
                            if w11 + w21 == 0:
                                if w105 + w205 != 0: # else: do nothing, no window wins
                                    if w105 >= w205:
                                        # w1 wins
                                        vote = votePartialWindow(read1, read2, w1Start, w1End)
                                        windows[i] = winVote(windows[i], vote)
                                    else:
                                        # w2 wins
                                        vote = votePartialWindow(read1, read2, w2Start, w2End)
                                        windows[i] = winVote(windows[i], vote)
                                else:
                                    multiWinNoVote += 1 # read across multiple windows without any window "winning" the vote to count as the one receiving the vote
                            elif w11 >= w21:
                                # w1 wins
                                vote = votePartialWindow(read1, read2, w1Start, w1End)
                                windows[i] = winVote(windows[i], vote)
                            elif read2[1][-1][0] > w2End:
                                vOutOfB += 1
                            else:
                                # w2 wins
                                vote = votePartialWindow(read1, read2, w2Start, w2End)
                                windows[i+1] = winVote(windows[i+1], vote)
                    elif len(read1[1]) > 0: # only fails the first if when read1 or 2 have no SNPs, so it is possible that only one of them have SNPs.
                        vote = voteFullPair(read1,read2)
                        v1Read += 1
                        windows[i] = winVote(windows[i], vote)
                    elif len(read2[1]) > 0:
                        vote = voteFullPair(read1,read2)
                        v1Read += 1
                        if read2[1][0][0] <= w1End or i == (len(windows)-1):
                            windows[i] = winVote(windows[i], vote)
                        elif read2[1][-1][0] <= w2End:
                            windows[i+1] = winVote(windows[i+1], vote)
                    else: # neither read1 or read2 has any SNPs
                        noSNP_L2 += 1 # the other noSNP is inferred during the creating of the dictionaries, not when summarizing windows
                        vote = 'NoSNP'
                        if (read1[0] + read2[0])/2 <= w1End or i == (len(windows)-1):
                            windows[i] = winVote(windows[i], vote)
                        else:
                            windows[i+1] = winVote(windows[i+1], vote)
                else: # read1 or read2 == 0
                    missingPairL += 1
        if wFoundCheck == 0:
            noWin += 1
            print('readpair outside any window (was the last win analyzed?):')
            print(read1)
            print(read2)
            quit()


    print('\n')
    print('Missing Pair L')
    print(missingPairL)
    print('noSNP_L2:')
    print(noSNP_L2)

    print('Read pair votes cast: Fully on Win1, Fully on Win2, using only 1 read, OutOfBounds, CheckWinPair:')
    print(vFull1)
    print(vFull2)
    print(v1Read)
    print(vOutOfB)
    print(vWinPair)

    print('multiWin No Vote:')
    print(multiWinNoVote)

    print('No Window')
    print(noWin)

    return windows
