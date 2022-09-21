# script to split chrm 2 and 3 in 3 parts for r/qtl
# tribeiro@wisc.edu

# Arg 1: input file

import re, sys
import numpy as np

out_filename = sys.argv[1][:-4]+'_7armsExtended.csv'

fullfile = []
with open(sys.argv[1]) as inpfile:
    for row in inpfile:
        row = re.split(',',row[:-1])
        fullfile.append(row)

countx = fullfile[1].count('1')
count2 = fullfile[1].count('2')
tracks2 = count2//3
count3 = fullfile[1].count('5')
tracks3 = count3//3

#print([countx,count2,count3,tracks2,tracks3])

i_3start = 2+countx+tracks2 # 2 empty, entire X, the track remaining as 2
i_3end = i_3start + tracks2 # whole track as 3

i_4start = i_3end 
i_4end = 2+countx+count2

i_6start = i_4end+tracks3
i_6end = i_6start + tracks3

i_7start = i_6end

#print([i_3start, i_3end, i_4start, i_4end, i_6start, i_6end, i_7start])
# just to check the indexes are right visually
if 0:
    print(fullfile[1][i_3start]) # 2
    print(fullfile[1][i_3end]) # 2
    print(fullfile[1][i_4start]) # 2
    print(fullfile[1][i_4end-1]) # 2
    print(fullfile[1][i_4end]) # 3
    print(fullfile[1][i_6start]) # 3
    print(fullfile[1][i_6end]) # 3
    print(fullfile[1][i_7start]) # 3

# show marker in the boudary regions - to manually check their rec_rate. Shown here for FR avoidance, minor changes depending on dateset based on merging neighboring windows variations because of different group of phenotyped RILs.
if 0:
    print(fullfile[0][i_3start]) # 2L_11082577_11097392. c (cM/Mb) = 3.987352355437420 (arm avg: 2.395)
    print(fullfile[0][i_3end]) # 2R_7878995_7893629. c (cM/Mb) = 2.928532960423250 (arm avg: 2.662)
    print(fullfile[0][i_6start]) # 3L_13391777_13410079. c (cM/Mb) = 1.522968800129060 (arm avg: 1.792)
    print(fullfile[0][i_6end]) # 3R_13667834_13679222. c (cM/Mb) = 0.544070697079204 (arm avg: 1.963) (this comeron window is surrounded by high c - 3.26 and 1.96)

for i in range(i_3start,i_3end):
    fullfile[1][i] = '3'
for i in range(i_4start,i_4end):
    fullfile[1][i] = '4'
for i in range(i_6start,i_6end):
    fullfile[1][i] = '6'
for i in range(i_7start,len(fullfile[1])):
    fullfile[1][i] = '7'

# check visually the boundaries of the chrm changes and cMorg changes
if 0:
    print(fullfile[1][i_3start-1:i_3start+1])
    print(fullfile[1][i_3end-1:i_3end+1])
    print(fullfile[1][i_4start-1:i_4start+1])
    print(fullfile[1][i_4end-1:i_4end+1])
    print(fullfile[1][i_6start-1:i_6start+1])
    print(fullfile[1][i_6end-1:i_6end+1])

output = open(out_filename,'w')
incr=50 # number of windows to be added expanding the chrm 2 and 3 subsets
row_count = 0
for row in fullfile:
    row_count += 1
    sub1 = row[i_3start:i_3start+incr] # add to end of 2
    sub2 = row[i_3start-incr:i_3start] # add to start of 3
    sub3 = row[i_3end:i_3end+incr] # add to end of 3
    sub4 = row[i_3end-incr:i_3end] # add to start of 4

    sub5 = row[i_6start:i_6start+incr] # add to end of 5
    sub6 = row[i_6start-incr:i_6start] # add to start of 6
    sub7 = row[i_6end:i_6end+incr] # add to end of 6
    sub8 = row[i_6end-incr:i_6end] # add to start of

    if row_count == 1:
        sub1 = list(map(lambda x: x+'_ext', sub1))
        sub2 = list(map(lambda x: x+'_ext', sub2))
        sub3 = list(map(lambda x: x+'_ext', sub3))
        sub4 = list(map(lambda x: x+'_ext', sub4))
        sub5 = list(map(lambda x: x+'_ext', sub5))
        sub6 = list(map(lambda x: x+'_ext', sub6))
        sub7 = list(map(lambda x: x+'_ext', sub7))
        sub8 = list(map(lambda x: x+'_ext', sub8))


    elif row_count == 2: # change the chrm label for the subsets to be added to the extensions, to make them count as part of the right chrm
        # check if the correct size was sliced into a new subset AND if right chrm identities will be modified
        if sub1.count('3') != incr:
            print('Pre Extension increment sub1 wrong. Breaking script!!!')
            quit()
        if sub2.count('2') != incr:
            print('Pre Extension increment sub2 wrong. Breaking script!!!')
            quit()
        if sub3.count('4') != incr:
            print('Pre Extension increment sub3 wrong. Breaking script!!!')
            quit()
        if sub4.count('3') != incr:
            print('Pre Extension increment sub4 wrong. Breaking script!!!')
            quit()
        if sub5.count('6') != incr:
            print('Pre Extension increment sub5 wrong. Breaking script!!!')
            quit()
        if sub6.count('5') != incr:
            print('Pre Extension increment sub6 wrong. Breaking script!!!')
            quit()
        if sub7.count('7') != incr:
            print('Pre Extension increment sub7 wrong. Breaking script!!!')
            quit()
        if sub8.count('6') != incr:
            print('Pre Extension increment sub8 wrong. Breaking script!!!')
            quit()

        sub1 = list(map(lambda x: x.replace('3', '2'), sub1))
        sub2 = list(map(lambda x: x.replace('2', '3'), sub2))
        sub3 = list(map(lambda x: x.replace('4', '3'), sub3))
        sub4 = list(map(lambda x: x.replace('3', '4'), sub4))
        sub5 = list(map(lambda x: x.replace('6', '5'), sub5))
        sub6 = list(map(lambda x: x.replace('5', '6'), sub6))
        sub7 = list(map(lambda x: x.replace('7', '6'), sub7))
        sub8 = list(map(lambda x: x.replace('6', '7'), sub8))

        # check if the correct size was added AND if right chrm identities were modified
        if sub1.count('2') != incr:
            print('Extension increment sub1 wrong. Breaking script!!!')
            quit()
        if sub2.count('3') != incr:
            print('Extension increment sub2 wrong. Breaking script!!!')
            quit()
        if sub3.count('3') != incr:
            print('Extension increment sub3 wrong. Breaking script!!!')
            quit()
        if sub4.count('4') != incr:
            print('Extension increment sub4 wrong. Breaking script!!!')
            quit()
        if sub5.count('5') != incr:
            print('Extension increment sub5 wrong. Breaking script!!!')
            quit()
        if sub6.count('6') != incr:
            print('Extension increment sub6 wrong. Breaking script!!!')
            quit()
        if sub7.count('6') != incr:
            print('Extension increment sub7 wrong. Breaking script!!!')
            quit()
        if sub8.count('7') != incr:
            print('Extension increment sub8 wrong. Breaking script!!!')
            quit()
        print('Subset for extensions made successfully.')

    n_row = row[:i_3start] + sub1 + sub2 + row[i_3start:i_3end] + sub3 + sub4 + row[i_3end:i_6start] + sub5 + sub6 + row[i_6start:i_6end] + sub7 + sub8 + row[i_6end:]

    output.write(r','.join(n_row)+'\n')
output.close()

