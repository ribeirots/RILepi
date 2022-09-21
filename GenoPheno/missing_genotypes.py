# script to replace -999 genotypes on ancestry calls with an inferred genotyped.
# tribeiro@wisc.edu


# arg1 = input file
# arg2 = output file
import re, sys



gap_size = 10

windows = []
with open(sys.argv[1]) as genotypes:
    for row in genotypes:
        row = re.split('\t',row[:-1])
        windows.append(row)

for l in range(3,len(windows[1])): # loop through RIL lines
    for w in range(1,len(windows)): # loop through windows
        if '-999' in windows[w][l]: # found the missing data window
            if windows[w-1][l] == windows[w+1][l]: # if surrounded by same genotype
                windows[w][l] = windows[w-1][l]
            elif windows[w-1][0] != windows[w+1][0]: # if surrounded by windows in different chromosomes, use the one in which the missing data is
                if windows[w-1][0] == windows[w][0]:
                    windows[w][l] = windows[w-1][l]
                else:
                    windows[w][l] = windows[w+1][l]
            elif windows[w+1][l] == "-999":
                large_gap = 1
                while windows[w+large_gap][l] == "-999" and windows[w+large_gap][0] == windows[w][0]:
                    large_gap = large_gap + 1
                # large_gap = the window after the sequence of -999s
                if windows[w+large_gap][l] == windows[w-1][l] or windows[w+large_gap][0] != windows[w][0]: # if surrounded by same genotype or ending at different arm
                    if large_gap <= gap_size:
                        for fill_gap in range(0,large_gap):
                            windows[w+fill_gap][l] = windows[w-1][l]


output = open(sys.argv[2],'w')
for row in windows:
    output.write('\t'.join(row)+'\n')
output.close()
