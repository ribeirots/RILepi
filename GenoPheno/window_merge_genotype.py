# Script to merge nearby windows with they have the same genotypes for all individuals
# tribeiro@wisc.edu

print('merging windows.')

import re, sys #arg1: input
import numpy as np

chr_genotype = []
header = 'empty'

gtyps = []
with open(sys.argv[1]) as genotypes:
    for win in genotypes:
        win = re.split(',',win[:-1])
        gtyps.append(win)

gt_array = np.array(gtyps)
transpose = gt_array.T

gt_list = transpose.tolist()

genotypes = []
for win in gt_list:
    w0 = re.split('_',win[0])
    new_w = w0+win[1:]
    genotypes.append(new_w)

header = []
for win in genotypes:
    if win[0] != 'X' and win[0] != '2L' and win[0] != '2R' and win[0] != '3L' and win[0] != '3R':
        header.append(win)
    else:
        chr_genotype.append(win)

win_start = chr_genotype[0][1]
win_end = chr_genotype[0][2]
merged = []
pos_check = 0
for i in range(0,len(chr_genotype)-1):
    if chr_genotype[i][3:] == chr_genotype[i+1][3:]:
        win_end = chr_genotype[i+1][2]
    else:
        merged_windows = [chr_genotype[i][0], win_start, win_end] + chr_genotype[i][3:]
        merged.append(merged_windows)
        win_start = chr_genotype[i+1][1]
        win_end = chr_genotype[i+1][2]

if chr_genotype[-2][3:] == chr_genotype[-1][3:]:
    merged_windows = [chr_genotype[-1][0], win_start, win_end] + chr_genotype[-1][3:]
    merged.append(merged_windows)

print(len(chr_genotype))
print(len(merged))


output = []

for row in header:
    output.append(row)
for row in merged:
    n_row = ['_'.join(row[:3])] + row[3:]
    row = n_row
    output.append(row)

out_array = np.array(output)
transpose = out_array.T
outputting = transpose.tolist()

output = open(sys.argv[1][:-4]+'_merged.csv','w')
for row in outputting:
    output.write(','.join(row)+'\n')
output.close()
    
