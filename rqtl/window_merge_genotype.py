# Script to merge nearby windows with they have the same genotypes for all individuals
# tribeiro@wisc.edu

print('merging windows.')

import re, sys #arg1: input
import numpy as np

chrX_genotype = []
chr2L_genotype = []
chr2R_genotype = []
chr3L_genotype = []
chr3R_genotype = []

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
    if win[0] == "X":
        chrX_genotype.append(win)
    elif win[0] == "2L":
        chr2L_genotype.append(win)
    elif win[0] == "2R":
        chr2R_genotype.append(win)
    elif win[0] == "3L":
        chr3L_genotype.append(win)
    elif win[0] == "3R":
        chr3R_genotype.append(win)
    else:
        header.append(win)

win_start = chrX_genotype[0][1]
win_end = chrX_genotype[0][2]
merged_X = []
pos_check = 0
win_positions = [float(chrX_genotype[0][4])]
for i in range(0,len(chrX_genotype)-1):
    if chrX_genotype[i][5:] == chrX_genotype[i+1][5:]:
        win_end = chrX_genotype[i+1][2]
        win_positions.append(float(chrX_genotype[i+1][4]))
    else:
        merged_windows = [chrX_genotype[i][0], win_start, win_end] + [chrX_genotype[i][3]] + [str(sum(win_positions)/len(win_positions))] + chrX_genotype[i][5:]
        merged_X.append(merged_windows)
        win_start = chrX_genotype[i+1][1]
        win_end = chrX_genotype[i+1][2]
        win_positions = [float(chrX_genotype[i+1][4])]

if chrX_genotype[-2][5:] == chrX_genotype[-1][5:]:
    merged_windows = [chrX_genotype[-1][0], win_start, win_end] + [chrX_genotype[-1][3]] +[str(sum(win_positions)/len(win_positions))] + chrX_genotype[-1][5:]
    merged_X.append(merged_windows)

print(len(chrX_genotype))
print(len(merged_X))



win_start = chr2L_genotype[0][1]
win_end = chr2L_genotype[0][2]
merged_2L = []
win_positions = [float(chr2L_genotype[0][4])]
for i in range(0,len(chr2L_genotype)-1):
    if chr2L_genotype[i][5:] == chr2L_genotype[i+1][5:]:
        win_end = chr2L_genotype[i+1][2]
        win_positions.append(float(chr2L_genotype[i+1][4]))
    else:
        merged_windows = [chr2L_genotype[i][0], win_start, win_end] + [chr2L_genotype[i][3]] +[str(sum(win_positions)/len(win_positions))] + chr2L_genotype[i][5:]
        merged_2L.append(merged_windows)
        win_start = chr2L_genotype[i+1][1]
        win_end = chr2L_genotype[i+1][2]
        win_positions = [float(chr2L_genotype[i+1][4])]

if chr2L_genotype[-2][5:] == chr2L_genotype[-1][5:]:
    merged_windows = [chr2L_genotype[-1][0], win_start, win_end] + [chr2L_genotype[-1][3]] +[str(sum(win_positions)/len(win_positions))] + chr2L_genotype[-1][5:]
    merged_2L.append(merged_windows)

print(len(chr2L_genotype))
print(len(merged_2L))

win_start = chr2R_genotype[0][1]
win_end = chr2R_genotype[0][2]
merged_2R = []
win_positions = [float(chr2R_genotype[0][4])]
for i in range(0,len(chr2R_genotype)-1):
    if chr2R_genotype[i][5:] == chr2R_genotype[i+1][5:]:
        win_end = chr2R_genotype[i+1][2]
        win_positions.append(float(chr2R_genotype[i+1][4]))
    else:
        merged_windows = [chr2R_genotype[i][0], win_start, win_end] + [chr2R_genotype[i][3]] +[str(sum(win_positions)/len(win_positions))] + chr2R_genotype[i][5:]
        merged_2R.append(merged_windows)
        win_start = chr2R_genotype[i+1][1]
        win_end = chr2R_genotype[i+1][2]
        win_positions = [float(chr2R_genotype[i+1][4])]


if chr2R_genotype[-2][5:] == chr2R_genotype[-1][5:]:
    merged_windows = [chr2R_genotype[-1][0], win_start, win_end] + [chr2R_genotype[-1][3]] +[str(sum(win_positions)/len(win_positions))] + chr2R_genotype[-1][5:]
    merged_2R.append(merged_windows)

print(len(chr2R_genotype))
print(len(merged_2R))


win_start = chr3L_genotype[0][1]
win_end = chr3L_genotype[0][2]
merged_3L = []
win_positions = [float(chr3L_genotype[0][4])]
for i in range(0,len(chr3L_genotype)-1):
    if chr3L_genotype[i][5:] == chr3L_genotype[i+1][5:]:
        win_end = chr3L_genotype[i+1][2]
        win_positions.append(float(chr3L_genotype[i+1][4]))
    else:
        merged_windows = [chr3L_genotype[i][0], win_start, win_end] + [chr3L_genotype[i][3]] +[str(sum(win_positions)/len(win_positions))] + chr3L_genotype[i][5:]
        merged_3L.append(merged_windows)
        win_start = chr3L_genotype[i+1][1]
        win_end = chr3L_genotype[i+1][2]
        win_positions = [float(chr3L_genotype[i+1][4])]


if chr3L_genotype[-2][5:] == chr3L_genotype[-1][5:]:
    merged_windows = [chr3L_genotype[-1][0], win_start, win_end] + [chr3L_genotype[-1][3]] +[str(sum(win_positions)/len(win_positions))] + chr3L_genotype[-1][5:]
    merged_3L.append(merged_windows)

print(len(chr3L_genotype))
print(len(merged_3L))

win_start = chr3R_genotype[0][1]
win_end = chr3R_genotype[0][2]
merged_3R = []
win_positions = [float(chr3R_genotype[0][4])]
for i in range(0,len(chr3R_genotype)-1):
    if chr3R_genotype[i][5:] == chr3R_genotype[i+1][5:]:
        win_end = chr3R_genotype[i+1][2]
        win_positions.append(float(chr3R_genotype[i+1][4]))
    else:
        merged_windows = [chr3R_genotype[i][0], win_start, win_end] + [chr3R_genotype[i][3]] +[str(sum(win_positions)/len(win_positions))] + chr3R_genotype[i][5:]
        merged_3R.append(merged_windows)
        win_start = chr3R_genotype[i+1][1]
        win_end = chr3R_genotype[i+1][2]
        win_positions = [float(chr3R_genotype[i+1][4])]

if chr3R_genotype[-2][5:] == chr3R_genotype[-1][5:]:
    merged_windows = [chr3R_genotype[-1][0], win_start, win_end] + [chr3R_genotype[-1][3]] +[str(sum(win_positions)/len(win_positions))] + chr3R_genotype[-1][5:]
    merged_3R.append(merged_windows)

print(len(chr3R_genotype))
print(len(merged_3R))

output = []
out3R = []
for row in header:
    output.append(row)
    out3R.append(row)
for row in merged_X:
    n_row = ['_'.join(row[:3])] + row[3:]
    row = n_row
    output.append(row)
for row in merged_2L:
    n_row = ['_'.join(row[:3])] + row[3:]
    row = n_row
    output.append(row)
for row in merged_2R:
    n_row = ['_'.join(row[:3])] + row[3:]
    row = n_row
    output.append(row)
for row in merged_3L:
    n_row = ['_'.join(row[:3])] + row[3:]
    row = n_row
    output.append(row)
for row in merged_3R:
    n_row = ['_'.join(row[:3])] + row[3:]
    row = n_row
    output.append(row)
    out3R.append(row)

out_array = np.array(output)
transpose = out_array.T
outputting = transpose.tolist()

out3R_array = np.array(out3R)
transpose = out3R_array.T
outputting3R = transpose.tolist()

output = open(sys.argv[1][:-4]+'_merged.csv','w')
for row in outputting:
    output.write(','.join(row)+'\n')
output.close()
    
output = open(sys.argv[1][:-4]+'_merged_3R.csv','w')
for row in outputting3R:
    output.write(','.join(row)+'\n')
output.close()
