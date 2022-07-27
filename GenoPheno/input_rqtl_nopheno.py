# prepare data for R/qtl
# tribeiro@wisc.edu
print('rqtl input - no phenotypes.')
import re, sys # 1: input

with open(sys.argv[1]) as original:
    genotype = [x.split() for x in original]

z = []
for x in zip(*genotype):
    new_row = []
    for y in x:
        new_row.append(y)
    z.append(new_row)

chrom = ['']
for x in z[0]:
    if x == 'X':
        chrom.append('3') # 3 = X
    elif x == '2L':
        chrom.append('4') # 4 = 2L
    elif x == '2R':
        chrom.append('4') # 4 = 2R
    elif x == '3L':
        chrom.append('5') # 5 = 3L
    elif x == '3R':
        chrom.append('5') # 5 = 3R
#print(chrom[0:10])

labels = ['RIL']

for i in range(1,len(z[0])):
    labels.append(str(z[0][i])+"_" +str(z[1][i])+"_" +str(z[2][i]))


new_file = []
new_file.append(labels)
new_file.append(chrom)

for i in range(3,len(z)):
    new_file.append(z[i])


for i in range(2,len(new_file)):
    for j in range(1,len(new_file[i])):
        if new_file[i][j] == '2':
            new_file[i][j] = 'ZZ'
        elif new_file[i][j] == '1':
            new_file[i][j] = 'FZ'
        elif new_file[i][j] == '0':
            new_file[i][j] = 'FF'
        else:
            new_file[i][j] = 'NA'

output = open(sys.argv[1][:-4]+'_rqtl.csv','w')
for line in new_file:
    output.write(r','.join(line)+'\n')
output.close()

        
