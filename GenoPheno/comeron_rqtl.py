# script to include recombination distances in the rqtl input file
# tribeiro@wisc.edu

import re, math, sys # 1: input
print('rqtl - comeron cM.')
with open(sys.argv[1]) as file1:
    markers = file1.readline()
    markers = re.split(',',markers)

for i in range(2,len(markers)):
    markers[i] = re.split('_',markers[i])
    markers[i][1] = int(markers[i][1])
    markers[i][2] = int(markers[i][2])
    markers[i].append((markers[i][2]-markers[i][1] )/2 +markers[i][1] )

comeron = {}
with open('recomb_Comeron_2012_0001_FR_cM.csv') as com:
    next(com)
    for r in com:
        r = re.split(',',r[:-1])
        r[0] = r[0][3:] # chrm
        r[-1] = float(r[-1]) # increment perbase to position
        r[-4] = float(r[-4]) # position
        comeron[r[0]+'_'+r[1]] = [r[-1],r[-4]] #eg. key="X_1", value = [perbaseIncrement, position]

## Loop to find the position of the empirical window's midpoint
cM_output = ['']
for i in range(1,len(markers)):
    pos_bp = float(markers[i][-1])
    chrm = markers[i][0]
    pos1 = (math.floor(pos_bp/100000)+1)*100000-99999
    pos_key = chrm+'_'+str(pos1)
    perbase_incr = comeron[pos_key][0]
    pos_cM_incr = perbase_incr * (pos_bp - pos1)
    pos_cM = pos_cM_incr + comeron[pos_key][1]
    cM_output.append(str(pos_cM))

cM_output = ','.join(cM_output)+'\n'

output = open(sys.argv[1][:-4]+'_cMorg.csv','w')
row_check = 0
with open(sys.argv[1]) as file1:
    for r in file1:
        row_check += 1
        if row_check == 3:
            output.write(cM_output)
        output.write(r)

output.close()
