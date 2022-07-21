# the goal of this script is to check how many windows with missing data are present for each RIL.
# tribeiro@wisc.edu

# arg1: input
# arg2: output

import re, sys

rils=[]
ms_data = []

first_row = 1
total_row_count = 0
with open(sys.argv[1]) as rdata:
    for r1 in rdata:
        total_row_count += 1
        r = re.split('\t',r1[:-1])[3:]
        if first_row == 1:
            first_row = 0
            for ril in r:
                rils.append(ril)
                ms_data.append(0)
        else:
            for ril_i in range(0,len(rils)):
                if r[ril_i] == '-999':
                    ms_data[ril_i] += 1

output = open(sys.argv[2],'w')
for ril_i in range(0,len(rils)):
    outting = [rils[ril_i], str(ms_data[ril_i]), str(ms_data[ril_i]/total_row_count*100)+'%\n']
    output.write('\t'.join(outting))
    if ms_data[ril_i]/total_row_count > 0.1:
        print(rils[ril_i])
output.close()

