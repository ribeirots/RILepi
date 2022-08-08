# script to calculate cM positions based on comeron file for RILs - using classes for rqtl windows
# argv 1: input file
# tribeiro@wisc.edu

import re, sys
comeron_file = 'recomb_Comeron_2012_for_rils.csv'


# read comeron file and create a dictionary with Chr_Start_End as key, cM_F12 and perbase_increment as values
comeron = {}
with open(comeron_file) as com:
    next(com)
    for r in com:
        r = re.split(',',r[:-1])
        comkey = r[0]+'_'+str(int(r[1]) - 1)+'_'+str(int(r[2])-1) # key example: ChrX_0_99999, I made it 0-based 
        comvalue = [float(r[7]), float(r[8])] # value example: [0, 1.1E-09]
        comeron[comkey] = comvalue

# function to find the cM position of the midpoint of a window
def cM_midpoint(chrm, start, end): # Rqtl input is zero-based, comeron was not but I modified it when creating the dict above.

    # find start position
    key_start = start // 100000 * 100000  # comeron dict is not 0 based, start at 1, 100001, 200001, etc.
    key_end = key_start  + 99999 # final position of the comeron window
    rqtl_key = 'Chr'+chrm + '_' + str(key_start) + '_' + str(key_end)
    bp_dist = start - key_start # distance from the start of the comeron window
    start_pos_cM = comeron[rqtl_key][0] + comeron[rqtl_key][1] * bp_dist

    # find end position
    key_start = end // 100000 * 100000  # comeron dict is not 0 based, start at 1, 100001, 200001, etc.
    key_end = key_start  + 99999 # final position of the comeron window
    rqtl_key = 'Chr'+chrm + '_' + str(key_start) + '_' + str(key_end)
    bp_dist = end - key_start # distance from the start of the comeron window
    end_pos_cM = comeron[rqtl_key][0] + comeron[rqtl_key][1] * bp_dist

    len_cM = end_pos_cM - start_pos_cM
    half_len_cM = len_cM / 2
    mid_point_cM = start_pos_cM + half_len_cM

    return mid_point_cM

# read the list of markers
with open(sys.argv[1]) as file1:
    markers = file1.readline()
    markers = re.split(',',markers)

# split markers in Chr#, Start, End
for i in range(2,len(markers)):
    markers[i] = re.split('_',markers[i])
    markers[i][1] = int(markers[i][1])
    markers[i][2] = int(markers[i][2])

# loop to create the row with cM positions in the rqtl input file
cM_output = ['','']
order_check = []
lastChrm = 'X'
for marker in markers[2:]:
    midcM = cM_midpoint(marker[0], marker[1], marker[2])
    cM_output.append(str(midcM))

    # this loop is present to check if the midpoint positions for all the windows in the same chrm are properly ordered
    if marker[0][0] != lastChrm:
        lastChrm = marker[0][0] # indexing the marker[0] to only get the [0] to focus on chrm 2L and 2R as a single chrm.
        print(lastChrm)
        for i in range(0,len(order_check)-1):
            if order_check[i] > order_check[i+1]:
                print('Order check Error.')
                quit()
        order_check = []
    order_check.append(midcM)

cM_output = ','.join(cM_output)+'\n'

# output file including the cM positions on the 3rd row
output = open(sys.argv[1][:-4]+'_cMorg.csv','w')
row_check = 0
with open(sys.argv[1]) as file1:
    for r in file1:
        row_check += 1
        if row_check == 3:
            output.write(cM_output)
        output.write(r)
output.close()
