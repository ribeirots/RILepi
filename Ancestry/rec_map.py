# script to calculate new recombination map for panel files - ancestry hmm pipeline
# it only runs with one chrm arm
# tribeiro@wisc.edu

# Arg1: recombination map, tab delimited, with chrm\tstart\tend\trec_perbase eg: X 1 100000 0.001
# Arg2: panel file output from populate perl script

# the output will use the same name of the panel file, but with 2 added to the extension.

import re, sys, math

# First, make a dictionary for the chrm arm being analyzed
# here, if recombination rate in a window is 0, it is replaced with 1.0E-09 to ensure that two sites are not located "on top of each other".
recmap_dict = {}
with open(sys.argv[1]) as reccomeron:
    for r in reccomeron:
        r = re.split('\t',r[:-1])
        if r[3] == '0':
            r[3] = '1.0E-010'
        recmap_dict[r[2]] = r[3]

# this function will receive 2 positions and a recombination map and 
# it returns the distance between them. Rec map set as default.
def get_distance(site1, site2, recomb_map=recmap_dict):
    window_s1 = str(math.ceil(int(site1)/100000)*100000)
    window_s2 = str(math.ceil(int(site2)/100000)*100000)
    if window_s1 == window_s2: # are the two sites in the same "100kbp recombination rate window"?
        recrate = float(recmap_dict[window_s1])
        distance = site2 - site1
        return recrate*distance
    elif int(window_s1) == (int(window_s2) - 100000): # are the two sites in windows right next to each other?
        dist1 = int(window_s1) - int(site1) + 1 # adding 1, it is the distance from this window to the next.
        dist2 = int(site2) - ( int(window_s1) + 1)

        recdist1 = dist1 * float(recmap_dict[window_s1])
        recdist2 = dist2 * float(recmap_dict[window_s2])
        return recdist1 + recdist2
    else: # if there are windows between the ones in which the sites are located, they need to be accounted
        dist1 = int(window_s1) - int(site1) + 1
        dist2 = int(site2) - (int(window_s2) - 99999)
        dist_inter = 0

        recdist1 = dist1 * float(recmap_dict[window_s1])
        recdist2 = dist2 * float(recmap_dict[window_s2])

        window_inter = int(window_s1) + 100000
        while window_inter < int(window_s2):
            dist_inter = 100000 * float(recmap_dict[str(window_inter )]) + dist_inter
            window_inter = window_inter + 100000
        
        return recdist1 + recdist2 + dist_inter


######### Read through the file to be updated

output = open(sys.argv[2]+'2','w')

previous_site = 1 # the first recombination distance is not used by AncHMM, so the initial previous site doesnt matter
with open(sys.argv[2]) as panelhmm:
    for r in panelhmm:
        r = re.split('\t',r[:-1])
        r[6] = str(get_distance(previous_site,int(r[1])))
        previous_site = int(r[1])
        output.write('\t'.join(r)+'\n')
output.close()
