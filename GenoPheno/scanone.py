# Script to calculate scantwo() - designed for parallelizations, one window1 at a time
# Author: tiaaagosr@gmail.com
# arg1: input rqtl data
# arg2: pheno column
# arg3: column of window1 with genotype
# arg4: number of permutations (0 for no permutation and an analysis of the full dataset)
# arg5 (optional): random seed
# e.g. command use: python3 scantwo_full.py FR_avoid_gap10_merged_cMorg_7armsExtended.csv 1 2 2

print('On @marula, conda activate sklearnpy\n')
print('This script has been modified to run the dominance model, apply changes to the "scanone" function called - and not the function itself - if you would like to undo it.\nRemember that the function is called twice, once in the block for permutations and once in the block for the empirical data.\n')

import re, sys
import numpy as np
from sklearn.linear_model import LinearRegression
import statsmodels.formula.api as smf
import pandas as pd
import random, math
from scipy import stats

def scanone(pheno, l1, l2, method=0):
    boxcox_modeling = method # 0: pheno as inputed. 1: output boxcox correction
    if method == 'dominance':
        # next block converts to dominance model, Ht == 2nd parental, 1 = 2
        ancestry1 = 0
        ancestry0 = 0
        for v in range(0,len(l1)):
            if l1[v] == '1':
                ancestry1 += 1
                l1[v] = '2'
                l2[v] = 2
            elif l1[v] == '0':
                ancestry0 += 1
        # end of dominance code

    if boxcox_modeling:
        fitted_data, fitted_lambda = stats.boxcox(pheno)
        d = {'y': pheno,'yboxcox': fitted_data, 'locus1': l1, 'locus2':l2}
    else:
        d = {'y': pheno, 'locus1': l1, 'locus2':l2}
    df = pd.DataFrame(data=d).dropna() # will drop rows (RILs) with missing data in at least one of the windows being compared
    if method == 'dominance':
        zAncestry = (ancestry0*2 + ancestry1 )/(len(df.index)*2) # needed for dominance model
    else:
        zAncestry = ((df.locus2.value_counts()[0]*2) + (df.locus2.value_counts()[1]) )/(len(df.index)*2)

    if len(df) < (len(pheno)/2): # if more than half the RILs were dropped, skips this window comparison
        print(len(df))
        return ['too many NaN']
    elif zAncestry >= 0.9:
        return [-0.25, -0.25, -0.25, -0.25, -0.25, -0.25]

    
    if boxcox_modeling:
        mod_f = smf.ols(formula='yboxcox ~ locus1', data=df) # Full Model
        res_f = mod_f.fit()

        mod_f_add = smf.ols(formula='yboxcox ~ locus2', data=df) # Full Model
        res_f_add = mod_f_add.fit()
        
        mod_0 = smf.ols(formula='yboxcox ~ 1', data=df) # Null Model
        res_0 = mod_0.fit()

        lod = (len(df.index)/2)*math.log((res_0.ssr/res_f.ssr),10)
        lod_add =  (len(l1)/2)*math.log((res_0.ssr/res_f_add.ssr),10)
        
        f_f = res_f.fvalue # F-stats
        f_f_add = res_f_add.fvalue # F-stats
        
        zz_prop = l1.count('0')/len(df.index)
    else:
        # categorical
        mod_f = smf.ols(formula='y ~ locus1', data=df) # Full Model
        res_f = mod_f.fit()

        mod_f_add = smf.ols(formula='y ~ locus2', data=df) # Full Model
        res_f_add = mod_f_add.fit()
        
        mod_0 = smf.ols(formula='y ~ 1', data=df) # Null Model
        res_0 = mod_0.fit()

        lod = (len(df.index)/2)*math.log((res_0.ssr/res_f.ssr),10)
        lod_add =  (len(l1)/2)*math.log((res_0.ssr/res_f_add.ssr),10)
        
        f_f = res_f.fvalue # F-stats
        f_f_add = res_f_add.fvalue # F-stats
        
        zz_prop = df.locus2.value_counts()[0]/len(df.index)


    return [lod, lod_add, f_f, f_f_add, zz_prop, len(df.index)]

perm = int(sys.argv[4]) # input 0 for False (running no permutation)

geno1 = 'ZZ' ## coded as 0
geno2 = 'FF' ## coded as 2
ht_geno = 'FZ' ## coded as 1. Missing data coded as NaN

with open(sys.argv[1]) as ril_data:
    markers = ril_data.readline()
    markers = re.split(r',',markers[:-1])
    next(ril_data)
    next(ril_data)
    p = []
    for line in ril_data:
        line = re.split(r',',line[:-1])
        p.append(float(line[int(sys.argv[2])]))
  
first_row_to_analyze = int(sys.argv[3])
#last_row_to_analyze = first_row_to_analyze + 1
#last_row_to_analyze = 6
last_row_to_analyze = len(markers)

nseed='all'
if perm > 0:
    if len(sys.argv) == 6:
        nseed = int(sys.argv[5])
        random.seed(251654+nseed)
    outfilename = sys.argv[1][:-4]+"_scanone_perm_"+str(nseed)+"_"+str(first_row_to_analyze)+sys.argv[1][-4:]
    perm_scanone = open(outfilename,"w")

    header = ["Marker1","LOD", "LOD_num", "f_f","f_f_num" ,"ZZ_proportion", "n","ZZ%_num","n_num", "Marker1_num"]  # add = num ZZ -> FZ -> FF with numbers instead of categories
    header = "\t".join(header) + "\n"
    perm_scanone.write(header)

    for permuting in range(0,perm):
        random.shuffle(p)
        if len(sys.argv) != 5 and len(sys.argv) != 6:
            print('Wrong input command.')
            quit()
        max_lods = []

        lod_list = []
        for i in range(first_row_to_analyze,last_row_to_analyze):
            for j in [0]:
                locus1 = []
                locus2 = []
                with open(sys.argv[1]) as ril_data:
                    next(ril_data)
                    next(ril_data)
                    next(ril_data)
                    for line in ril_data:
                        line = re.split(r',',line[:-1])
                        
                        if line[i] == geno1:
                            locus1.append('0')
                        elif line[i] == geno2:
                            locus1.append('2')
                        elif line[i] == ht_geno:
                            locus1.append('1')
                        else:
                            locus1.append(np.nan)
                                
                        if line[i] == geno1:
                            locus2.append(0)
                        elif line[i] == geno2:
                            locus2.append(2)
                        elif line[i] == ht_geno:
                            locus2.append(1)
                        else:
                            locus2.append(np.nan)

                scan_results = scanone(p,locus1,locus2, method="dominance")
                if len(scan_results) > 1: # if == 1, too many NA, skip the window comparison
                    one_window = [i]+scan_results # storing i and j instead of marker name for outputting purposes
                    lod_list.append(one_window)
        lod_perm = pd.DataFrame(data=lod_list)
        maxperm = lod_perm.max().tolist()
        lod_maxid = lod_perm.idxmax()[1]
        lod_add_maxid = lod_perm.idxmax()[2]
        maxperm[0] = markers[lod_maxid] # window1 from max lodf
        maxperm[5] = lod_list[lod_maxid][5]
        maxperm[6] = lod_list[lod_maxid][6]
        maxperm.append(lod_list[lod_add_maxid][5]) # zz%
        maxperm.append(lod_list[lod_add_maxid][6]) # n
        maxperm.append(markers[lod_add_maxid]) # add window1 from max glodf to the end
        maxperm = list(map(str,maxperm))
        perm_scanone.write('\t'.join(maxperm)+'\n')

    perm_scanone.close()


### code for no permutation, full data set analysis
if perm == 0:
    print('Running script for the full data set, without permutations.')
    outfilename = sys.argv[1][:-4]+"_scanone_"+str(first_row_to_analyze)+sys.argv[1][-4:]
    full_scantwo = open(outfilename,"w")

    #g = GLM
    header = ["Marker1","LOD", "LOD_num", "f_f","f_f_num" ,"ZZ_proportion", "n"]  # add = num ZZ -> FZ -> FF with numbers instead of categories
    header = "\t".join(header) + "\n"
    full_scantwo.write(header)

    lod_list = []
    for i in range(first_row_to_analyze,last_row_to_analyze):
        for j in [0]:
            locus1 = []
            locus2 = []
            marker1 = markers[i]
            with open(sys.argv[1]) as ril_data:
                next(ril_data)
                next(ril_data)
                next(ril_data)
                for line in ril_data:
                    line = re.split(r',',line[:-1])
                    
                    if line[i] == geno1:
                        locus1.append('0')
                    elif line[i] == geno2:
                        locus1.append('2')
                    elif line[i] == ht_geno:
                        locus1.append('1')
                    else:
                        locus1.append(np.nan)
                            
                    if line[i] == geno1:
                        locus2.append(0)
                    elif line[i] == geno2:
                        locus2.append(2)
                    elif line[i] == ht_geno:
                        locus2.append(1)
                    else:
                        locus2.append(np.nan)
            scan_results = scanone(p,locus1,locus2, method="dominance")
            if len(scan_results) > 1: # if == 1, too many NA, skip the window comparison
                one_window = [marker1]+scan_results
                lod_list.append(scan_results)
                one_window = list(map(str, one_window))
                one_window = "\t".join(one_window) + "\n"
                full_scantwo.write(one_window)


    full_scantwo.close()


if 0: # change to 1 to get screen average outputs
    lod_avgs = np.asarray(lod_list)
    print('mean LODf')
    print(sum(lod_avgs[:,0])/len(lod_avgs[:,0]))

    print('mean log-lik_f')
    print(sum(lod_avgs[:,4])/len(lod_avgs[:,4]))

    print('mean gLODf')
    print(sum(lod_avgs[:,13])/len(lod_avgs[:,13]))

    print('mean glog-lik_f')
    print(sum(lod_avgs[:,17])/len(lod_avgs[:,17]))

