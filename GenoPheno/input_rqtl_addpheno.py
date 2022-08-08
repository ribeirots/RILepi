# script to read rqtl input without phenotype and add the phenotype column.
# Next step is to add the recombination row
# tribeiro@wisc.eud

# Arg1: input geno
# Arg2: input pheno
# Arg3: pheno-name: one word
# Arg4: output
import re, sys 
print('rqtl input - add pheno')

geno = []
with open(sys.argv[1]) as g:
    for r in g:
        geno.append(re.split(',',r[:-1]))

pheno = []
with open(sys.argv[2]) as p:
    for r in p:
        r = re.split(',',r[:-1])
#        r[0] = r[0][1:-1] # remove the "F_R" from the ril name
        pheno.append(r)


output = open(sys.argv[4],'w')

# Markers
new_row1 = [geno[0][0]] + [sys.argv[3]] + geno[0][1:] # header

# Chrm. cM
new_row2 = [geno[1][0]] + [''] + geno[1][1:] # chrm
#new_row3 = [geno[2][0]] + [''] + geno[2][1:] # distance


print(len(new_row1))
output.write(','.join(new_row1)+'\n')
output.write(','.join(new_row2)+'\n')
#output.write(','.join(new_row3)+'\n')

for g in geno:
    new_row = []
    ril_check = 0
    for p in pheno:
        if g[0][1:-1] == p[0]: # RIL on g has F_R but on p it is just the number. This will likely need to be modified file by file.
            ril_check = 1
            new_row = [g[0]] + [p[1]] + g[1:]
            output.write(','.join(new_row)+'\n')
    if ril_check == 0:
        print(g[0])

output.close()
