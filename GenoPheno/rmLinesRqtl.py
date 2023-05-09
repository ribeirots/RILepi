# script to remove selected lines from rqtl input
# tribeiro@wisc.edu

# arg1: input file
# arg2: output file
import re, sys

lines = ["F198R", "F310R", "F385R", "F333R", "F316R"] # removed due to high mismatch and switches

#lines = ["F120R", "F134R", "F159R","F161R", "F165R", "F175R", "F181R", "F213R", "F220R", "F221R", "F284R", "F294R", "F312R", "F329R", "F351R", "F352R", "F371R", "F54R", "F71R", "F95R", "F245R"] # removed due to not having a match phenotype for Thermal Prefence (Marco Gallio's Lab)

output = open(sys.argv[2],'w')
with open(sys.argv[1]) as originalInput:
    for r in originalInput:
        rsplit = re.split(',',r)
        if rsplit[0] not in lines:
            output.write(r)
output.close()
