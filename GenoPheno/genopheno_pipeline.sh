#!/bin/bash

# arguments:
# 1: input file from ancestry calls per window following AHMM, each row is a window each columns is a RIL
# 2: phenotype data to be used
# 3: label for output filename

# Example cmd usage: bash ./genopheno_pipeline.sh ancprop_nofix_FR_RIL_genotypes_winZI1k.txt phenoAvoidanceMay2023.csv FR_thermal

python3 input_rqtl_nopheno.py $1

rqtlF=${1%????}

python3 rmLinesRqtl.py ${rqtlF}_rqtl.csv genotypes_postRILremoval.csv

python3 input_rqtl_addpheno.py genotypes_postRILremoval.csv $2 pheno ${3}.csv

python3 missing_genotypes.py ${3}.csv ${3}_gap10.csv

python3 window_merge_genotype.py ${3}_gap10.csv

python3 comeron_rqtl.py ${3}_gap10_merged.csv

# final output file: ${3}_gap10_merged_cMorg.csv
# also include a final file called "genotypes_postRILremoval.csv" with the genotypes used in the current analysis.

rm ${3}.csv
rm *rqtl.csv
rm *gap10.csv
rm *merged.csv
