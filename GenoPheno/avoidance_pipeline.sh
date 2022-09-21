#!/bin/bash

python3 input_rqtl_nopheno.py FR_RIL_genotypes_winZI1k.txt

python3 input_rqtl_addpheno.py FR_RIL_genotypes_winZI1k_rqtl.csv avoidance_index_pheno.csv avoid FR_avoid.csv

# python3 input_rqtl_addpheno.py FR_RIL_genotypes_winZI1k_rqtl.csv [phenofile] [phenoword] [output]

python3 missing_genotypes.py FR_avoid.csv FR_avoid_gap10.csv

python3 window_merge_genotype.py FR_avoid_gap10.csv # not rm anything anymore

python3 comeron_rqtl.py FR_avoid_gap10_merged.csv

python3 split_in_7.py FR_avoid_gap10_merged_cMorg.csv

# final output file: FR_avoid_gap10_merged_cMorg_7armsExtended.csv
