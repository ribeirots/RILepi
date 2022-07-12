#!/bin/bash
# Script to process block bam files after they were mapped on CHTC.
# The argument this script takes is the sample ID, everything before "_Block*.bam" for each sample.

samtools merge ${1}_merged.bam ${1}_Block*.bam

samtools markdup -r -d 2500 -f stats.txt ${1}_merged.bam ${1}_rmdup.bam

rm ${1}_Block*.bam
rm ${1}_merged.bam
