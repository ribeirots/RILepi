#!/bin/bash

#untar installations
tar -xzf bwa.tar.gz
#tar -xzf samtools_1.12.tar.gz
tar -xzf DmelRef.fasta.tar.gz

#run script

gunzip $1_R1.fastq.gz
gunzip $1_R2.fastq.gz 


export PATH="$(pwd)/bwa:$PATH"
#export PATH="$(pwd)/samtools_1.12/bin:$PATH"

bwa index DmelRef.fasta
bwa mem DmelRef.fasta $1_R1.fastq $1_R2.fastq > ${1}_align_indexed.bam

### Remove duplicate sequences
#samtools collate -o bam_collate.bam ${1}_align_indexed.bam
#samtools fixmate -m bam_collate.bam bam_fixmate.bam
#samtools sort -o bam_position.bam bam_fixmate.bam
#samtools markdup -r -d 2500 -f stats.txt bam_position.bam ${1}_rmdup.bam

#rm bam_collate.bam
#rm bam_fixmate.bam
#rm bam_position.bam


rm $1_R1.fastq
rm $1_R2.fastq
rm DmelRef.fast*

