#!/bin/bash
# Script for the RILs.
# Script starting with bam files with duplicated removed. The file has to be labeled FILEID_rmdup.bam
# FILEID is the input argument that the script takes
# The final output is a bam file with indels realigned

# It requires picard and GATK, specific versions can be seen in the commands below.
# It requires Dmel reference fasta

# creates tmp/ directory if one does not already exists.
if [ ! -d tmp/ ]
then
	mkdir tmp
fi

java -Xmx5g -Djava.io.tmpdir=./tmp -jar ../alignment_software/picard-tools-1.79/picard-tools-1.79/AddOrReplaceReadGroups.jar RGLB=Lane1 RGPL=Illumina RGPU=TTAGGC RGSM=$1 INPUT=${1}_rmdup.bam OUTPUT=${1}_header.bam

java -Xmx5g -Djava.io.tmpdir=./tmp -jar ../alignment_software/picard-tools-1.79/picard-tools-1.79/BuildBamIndex.jar INPUT=${1}_header.bam

java -Xmx5g -Djava.io.tmpdir=./tmp -jar ../alignment_software/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ../alignment_software/dmel_ref/DmelRef.fasta -I ${1}_header.bam -o ${1}.intervals 2> ${1}_realigntarget.log

java -Xmx5g -Djava.io.tmpdir=./tmp -jar ../alignment_software/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T IndelRealigner -R ../alignment_software/dmel_ref/DmelRef.fasta -targetIntervals ${1}.intervals -I ${1}_header.bam -o ${1}_realign.bam 2> ${1}_indelrealign.log

rm ${1}_header.bam
rm ${1}.intervals
